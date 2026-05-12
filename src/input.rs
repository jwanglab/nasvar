//! Unified alignment input that supports both BAM and CRAM files via noodles.
//!
//! Supports transparent multi-BAM input: if the path ends in `.txt` or `.list`,
//! it is read as a file-of-filenames (one BAM/CRAM path per line). All readers
//! share the same reference sequences; headers are validated on open.

use std::fs::File;
use std::io::BufReader;
use anyhow::{Result, bail};
use log::info;

use noodles::bam;
use noodles::bgzf;
#[cfg(feature = "cram")]
use noodles::cram;
use noodles::csi::binning_index::{
    Indexer as CsiIndexer,
    index::reference_sequence::{bin::Chunk, index::LinearIndex},
};
use noodles::fasta;
pub use noodles::sam;
use noodles::core::Region;

use crate::bam::ContigMapper;

pub use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

/// Simplified header info extracted from noodles sam::Header.
#[derive(Debug, Clone)]
pub struct AlignmentHeader {
    pub text: String,
    pub refs: Vec<String>,
    pub lengths: Vec<i32>,
}

impl AlignmentHeader {
    pub fn from_sam_header(header: &sam::Header) -> Self {
        let mut refs = Vec::new();
        let mut lengths = Vec::new();
        let ref_seqs = header.reference_sequences();
        for (name, map) in ref_seqs.iter() {
            refs.push(String::from_utf8_lossy(name).to_string());
            lengths.push(map.length().get() as i32);
        }
        // Build SAM header text manually
        let mut text = String::new();
        // @HD line
        if let Some(hd) = header.header() {
            text.push_str("@HD\tVN:");
            text.push_str(&hd.version().to_string());
            text.push('\n');
        }
        // @SQ lines
        for (name, map) in ref_seqs.iter() {
            text.push_str("@SQ\tSN:");
            text.push_str(&String::from_utf8_lossy(name));
            text.push_str("\tLN:");
            text.push_str(&map.length().to_string());
            text.push('\n');
        }
        // @RG lines - preserve all fields
        for (id, rg) in header.read_groups().iter() {
            text.push_str("@RG\tID:");
            text.push_str(&String::from_utf8_lossy(id));
            // Add other fields from the read group
            for (tag, value) in rg.other_fields().iter() {
                text.push('\t');
                text.push_str(&String::from_utf8_lossy(tag.as_ref()));
                text.push(':');
                text.push_str(&String::from_utf8_lossy(value));
            }
            text.push('\n');
        }
        // @CO lines
        for comment in header.comments() {
            text.push_str("@CO\t");
            text.push_str(&String::from_utf8_lossy(comment));
            text.push('\n');
        }
        AlignmentHeader { text, refs, lengths }
    }
}

/// A record wrapper that provides a uniform interface over noodles BAM/CRAM records.
#[derive(Debug, Clone)]
pub struct AlignmentRecord {
    pub name: Option<String>,
    pub ref_id: i32,
    pub pos: i32, // 0-based position
    pub flag: u16,
    pub mapq: u8,
    pub seq: Vec<u8>, // ASCII bases
    pub qual: Vec<u8>,
    pub cigar: Vec<(CigarKind, usize)>, // decoded cigar ops
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub tlen: i32,
}

impl Default for AlignmentRecord {
    fn default() -> Self {
        AlignmentRecord {
            name: None,
            ref_id: -1,
            pos: -1,
            flag: 0,
            mapq: 0,
            seq: Vec::new(),
            qual: Vec::new(),
            cigar: Vec::new(),
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }
    }
}

impl AlignmentRecord {
    /// Returns the read name, if present.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Returns the 1-based alignment start position, or None if unmapped.
    pub fn alignment_start(&self) -> Option<usize> {
        if self.pos >= 0 {
            Some((self.pos + 1) as usize)
        } else {
            None
        }
    }

    /// Returns the alignment span on the reference (sum of M/D/N/=/X ops).
    pub fn alignment_span(&self) -> usize {
        let mut span = 0;
        for &(kind, len) in &self.cigar {
            match kind {
                CigarKind::Match
                | CigarKind::Deletion
                | CigarKind::Skip
                | CigarKind::SequenceMatch
                | CigarKind::SequenceMismatch => span += len,
                _ => {}
            }
        }
        span
    }

    /// Returns the raw flag bits.
    pub fn flags(&self) -> u16 {
        self.flag
    }

    /// Returns the decoded sequence as a byte slice.
    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    /// Iterates over CIGAR operations as (Kind, length) pairs.
    pub fn cigar_ops(&self) -> &[(CigarKind, usize)] {
        &self.cigar
    }

    /// Returns the raw encoded CIGAR as u32 values (for compatibility with fix_sam_coords).
    pub fn cigar_raw(&self) -> Vec<u32> {
        self.cigar.iter().map(|&(kind, len)| {
            let op_code: u32 = match kind {
                CigarKind::Match => 0,
                CigarKind::Insertion => 1,
                CigarKind::Deletion => 2,
                CigarKind::Skip => 3,
                CigarKind::SoftClip => 4,
                CigarKind::HardClip => 5,
                CigarKind::Pad => 6,
                CigarKind::SequenceMatch => 7,
                CigarKind::SequenceMismatch => 8,
            };
            ((len as u32) << 4) | op_code
        }).collect()
    }
}

/// Decode a noodles BAM record into our AlignmentRecord.
fn decode_bam_record(rec: &bam::Record, _header: &sam::Header) -> Result<AlignmentRecord> {
    let name = rec.name().map(|n| {
        String::from_utf8_lossy(n.as_ref()).to_string()
    });

    let ref_id = match rec.reference_sequence_id() {
        Some(Ok(id)) => id as i32,
        Some(Err(e)) => return Err(e.into()),
        None => -1,
    };

    let pos = match rec.alignment_start() {
        Some(Ok(p)) => (p.get() as i32) - 1, // convert 1-based to 0-based
        Some(Err(e)) => return Err(e.into()),
        None => -1,
    };

    let flag = rec.flags().bits();
    let mapq = match rec.mapping_quality() {
        Some(q) => q.get(),
        None => 255,
    };

    // Decode sequence
    let seq: Vec<u8> = rec.sequence().iter().collect();

    // Decode quality
    let qual: Vec<u8> = rec.quality_scores().iter().collect();

    // Decode CIGAR
    let mut cigar = Vec::new();
    for op_result in rec.cigar().iter() {
        let op = op_result?;
        cigar.push((op.kind(), op.len()));
    }

    let next_ref_id = match rec.mate_reference_sequence_id() {
        Some(Ok(id)) => id as i32,
        _ => -1,
    };
    let next_pos = match rec.mate_alignment_start() {
        Some(Ok(p)) => (p.get() as i32) - 1,
        _ => -1,
    };
    let tlen = rec.template_length();

    Ok(AlignmentRecord {
        name,
        ref_id,
        pos,
        flag,
        mapq,
        seq,
        qual,
        cigar,
        next_ref_id,
        next_pos,
        tlen,
    })
}

/// Decode any type implementing the noodles alignment Record trait into AlignmentRecord.
pub fn decode_alignment_record(rec: &dyn sam::alignment::Record, header: &sam::Header) -> Result<AlignmentRecord> {
    let name = rec.name().map(|n| {
        String::from_utf8_lossy(n.as_ref()).to_string()
    });

    let ref_id = match rec.reference_sequence_id(header) {
        Some(Ok(id)) => id as i32,
        Some(Err(e)) => return Err(e.into()),
        None => -1,
    };

    let pos = match rec.alignment_start() {
        Some(Ok(p)) => (p.get() as i32) - 1,
        Some(Err(e)) => return Err(e.into()),
        None => -1,
    };

    let flag = rec.flags()?.bits();
    let mapq = match rec.mapping_quality() {
        Some(Ok(q)) => q.get(),
        Some(Err(e)) => return Err(e.into()),
        None => 255,
    };

    let seq: Vec<u8> = rec.sequence().iter().collect();
    let qual: Vec<u8> = rec.quality_scores().iter()
        .map(|r| r.unwrap_or(255))
        .collect();

    let mut cigar = Vec::new();
    for op_result in rec.cigar().iter() {
        let op = op_result?;
        cigar.push((op.kind(), op.len()));
    }

    let next_ref_id = match rec.mate_reference_sequence_id(header) {
        Some(Ok(id)) => id as i32,
        _ => -1,
    };
    let next_pos = match rec.mate_alignment_start() {
        Some(Ok(p)) => (p.get() as i32) - 1,
        _ => -1,
    };
    let tlen = rec.template_length()?;

    Ok(AlignmentRecord {
        name,
        ref_id,
        pos,
        flag,
        mapq,
        seq,
        qual,
        cigar,
        next_ref_id,
        next_pos,
        tlen,
    })
}

/// Feed one alignment into the inline indexer. Also enforces coordinate
/// sortedness: if the BAM is out of (ref_id, pos) order we bail, since a
/// BAI built from an unsorted BAM is useless.
fn update_inline_index(
    indexer: &mut CsiIndexer<LinearIndex>,
    chunk: Chunk,
    record: &AlignmentRecord,
    last_seen: &mut Option<(i32, i32)>,
    n_refs_seen: &mut usize,
) -> Result<()> {
    use noodles::core::Position;

    // Unmapped / no reference: passed to indexer as None.
    if record.ref_id < 0 || record.pos < 0 {
        indexer.add_record(None, chunk).map_err(|e| anyhow::anyhow!("{}", e))?;
        return Ok(());
    }

    // Enforce (ref_id, pos) non-decreasing. The csi indexer itself also
    // rejects ref_id regressions, but surfacing it here gives a clearer
    // "BAM must be sorted" message.
    if let Some((prev_ref, prev_pos)) = *last_seen
        && (record.ref_id < prev_ref
            || (record.ref_id == prev_ref && record.pos < prev_pos)) {
        bail!(
            "BAM is not coordinate-sorted (saw {}:{} after {}:{}); \
             cannot build an index on the fly. Pre-sort the BAM or \
             provide a .bai.",
            record.ref_id, record.pos, prev_ref, prev_pos,
        );
    }
    *last_seen = Some((record.ref_id, record.pos));

    // Coordinates for BAI/CSI are 1-based inclusive Position.
    let start_1b = (record.pos as usize) + 1;
    let span = record.alignment_span().max(1);
    let end_1b = start_1b + span - 1;
    let start = Position::try_from(start_1b)
        .map_err(|e| anyhow::anyhow!("invalid start position: {}", e))?;
    let end = Position::try_from(end_1b)
        .map_err(|e| anyhow::anyhow!("invalid end position: {}", e))?;

    let is_mapped = (record.flag & 0x4) == 0;
    let ref_id = record.ref_id as usize;
    *n_refs_seen = (*n_refs_seen).max(ref_id + 1);

    indexer
        .add_record(Some((ref_id, start, end, is_mapped)), chunk)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

/// Inner reader enum
enum Inner {
    Bam(bam::io::IndexedReader<bgzf::io::Reader<File>>),
    BamNoIndex(bam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    #[cfg(feature = "cram")]
    Cram(cram::io::IndexedReader<File>),
    #[cfg(feature = "cram")]
    CramNoIndex(cram::io::Reader<BufReader<File>>),
    Sam(sam::io::Reader<BufReader<File>>),
}

/// Per-file reader state for multi-BAM support.
struct ReaderState {
    inner: Inner,
    file_path: String,
    #[cfg_attr(not(feature = "cram"), allow(dead_code))]
    fasta_repo: fasta::Repository,
    #[cfg(feature = "cram")]
    cram_record_buf: std::collections::VecDeque<AlignmentRecord>,
    /// File size in bytes (cached on open for progress tracking)
    file_size: u64,
    /// When Some, we're accumulating a BAI in memory during sequential reads
    /// of a sorted BAM that lacked an index file. After the first pass,
    /// `finalize_inline_index` converts this into an IndexedReader so
    /// subsequent `query()` calls have random access.
    inline_indexer: Option<CsiIndexer<LinearIndex>>,
    /// (ref_id, pos) of the most-recent record we indexed, used to detect
    /// sortedness violations while building the index.
    inline_last_seen: Option<(i32, i32)>,
    /// Total aligned record count (across all reference sequences); used to
    /// size the indexer's reference sequence table at build time.
    inline_n_refs_seen: usize,
    /// Cached copy of the finalized inline BAI. When Some, seek_to_start
    /// uses it to re-open as IndexedReader so the index survives any
    /// `seek(0)` call (e.g. fusion pass 2's bam.seek(bam.start_pos)).
    inline_built_index: Option<bam::bai::Index>,
}

impl ReaderState {
    /// Enable on-the-fly BAI building during the next sequential scan.
    /// Only meaningful for a BAM that was opened without an index (i.e.
    /// currently `Inner::BamNoIndex`).
    fn enable_inline_index(&mut self) {
        if matches!(self.inner, Inner::BamNoIndex(_)) && self.inline_indexer.is_none() {
            // LinearIndex with default (min_shift=14, depth=5) is the BAI
            // parameter set, so the built index is interchangeable with a
            // real `.bai`.
            self.inline_indexer = Some(CsiIndexer::<LinearIndex>::default());
            self.inline_last_seen = None;
            self.inline_n_refs_seen = 0;
        }
    }

    /// Finalize the in-progress inline BAI (if any). If one was being built,
    /// reopen the underlying BAM wrapped in an IndexedReader that uses it,
    /// so subsequent `query_inner` calls have random access. Idempotent:
    /// no-op if no index was being built.
    fn finalize_inline_index(&mut self, header_n_refs: usize) -> Result<bool> {
        let Some(indexer) = self.inline_indexer.take() else {
            return Ok(false);
        };
        if !matches!(self.inner, Inner::BamNoIndex(_)) {
            // Shouldn't happen — enable_inline_index guards against this.
            return Ok(false);
        }
        let n_refs = header_n_refs.max(self.inline_n_refs_seen).max(1);
        let index: bam::bai::Index = indexer.build(n_refs);
        info!(
            "Built inline BAI for {} ({} reference sequences indexed; pass 1 saw reads on {})",
            self.file_path, n_refs, self.inline_n_refs_seen
        );

        // Cache the index so subsequent seek_to_start calls can re-use it
        // (otherwise they'd reopen the BAM without an index, losing it).
        self.inline_built_index = Some(index.clone());

        // Also persist the BAI alongside the BAM so external consumers
        // (e.g. an embedded genome browser in the hosting page) can pick
        // it up from the VFS without having to re-scan the BAM. Best-effort
        // — a write failure shouldn't kill the pipeline.
        let bai_path = format!("{}.bai", self.file_path);
        match bam::bai::fs::write(&bai_path, &index) {
            Ok(()) => info!("Wrote inline BAI to {}", bai_path),
            Err(e) => log::warn!("Failed to write inline BAI to {}: {}", bai_path, e),
        }

        // Reopen the BAM path with an IndexedReader backed by the in-memory
        // index. The underlying VFS (native fs, OPFS SyncAccessHandle, or
        // the SAB-bridge Fd in the web worker) already supports random
        // access via fd_seek, so the query() path just works.
        let file = File::open(&self.file_path)
            .map_err(|e| anyhow::anyhow!("Failed to reopen BAM {}: {}", self.file_path, e))?;
        let mut reader = bam::io::IndexedReader::new(file, index);
        reader.read_header()?;
        self.inner = Inner::Bam(reader);
        info!("Reader for {} upgraded to indexed mode", self.file_path);
        Ok(true)
    }

    /// Read the next record from this reader.
    fn read_record_inner(&mut self, sam_header: &sam::Header) -> Result<Option<AlignmentRecord>> {
        // For CRAM, use the record buffer
        #[cfg(feature = "cram")]
        if matches!(self.inner, Inner::Cram(_) | Inner::CramNoIndex(_)) {
            if self.cram_record_buf.is_empty() {
                self.fill_cram_buffer(sam_header)?;
                if self.cram_record_buf.is_empty() {
                    return Ok(None); // EOF
                }
            }
            return Ok(self.cram_record_buf.pop_front());
        }

        match &mut self.inner {
            Inner::Bam(r) => {
                let mut buf = bam::Record::default();
                match r.read_record(&mut buf) {
                    Ok(0) => Ok(None),
                    Ok(_) => Ok(Some(decode_bam_record(&buf, sam_header)?)),
                    Err(e) => Err(e.into()),
                }
            }
            Inner::BamNoIndex(r) => {
                let chunk_start = self.inline_indexer.as_ref()
                    .map(|_| r.get_ref().virtual_position());
                let mut buf = bam::Record::default();
                match r.read_record(&mut buf) {
                    Ok(0) => Ok(None),
                    Ok(_) => {
                        let record = decode_bam_record(&buf, sam_header)?;
                        // Update the in-memory index with this record's chunk,
                        // if inline-index-build mode is active.
                        if let Some(indexer) = self.inline_indexer.as_mut() {
                            let chunk_end = r.get_ref().virtual_position();
                            let chunk = Chunk::new(chunk_start.unwrap(), chunk_end);
                            update_inline_index(
                                indexer, chunk, &record,
                                &mut self.inline_last_seen,
                                &mut self.inline_n_refs_seen,
                            )?;
                        }
                        Ok(Some(record))
                    }
                    Err(e) => Err(e.into()),
                }
            }
            Inner::Sam(r) => {
                let mut buf = sam::Record::default();
                match r.read_record(&mut buf) {
                    Ok(0) => Ok(None),
                    Ok(_) => Ok(Some(decode_alignment_record(&buf, sam_header)?)),
                    Err(e) => Err(e.into()),
                }
            }
            #[cfg(feature = "cram")]
            _ => unreachable!(),
        }
    }

    /// Fill the CRAM record buffer by reading and decoding the next container.
    #[cfg(feature = "cram")]
    fn fill_cram_buffer(&mut self, sam_header: &sam::Header) -> Result<()> {
        let mut container = cram::io::reader::Container::default();

        let bytes_read = match &mut self.inner {
            Inner::Cram(r) => r.read_container(&mut container)?,
            Inner::CramNoIndex(r) => r.read_container(&mut container)?,
            _ => unreachable!(),
        };

        if bytes_read == 0 {
            return Ok(()); // EOF
        }

        let compression_header = container.compression_header()?;

        for slice_result in container.slices() {
            let slice = slice_result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let cram_records = slice.records(
                self.fasta_repo.clone(),
                sam_header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;
            for cram_rec in &cram_records {
                let ar = decode_alignment_record(cram_rec, sam_header)?;
                self.cram_record_buf.push_back(ar);
            }
        }

        Ok(())
    }

    /// Query records overlapping a region from this reader.
    fn query_inner(&mut self, sam_header: &sam::Header, region: &Region) -> Result<Vec<AlignmentRecord>> {
        match &mut self.inner {
            Inner::Bam(r) => {
                r.query(sam_header, region)?
                    .records()
                    .map(|result| -> Result<AlignmentRecord> {
                        let rec = result?;
                        decode_bam_record(&rec, sam_header)
                    })
                    .collect::<Result<Vec<_>>>()
            }
            #[cfg(feature = "cram")]
            Inner::Cram(r) => {
                r.query(sam_header, region)?
                    .map(|result| -> Result<AlignmentRecord> {
                        let rec = result?;
                        decode_alignment_record(&rec, sam_header)
                    })
                    .collect::<Result<Vec<_>>>()
            }
            _ => bail!("Region queries require an indexed BAM/CRAM file"),
        }
    }

    /// Get current virtual file offset for progress tracking.
    fn tell_inner(&self) -> u64 {
        match &self.inner {
            Inner::Bam(r) => r.get_ref().virtual_position().compressed(),
            Inner::BamNoIndex(r) => r.get_ref().virtual_position().compressed(),
            #[cfg(feature = "cram")]
            Inner::Cram(_) | Inner::CramNoIndex(_) => 0,
            Inner::Sam(_) => 0,
        }
    }

    fn is_sam(&self) -> bool {
        matches!(self.inner, Inner::Sam(_))
    }

    /// Reset this reader to the beginning of the file.
    fn seek_to_start(&mut self) -> Result<()> {
        let path = self.file_path.clone();
        if self.is_sam() {
            let file = File::open(&path)
                .map_err(|e| anyhow::anyhow!("Failed to re-open SAM {}: {}", path, e))?;
            let mut reader = sam::io::Reader::new(BufReader::new(file));
            reader.read_header()?;
            self.inner = Inner::Sam(reader);
            return Ok(());
        }
        #[cfg(feature = "cram")]
        if self.is_cram() {
            self.cram_record_buf.clear();
            let repo = self.fasta_repo.clone();
            let index_path = format!("{}.crai", path);
            if std::path::Path::new(&index_path).exists() {
                let mut reader = cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_path(&path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open CRAM {}: {}", path, e))?;
                reader.read_header()?;
                self.inner = Inner::Cram(reader);
            } else {
                let file = File::open(&path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open CRAM {}: {}", path, e))?;
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_reader(BufReader::new(file));
                reader.read_header()?;
                self.inner = Inner::CramNoIndex(reader);
            }
            return Ok(());
        }
        let index_path = format!("{}.bai", path);
        if std::path::Path::new(&index_path).exists() {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(&path)
                .map_err(|e| anyhow::anyhow!("Failed to re-open BAM {}: {}", path, e))?;
            reader.read_header()?;
            self.inner = Inner::Bam(reader);
        } else if let Some(cached) = self.inline_built_index.as_ref() {
            // An inline BAI was built during a previous pass. Re-open
            // with it so seek(0) doesn't downgrade us back to BamNoIndex
            // and break subsequent `query()` calls.
            let file = File::open(&path)
                .map_err(|e| anyhow::anyhow!("Failed to re-open BAM {}: {}", path, e))?;
            let mut reader = bam::io::IndexedReader::new(file, cached.clone());
            reader.read_header()?;
            self.inner = Inner::Bam(reader);
        } else {
            let file = File::open(&path)
                .map_err(|e| anyhow::anyhow!("Failed to re-open BAM {}: {}", path, e))?;
            let mut reader = bam::io::Reader::new(BufReader::new(file));
            reader.read_header()?;
            self.inner = Inner::BamNoIndex(reader);
        }
        Ok(())
    }

    #[cfg(feature = "cram")]
    fn is_cram(&self) -> bool {
        matches!(self.inner, Inner::Cram(_) | Inner::CramNoIndex(_))
    }

    fn has_index(&self) -> bool {
        #[cfg(feature = "cram")]
        {
            matches!(self.inner, Inner::Bam(_) | Inner::Cram(_))
        }
        #[cfg(not(feature = "cram"))]
        {
            matches!(self.inner, Inner::Bam(_))
        }
    }
}

/// Unified alignment input wrapping one or more BAM/CRAM files via noodles.
///
/// For multi-BAM input, pass a `.txt` or `.list` file containing one BAM/CRAM
/// path per line. All files must share identical reference sequences.
pub struct AlignmentInput {
    readers: Vec<ReaderState>,
    current_reader_idx: usize,
    sam_header: sam::Header,
    pub header: AlignmentHeader,
    pub contig_mapper: ContigMapper,
    pub start_pos: u64,
    /// Original path passed to open() (list file or single BAM)
    input_path: String,
    /// Count of malformed records skipped (SAM only)
    skipped_records: u64,
}

impl AlignmentInput {
    /// Open an alignment file (BAM or CRAM), auto-detected by magic bytes
    /// (falling back to file extension).
    ///
    /// If `path` ends in `.txt` or `.list`, it is treated as a file-of-filenames
    /// containing one BAM/CRAM path per line. All files must have identical
    /// reference sequences.
    pub fn open(path: &str, ref_path: Option<&str>) -> Result<Self> {
        if path.ends_with(".txt") || path.ends_with(".list") {
            Self::open_multi(path, ref_path)
        } else {
            Self::open_single(path, ref_path)
        }
    }

    /// Open a single alignment file.
    fn open_single(path: &str, ref_path: Option<&str>) -> Result<Self> {
        let fasta_repo = Self::build_fasta_repo(ref_path)?;

        let (reader_state, sam_header, header) = Self::open_any_reader(path, fasta_repo)?;

        let contig_mapper = ContigMapper::from_refs(&header.refs);

        Ok(AlignmentInput {
            readers: vec![reader_state],
            current_reader_idx: 0,
            sam_header,
            header,
            contig_mapper,
            start_pos: 0,
            input_path: path.to_string(),
            skipped_records: 0,
        })
    }

    /// Open multiple alignment files from a file-of-filenames.
    fn open_multi(list_path: &str, ref_path: Option<&str>) -> Result<Self> {
        let content = std::fs::read_to_string(list_path)
            .map_err(|e| anyhow::anyhow!("Failed to read file list '{}': {}", list_path, e))?;

        let paths: Vec<&str> = content.lines()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect();

        if paths.is_empty() {
            bail!("File list '{}' contains no paths", list_path);
        }

        let fasta_repo = Self::build_fasta_repo(ref_path)?;

        // Open the first file to get the canonical header
        let (first_reader, sam_header, canonical_header) = Self::open_any_reader(paths[0], fasta_repo.clone())?;

        let mut readers = vec![first_reader];

        // Open remaining files and validate headers
        for &p in &paths[1..] {
            let (reader_state, _, other_header) = Self::open_any_reader(p, fasta_repo.clone())?;

            Self::validate_headers_match(&canonical_header, &other_header, p)?;
            readers.push(reader_state);
        }

        let contig_mapper = ContigMapper::from_refs(&canonical_header.refs);

        info!("Multi-BAM: opened {} files from {}", readers.len(), list_path);

        Ok(AlignmentInput {
            readers,
            current_reader_idx: 0,
            sam_header,
            header: canonical_header,
            contig_mapper,
            start_pos: 0,
            input_path: list_path.to_string(),
            skipped_records: 0,
        })
    }

    /// Build a FASTA repository from an optional reference path.
    fn build_fasta_repo(ref_path: Option<&str>) -> Result<Option<fasta::Repository>> {
        let Some(rp) = ref_path else { return Ok(None); };

        if !std::path::Path::new(rp).exists() {
            bail!("FASTA reference file not found: {}", rp);
        }
        let fai_path = format!("{}.fai", rp);
        if !std::path::Path::new(&fai_path).exists() {
            bail!(
                "FASTA index not found: {} (FASTA itself exists at {}). \
                 Create one with: samtools faidx {}",
                fai_path, rp, rp
            );
        }

        let indexed_reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(rp)
            .map_err(|e| anyhow::anyhow!("Failed to read FASTA reference {} (index {}): {}", rp, fai_path, e))?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        Ok(Some(fasta::Repository::new(adapter)))
    }

    /// Validate that two headers have identical reference sequences.
    pub fn validate_headers_match(canonical: &AlignmentHeader, other: &AlignmentHeader, other_path: &str) -> Result<()> {
        if canonical.refs.len() != other.refs.len() {
            bail!(
                "BAM {} has {} reference sequences, expected {} (must match first BAM)",
                other_path, other.refs.len(), canonical.refs.len()
            );
        }
        for (i, (a, b)) in canonical.refs.iter().zip(other.refs.iter()).enumerate() {
            if a != b {
                bail!(
                    "BAM {} reference mismatch at index {}: '{}' vs '{}'",
                    other_path, i, a, b
                );
            }
        }
        Ok(())
    }

    /// Detect whether a file is SAM (plain text) by extension.
    fn is_sam_file(path: &str) -> bool {
        path.ends_with(".sam")
    }

    /// Detect whether a file is CRAM by reading the first 4 magic bytes ("CRAM"),
    /// falling back to file extension if the file can't be read.
    #[cfg(feature = "cram")]
    fn is_cram_file(path: &str) -> bool {
        if let Ok(mut f) = File::open(path) {
            let mut magic = [0u8; 4];
            if std::io::Read::read_exact(&mut f, &mut magic).is_ok() {
                return &magic == b"CRAM";
            }
        }
        path.ends_with(".cram")
    }

    /// Open any alignment file (SAM, BAM, or CRAM) based on auto-detection.
    fn open_any_reader(path: &str, fasta_repo: Option<fasta::Repository>) -> Result<(ReaderState, sam::Header, AlignmentHeader)> {
        if Self::is_sam_file(path) {
            return Self::open_sam_reader(path);
        }
        #[cfg(feature = "cram")]
        if Self::is_cram_file(path) {
            return Self::open_cram_reader(path, fasta_repo);
        }
        Self::open_bam_reader(path, fasta_repo)
    }

    /// Open a single SAM file, returning the ReaderState, sam::Header, and AlignmentHeader.
    fn open_sam_reader(path: &str) -> Result<(ReaderState, sam::Header, AlignmentHeader)> {
        let file_size = std::fs::metadata(path)
            .map_err(|e| anyhow::anyhow!("Failed to stat SAM {}: {}", path, e))?
            .len();

        let file = File::open(path)
            .map_err(|e| anyhow::anyhow!("Failed to open SAM {}: {}", path, e))?;
        let mut reader = sam::io::Reader::new(BufReader::new(file));
        let sam_header = reader.read_header()?;
        let header = AlignmentHeader::from_sam_header(&sam_header);
        let state = ReaderState {
            inner: Inner::Sam(reader),
            file_path: path.to_string(),
            fasta_repo: fasta::Repository::default(),
            #[cfg(feature = "cram")]
            cram_record_buf: std::collections::VecDeque::new(),
            file_size,
            inline_indexer: None,
            inline_last_seen: None,
            inline_n_refs_seen: 0,
            inline_built_index: None,
        };
        Ok((state, sam_header, header))
    }

    /// Open a single BAM file, returning the ReaderState, sam::Header, and AlignmentHeader.
    fn open_bam_reader(path: &str, _fasta_repo: Option<fasta::Repository>) -> Result<(ReaderState, sam::Header, AlignmentHeader)> {
        let file_size = std::fs::metadata(path)
            .map_err(|e| anyhow::anyhow!("Failed to stat BAM {}: {}", path, e))?
            .len();

        let index_path = format!("{}.bai", path);
        let has_index = std::path::Path::new(&index_path).exists();

        if has_index {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| anyhow::anyhow!("Failed to open BAM {}: {}", path, e))?;
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let state = ReaderState {
                inner: Inner::Bam(reader),
                file_path: path.to_string(),
                fasta_repo: fasta::Repository::default(),
                #[cfg(feature = "cram")]
                cram_record_buf: std::collections::VecDeque::new(),
                file_size,
                inline_indexer: None,
                inline_last_seen: None,
                inline_n_refs_seen: 0,
                inline_built_index: None,
            };
            Ok((state, sam_header, header))
        } else {
            let file = File::open(path)
                .map_err(|e| anyhow::anyhow!("Failed to open BAM {}: {}", path, e))?;
            let mut reader = bam::io::Reader::new(BufReader::new(file));
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let state = ReaderState {
                inner: Inner::BamNoIndex(reader),
                file_path: path.to_string(),
                fasta_repo: fasta::Repository::default(),
                #[cfg(feature = "cram")]
                cram_record_buf: std::collections::VecDeque::new(),
                file_size,
                inline_indexer: None,
                inline_last_seen: None,
                inline_n_refs_seen: 0,
                inline_built_index: None,
            };
            Ok((state, sam_header, header))
        }
    }

    /// Open a single CRAM file, returning the ReaderState, sam::Header, and AlignmentHeader.
    #[cfg(feature = "cram")]
    fn open_cram_reader(path: &str, fasta_repo: Option<fasta::Repository>) -> Result<(ReaderState, sam::Header, AlignmentHeader)> {
        let file_size = std::fs::metadata(path)
            .map_err(|e| anyhow::anyhow!("Failed to stat CRAM {}: {}", path, e))?
            .len();

        let repo = fasta_repo.unwrap_or_default();
        let index_path = format!("{}.crai", path);
        let has_index = std::path::Path::new(&index_path).exists();

        if has_index {
            let mut reader = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_path(path)
                .map_err(|e| {
                    let msg = e.to_string();
                    if msg.contains("invalid gzip header") || msg.contains("gzip") {
                        anyhow::anyhow!(
                            "Failed to open CRAM index {}.crai: the index file does not appear \
                             to be gzip-compressed. noodles requires gzip-compressed .crai files. \
                             Re-create the index with: samtools index {}",
                            path, path
                        )
                    } else {
                        anyhow::anyhow!("Failed to open CRAM {}: {}", path, e)
                    }
                })?;
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let state = ReaderState {
                inner: Inner::Cram(reader),
                file_path: path.to_string(),
                fasta_repo: repo,
                cram_record_buf: std::collections::VecDeque::new(),
                file_size,
                inline_indexer: None,
                inline_last_seen: None,
                inline_n_refs_seen: 0,
                inline_built_index: None,
            };
            Ok((state, sam_header, header))
        } else {
            let file = File::open(path)
                .map_err(|e| anyhow::anyhow!("Failed to open CRAM {}: {}", path, e))?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_reader(BufReader::new(file));
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let state = ReaderState {
                inner: Inner::CramNoIndex(reader),
                file_path: path.to_string(),
                fasta_repo: repo,
                cram_record_buf: std::collections::VecDeque::new(),
                file_size,
                inline_indexer: None,
                inline_last_seen: None,
                inline_n_refs_seen: 0,
                inline_built_index: None,
            };
            Ok((state, sam_header, header))
        }
    }

    /// Read the next alignment record, chaining across readers for multi-BAM.
    ///
    /// For SAM readers, malformed records are skipped with a warning rather than
    /// aborting (truncated lines are common in SAM files from interrupted writes).
    pub fn read_record(&mut self) -> Result<Option<AlignmentRecord>> {
        loop {
            if self.current_reader_idx >= self.readers.len() {
                if self.skipped_records > 0 {
                    log::warn!(
                        "Skipped {} malformed record(s) in {}",
                        self.skipped_records, self.input_path
                    );
                }
                return Ok(None);
            }
            let sam_header = &self.sam_header;
            let reader = &mut self.readers[self.current_reader_idx];
            match reader.read_record_inner(sam_header) {
                Ok(Some(rec)) => return Ok(Some(rec)),
                Ok(None) => {
                    // Current reader exhausted, advance to next
                    self.current_reader_idx += 1;
                }
                Err(e) if reader.is_sam() => {
                    self.skipped_records += 1;
                    if self.skipped_records <= 5 {
                        log::warn!("Skipping malformed SAM record: {}", e);
                    }
                    if self.skipped_records == 5 {
                        log::warn!("(further malformed-record warnings suppressed)");
                    }
                    // continue to try next record
                }
                Err(e) => return Err(e),
            }
        }
    }

    /// Get current virtual file offset (for progress tracking).
    /// For multi-BAM, returns a composite position across all readers.
    pub fn tell(&self) -> u64 {
        let mut total = 0u64;
        for (i, reader) in self.readers.iter().enumerate() {
            if i < self.current_reader_idx {
                // Completed reader: use its full file size
                total += reader.file_size;
            } else if i == self.current_reader_idx {
                total += reader.tell_inner();
            }
            // Future readers contribute 0
        }
        total
    }

    /// Total file size in BGZF block units across all readers.
    /// Use as the denominator for progress percentage: `bam.tell() >> 16` * 100 / `bam.total_file_blocks()`.
    pub fn total_file_blocks(&self) -> u64 {
        let total_bytes: u64 = self.readers.iter().map(|r| r.file_size).sum();
        (total_bytes >> 16) + 1
    }

    /// Seek to a virtual file offset. Only seek(0) is supported
    /// (resets all readers to the beginning).
    pub fn seek(&mut self, pos: u64) -> Result<()> {
        if pos != 0 {
            bail!("Only seek(0) is supported");
        }
        for reader in &mut self.readers {
            reader.seek_to_start()?;
        }
        self.current_reader_idx = 0;
        self.skipped_records = 0;
        Ok(())
    }

    /// Query records overlapping a region, returning an iterator.
    /// For multi-BAM, queries all readers and merges results sorted by (ref_id, pos).
    pub fn query(&mut self, region: &str) -> Result<RegionIterator> {
        let parsed_region: Region = region.parse()
            .map_err(|e| anyhow::anyhow!("Invalid region '{}': {}", region, e))?;

        let mut all_records: Vec<AlignmentRecord> = Vec::new();

        for reader in &mut self.readers {
            let records = reader.query_inner(&self.sam_header, &parsed_region)?;
            all_records.extend(records);
        }

        // Sort by (ref_id, pos) for multi-BAM coordinate ordering
        if self.readers.len() > 1 {
            all_records.sort_by_key(|r| (r.ref_id, r.pos));
        }

        Ok(RegionIterator { records: all_records, index: 0 })
    }

    /// Returns the noodles sam::Header.
    pub fn sam_header(&self) -> &sam::Header {
        &self.sam_header
    }

    /// Returns true if any reader is CRAM.
    pub fn is_cram(&self) -> bool {
        #[cfg(feature = "cram")]
        {
            self.readers.iter().any(|r| r.is_cram())
        }
        #[cfg(not(feature = "cram"))]
        {
            false
        }
    }

    /// Returns true if all readers have an index file.
    pub fn has_index(&self) -> bool {
        self.readers.iter().all(|r| r.has_index())
    }

    /// Turn on inline-BAI-building for every currently-non-indexed BAM
    /// reader. No effect on readers that already have an on-disk index or
    /// are SAM/CRAM.
    pub fn enable_inline_index(&mut self) {
        for r in &mut self.readers {
            r.enable_inline_index();
        }
    }

    /// Finalize any inline indexes that were built during the most recent
    /// sequential scan. After this returns, `query()` works on readers that
    /// didn't ship with a .bai. Idempotent.
    pub fn finalize_inline_index(&mut self) -> Result<()> {
        // Use the full reference count from the BAM header so the built
        // index covers every contig the pipeline may query — even ones
        // that happened to have no reads in pass 1.
        let n_refs = self.header.refs.len();
        for r in &mut self.readers {
            r.finalize_inline_index(n_refs)?;
        }
        Ok(())
    }

    /// Check that all readers have an index, returning a clear error if not.
    pub fn require_index(&self, path: &str) -> Result<()> {
        if self.readers.len() == 1 {
            // Single-file mode: original error message
            if !self.readers[0].has_index() {
                #[cfg(feature = "cram")]
                let expected = if self.readers[0].is_cram() {
                    format!("{}.crai", self.readers[0].file_path)
                } else {
                    format!("{}.bai", self.readers[0].file_path)
                };
                #[cfg(not(feature = "cram"))]
                let expected = format!("{}.bai", self.readers[0].file_path);
                bail!(
                    "Index file not found for '{}'. Expected '{}'. \
                     Create one with 'samtools index'.",
                    path, expected
                );
            }
        } else {
            // Multi-file mode: report all missing indices
            let missing: Vec<&str> = self.readers.iter()
                .filter(|r| !r.has_index())
                .map(|r| r.file_path.as_str())
                .collect();
            if !missing.is_empty() {
                bail!(
                    "{} of {} BAM files from '{}' are missing index files: {}. \
                     Create indices with 'samtools index'.",
                    missing.len(), self.readers.len(), path,
                    missing.join(", ")
                );
            }
        }
        Ok(())
    }

    /// Returns the number of alignment files being read.
    pub fn num_files(&self) -> usize {
        self.readers.len()
    }

    /// Returns the original path passed to open() (list file or single BAM path).
    pub fn input_path(&self) -> &str {
        &self.input_path
    }
}

/// Iterator over alignment records in a genomic region.
pub struct RegionIterator {
    records: Vec<AlignmentRecord>,
    index: usize,
}

impl Iterator for RegionIterator {
    type Item = Result<AlignmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.records.len() {
            let rec = self.records[self.index].clone();
            self.index += 1;
            Some(Ok(rec))
        } else {
            None
        }
    }
}
