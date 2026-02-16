//! Unified alignment input that supports both BAM and CRAM files via noodles.

use std::fs::File;
use std::io::BufReader;
use anyhow::{Result, bail};

use noodles::bam;
use noodles::bgzf;
use noodles::cram;
use noodles::fasta;
use noodles::sam;
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
    fn from_sam_header(header: &sam::Header) -> Self {
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
fn decode_alignment_record(rec: &dyn sam::alignment::Record, header: &sam::Header) -> Result<AlignmentRecord> {
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

/// Inner reader enum
enum Inner {
    Bam(bam::io::IndexedReader<bgzf::io::Reader<File>>),
    BamNoIndex(bam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    Cram(cram::io::IndexedReader<File>),
    CramNoIndex(cram::io::Reader<BufReader<File>>),
}

/// Unified alignment input wrapping either BAM or CRAM via noodles.
pub struct AlignmentInput {
    inner: Inner,
    sam_header: sam::Header,
    pub header: AlignmentHeader,
    pub contig_mapper: ContigMapper,
    pub start_pos: u64,
    /// Buffer for CRAM records (one container's worth at a time)
    cram_record_buf: std::collections::VecDeque<AlignmentRecord>,
    /// FASTA reference repository for CRAM decoding
    fasta_repo: fasta::Repository,
    /// Path to the alignment file (for re-opening on seek)
    file_path: String,
}

impl AlignmentInput {
    /// Open an alignment file (BAM or CRAM), auto-detected by magic bytes
    /// (falling back to file extension).
    pub fn open(path: &str, ref_path: Option<&str>) -> Result<Self> {
        let fasta_repo = if let Some(rp) = ref_path {
            let indexed_reader = fasta::io::indexed_reader::Builder::default()
                .build_from_path(rp)
                .map_err(|e| anyhow::anyhow!("Failed to open FASTA reference {}: {}", rp, e))?;
            let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
            Some(fasta::Repository::new(adapter))
        } else {
            None
        };

        if Self::is_cram_file(path) {
            Self::open_cram(path, fasta_repo)
        } else {
            Self::open_bam(path, fasta_repo)
        }
    }

    /// Detect whether a file is CRAM by reading the first 4 magic bytes ("CRAM"),
    /// falling back to file extension if the file can't be read.
    fn is_cram_file(path: &str) -> bool {
        if let Ok(mut f) = File::open(path) {
            let mut magic = [0u8; 4];
            if std::io::Read::read_exact(&mut f, &mut magic).is_ok() {
                return &magic == b"CRAM";
            }
        }
        path.ends_with(".cram")
    }

    fn open_bam(path: &str, _fasta_repo: Option<fasta::Repository>) -> Result<Self> {
        let index_path = format!("{}.bai", path);
        let has_index = std::path::Path::new(&index_path).exists();

        if has_index {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| anyhow::anyhow!("Failed to open BAM {}: {}", path, e))?;
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let contig_mapper = ContigMapper::from_refs(&header.refs);
            Ok(AlignmentInput {
                inner: Inner::Bam(reader),
                sam_header,
                header,
                contig_mapper,
                start_pos: 0,
                cram_record_buf: std::collections::VecDeque::new(),
                fasta_repo: fasta::Repository::default(),
                file_path: path.to_string(),
            })
        } else {
            let file = File::open(path)
                .map_err(|e| anyhow::anyhow!("Failed to open BAM {}: {}", path, e))?;
            let mut reader = bam::io::Reader::new(BufReader::new(file));
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let contig_mapper = ContigMapper::from_refs(&header.refs);
            Ok(AlignmentInput {
                inner: Inner::BamNoIndex(reader),
                sam_header,
                header,
                contig_mapper,
                start_pos: 0,
                cram_record_buf: std::collections::VecDeque::new(),
                fasta_repo: fasta::Repository::default(),
                file_path: path.to_string(),
            })
        }
    }

    fn open_cram(path: &str, fasta_repo: Option<fasta::Repository>) -> Result<Self> {
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
            let contig_mapper = ContigMapper::from_refs(&header.refs);
            Ok(AlignmentInput {
                inner: Inner::Cram(reader),
                sam_header,
                header,
                contig_mapper,
                start_pos: 0,
                cram_record_buf: std::collections::VecDeque::new(),
                fasta_repo: repo,
                file_path: path.to_string(),
            })
        } else {
            let file = File::open(path)
                .map_err(|e| anyhow::anyhow!("Failed to open CRAM {}: {}", path, e))?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_reader(BufReader::new(file));
            let sam_header = reader.read_header()?;
            let header = AlignmentHeader::from_sam_header(&sam_header);
            let contig_mapper = ContigMapper::from_refs(&header.refs);
            Ok(AlignmentInput {
                inner: Inner::CramNoIndex(reader),
                sam_header,
                header,
                contig_mapper,
                start_pos: 0,
                cram_record_buf: std::collections::VecDeque::new(),
                fasta_repo: repo,
                file_path: path.to_string(),
            })
        }
    }

    /// Fill the CRAM record buffer by reading and decoding the next container.
    fn fill_cram_buffer(&mut self) -> Result<()> {
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
                &self.sam_header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;
            for cram_rec in &cram_records {
                let ar = decode_alignment_record(cram_rec, &self.sam_header)?;
                self.cram_record_buf.push_back(ar);
            }
        }

        Ok(())
    }

    /// Read the next alignment record.
    pub fn read_record(&mut self) -> Result<Option<AlignmentRecord>> {
        // For CRAM, use the record buffer (filled on first call)
        if matches!(self.inner, Inner::Cram(_) | Inner::CramNoIndex(_)) {
            if self.cram_record_buf.is_empty() {
                self.fill_cram_buffer()?;
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
                    Ok(_) => Ok(Some(decode_bam_record(&buf, &self.sam_header)?)),
                    Err(e) => Err(e.into()),
                }
            }
            Inner::BamNoIndex(r) => {
                let mut buf = bam::Record::default();
                match r.read_record(&mut buf) {
                    Ok(0) => Ok(None),
                    Ok(_) => Ok(Some(decode_bam_record(&buf, &self.sam_header)?)),
                    Err(e) => Err(e.into()),
                }
            }
            _ => unreachable!(),
        }
    }

    /// Get current virtual file offset (for progress tracking).
    /// For BAM files, returns the compressed byte offset in the BGZF file.
    /// For CRAM files, returns 0 (position tracking not implemented).
    pub fn tell(&self) -> u64 {
        match &self.inner {
            Inner::Bam(r) => r.get_ref().virtual_position().compressed(),
            Inner::BamNoIndex(r) => r.get_ref().virtual_position().compressed(),
            // CRAM position tracking not implemented
            Inner::Cram(_) | Inner::CramNoIndex(_) => 0,
        }
    }

    /// Seek to a virtual file offset. Only seek(0) is supported
    /// (re-opens the reader to reset to the beginning).
    pub fn seek(&mut self, pos: u64) -> Result<()> {
        if pos != 0 {
            bail!("Only seek(0) is supported");
        }
        if self.is_cram() {
            self.cram_record_buf.clear();
            let repo = self.fasta_repo.clone();
            let path = &self.file_path;
            let index_path = format!("{}.crai", path);
            if std::path::Path::new(&index_path).exists() {
                let mut reader = cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_path(path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open CRAM {}: {}", path, e))?;
                reader.read_header()?;
                self.inner = Inner::Cram(reader);
            } else {
                let file = File::open(path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open CRAM {}: {}", path, e))?;
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_reader(BufReader::new(file));
                reader.read_header()?;
                self.inner = Inner::CramNoIndex(reader);
            }
        } else {
            // BAM: Re-open file to reset to beginning
            let path = self.file_path.clone();
            let index_path = format!("{}.bai", path);
            if std::path::Path::new(&index_path).exists() {
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(&path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open BAM {}: {}", path, e))?;
                reader.read_header()?;
                self.inner = Inner::Bam(reader);
            } else {
                let file = File::open(&path)
                    .map_err(|e| anyhow::anyhow!("Failed to re-open BAM {}: {}", path, e))?;
                let mut reader = bam::io::Reader::new(BufReader::new(file));
                reader.read_header()?;
                self.inner = Inner::BamNoIndex(reader);
            }
        }
        Ok(())
    }

    /// Query records overlapping a region, returning an iterator.
    pub fn query(&mut self, region: &str) -> Result<RegionIterator> {
        let parsed_region: Region = region.parse()
            .map_err(|e| anyhow::anyhow!("Invalid region '{}': {}", region, e))?;

        match &mut self.inner {
            Inner::Bam(r) => {
                let header = &self.sam_header;
                let records: Vec<AlignmentRecord> = r.query(header, &parsed_region)?
                    .records()
                    .map(|result| -> Result<AlignmentRecord> {
                        let rec = result?;
                        decode_bam_record(&rec, header)
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(RegionIterator { records, index: 0 })
            }
            Inner::Cram(r) => {
                let header = &self.sam_header;
                let records: Vec<AlignmentRecord> = r.query(header, &parsed_region)?
                    .map(|result| -> Result<AlignmentRecord> {
                        let rec = result?;
                        decode_alignment_record(&rec, header)
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(RegionIterator { records, index: 0 })
            }
            _ => bail!("Region queries require an indexed BAM/CRAM file"),
        }
    }

    /// Returns the noodles sam::Header.
    pub fn sam_header(&self) -> &sam::Header {
        &self.sam_header
    }

    /// Returns true if this is a CRAM input.
    pub fn is_cram(&self) -> bool {
        matches!(self.inner, Inner::Cram(_) | Inner::CramNoIndex(_))
    }

    /// Returns true if an index file was found for this reader.
    pub fn has_index(&self) -> bool {
        matches!(self.inner, Inner::Bam(_) | Inner::Cram(_))
    }

    /// Check that an index exists, returning a clear error if not.
    pub fn require_index(&self, path: &str) -> Result<()> {
        if !self.has_index() {
            let expected = if self.is_cram() {
                format!("{}.crai", path)
            } else {
                format!("{}.bai", path)
            };
            bail!(
                "Index file not found for '{}'. Expected '{}'. \
                 Create one with 'samtools index'.",
                path, expected
            );
        }
        Ok(())
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
