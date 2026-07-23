use bitvec::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use crate::bam::ContigMapper;
use crate::config::Contig;
use crate::input::{AlignmentInput, AlignmentHeader, AlignmentRecord};
use crate::utils::bed::BedRegion;
use log::{info, debug, warn};

/// Pre-computed per-bin GC content loaded from a TSV (chrom, bin_start,
/// bin_end, gc_content). Used to bypass FASTA-based `compute_gc_content`
/// when a complete reference isn't shipped (e.g. the browser bundle uses a
/// sparse FASTA). Keyed on (chrom, bin_index) where bin_index = start / bin_size.
pub type GcBinMap = HashMap<(String, usize), f32>;

/// Load a pre-computed GC-per-bin TSV. File format (tab-separated):
///   chromosome<TAB>bin_start<TAB>bin_end<TAB>gc_content
/// Header line starting with "chromosome" is skipped. NaN / "NA" entries
/// are accepted and stored as NaN so the downstream correction knows to
/// skip them.
pub fn load_gc_bins(path: &str, bin_size: usize) -> Result<GcBinMap, Box<dyn std::error::Error>> {
    let file = File::open(path)
        .map_err(|e| std::io::Error::other(format!("Error opening GC bins file {}: {}", path, e)))?;
    let reader = BufReader::new(file);
    let mut map: GcBinMap = HashMap::new();
    for (lineno, line) in reader.lines().enumerate() {
        let l = line?;
        if l.trim().is_empty() || l.starts_with('#') || l.starts_with("chromosome") {
            continue;
        }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() < 4 {
            warn!("gc_bins line {}: expected 4+ columns, got {}; skipping", lineno + 1, parts.len());
            continue;
        }
        let chrom = parts[0].to_string();
        let start: usize = parts[1].parse()
            .map_err(|e| format!("gc_bins line {}: bad start: {}", lineno + 1, e))?;
        let gc: f32 = if parts[3] == "NA" || parts[3] == "NaN" {
            f32::NAN
        } else {
            parts[3].parse()
                .map_err(|e| format!("gc_bins line {}: bad gc value: {}", lineno + 1, e))?
        };
        let bin_idx = start / bin_size;
        map.insert((chrom, bin_idx), gc);
    }
    info!("Loaded {} pre-computed GC bins from {}", map.len(), path);
    Ok(map)
}

// ==================== Shared CoverageAccumulator ====================
// Used by both standalone read_depth() and pipeline mode

struct ChromData {
    mapped_name: String,
    len: usize,
    mask: BitVec,
    bin_read_counts: Vec<u32>,
    bin_unmasked_counts: Vec<u32>,
    bin_gc_content: Vec<f32>,
}

/// Shared coverage accumulator - used by both CLI and pipeline.
/// Processes records one at a time, accumulating read counts per 1Mb bin.
pub struct CoverageAccumulator {
    chroms: Vec<Option<ChromData>>,
    bin_size: usize,
    reads_aligned: u64,
}

fn lookup_gc_bins(gc: &GcBinMap, mapped_name: &str, n_bins: usize) -> Vec<f32> {
    let mut out = Vec::with_capacity(n_bins);
    let mut hits = 0usize;
    for i in 0..n_bins {
        match gc.get(&(mapped_name.to_string(), i)) {
            Some(&v) => { out.push(v); hits += 1; }
            None => out.push(f32::NAN),
        }
    }
    debug!("gc_bins: {} has {}/{} bins with pre-computed GC", mapped_name, hits, n_bins);
    out
}

impl CoverageAccumulator {
    /// Create new accumulator with repeat masks and optional GC content.
    pub fn new(header: &AlignmentHeader, repeats: &[BedRegion], ref_path: Option<&str>, contigs: &[Contig]) -> Self {
        Self::new_with_gc(header, repeats, ref_path, None, contigs)
    }

    /// Like `new`, but if `gc_bins` is Some, pre-computed GC values are used
    /// instead of scanning the FASTA. Keys in `gc_bins` are (mapped_name,
    /// bin_index) where mapped_name is the chr-style name (e.g. "chr1") and
    /// bin_index = start_bp / bin_size. Used by the browser-bundle path
    /// where the shipped FASTA is sparse (mostly N).
    pub fn new_with_gc(
        header: &AlignmentHeader,
        repeats: &[BedRegion],
        ref_path: Option<&str>,
        gc_bins: Option<&GcBinMap>,
        contigs: &[Contig],
    ) -> Self {
        info!("Initializing coverage accumulator...");
        let bin_size = 1_000_000;
        let mut repeats_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

        for r in repeats {
            repeats_map
                .entry(r.segment.clone())
                .or_default()
                .push((r.start as usize, r.end as usize));
        }

        // Build set of chromosomes that have repeats (we only need to mask those)
        let chroms_with_repeats: HashSet<&String> = repeats_map.keys().collect();

        // Open indexed FASTA reader for GC content calculation
        let mut fasta_reader = ref_path.and_then(|rp| {
            noodles::fasta::io::indexed_reader::Builder::default()
                .build_from_path(rp)
                .ok()
        });

        // Detect FASTA naming convention from the .fai index
        let fasta_mapper = ref_path.and_then(|rp| {
            let fai_path = format!("{}.fai", rp);
            ContigMapper::from_fai(&fai_path, contigs).ok()
        });

        let mut chroms = Vec::with_capacity(header.refs.len());
        let mapper = ContigMapper::from_contigs(contigs);

        for (idx, (name, len)) in header.refs.iter().zip(header.lengths.iter()).enumerate() {
            let len = *len as usize;
            if len == 0 {
                chroms.push(None);
                continue;
            }

            let mapped_name = mapper.to_chr_name(name);

            // Only create full mask for chromosomes with repeats
            // Others get a lightweight entry (no mask, assume all unmasked)
            let has_repeats = chroms_with_repeats.contains(&mapped_name);

            let n_bins = len.div_ceil(bin_size);

            if has_repeats {
                debug!("Processing chromosome {} ({}/{})...", name, idx + 1, header.refs.len());

                let mut mask = bitvec![0; len + 1];

                // Apply repeats
                if let Some(reps) = repeats_map.get(&mapped_name) {
                    for &(s, e) in reps {
                        let start = s.min(len);
                        let end = e.min(len);
                        for mut item in mask.iter_mut().take(end).skip(start) {
                            *item = true;
                        }
                    }
                }

                let mut bin_unmasked_counts = vec![0u32; n_bins];
                for (i, bin_count) in bin_unmasked_counts.iter_mut().enumerate() {
                    let start = i * bin_size;
                    let end = std::cmp::min((i + 1) * bin_size, len);
                    let mut count = 0;
                    for item in mask.iter().take(end).skip(start) {
                        if !item {
                            count += 1;
                        }
                    }
                    *bin_count = count;
                }

                // Per-bin GC content: prefer pre-computed table if supplied
                // (exact bin-level GC from a full reference), otherwise scan
                // the FASTA. The table is keyed by (mapped_name, bin_index).
                let bin_gc_content = if let Some(gc) = gc_bins {
                    lookup_gc_bins(gc, &mapped_name, n_bins)
                } else {
                    Self::compute_gc_content(
                        &mut fasta_reader, name, fasta_mapper.as_ref(), len, n_bins, bin_size, &mask,
                    )
                };

                chroms.push(Some(ChromData {
                    mapped_name,
                    len,
                    mask,
                    bin_read_counts: vec![0u32; n_bins],
                    bin_unmasked_counts,
                    bin_gc_content,
                }));
            } else {
                // No repeats for this chromosome - use lightweight initialization
                // All positions are unmasked, bin_unmasked_counts = bin_size for each bin
                let bin_unmasked_counts: Vec<u32> = (0..n_bins)
                    .map(|i| {
                        let start = i * bin_size;
                        let end = std::cmp::min((i + 1) * bin_size, len);
                        (end - start) as u32
                    })
                    .collect();

                // GC: table lookup if supplied; otherwise NaN (no repeats means
                // we skipped GC-from-FASTA for speed — same as before).
                let bin_gc_content = if let Some(gc) = gc_bins {
                    lookup_gc_bins(gc, &mapped_name, n_bins)
                } else {
                    vec![f32::NAN; n_bins]
                };

                chroms.push(Some(ChromData {
                    mapped_name,
                    len,
                    mask: bitvec![0; 1], // minimal mask - we'll handle specially in process()
                    bin_read_counts: vec![0u32; n_bins],
                    bin_unmasked_counts,
                    bin_gc_content,
                }));
            }
        }

        info!("Coverage accumulator initialized.");
        Self {
            chroms,
            bin_size,
            reads_aligned: 0,
        }
    }

    fn compute_gc_content<R: std::io::BufRead + std::io::Seek>(
        fasta_reader: &mut Option<noodles::fasta::io::IndexedReader<R>>,
        name: &str,
        fasta_mapper: Option<&ContigMapper>,
        len: usize,
        n_bins: usize,
        bin_size: usize,
        mask: &BitVec,
    ) -> Vec<f32> {
        let reader = match fasta_reader.as_mut() {
            Some(r) => r,
            None => return vec![f32::NAN; n_bins],
        };

        // Translate BAM chromosome name to FASTA naming convention
        let fasta_name = if let Some(fm) = fasta_mapper {
            fm.to_bam_name(name)
        } else {
            name.to_string()
        };

        // Query the full chromosome sequence (1-based, inclusive)
        let region_str = format!("{}:{}-{}", fasta_name, 1, len);
        let region: noodles::core::Region = match region_str.parse() {
            Ok(r) => r,
            Err(_) => {
                debug!("Could not parse FASTA region for {} (fasta: {}), skipping GC", name, fasta_name);
                return vec![f32::NAN; n_bins];
            }
        };
        let record: noodles::fasta::Record = match reader.query(&region) {
            Ok(r) => r,
            Err(_) => {
                debug!("Could not query FASTA for {} (fasta: {}), skipping GC", name, fasta_name);
                return vec![f32::NAN; n_bins];
            }
        };
        let seq = record.sequence().as_ref();

        let mut gc_content = Vec::with_capacity(n_bins);
        for i in 0..n_bins {
            let start = i * bin_size;
            let end = std::cmp::min((i + 1) * bin_size, len);
            let mut gc = 0u32;
            let mut at = 0u32;
            for pos in start..end {
                if mask[pos] {
                    continue;
                }
                if pos < seq.len() {
                    match seq[pos] | 0x20 {  // lowercase
                        b'g' | b'c' => gc += 1,
                        b'a' | b't' => at += 1,
                        _ => {} // N or other ambiguous bases
                    }
                }
            }
            let total = gc + at;
            if total > 0 {
                gc_content.push(gc as f32 / total as f32);
            } else {
                gc_content.push(f32::NAN);
            }
        }
        gc_content
    }

    /// Process a single record - called per-record during BAM scan.
    pub fn process(&mut self, record: &AlignmentRecord) {
        // Count primary aligned reads (not unmapped, not secondary, not supplementary)
        let flag = record.flags();
        if (flag & 0x904) == 0 {
            self.reads_aligned += 1;
        }

        if record.ref_id < 0 {
            return;
        }
        let id = record.ref_id as usize;

        if id < self.chroms.len()
            && let Some(data) = &mut self.chroms[id]
                && let Some(pos) = record.alignment_start() {
                    // Check if position is masked (for lightweight masks, all positions are unmasked)
                    let is_masked = if data.mask.len() > pos {
                        data.mask[pos]
                    } else {
                        false // Lightweight mask - position is unmasked
                    };

                    if pos < data.len && !is_masked {
                        let bin = pos / self.bin_size;
                        if bin < data.bin_read_counts.len() {
                            data.bin_read_counts[bin] += 1;
                        }
                    }
                }
    }

    /// Get the total count of aligned reads (primary alignments only).
    pub fn reads_aligned(&self) -> u64 {
        self.reads_aligned
    }

    /// Return the accumulator's bins in memory in the same shape and order
    /// the TSV would emit. Any bin whose unmasked fraction is at or below
    /// 25 % of `bin_size` is dropped (same threshold as `write_output`), so
    /// the resulting Vec is the exact input the karyotype call expects.
    ///
    /// Callers that need to also render / carry gc-corrected values fill
    /// `coverage_gc_adj` themselves after running the fit; leave it `None`
    /// when only the raw pass exists.
    pub fn to_bins(&self) -> Vec<CoverageBin> {
        let mut out = Vec::new();
        for data in self.chroms.iter().flatten() {
            for (i, (&unmasked, &reads)) in data
                .bin_unmasked_counts
                .iter()
                .zip(data.bin_read_counts.iter())
                .enumerate()
            {
                let threshold = (self.bin_size as f64 * 0.25) as u32;
                if unmasked <= threshold { continue; }

                let start = (i * self.bin_size) as u32;
                let end = std::cmp::min((i + 1) * self.bin_size, data.len) as u32;
                let mut value = reads as f64;
                value *= self.bin_size as f64 / unmasked as f64;

                let gc = data.bin_gc_content[i] as f64;
                let gc_opt = if gc.is_nan() { None } else { Some(gc) };

                out.push(CoverageBin {
                    chrom: data.mapped_name.clone(),
                    start,
                    end,
                    coverage: value,
                    coverage_gc_adj: None,
                    maf_data: None,
                    gc_content: gc_opt,
                });
            }
        }
        out
    }

    /// Serialize a set of bins to disk in the unified 6-column format used
    /// by the whole pipeline. Column 6 (`n_reads_gc_adj`) is populated when
    /// the bin's `coverage_gc_adj` is `Some`; otherwise emitted as `NA`.
    ///
    /// The single-writer path replaces the historical split between
    /// `.coverage.tsv` (5-col, raw) and `.coverage.gc_adjusted.tsv` (5-col,
    /// corrected). Downstream readers pick a column based on what they need.
    pub fn write_bins_tsv(path: &str, bins: &[CoverageBin]) -> std::io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "chromosome\tbin_start\tbin_end\tn_reads\tgc_content\tn_reads_gc_adj")?;
        for b in bins {
            let gc = match b.gc_content {
                Some(v) => format!("{:.4}", v),
                None => "NA".to_string(),
            };
            let adj = match b.coverage_gc_adj {
                Some(v) => format!("{:.4}", v),
                None => "NA".to_string(),
            };
            writeln!(
                file, "{}\t{}\t{}\t{:.0}\t{}\t{}",
                b.chrom, b.start, b.end, b.coverage, gc, adj,
            )?;
        }
        Ok(())
    }

    /// Legacy: write the raw 5-column TSV. Kept only for the standalone
    /// `nasvar coverage` subcommand and its downstream (which may not have
    /// a karyotype step to produce gc-adjusted values). Prefer `to_bins()`
    /// + `write_bins_tsv()` in new code so we stay on the unified schema.
    pub fn write_output(&self, path: &str, include_gc: bool) -> std::io::Result<()> {
        let bins = self.to_bins();
        // include_gc is retained for signature compatibility; the unified
        // writer always emits GC (as `NA` when unknown), so a caller passing
        // false effectively gets the same 6-col file. Callers that truly
        // need a 4-col format should filter after read.
        let _ = include_gc;
        Self::write_bins_tsv(path, &bins)
    }
}

// ============================================================================
// CoverageBin — the pipeline's per-bin coverage row, shared between the
// accumulator, the karyotype call, CNV, and the TSV reader/writer. Kept in
// this module because the accumulator is the canonical producer; karyotype
// re-exports it so existing `use crate::karyotype::CoverageBin` imports
// still resolve.
// ============================================================================

#[derive(Debug, Clone)]
pub struct CoverageBin {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    /// Mask-adjusted coverage (raw, `reads × bin_size / unmasked`).
    pub coverage: f64,
    /// GC-corrected coverage (raw × GC-fit scale). `None` until the
    /// karyotype step's fit has been applied.
    pub coverage_gc_adj: Option<f64>,
    /// Median MAF at bin locations (populated by the karyotype step for
    /// BAF-heavy plots; unused by the accumulator itself).
    pub maf_data: Option<f64>,
    /// GC fraction across the (unmasked) bases in this bin. `None` when
    /// the FASTA wasn't available so GC couldn't be computed.
    pub gc_content: Option<f64>,
}

/// Read the unified 6-column TSV (or a legacy 5-column raw file) back into
/// `CoverageBin`s. Column 6 (`n_reads_gc_adj`) is optional; when absent or
/// `"NA"`, the returned bin's `coverage_gc_adj` is `None`.
pub fn parse_coverage_tsv(cov_path: &str) -> std::io::Result<Vec<CoverageBin>> {
    let file = File::open(cov_path)?;
    let reader = BufReader::new(file);
    let mut bins = Vec::new();
    for line in reader.lines() {
        let l = line?;
        if l.starts_with("chromosome") { continue; }
        let p: Vec<&str> = l.split('\t').collect();
        if p.len() < 4 { continue; }
        let chrom = p[0].to_string();
        let start: u32 = p[1].parse().unwrap_or(0);
        let end:   u32 = p[2].parse().unwrap_or(0);
        let cov:   f64 = p[3].parse().unwrap_or(0.0);
        let gc = if p.len() >= 5 {
            p[4].parse::<f64>().ok().filter(|v| v.is_finite())
        } else { None };
        let gc_adj = if p.len() >= 6 {
            p[5].parse::<f64>().ok().filter(|v| v.is_finite())
        } else { None };
        bins.push(CoverageBin {
            chrom, start, end,
            coverage: cov,
            coverage_gc_adj: gc_adj,
            maf_data: None,
            gc_content: gc,
        });
    }
    Ok(bins)
}

// ==================== End CoverageAccumulator ====================

/// Calculate read depth per 1Mb bin and write to TSV file.
/// Returns the total count of aligned reads (primary alignments only).
///
/// Uses the shared CoverageAccumulator for consistency with pipeline mode.
pub fn read_depth(
    bam: &mut AlignmentInput,
    repeats: &[BedRegion],
    out_prefix: &str,
    as_alignments: Option<&str>,
    ref_path: Option<&str>,
) -> Result<u64, Box<dyn std::error::Error>>
{
    info!("Calculating coverage...");

    let out_file = format!("{}.coverage.tsv", out_prefix);
    let include_gc = ref_path.is_some();

    if let Some(as_path) = as_alignments {
        // Use adaptive sampling alignments exclusively for coverage
        let accumulator = scan_as_alignments(as_path, &bam.header, repeats, ref_path, &bam.contigs)?;
        accumulator.write_output(&out_file, include_gc)?;
        let reads_aligned = accumulator.reads_aligned();
        info!("Total aligned reads (primary, from AS): {}", reads_aligned);
        Ok(reads_aligned)
    } else {
        // Standard: scan main BAM for coverage
        let mut accumulator = CoverageAccumulator::new(&bam.header, repeats, ref_path, &bam.contigs);
        bam.seek(bam.start_pos)?;
        while let Some(record) = bam.read_record()? {
            accumulator.process(&record);
        }
        accumulator.write_output(&out_file, include_gc)?;
        let reads_aligned = accumulator.reads_aligned();
        info!("Total aligned reads (primary): {}", reads_aligned);
        Ok(reads_aligned)
    }
}

/// Scan an adaptive sampling alignment file (SAM/BAM/CRAM) for coverage,
/// creating and returning a new CoverageAccumulator.
///
/// Used when --as-alignments is provided to replace the main BAM for
/// coverage/karyotyping. The file does not need to be sorted or indexed.
/// Its reference sequences must match the primary BAM header.
pub fn scan_as_alignments(
    path: &str,
    bam_header: &AlignmentHeader,
    repeats: &[BedRegion],
    ref_path: Option<&str>,
    contigs: &[Contig],
) -> Result<CoverageAccumulator, Box<dyn std::error::Error>>
{
    scan_as_alignments_with_gc(path, bam_header, repeats, ref_path, None, contigs, None)
}

/// Scan an adaptive-sampling BAM for coverage. Returns the mask-adjusted
/// 1 Mb accumulator (drives karyotype + CNV, as before). When AS is
/// providing coverage for a run, the caller passes a `BinaryCoverage100k`
/// via `cov100k` so the 100 kb whole-genome track is built from the same
/// reads — without this, the viewer's coverage panel showed the *main*
/// BAM's read distribution while every karyotype/CNV number was computed
/// off the AS BAM.
///
/// # Read-name dedup
///
/// MinKNOW's live adaptive-sampling BAM stream re-emits the *same read*
/// multiple times, each time with FLAG 0 (a spec-violating "second
/// primary"). Counting every occurrence inflates coverage 2–5× on
/// re-emitted reads. We keep only the first sighting of each read name;
/// later sightings are skipped for every downstream accumulator.
///
/// The dedup table is a `HashSet<String>` sized to the AS BAM's unique
/// read count. On a typical AS run (~100 k unique reads, ONT UUID names
/// ~36 B each) that's a few MB — negligible next to the pipeline's other
/// scans. Memory scales linearly with unique reads.
pub fn scan_as_alignments_with_gc(
    path: &str,
    bam_header: &AlignmentHeader,
    repeats: &[BedRegion],
    ref_path: Option<&str>,
    gc_bins: Option<&GcBinMap>,
    contigs: &[Contig],
    mut cov100k: Option<&mut crate::var::coverage_binary::BinaryCoverage100k>,
) -> Result<CoverageAccumulator, Box<dyn std::error::Error>>
{
    info!("Scanning adaptive sampling alignments: {}", path);

    let mut input = AlignmentInput::open(path, None, contigs)?;

    // Validate reference sequences match the primary input
    AlignmentInput::validate_headers_match(bam_header, &input.header, path)?;

    let mut accumulator = CoverageAccumulator::new_with_gc(bam_header, repeats, ref_path, gc_bins, contigs);

    // First-sighting-wins dedup table for AS re-emit. Records without a
    // read name pass through (nothing to key on); this is rare enough in
    // practice that inflating them by a factor of ~1 is preferable to
    // dropping them entirely.
    let mut seen: HashSet<String> = HashSet::new();
    let mut count = 0u64;
    let mut duplicates = 0u64;
    while let Some(record) = input.read_record()? {
        count += 1;
        if count.is_multiple_of(100_000) {
            eprint!("\rAS alignments: {} records...", count);
            std::io::Write::flush(&mut std::io::stderr())?;
        }

        // First-sighting-wins: skip re-emitted reads. `insert` returns
        // false when the name was already present.
        if let Some(name) = record.name()
            && !seen.insert(name.to_string())
        {
            duplicates += 1;
            continue;
        }

        accumulator.process(&record);
        if let Some(c) = cov100k.as_deref_mut() { c.process(&record); }
    }
    if count >= 100_000 {
        eprintln!();
    }

    let unique = count - duplicates;
    if duplicates > 0 {
        info!(
            "AS alignments: processed {} records from {} ({} unique, {} duplicate names skipped)",
            count, path, unique, duplicates,
        );
    } else {
        info!("AS alignments: processed {} records from {}", count, path);
    }
    Ok(accumulator)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch_path(tag: &str) -> String {
        std::env::temp_dir()
            .join(format!("nasvar_covtsv_{}_{}.tsv", tag, std::process::id()))
            .to_string_lossy()
            .into_owned()
    }

    #[test]
    fn unified_tsv_roundtrips_with_and_without_gc_adj() {
        let bins = vec![
            CoverageBin {
                chrom: "chr1".to_string(), start: 0, end: 1_000_000,
                coverage: 1234.0,
                coverage_gc_adj: Some(1200.5),
                maf_data: None,
                gc_content: Some(0.412),
            },
            CoverageBin {
                chrom: "chr1".to_string(), start: 1_000_000, end: 2_000_000,
                coverage: 987.0,
                coverage_gc_adj: None,     // pass-1-only bin
                maf_data: None,
                gc_content: None,          // GC unknown
            },
        ];
        let path = scratch_path("both");
        CoverageAccumulator::write_bins_tsv(&path, &bins).unwrap();
        let back = parse_coverage_tsv(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].chrom, "chr1");
        assert_eq!(back[0].start, 0);
        assert_eq!(back[0].coverage, 1234.0);
        assert_eq!(back[0].coverage_gc_adj, Some(1200.5));
        assert!((back[0].gc_content.unwrap() - 0.412).abs() < 1e-4);
        assert_eq!(back[1].coverage, 987.0);
        assert!(back[1].coverage_gc_adj.is_none());
        assert!(back[1].gc_content.is_none());
    }

    #[test]
    fn parse_accepts_legacy_five_column_files() {
        // Older packages predate the unified schema and only ship 5
        // columns (no gc_adj). Reader should still produce sane bins.
        let path = scratch_path("legacy5");
        std::fs::write(
            &path,
            b"chromosome\tbin_start\tbin_end\tn_reads\tgc_content\n\
              chr1\t0\t1000000\t1234\t0.412\n\
              chr2\t0\t500000\t42\tNA\n",
        ).unwrap();
        let back = parse_coverage_tsv(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].coverage, 1234.0);
        assert!(back[0].coverage_gc_adj.is_none()); // column 6 absent
        assert!((back[0].gc_content.unwrap() - 0.412).abs() < 1e-4);
        assert_eq!(back[1].coverage, 42.0);
        assert!(back[1].gc_content.is_none());     // "NA" → None
    }
}

