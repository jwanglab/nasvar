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

/// Per-bin GC fraction over `seq` (uppercase or lowercase bases). When `mask`
/// is `Some`, positions where the bit is set are excluded (repeat-masked);
/// `None` counts every position. Bins with no A/T/G/C bases (all-N or fully
/// masked) yield `NaN`. Kept free of FASTA I/O so it can be unit-tested on a
/// synthetic sequence — this is the GC math that must run regardless of
/// whether repeat data is present.
fn gc_fraction_per_bin(
    seq: &[u8],
    len: usize,
    n_bins: usize,
    bin_size: usize,
    mask: Option<&BitVec>,
) -> Vec<f32> {
    let mut gc_content = Vec::with_capacity(n_bins);
    for i in 0..n_bins {
        let start = i * bin_size;
        let end = std::cmp::min((i + 1) * bin_size, len);
        let mut gc = 0u32;
        let mut at = 0u32;
        for pos in start..end {
            // A lightweight/absent mask leaves every position unmasked; a
            // position past the mask's end is treated as unmasked too.
            if let Some(m) = mask
                && pos < m.len()
                && m[pos] {
                    continue;
                }
            if pos < seq.len() {
                match seq[pos] | 0x20 {
                    // lowercase
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

/// Result of mapping a repeats BED onto the reference's contig naming.
struct RepeatsMapStats {
    /// Repeat intervals keyed by chr-normalized contig name.
    map: HashMap<String, Vec<(usize, usize)>>,
    /// Total number of regions in the input BED.
    total_regions: usize,
    /// Number of regions whose normalized name matches a reference contig.
    mapped_regions: usize,
}

/// Build a repeats map keyed by chr-normalized contig name (so a BED using
/// `1`, `chr1`, or an accession all resolve to the same key the coverage
/// lookup uses). `known_contigs` holds the normalized reference contig names
/// so callers can detect a BED whose names don't match the reference.
fn build_normalized_repeats_map(
    repeats: &[BedRegion],
    mapper: &ContigMapper,
    known_contigs: &HashSet<String>,
) -> RepeatsMapStats {
    let mut map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut mapped_regions = 0;
    for r in repeats {
        let key = mapper.to_chr_name(&r.segment);
        if known_contigs.contains(&key) {
            mapped_regions += 1;
        }
        map.entry(key)
            .or_default()
            .push((r.start as usize, r.end as usize));
    }
    RepeatsMapStats { map, total_regions: repeats.len(), mapped_regions }
}

/// Verdict on a FASTA-based GC scan, used to surface silent failures where
/// the FASTA is present but produces no usable GC (a contig-naming mismatch
/// between the alignment/reference config and the FASTA index, or an index
/// that couldn't be opened).
#[derive(Debug, PartialEq, Eq)]
enum GcScanOutcome {
    /// GC either wasn't scanned from the FASTA, or produced finite values.
    Ok,
    /// A FASTA was given but its index couldn't be opened.
    ReaderUnavailable,
    /// The FASTA opened, bins exist, but no contig yielded any finite GC —
    /// almost always a `.fai`-vs-alignment contig-naming mismatch.
    NoContigResolved,
}

/// Classify the result of GC computation so the caller can warn loudly when
/// a supplied FASTA silently contributed nothing. Only meaningful when GC was
/// meant to come from the FASTA (`ref_path_given && !used_gc_bins`).
fn assess_fasta_gc_scan(
    ref_path_given: bool,
    used_gc_bins: bool,
    reader_opened: bool,
    total_bins: usize,
    finite_gc_bins: usize,
) -> GcScanOutcome {
    // GC wasn't sourced from the FASTA (no reference, or a pre-computed table
    // was used) — nothing to diagnose here.
    if !ref_path_given || used_gc_bins {
        return GcScanOutcome::Ok;
    }
    if !reader_opened {
        return GcScanOutcome::ReaderUnavailable;
    }
    // Opened and there were bins to fill, yet none resolved to a finite value.
    if total_bins > 0 && finite_gc_bins == 0 {
        return GcScanOutcome::NoContigResolved;
    }
    GcScanOutcome::Ok
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
        let mapper = ContigMapper::from_contigs(contigs);

        // Normalize repeats-BED contig names into the same chr-style space the
        // per-chromosome lookup uses, so a BED written with `1`, `chr1`, or an
        // accession all match. `known_contigs` are the reference contigs
        // present in this BAM header (normalized), used to flag a BED whose
        // names don't line up with the reference at all.
        let known_contigs: HashSet<String> =
            header.refs.iter().map(|n| mapper.to_chr_name(n)).collect();
        let repeats_stats = build_normalized_repeats_map(repeats, &mapper, &known_contigs);

        // A missing or mis-named repeats BED must not fail silently: it only
        // disables repeat masking (coverage/karyotype still run, and GC
        // correction now proceeds regardless — see the GC step below), so we
        // warn loudly rather than error.
        if repeats_stats.total_regions == 0 {
            warn!(
                "[coverage] repeats BED contained no usable regions (empty or comment-only); \
                 repeat masking is disabled"
            );
        } else if repeats_stats.mapped_regions == 0 {
            let example_acc = contigs
                .first()
                .map(|c| c.accession.as_str())
                .unwrap_or("NC_000001.11");
            warn!(
                "[coverage] none of {} repeats-BED regions matched any reference contig — \
                 check the BED's contig naming (e.g. 'chr1' vs '1' vs accession '{}'); \
                 repeat masking is disabled",
                repeats_stats.total_regions, example_acc
            );
        } else if repeats_stats.mapped_regions < repeats_stats.total_regions {
            warn!(
                "[coverage] {} of {} repeats-BED regions did not match any reference contig \
                 and will be ignored (check contig naming)",
                repeats_stats.total_regions - repeats_stats.mapped_regions,
                repeats_stats.total_regions
            );
        }

        let repeats_map = repeats_stats.map;
        // Chromosomes we build a repeat mask for. Unmatched BED keys can't
        // equal any normalized reference name, so they're harmless here.
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

        // Track whether FASTA-based GC actually produced anything, so a
        // reference that silently contributes no GC (unopenable index, or a
        // contig-naming mismatch) surfaces as a warning rather than only a
        // downstream "insufficient data points" failure.
        let reader_opened = fasta_reader.is_some();
        let mut total_gc_bins = 0usize;
        let mut finite_gc_bins = 0usize;

        let mut chroms = Vec::with_capacity(header.refs.len());

        for (idx, (name, len)) in header.refs.iter().zip(header.lengths.iter()).enumerate() {
            let len = *len as usize;
            if len == 0 {
                chroms.push(None);
                continue;
            }

            let mapped_name = mapper.to_chr_name(name);

            // A full repeat mask is built only for chromosomes present in the
            // repeats BED; others stay unmasked (lightweight minimal mask that
            // process() treats as all-unmasked).
            let has_repeats = chroms_with_repeats.contains(&mapped_name);

            let n_bins = len.div_ceil(bin_size);

            let (mask, bin_unmasked_counts): (BitVec, Vec<u32>) = if has_repeats {
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

                let mut counts = vec![0u32; n_bins];
                for (i, bin_count) in counts.iter_mut().enumerate() {
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
                (mask, counts)
            } else {
                // No repeats: all positions unmasked, so each bin's unmasked
                // count is simply its base span.
                let counts: Vec<u32> = (0..n_bins)
                    .map(|i| {
                        let start = i * bin_size;
                        let end = std::cmp::min((i + 1) * bin_size, len);
                        (end - start) as u32
                    })
                    .collect();
                // Minimal mask - handled specially in process().
                (bitvec![0; 1], counts)
            };

            // GC content is computed for EVERY chromosome, independent of
            // whether repeat data is present: prefer the pre-computed table,
            // otherwise scan the FASTA. Repeat-masked positions are excluded
            // from GC only when a real mask exists (has_repeats). This is what
            // lets GC correction proceed even when the repeats BED is missing
            // or mis-named. The table is keyed by (mapped_name, bin_index).
            let bin_gc_content = if let Some(gc) = gc_bins {
                lookup_gc_bins(gc, &mapped_name, n_bins)
            } else {
                let mask_opt = if has_repeats { Some(&mask) } else { None };
                Self::compute_gc_content(
                    &mut fasta_reader, name, fasta_mapper.as_ref(), len, n_bins, bin_size, mask_opt,
                )
            };

            total_gc_bins += bin_gc_content.len();
            finite_gc_bins += bin_gc_content.iter().filter(|v| !v.is_nan()).count();

            chroms.push(Some(ChromData {
                mapped_name,
                len,
                mask,
                bin_read_counts: vec![0u32; n_bins],
                bin_unmasked_counts,
                bin_gc_content,
            }));
        }

        // Diagnose a FASTA that was supplied but contributed no GC.
        match assess_fasta_gc_scan(
            ref_path.is_some(),
            gc_bins.is_some(),
            reader_opened,
            total_gc_bins,
            finite_gc_bins,
        ) {
            GcScanOutcome::ReaderUnavailable => warn!(
                "[coverage] could not open FASTA index for GC computation (reference '{}', \
                 expected index '{}.fai'); GC correction will be skipped",
                ref_path.unwrap_or("<none>"),
                ref_path.unwrap_or("<none>"),
            ),
            GcScanOutcome::NoContigResolved => warn!(
                "[coverage] FASTA '{}' opened but no contig produced GC values — almost \
                 certainly a contig-naming mismatch between the alignment/reference config \
                 and the FASTA index (check that '{}.fai' sequence names match the BAM/config, \
                 e.g. 'chr1' vs '1' vs accession); GC correction will be skipped",
                ref_path.unwrap_or("<none>"),
                ref_path.unwrap_or("<none>"),
            ),
            GcScanOutcome::Ok => {}
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
        mask: Option<&BitVec>,
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
        gc_fraction_per_bin(seq, len, n_bins, bin_size, mask)
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

    // ---- gc_fraction_per_bin: GC math must work with or without a mask ----

    #[test]
    fn gc_fraction_computes_per_bin_without_mask() {
        // "GGCC" -> 100% GC, "AATT" -> 0% GC
        let out = gc_fraction_per_bin(b"GGCCAATT", 8, 2, 4, None);
        assert_eq!(out.len(), 2);
        assert!((out[0] - 1.0).abs() < 1e-6);
        assert!((out[1] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn gc_fraction_excludes_masked_positions() {
        // seq "GATC": G,A,T,C. Mask positions 2 (T) and 3 (C); remaining
        // counted bases are G and A -> gc=1, at=1 -> 0.5.
        let mut mask = bitvec![0; 4];
        mask.set(2, true);
        mask.set(3, true);
        let out = gc_fraction_per_bin(b"GATC", 4, 1, 4, Some(&mask));
        assert_eq!(out.len(), 1);
        assert!((out[0] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn gc_fraction_all_ambiguous_bin_is_nan() {
        let out = gc_fraction_per_bin(b"NNNN", 4, 1, 4, None);
        assert_eq!(out.len(), 1);
        assert!(out[0].is_nan());
    }

    // ---- build_normalized_repeats_map: name normalization + match stats ----

    fn test_contigs() -> Vec<Contig> {
        vec![
            Contig { name: "chr1".into(), accession: "NC_000001.11".into() },
            Contig { name: "chr2".into(), accession: "NC_000002.12".into() },
        ]
    }

    fn region(segment: &str) -> BedRegion {
        BedRegion { segment: segment.into(), start: 0, end: 100, name: segment.into() }
    }

    #[test]
    fn repeats_map_normalizes_short_and_accession_names() {
        let contigs = test_contigs();
        let mapper = ContigMapper::from_contigs(&contigs);
        let known: HashSet<String> = contigs.iter().map(|c| c.name.clone()).collect();
        let repeats = vec![region("1"), region("NC_000002.12")];
        let stats = build_normalized_repeats_map(&repeats, &mapper, &known);
        assert!(stats.map.contains_key("chr1")); // short "1" -> "chr1"
        assert!(stats.map.contains_key("chr2")); // accession -> "chr2"
        assert_eq!(stats.total_regions, 2);
        assert_eq!(stats.mapped_regions, 2);
    }

    #[test]
    fn repeats_map_flags_unmapped_names() {
        let contigs = test_contigs();
        let mapper = ContigMapper::from_contigs(&contigs);
        let known: HashSet<String> = contigs.iter().map(|c| c.name.clone()).collect();
        let stats = build_normalized_repeats_map(&[region("contig_xyz")], &mapper, &known);
        assert_eq!(stats.total_regions, 1);
        assert_eq!(stats.mapped_regions, 0); // name matches no reference contig
    }

    #[test]
    fn repeats_map_empty_input_has_zero_regions() {
        let contigs = test_contigs();
        let mapper = ContigMapper::from_contigs(&contigs);
        let known: HashSet<String> = contigs.iter().map(|c| c.name.clone()).collect();
        let stats = build_normalized_repeats_map(&[], &mapper, &known);
        assert_eq!(stats.total_regions, 0);
        assert_eq!(stats.mapped_regions, 0);
        assert!(stats.map.is_empty());
    }

    // ---- assess_fasta_gc_scan: surface silent FASTA-GC failures ----

    #[test]
    fn gc_scan_ok_when_finite_values_produced() {
        // FASTA scanned, some bins resolved.
        let out = assess_fasta_gc_scan(true, false, true, 10, 7);
        assert_eq!(out, GcScanOutcome::Ok);
    }

    #[test]
    fn gc_scan_flags_reader_unavailable() {
        // Reference given but the index couldn't be opened.
        let out = assess_fasta_gc_scan(true, false, false, 0, 0);
        assert_eq!(out, GcScanOutcome::ReaderUnavailable);
    }

    #[test]
    fn gc_scan_flags_no_contig_resolved() {
        // FASTA opened, bins existed, but nothing resolved — naming mismatch.
        let out = assess_fasta_gc_scan(true, false, true, 10, 0);
        assert_eq!(out, GcScanOutcome::NoContigResolved);
    }

    #[test]
    fn gc_scan_ok_when_gc_bins_used() {
        // Pre-computed table path is not diagnosed here even if empty.
        let out = assess_fasta_gc_scan(true, true, true, 10, 0);
        assert_eq!(out, GcScanOutcome::Ok);
    }

    #[test]
    fn gc_scan_ok_when_no_reference() {
        let out = assess_fasta_gc_scan(false, false, false, 0, 0);
        assert_eq!(out, GcScanOutcome::Ok);
    }

    #[test]
    fn gc_scan_ok_when_no_bins_to_fill() {
        // Reader opened but header had no usable contigs — nothing to warn about.
        let out = assess_fasta_gc_scan(true, false, true, 0, 0);
        assert_eq!(out, GcScanOutcome::Ok);
    }

    // ---- end-to-end: GC correction proceeds with NO repeats provided ----

    #[test]
    fn gc_content_populated_when_no_repeats_provided() {
        // Regression guard for the decoupling fix: previously a chromosome
        // absent from the repeats BED got NaN GC (→ None), starving GC
        // correction. Now GC is scanned from the FASTA regardless of repeats.
        //
        // The chromosome must exceed to_bins()'s 25%-of-bin_size unmasked
        // threshold (250 kb) to survive, so use a 300 kb all-GC sequence.
        let seq_len = 300_000usize;
        let seq: String = std::iter::repeat_n("GC", seq_len / 2).collect();

        let base = std::env::temp_dir()
            .join(format!("nasvar_gc_norep_{}.fa", std::process::id()));
        let fa = base.to_string_lossy().into_owned();
        let fai = format!("{}.fai", fa);
        // Single-line FASTA: sequence starts after ">chr1\n" (6 bytes).
        std::fs::write(&fa, format!(">chr1\n{}\n", seq)).unwrap();
        std::fs::write(
            &fai,
            format!("chr1\t{}\t6\t{}\t{}\n", seq_len, seq_len, seq_len + 1),
        )
        .unwrap();

        let header = AlignmentHeader {
            text: String::new(),
            refs: vec!["chr1".to_string()],
            lengths: vec![seq_len as i32],
        };
        let contigs = vec![Contig {
            name: "chr1".into(),
            accession: "NC_000001.11".into(),
        }];

        // No repeats at all.
        let acc = CoverageAccumulator::new(&header, &[], Some(&fa), &contigs);
        let bins = acc.to_bins();

        std::fs::remove_file(&fa).ok();
        std::fs::remove_file(&fai).ok();

        assert_eq!(bins.len(), 1, "300 kb bin should survive the unmasked-fraction filter");
        let gc = bins[0].gc_content.expect("GC must be computed even without repeats");
        assert!((gc - 1.0).abs() < 1e-6, "all-GC sequence should read as 1.0, got {gc}");
    }
}

