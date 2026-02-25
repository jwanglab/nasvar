use bitvec::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

use crate::bam::ContigMapper;
use crate::input::{AlignmentInput, AlignmentHeader, AlignmentRecord};
use crate::utils::bed::BedRegion;
use log::{info, debug};

// ==================== Shared CoverageAccumulator ====================
// Used by both standalone read_depth() and pipeline mode

struct ChromData {
    name: String,
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

impl CoverageAccumulator {
    /// Create new accumulator with repeat masks and optional GC content.
    pub fn new(header: &AlignmentHeader, repeats: &[BedRegion], ref_path: Option<&str>) -> Self {
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
            ContigMapper::from_fai(&fai_path).ok()
        });

        let mut chroms = Vec::with_capacity(header.refs.len());
        let mapper = ContigMapper::new();

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

                // Compute per-bin GC content from unmasked bases
                let bin_gc_content = Self::compute_gc_content(
                    &mut fasta_reader, name, fasta_mapper.as_ref(), len, n_bins, bin_size, &mask,
                );

                chroms.push(Some(ChromData {
                    name: name.clone(),
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

                // Skip GC content for chromosomes without repeats (much faster)
                // GC will be NA for these bins
                let bin_gc_content = vec![f32::NAN; n_bins];

                chroms.push(Some(ChromData {
                    name: name.clone(),
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

    /// Write coverage output to file.
    /// If include_gc is true, writes GC content column (pipeline format).
    /// If include_gc is false, writes 4-column format (standalone format).
    pub fn write_output(&self, path: &str, include_gc: bool) -> std::io::Result<()> {
        let mut file = File::create(path)?;

        if include_gc {
            writeln!(file, "chromosome\tbin_start\tbin_end\tn_reads\tgc_content")?;
        } else {
            writeln!(file, "chromosome\tbin_start\tbin_end\tn_reads")?;
        }

        for data in self.chroms.iter().flatten() {
            debug!("Writing coverage for {}...", data.name);

            for (i, (&unmasked, &reads)) in data
                .bin_unmasked_counts
                .iter()
                .zip(data.bin_read_counts.iter())
                .enumerate()
            {
                let start = i * self.bin_size;
                let end = std::cmp::min((i + 1) * self.bin_size, data.len);

                let threshold = (self.bin_size as f64 * 0.25) as u32;
                if unmasked > threshold {
                    let mut value = reads as f64;
                    value *= self.bin_size as f64 / unmasked as f64;

                    if include_gc {
                        let gc = data.bin_gc_content[i];
                        if gc.is_nan() {
                            writeln!(
                                file,
                                "{}\t{}\t{}\t{:.0}\tNA",
                                data.mapped_name, start, end, value
                            )?;
                        } else {
                            writeln!(
                                file,
                                "{}\t{}\t{}\t{:.0}\t{:.4}",
                                data.mapped_name, start, end, value, gc
                            )?;
                        }
                    } else {
                        // Truncate depth to integer for output
                        let value_int = value as u64;
                        writeln!(file, "{}\t{}\t{}\t{}", data.mapped_name, start, end, value_int)?;
                    }
                }
            }
        }
        Ok(())
    }
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
) -> Result<u64, Box<dyn std::error::Error>>
{
    info!("Calculating coverage...");

    // Create accumulator (no reference = no GC content)
    let mut accumulator = CoverageAccumulator::new(&bam.header, repeats, None);

    // Single-pass BAM scan
    bam.seek(bam.start_pos)?;

    while let Some(record) = bam.read_record()? {
        accumulator.process(&record);
    }

    // Write output (standalone format without GC column)
    let out_file = format!("{}.coverage.tsv", out_prefix);
    accumulator.write_output(&out_file, false)?;

    let reads_aligned = accumulator.reads_aligned();
    info!("Total aligned reads (primary): {}", reads_aligned);
    Ok(reads_aligned)
}

