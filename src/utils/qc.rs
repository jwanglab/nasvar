//! QC data structures for pipeline metrics

use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use crate::config::Contig;
use crate::input::{AlignmentHeader, AlignmentRecord};
use crate::utils::bed::BedRegion;
use crate::utils::contig::ContigMapper;

/// QC metrics accumulated during pipeline execution
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PipelineQcData {
    pub nt_on_target: f64,
    pub reads_on_target: f64,
    pub target_regions_nt: f64,
    /// Median per-read quality (mean-error-probability definition, matching
    /// Dorado/MinKNOW) of primary reads overlapping enriched regions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub median_read_quality_enriched: Option<f64>,
}

/// Summary of per-read quality over a read set.
#[derive(Clone, Copy, Debug)]
pub struct ReadQStats {
    pub median: f64,
    pub reads: u64,
}

/// Per-read quality bins: 0.1-Q resolution over Q0..=Q60 (values above are
/// clamped into the top bin). 601 u64 counters, so memory is constant.
const Q_HIST_BINS: usize = 601;

/// Accumulates per-read quality scores using the mean-error-probability
/// definition: `Q_read = -10 * log10(mean(10^(-q_i/10)))` over the read's
/// base qualities. This matches the read Q reported by Dorado, MinKNOW,
/// NanoPlot, and cramino (arithmetic mean of Phred scores would bias high).
pub struct ReadQHist {
    /// 10^(-q/10) for q in 0..128 (BAM caps base quality at 93)
    err_lut: [f64; 128],
    hist: [u64; Q_HIST_BINS],
    count: u64,
}

impl Default for ReadQHist {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadQHist {
    pub fn new() -> Self {
        let mut err_lut = [0f64; 128];
        for (q, e) in err_lut.iter_mut().enumerate() {
            *e = 10f64.powf(-(q as f64) / 10.0);
        }
        Self {
            err_lut,
            hist: [0u64; Q_HIST_BINS],
            count: 0,
        }
    }

    /// Compute one read's quality from its base qualities and bin it.
    /// Reads with no quality string are skipped: BAM fills QUAL with 0xFF
    /// when absent, and our SAM decode path maps unparseable scores to 255.
    pub fn add(&mut self, qual: &[u8]) {
        if qual.is_empty() || qual[0] == 0xFF {
            return;
        }
        let mut err_sum = 0f64;
        for &b in qual {
            err_sum += self.err_lut[(b & 0x7F) as usize];
        }
        let read_q = -10.0 * (err_sum / qual.len() as f64).log10();
        let bin = ((read_q * 10.0) as usize).min(Q_HIST_BINS - 1);
        self.hist[bin] += 1;
        self.count += 1;
    }

    /// Median of per-read quality, or None if no reads were added.
    /// Exact to the 0.1-Q bin resolution.
    pub fn stats(&self) -> Option<ReadQStats> {
        if self.count == 0 {
            return None;
        }
        let half = self.count.div_ceil(2);
        let mut cum = 0u64;
        let mut median = 0f64;
        for (i, &c) in self.hist.iter().enumerate() {
            cum += c;
            if cum >= half {
                median = i as f64 * 0.1;
                break;
            }
        }
        Some(ReadQStats {
            median,
            reads: self.count,
        })
    }
}

/// Accumulates whole-read quality for primary reads whose alignment overlaps
/// any enriched region. Membership is an interval test on the alignment span
/// (no CIGAR clipping): a read either contributes all of its base qualities
/// or none, since per-read Q is a property of the read.
pub struct EnrichedReadQ {
    /// ref_id -> merged, sorted (start, end) intervals, 0-based half-open
    by_ref: HashMap<usize, Vec<(i64, i64)>>,
    q: ReadQHist,
}

impl EnrichedReadQ {
    pub fn new(header: &AlignmentHeader, regions: &[BedRegion], contigs: &[Contig]) -> Self {
        let mapper = ContigMapper::from_contigs_and_refs(contigs, &header.refs);
        let mut by_ref: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
        for r in regions {
            let bam_chrom = mapper.to_bam_name(&r.segment);
            if let Some(ref_id) = header.refs.iter().position(|n| n == &bam_chrom) {
                by_ref
                    .entry(ref_id)
                    .or_default()
                    .push((r.start as i64, r.end as i64));
            }
        }
        // Sort and merge per ref: merged intervals have strictly increasing
        // starts AND ends, so the overlap lookup is a single binary search.
        for ivs in by_ref.values_mut() {
            ivs.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(ivs.len());
            for &(s, e) in ivs.iter() {
                match merged.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            *ivs = merged;
        }
        Self {
            by_ref,
            q: ReadQHist::new(),
        }
    }

    pub fn process(&mut self, record: &AlignmentRecord) {
        // Primary aligned reads only (not unmapped/secondary/supplementary),
        // matching QcAccumulator's reads_aligned semantics. Supplementary
        // records carry hard-clipped partial quals and secondaries usually
        // have none, so counting them would skew the distribution.
        if (record.flags() & 0x904) != 0 || record.ref_id < 0 || record.pos < 0 {
            return;
        }
        let Some(intervals) = self.by_ref.get(&(record.ref_id as usize)) else {
            return;
        };
        let start = record.pos as i64;
        let end = start + record.alignment_span().max(1) as i64;
        // Only the last interval starting before `end` can overlap.
        let idx = intervals.partition_point(|&(s, _)| s < end);
        if idx == 0 || intervals[idx - 1].1 <= start {
            return;
        }
        self.q.add(&record.qual);
    }

    pub fn stats(&self) -> Option<ReadQStats> {
        self.q.stats()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::CigarKind;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn test_readq_uniform() {
        // All bases Q20 -> read Q exactly 20
        let mut h = ReadQHist::new();
        h.add(&vec![20u8; 1000]);
        let s = h.stats().unwrap();
        assert_eq!(s.reads, 1);
        assert!(approx(s.median, 20.0, 0.1 + 1e-9), "median {}", s.median);
    }

    #[test]
    fn test_readq_mean_error_definition() {
        // Half Q10, half Q30: mean err = (0.1 + 0.001)/2 = 0.0505
        // -> Q = -10*log10(0.0505) = 12.9666... -> 0.1-bin 12.9
        // (arithmetic mean of scores would give 20 — the definitions differ)
        let mut qual = vec![10u8; 500];
        qual.extend(vec![30u8; 500]);
        let mut h = ReadQHist::new();
        h.add(&qual);
        let s = h.stats().unwrap();
        assert!(approx(s.median, 12.9, 1e-9), "median {}", s.median);
    }

    #[test]
    fn test_readq_median() {
        let mut h = ReadQHist::new();
        h.add(&vec![10u8; 100]);
        h.add(&vec![15u8; 100]);
        h.add(&vec![40u8; 100]);
        let s = h.stats().unwrap();
        assert_eq!(s.reads, 3);
        assert!(approx(s.median, 15.0, 0.1 + 1e-9), "median {}", s.median);
    }

    #[test]
    fn test_readq_skips_missing_qual() {
        let mut h = ReadQHist::new();
        h.add(&[]);
        h.add(&[0xFF; 50]); // BAM's "no quality" fill
        assert!(h.stats().is_none());
    }

    fn test_header() -> AlignmentHeader {
        AlignmentHeader {
            text: String::new(),
            refs: vec!["chr1".to_string(), "chr2".to_string()],
            lengths: vec![1_000_000, 1_000_000],
        }
    }

    fn test_contigs() -> Vec<Contig> {
        vec![
            Contig { name: "chr1".to_string(), accession: "NC_000001.11".to_string() },
            Contig { name: "chr2".to_string(), accession: "NC_000002.12".to_string() },
        ]
    }

    fn mapped_record(ref_id: i32, pos: i32, len: usize, flag: u16) -> AlignmentRecord {
        AlignmentRecord {
            ref_id,
            pos,
            flag,
            qual: vec![20u8; len],
            cigar: vec![(CigarKind::Match, len)],
            ..Default::default()
        }
    }

    fn enriched_regions() -> Vec<BedRegion> {
        vec![
            BedRegion { segment: "chr1".to_string(), start: 1000, end: 2000, name: "r1".to_string() },
            // Overlapping entries must merge
            BedRegion { segment: "chr1".to_string(), start: 1500, end: 3000, name: "r2".to_string() },
            BedRegion { segment: "chr2".to_string(), start: 500, end: 600, name: "r3".to_string() },
        ]
    }

    #[test]
    fn test_enriched_overlap() {
        let mut acc = EnrichedReadQ::new(&test_header(), &enriched_regions(), &test_contigs());

        acc.process(&mapped_record(0, 900, 200, 0)); // spans into [1000,3000): counted
        acc.process(&mapped_record(0, 2999, 100, 0)); // starts on last base: counted
        acc.process(&mapped_record(0, 3000, 100, 0)); // starts at merged end: not counted
        acc.process(&mapped_record(0, 0, 1000, 0)); // ends at region start: not counted
        acc.process(&mapped_record(1, 550, 10, 0)); // inside chr2 region: counted
        acc.process(&mapped_record(0, 100_000, 100, 0)); // far away: not counted

        assert_eq!(acc.stats().unwrap().reads, 3);
    }

    #[test]
    fn test_enriched_reads_quality_from_base_quals() {
        let mut acc = EnrichedReadQ::new(&test_header(), &enriched_regions(), &test_contigs());
        acc.process(&mapped_record(0, 1500, 100, 0)); // Q20 bases
        let s = acc.stats().unwrap();
        assert_eq!(s.reads, 1);
        assert!(approx(s.median, 20.0, 0.1 + 1e-9), "median {}", s.median);
    }

    #[test]
    fn test_enriched_skips_non_primary() {
        let mut acc = EnrichedReadQ::new(&test_header(), &enriched_regions(), &test_contigs());
        acc.process(&mapped_record(0, 1500, 100, 0x100)); // secondary
        acc.process(&mapped_record(0, 1500, 100, 0x800)); // supplementary
        acc.process(&mapped_record(0, 1500, 100, 0x4)); // unmapped
        acc.process(&mapped_record(-1, -1, 100, 0)); // no ref
        assert!(acc.stats().is_none());
    }
}
