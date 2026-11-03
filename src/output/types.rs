//! Output data structures for unified pipeline output
//!
//! This module defines the structures used to collect and serialize
//! all pipeline outputs into a unified JSON format.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::utils::metadata::SequencingMetaData;
use crate::utils::qc::PipelineQcData;

/// Top-level unified output structure containing all pipeline results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct UnifiedOutput {
    /// Pipeline version
    pub version: String,

    /// Timestamp of analysis (ISO 8601 format)
    pub timestamp: String,

    /// Sequencing metadata from BAM headers
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<SequencingMetaData>,

    /// QC metrics
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qc: Option<QcOutput>,

    /// Gene fusion results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fusions: Option<FusionsOutput>,

    /// Karyotype analysis results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub karyotype: Option<KaryotypeOutput>,

    /// Copy number variation results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cnv: Option<CnvOutput>,

    /// Single nucleotide variant results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snv: Option<SnvOutput>,

    /// Internal tandem duplication results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub itd: Option<ItdOutput>,

    /// Fusion breakpoint consensus sequences
    #[serde(skip_serializing_if = "Option::is_none")]
    pub breakpoint_consensus: Option<BreakpointConsensusOutput>,
}

// ============================================================================
// QC Output
// ============================================================================

/// QC metrics output - extends PipelineQcData with additional metrics
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct QcOutput {
    /// Total nucleotides aligned to target regions
    pub nt_on_target: f64,

    /// Number of reads overlapping target regions
    pub reads_on_target: f64,

    /// Total nucleotides in target region definitions
    pub target_regions_nt: f64,

    /// Total primary aligned reads genome-wide
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reads_aligned: Option<u64>,

    /// Mean coverage across targets
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_coverage: Option<f64>,

    /// Per-target average coverage (gene name -> mean depth)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub target_coverage: Option<HashMap<String, f64>>,

    /// Median per-read quality (mean-error-probability definition) of primary
    /// reads overlapping enriched regions
    #[serde(skip_serializing_if = "Option::is_none")]
    pub median_read_quality_enriched: Option<f64>,
}

impl From<PipelineQcData> for QcOutput {
    fn from(qc: PipelineQcData) -> Self {
        Self {
            nt_on_target: qc.nt_on_target,
            reads_on_target: qc.reads_on_target,
            target_regions_nt: qc.target_regions_nt,
            reads_aligned: None,
            mean_coverage: None,
            target_coverage: None,
            median_read_quality_enriched: qc.median_read_quality_enriched,
        }
    }
}

// ============================================================================
// Fusion Output
// ============================================================================

/// Container for fusion calling results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct FusionsOutput {
    /// Called gene fusions
    pub fusions: Vec<FusionEvent>,

    /// Spike-in control fusions (if detected)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spike_in: Option<Vec<FusionEvent>>,
}

/// A single fusion event between two genes
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct FusionEvent {
    /// First gene partner
    pub gene1: GeneInfo,

    /// Second gene partner
    pub gene2: GeneInfo,

    /// Number of reads supporting the fusion
    pub supporting_reads: usize,

    /// Number of reads in repetitive regions
    #[serde(default)]
    pub repetitive_reads: usize,

    /// Breakpoint locations for this fusion
    pub breakpoints: Vec<FusionBreakpoint>,
}

/// Gene location information
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct GeneInfo {
    /// Gene name (empty string if unknown/intergenic)
    pub name: String,

    /// Chromosome
    pub chr: String,

    /// Position on chromosome
    pub pos: u32,
}

/// Direction of read support relative to a breakpoint position.
///
/// `Left` means reads extend leftward from the breakpoint (breakpoint at template end).
/// `Right` means reads extend rightward from the breakpoint (breakpoint at template start).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum BreakpointDirection {
    Left,
    Right,
}

impl BreakpointDirection {
    /// Returns the arrow text representation for report display.
    pub fn arrow(&self) -> &'static str {
        match self {
            BreakpointDirection::Left => "<<",
            BreakpointDirection::Right => ">>",
        }
    }
}

/// A fusion breakpoint location
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct FusionBreakpoint {
    /// Gene 0 name
    pub gene0_name: String,

    /// Gene 0 chromosome
    pub gene0_chr: String,

    /// Gene 0 breakpoint position
    pub gene0_pos: u32,

    /// Gene 0 direction
    pub gene0_dir: BreakpointDirection,

    /// Gene 1 name
    pub gene1_name: String,

    /// Gene 1 chromosome
    pub gene1_chr: String,

    /// Gene 1 breakpoint position
    pub gene1_pos: u32,

    /// Gene 1 direction
    pub gene1_dir: BreakpointDirection,

    /// Number of reads supporting this breakpoint
    pub n_reads: usize,

    /// Median query-space overlap/gap at junction across reads.
    /// Positive = microhomology (shared bases). Negative = insertion (unaligned bases).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[schemars(skip)]
    pub overlap_gap: Option<i32>,
}

// ============================================================================
// Karyotype Output
// ============================================================================

/// Genome order of a karyotype segment name: 1p, 1q, 2p, ..., 13, ..., 22, X, Y,
/// then anything else (e.g. PAR) by name. An optional "chr" prefix is ignored.
fn segment_order(seg: &str) -> (u8, u8, &str) {
    let s = seg.strip_prefix("chr").unwrap_or(seg);
    let (chrom, arm) = match s.strip_suffix('p') {
        Some(c) => (c, 0),
        None => match s.strip_suffix('q') {
            Some(c) => (c, 1),
            None => (s, 0),
        },
    };
    let rank = match chrom {
        "X" => 23,
        "Y" => 24,
        c => match c.parse::<u8>() {
            Ok(n) if (1..=22).contains(&n) => n,
            _ => return (u8::MAX, 0, seg),
        },
    };
    (rank, arm, "")
}

/// Writes a segment-keyed map in genome order (HashMap order is arbitrary).
fn serialize_segments<S: serde::Serializer, V: Serialize>(
    map: &HashMap<String, V>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    use serde::ser::SerializeMap;
    let mut entries: Vec<_> = map.iter().collect();
    entries.sort_by(|a, b| segment_order(a.0).cmp(&segment_order(b.0)));
    let mut m = serializer.serialize_map(Some(entries.len()))?;
    for (k, v) in entries {
        m.serialize_entry(k, v)?;
    }
    m.end()
}

/// Writes an optional map with its keys sorted.
fn serialize_opt_sorted<S: serde::Serializer, K: Ord + Serialize, V: Serialize>(
    map: &Option<HashMap<K, V>>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    match map {
        Some(m) => serializer.collect_map(m.iter().collect::<std::collections::BTreeMap<_, _>>()),
        None => serializer.serialize_none(),
    }
}

/// Karyotype analysis results
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct KaryotypeOutput {
    /// Copy number per chromosome arm (e.g., "1p" -> 2), written in genome order
    #[serde(serialize_with = "serialize_segments")]
    pub karyotype: HashMap<String, usize>,

    /// Median coverage per chromosome arm, written in genome order
    #[serde(serialize_with = "serialize_segments")]
    pub medians: HashMap<String, f64>,

    /// Coverage levels detected (coverage_value, copy_number)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub levels_found: Option<Vec<(f64, usize)>>,

    /// Human-readable karyotype string (e.g., "46; XX")
    pub karyotype_string: String,

    /// ISCN nomenclature string
    pub iscn_string: String,

    /// Estimated blast/tumor ratio (0.0 - 1.0)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub blast_ratio: Option<f64>,

    /// Total aligned reads (for QC)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reads_aligned: Option<u64>,

    /// Within-segment spread metric
    #[serde(skip_serializing_if = "Option::is_none")]
    pub within_segment_spread: Option<f64>,

    /// Warning messages
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warnings: Option<String>,

    /// MAF peaks for copy number validation (CN -> MAF value), written by CN
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "serialize_opt_sorted")]
    pub maf_peaks: Option<HashMap<usize, f64>>,
}

// ============================================================================
// CNV Output
// ============================================================================

/// Copy number variation results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct CnvOutput {
    /// CNV results per gene
    pub genes: HashMap<String, CnvGeneResult>,

    /// How focal CN was estimated: "read_count_nb", or "coverage" when the NB
    /// could not be fitted
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_method: Option<String>,

    /// Blast ratio used for focal CN (1.0 when none was provided or estimated)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_blast_ratio: Option<f64>,

    /// Diploid baseline genes for focal CN
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_baseline_genes: Option<usize>,

    /// Negative binomial fit of the baseline genes' read counts
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_nb_fit: Option<FocalNbFit>,

    /// Warnings shown under "Focal copy-number variation" in the report
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warnings: Option<Vec<String>>,
}

/// Negative binomial fit of read counts: reads ~ NB(mu, mu + alpha mu^2),
/// mu = reads_per_base * (target length + length_offset)
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct FocalNbFit {
    pub reads_per_base: f64,
    pub length_offset: f64,
    pub alpha: f64,
}

/// CNV analysis result for a single gene
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct CnvGeneResult {
    /// Focal (gene-level) copy number estimate, adjusted for blast ratio,
    /// clipped at 0. From read counts, or coverage if the NB was not fitted.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal: Option<f64>,

    /// Copy number unadjusted for blast percentage (2 * reads / expected reads at 2n, or
    /// 2 * depth / baseline depth). Independent of blast ratio.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_ratio: Option<f64>,

    /// Two-sided negative binomial p-value of the read count under 2n
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_p_2n: Option<f64>,

    /// Primary reads with >= 1 aligned base in the gene's targets
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_reads: Option<u64>,

    /// Expected reads at 2n under the NB fit
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_expected_reads: Option<f64>,

    /// Coverage-based focal CN (the previous method), kept for comparison
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_coverage: Option<f64>,

    /// Set when the rounded focal CN is not 2 but p(2n) >= 0.05
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal_flag: Option<String>,

    /// Local (region-level) copy number estimate
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local: Option<f64>,

    /// Detected intragenic deletions
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deletions: Option<Vec<StructuralVariant>>,

    /// Detected intragenic duplications
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duplications: Option<Vec<StructuralVariant>>,
}

/// A structural variant (deletion or duplication)
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct StructuralVariant {
    /// Chromosome
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chrom: Option<String>,

    /// Start position (0-based)
    pub start: i64,

    /// End position (0-based)
    pub end: i64,

    /// Number of supporting reads
    pub reads: usize,
}

// ============================================================================
// SNV Output
// ============================================================================

/// Single nucleotide variant results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct SnvOutput {
    /// SNV results per gene
    pub genes: HashMap<String, SnvGeneResult>,
}

/// SNV analysis result for a single gene
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SnvGeneResult {
    /// Called genotype (e.g., "*1/*3A")
    pub genotype: String,

    /// Average sequencing depth across the gene
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coverage: Option<f64>,

    /// Detected mutations (position+change -> allele frequency)
    #[serde(default)]
    pub mutations: HashMap<String, String>,

    /// Amino acid changes (e.g., ["p.Ile105Thr"])
    #[serde(default)]
    pub aa_changes: Vec<String>,

    /// Phase information
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase: Option<serde_json::Value>,

    /// Unassigned mutations (couldn't be phased)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unassigned_muts: Option<Vec<String>>,
}

// ============================================================================
// ITD Output
// ============================================================================

/// Internal tandem duplication results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct ItdOutput {
    /// ITD events per gene
    pub genes: HashMap<String, Vec<ItdEvent>>,
}

/// A single ITD event
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ItdEvent {
    /// Chromosome name (from the ITD region config). Carried on the event
    /// so downstream consumers don't need the pipeline config to render a
    /// "view alignments" link.
    #[serde(default)]
    pub chrom: String,

    /// Insertion position on reference
    pub position: i64,

    /// Length of the insertion
    pub length: i64,

    /// Merged read count (supporting reads after clustering)
    pub merged: i64,

    /// Total coverage at the position
    pub coverage: i64,

    /// Amino acid position of the insertion (1-based), if a transcript was configured
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aa_position: Option<i64>,

    /// Reference amino acid at the insertion position (single-letter code)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aa_ref: Option<char>,
}

// ============================================================================
// Fusion Breakpoint Consensus Output
// ============================================================================

/// Container for fusion breakpoint consensus results
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct BreakpointConsensusOutput {
    /// Consensus results for each breakpoint
    pub breakpoints: Vec<FusionBreakpointConsensus>,
}

/// Consensus sequence reconstructed across a single fusion breakpoint
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct FusionBreakpointConsensus {
    /// Gene 0 name
    pub gene0_name: String,

    /// Gene 0 chromosome
    pub gene0_chr: String,

    /// Gene 0 breakpoint position
    pub gene0_pos: u32,

    /// Gene 1 name
    pub gene1_name: String,

    /// Gene 1 chromosome
    pub gene1_chr: String,

    /// Gene 1 breakpoint position
    pub gene1_pos: u32,

    /// Full consensus sequence spanning the breakpoint
    pub consensus_sequence: String,

    /// Primer-design-ready format: GENE0_SEQ[INSERTED]GENE1_SEQ
    pub bracketed_sequence: String,

    /// Inserted bases at the junction (empty string if none)
    pub inserted_bases: String,

    /// Number of reads used in consensus
    pub n_reads: usize,

    /// Mean per-position coverage across the consensus
    pub mean_coverage: f64,
}

#[cfg(test)]
mod segment_order_tests {
    use super::*;

    #[test]
    fn karyotype_maps_serialize_in_genome_order() {
        let segs = ["X", "22", "10q", "1q", "2p", "13", "Y", "1p", "10p", "PAR1", "21"];
        // Built from JSON so the test does not depend on the optional fields
        let karyo: KaryotypeOutput = serde_json::from_value(serde_json::json!({
            "karyotype": segs.iter().map(|s| (s.to_string(), 2)).collect::<HashMap<_, _>>(),
            "medians": segs.iter().map(|s| (s.to_string(), 1.0)).collect::<HashMap<_, _>>(),
            "karyotype_string": "",
            "iscn_string": "",
            "maf_peaks": {"3": 0.33, "1": 0.0, "2": 0.5},
        }))
        .unwrap();
        let json = serde_json::to_string(&karyo).unwrap();
        let order = ["1p", "1q", "2p", "10p", "10q", "13", "21", "22", "X", "Y", "PAR1"];
        for field in ["\"karyotype\":{", "\"medians\":{"] {
            let body = &json[json.find(field).unwrap()..];
            let pos: Vec<usize> = order.iter().map(|s| body.find(&format!("\"{}\":", s)).unwrap()).collect();
            assert!(pos.windows(2).all(|w| w[0] < w[1]), "{} order wrong: {}", field, body);
        }
        let peaks = &json[json.find("\"maf_peaks\":").unwrap()..];
        assert!(peaks.find("\"1\":").unwrap() < peaks.find("\"2\":").unwrap());
        assert!(peaks.find("\"2\":").unwrap() < peaks.find("\"3\":").unwrap());
        // Round trip is unaffected
        let back: KaryotypeOutput = serde_json::from_str(&json).unwrap();
        assert_eq!(back.karyotype.len(), segs.len());
    }
}
