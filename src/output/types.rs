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

/// Karyotype analysis results
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct KaryotypeOutput {
    /// Copy number per chromosome arm (e.g., "1p" -> 2)
    pub karyotype: HashMap<String, usize>,

    /// Median coverage per chromosome arm
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

    /// MAF peaks for copy number validation (CN -> MAF value)
    #[serde(skip_serializing_if = "Option::is_none")]
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
}

/// CNV analysis result for a single gene
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct CnvGeneResult {
    /// Focal (gene-level) copy number estimate
    #[serde(skip_serializing_if = "Option::is_none")]
    pub focal: Option<f64>,

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
    /// Insertion position on reference
    pub position: i64,

    /// Length of the insertion
    pub length: i64,

    /// Merged read count (supporting reads after clustering)
    pub merged: i64,

    /// Total coverage at the position
    pub coverage: i64,
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
