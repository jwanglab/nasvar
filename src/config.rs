//! Unified configuration system for nasvar pipeline.
//!
//! This module provides structs and loading functions for:
//! - Pipeline configuration (genes, thresholds)
//! - Reference genome data (centromeres, PAR regions)

use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

// ============================================================================
// Reference Genome Configuration
// ============================================================================

/// Reference genome data including centromeres and PAR regions
#[derive(Deserialize, Debug, Clone)]
pub struct ReferenceConfig {
    pub name: String,
    pub centromeres: HashMap<String, (u32, u32)>,
    pub par1: HashMap<String, (u32, u32)>,
    pub par2: HashMap<String, (u32, u32)>,
}

impl ReferenceConfig {
    /// Load reference configuration from a JSON file
    pub fn load(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let config: ReferenceConfig = serde_json::from_reader(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok(config)
    }
}

// ============================================================================
// Pipeline Configuration
// ============================================================================

/// ITD region configuration
#[derive(Deserialize, Debug, Clone)]
pub struct ItdRegion {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    #[serde(default = "default_min_length")]
    pub min_length: isize,
    #[serde(default)]
    pub min_frequency: f64,
}

fn default_min_length() -> isize { 3 }

/// CNV gene configuration
#[derive(Deserialize, Debug, Clone)]
pub struct CnvGeneConfig {
    pub focal_genes: Vec<String>,
    pub exclude_from_baseline: Vec<String>,
    pub deletion_genes: Vec<String>,
    pub duplication_genes: Vec<String>,
    /// Genes for which to display local (region-level) CN instead of focal CN in reports
    #[serde(default)]
    pub local_cn_genes: Vec<String>,
}

/// Spike-in control fusion configuration
#[derive(Deserialize, Debug, Clone)]
pub struct SpikeInConfig {
    pub gene1: String,
    pub gene2: String,
}

/// Fusion gene configuration
#[derive(Deserialize, Debug, Clone)]
pub struct FusionGeneConfig {
    #[serde(default)]
    pub special_margins: HashMap<String, u32>,
    #[serde(default)]
    pub skip_self_fusions: Vec<String>,
    #[serde(default)]
    pub blacklist_pairs: Vec<(String, String)>,
    /// Genes to report one-sided fusions for (promiscuous fusion partners)
    #[serde(default)]
    pub one_sided_genes: Vec<String>,
    /// Allowed fusion partners for specific one-sided genes (gene -> list of allowed partners)
    /// If a gene has partners configured, only fusions landing near those partners are reported
    #[serde(default)]
    pub one_sided_partners: HashMap<String, Vec<String>>,
    /// Genes that bypass the per-breakpoint read count filter (e.g., DUX4)
    /// Uses substring matching (like skip_self_fusions)
    #[serde(default)]
    pub bypass_breakpoint_filter: Vec<String>,
    /// Spike-in control fusion construct to separate from real fusions
    #[serde(default)]
    pub spike_in: Option<SpikeInConfig>,
}

/// SNV transcript configuration for a single gene
#[derive(Deserialize, Debug, Clone)]
pub struct SnvTranscriptConfig {
    /// Transcript ID (e.g., "rna-NM_001346817.1")
    #[serde(default)]
    pub transcript: Option<String>,
    /// Variant definitions: name -> list of CDS positions
    #[serde(default)]
    pub variants: HashMap<String, Vec<usize>>,
}

/// SNV gene configuration
#[derive(Deserialize, Debug, Clone)]
pub struct SnvGeneConfig {
    /// List of genes for pharmacogenomics reporting
    #[serde(default)]
    pub pharmacogenomics: Vec<String>,
    /// List of genes for pathogenic variant reporting
    #[serde(default)]
    pub pathogenic: Vec<String>,
    /// Transcript and variant definitions per gene
    #[serde(default)]
    pub transcripts: HashMap<String, SnvTranscriptConfig>,
}

/// All gene configurations
#[derive(Deserialize, Debug, Clone)]
pub struct GeneConfig {
    pub cnv: CnvGeneConfig,
    pub itd: HashMap<String, ItdRegion>,
    pub snv: SnvGeneConfig,
    pub fusions: FusionGeneConfig,
}

/// One-sided fusion specific thresholds
#[derive(Deserialize, Debug, Clone, Default)]
pub struct OneSidedThresholds {
    /// Min breakpoint reads for one-sided fusions (defaults to general min_breakpoint_reads if not set)
    #[serde(default)]
    pub min_breakpoint_reads: Option<i64>,
    /// Min fraction of supporting reads / average depth for one-sided fusions
    #[serde(default)]
    pub min_fraction: Option<f64>,
}

/// Fusion threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct FusionThresholds {
    #[serde(default = "default_fusion_margin")]
    pub default_margin: u32,
    #[serde(default = "default_min_anchor")]
    pub min_anchor: u32,
    #[serde(default = "default_min_supporting_reads")]
    pub min_supporting_reads: i64,
    #[serde(default = "default_min_breakpoint_reads")]
    pub min_breakpoint_reads: i64,
    /// One-sided fusion specific thresholds
    #[serde(default)]
    pub one_sided: OneSidedThresholds,
}

fn default_fusion_margin() -> u32 { 5000 }
fn default_min_anchor() -> u32 { 500 }
fn default_min_supporting_reads() -> i64 { 3 }
fn default_min_breakpoint_reads() -> i64 { 3 }

/// SNV threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct SnvThresholds {
    #[serde(default = "default_min_coverage")]
    pub min_coverage: f64,
    #[serde(default = "default_min_allele_freq")]
    pub min_allele_freq: f64,
    #[serde(default = "default_homozygosity_vaf")]
    pub homozygosity_vaf: f64,
    #[serde(default = "default_phasing_ratio")]
    pub phasing_ratio: i32,
}

fn default_min_coverage() -> f64 { 5.0 }
fn default_min_allele_freq() -> f64 { 0.3 }
fn default_homozygosity_vaf() -> f64 { 0.8 }
fn default_phasing_ratio() -> i32 { 10 }

/// CNV threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct CnvThresholds {
    #[serde(default = "default_gene_margin")]
    pub gene_margin: u32,
    #[serde(default = "default_cnv_min_supporting_reads")]
    pub min_supporting_reads: usize,
    #[serde(default = "default_breakpoint_margin")]
    pub breakpoint_margin: u32,
    #[serde(default = "default_cnv_min_anchor")]
    pub min_anchor: u32,
}

fn default_gene_margin() -> u32 { 50000 }
fn default_cnv_min_supporting_reads() -> usize { 3 }
fn default_breakpoint_margin() -> u32 { 50 }
fn default_cnv_min_anchor() -> u32 { 500 }

/// MAF threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct MafThresholds {
    #[serde(default = "default_min_depth")]
    pub min_depth: u32,
}

fn default_min_depth() -> u32 { 20 }

/// Coverage threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct CoverageThresholds {
    #[serde(default = "default_bin_size")]
    pub bin_size: usize,
}

fn default_bin_size() -> usize { 1_000_000 }

/// Karyotype threshold configuration
#[derive(Deserialize, Debug, Clone)]
pub struct KaryotypeThresholds {
    #[serde(default = "default_spread_warning")]
    pub spread_warning: f64,
    #[serde(default = "default_level_tolerance")]
    pub level_tolerance: f64,
    #[serde(default = "default_maf_peak_diploid")]
    pub maf_peak_diploid: f64,
    #[serde(default = "default_min_maf_sites")]
    pub min_maf_sites: usize,
}

fn default_spread_warning() -> f64 { 0.08 }
fn default_level_tolerance() -> f64 { 0.3 }
fn default_maf_peak_diploid() -> f64 { 0.4 }
fn default_min_maf_sites() -> usize { 50 }

/// All threshold configurations
#[derive(Deserialize, Debug, Clone)]
#[derive(Default)]
pub struct ThresholdConfig {
    #[serde(default)]
    pub fusions: FusionThresholds,
    #[serde(default)]
    pub snv: SnvThresholds,
    #[serde(default)]
    pub cnv: CnvThresholds,
    #[serde(default)]
    pub maf: MafThresholds,
    #[serde(default)]
    pub coverage: CoverageThresholds,
    #[serde(default)]
    pub karyotype: KaryotypeThresholds,
}

impl Default for FusionThresholds {
    fn default() -> Self {
        Self {
            default_margin: default_fusion_margin(),
            min_anchor: default_min_anchor(),
            min_supporting_reads: default_min_supporting_reads(),
            min_breakpoint_reads: default_min_breakpoint_reads(),
            one_sided: OneSidedThresholds::default(),
        }
    }
}

impl Default for SnvThresholds {
    fn default() -> Self {
        Self {
            min_coverage: default_min_coverage(),
            min_allele_freq: default_min_allele_freq(),
            homozygosity_vaf: default_homozygosity_vaf(),
            phasing_ratio: default_phasing_ratio(),
        }
    }
}

impl Default for CnvThresholds {
    fn default() -> Self {
        Self {
            gene_margin: default_gene_margin(),
            min_supporting_reads: default_cnv_min_supporting_reads(),
            breakpoint_margin: default_breakpoint_margin(),
            min_anchor: default_cnv_min_anchor(),
        }
    }
}

impl Default for MafThresholds {
    fn default() -> Self {
        Self {
            min_depth: default_min_depth(),
        }
    }
}

impl Default for CoverageThresholds {
    fn default() -> Self {
        Self {
            bin_size: default_bin_size(),
        }
    }
}

impl Default for KaryotypeThresholds {
    fn default() -> Self {
        Self {
            spread_warning: default_spread_warning(),
            level_tolerance: default_level_tolerance(),
            maf_peak_diploid: default_maf_peak_diploid(),
            min_maf_sites: default_min_maf_sites(),
        }
    }
}


/// Main pipeline configuration
#[derive(Deserialize, Debug, Clone)]
pub struct PipelineConfig {
    /// Path to reference genome configuration file
    #[serde(default)]
    pub reference_file: Option<String>,
    /// Gene configurations
    pub genes: GeneConfig,
    /// Threshold configurations
    #[serde(default)]
    pub thresholds: ThresholdConfig,
}

impl PipelineConfig {
    /// Load pipeline configuration from a JSON file
    pub fn load(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let config: PipelineConfig = serde_json::from_reader(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok(config)
    }

    /// Get gene margin for fusion calling, using special margin if configured
    pub fn get_fusion_margin(&self, gene_name: &str) -> u32 {
        self.genes.fusions.special_margins
            .get(gene_name)
            .copied()
            .unwrap_or(self.thresholds.fusions.default_margin)
    }

    /// Check if a gene should skip self-fusion filtering
    pub fn skip_self_fusion(&self, gene_name: &str) -> bool {
        self.genes.fusions.skip_self_fusions.iter().any(|g| gene_name.contains(g))
    }

    /// Check if a gene should bypass per-breakpoint read filtering
    pub fn bypass_breakpoint_filter(&self, gene_name: &str) -> bool {
        self.genes.fusions.bypass_breakpoint_filter.iter().any(|g| gene_name.contains(g))
    }

    /// Check if a fusion pair is blacklisted
    pub fn is_blacklisted_pair(&self, gene1: &str, gene2: &str) -> bool {
        for (a, b) in &self.genes.fusions.blacklist_pairs {
            if (gene1 == a && gene2 == b) || (gene1 == b && gene2 == a) {
                return true;
            }
        }
        false
    }

    /// Get the one-sided fusion genes as a HashSet
    pub fn get_one_sided_genes(&self) -> std::collections::HashSet<String> {
        self.genes.fusions.one_sided_genes.iter().cloned().collect()
    }

    /// Get allowed fusion partners for a one-sided gene, if configured
    pub fn get_one_sided_partners(&self, gene_name: &str) -> Option<&Vec<String>> {
        self.genes.fusions.one_sided_partners.get(gene_name)
    }

    /// Get min_breakpoint_reads threshold for one-sided fusions
    /// Falls back to general min_breakpoint_reads if not specifically configured
    pub fn get_one_sided_min_breakpoint_reads(&self) -> i64 {
        self.thresholds.fusions.one_sided.min_breakpoint_reads
            .unwrap_or(self.thresholds.fusions.min_breakpoint_reads)
    }

    /// Get min_fraction threshold for one-sided fusions
    /// Returns None if not configured (no filtering)
    pub fn get_one_sided_min_fraction(&self) -> Option<f64> {
        self.thresholds.fusions.one_sided.min_fraction
    }

    /// Check if a fusion matches the configured spike-in control
    pub fn is_spike_in_fusion(&self, g1: &str, g2: &str) -> bool {
        if let Some(si) = &self.genes.fusions.spike_in {
            (g1 == si.gene1 && g2 == si.gene2) || (g1 == si.gene2 && g2 == si.gene1)
        } else {
            false
        }
    }

    /// Collect all partner gene names that need coordinates loaded from GFF
    pub fn get_all_partner_genes(&self) -> std::collections::HashSet<String> {
        self.genes.fusions.one_sided_partners
            .values()
            .flatten()
            .cloned()
            .collect()
    }
}

// ============================================================================
// Aggregate Configuration (Separate from Pipeline)
// ============================================================================

/// Special gene handling for fusions in aggregate
#[derive(Deserialize, Debug, Clone, Default)]
pub struct FusionSpecialGene {
    #[serde(default)]
    pub skip_self_fusion: bool,
    #[serde(default)]
    pub bypass_breakpoint_filter: bool,
}

/// Aggregate report configuration
#[derive(Deserialize, Debug, Clone)]
pub struct AggregateConfig {
    /// ITD genes to report (from ITD caller)
    #[serde(default)]
    pub itd_genes: Vec<String>,
    /// CNV genes to report focal CN for
    #[serde(default)]
    pub cnv_genes: Vec<String>,
    /// Genes to report local (region-level) CN for
    #[serde(default)]
    pub local_cn_genes: Vec<String>,
    /// Genes with deletion reporting
    #[serde(default)]
    pub deletion_genes: Vec<String>,
    /// Genes with duplication reporting (CNV-detected intragenic duplications)
    #[serde(default)]
    pub duplication_genes: Vec<String>,
    /// SNV/pharmacogenomics genes to report
    #[serde(default)]
    pub snv_genes: Vec<String>,
    /// Fusion blacklist pairs
    #[serde(default)]
    pub fusion_blacklist: Vec<(String, String)>,
    /// Special fusion gene handling
    #[serde(default)]
    pub fusion_special_genes: HashMap<String, FusionSpecialGene>,
    /// Minimum supporting reads for fusion filtering
    #[serde(default = "default_min_supporting_reads")]
    pub min_supporting_reads: i64,
    /// Minimum breakpoint reads for fusion filtering
    #[serde(default = "default_min_breakpoint_reads")]
    pub min_breakpoint_reads: i64,
    /// Column order for output (optional - if empty, auto-generates from gene lists)
    #[serde(default)]
    pub columns: Vec<String>,
}

impl AggregateConfig {
    /// Load aggregate configuration from a JSON file
    pub fn load(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let config: AggregateConfig = serde_json::from_reader(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok(config)
    }

    /// Check if a fusion pair is blacklisted
    pub fn is_blacklisted_pair(&self, gene1: &str, gene2: &str) -> bool {
        for (a, b) in &self.fusion_blacklist {
            if (gene1 == a && gene2 == b) || (gene1 == b && gene2 == a) {
                return true;
            }
        }
        false
    }

    /// Check if a gene should skip self-fusion filtering
    pub fn skip_self_fusion(&self, gene_name: &str) -> bool {
        for (name, special) in &self.fusion_special_genes {
            if gene_name.contains(name) && special.skip_self_fusion {
                return true;
            }
        }
        false
    }

    /// Check if a gene should bypass breakpoint filter
    pub fn bypass_breakpoint_filter(&self, gene_name: &str) -> bool {
        for (name, special) in &self.fusion_special_genes {
            if gene_name.contains(name) && special.bypass_breakpoint_filter {
                return true;
            }
        }
        false
    }
}

impl Default for AggregateConfig {
    fn default() -> Self {
        Self {
            itd_genes: vec![],
            cnv_genes: vec![],
            local_cn_genes: vec![],
            deletion_genes: vec![],
            duplication_genes: vec![],
            snv_genes: vec![],
            fusion_blacklist: vec![
            ],
            fusion_special_genes: {
                HashMap::new()
            },
            min_supporting_reads: default_min_supporting_reads(),
            min_breakpoint_reads: default_min_breakpoint_reads(),
            columns: vec![],
        }
    }
}
