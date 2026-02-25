//! Output collector for aggregating pipeline results
//!
//! The `OutputCollector` provides a builder pattern for collecting
//! results from all pipeline modules into a unified output structure.

use std::collections::HashMap;
use std::fs::File;

use crate::utils::metadata::SequencingMetaData;
use crate::utils::qc::PipelineQcData;

use super::types::{
    BreakpointConsensusOutput, CnvOutput, FusionsOutput, ItdOutput, KaryotypeOutput, QcOutput,
    SnvOutput, UnifiedOutput,
};

/// Builder for collecting pipeline outputs into a unified structure
pub struct OutputCollector {
    output: UnifiedOutput,
}

impl OutputCollector {
    /// Create a new output collector with version and timestamp
    pub fn new() -> Self {
        Self {
            output: UnifiedOutput {
                version: env!("CARGO_PKG_VERSION").to_string(),
                timestamp: crate::utils::time::utc_now_iso8601(),
                ..Default::default()
            },
        }
    }

    /// Set sequencing metadata from BAM headers
    pub fn with_metadata(mut self, meta: SequencingMetaData) -> Self {
        self.output.metadata = Some(meta);
        self
    }

    /// Set QC data from pipeline metrics
    pub fn with_qc(mut self, qc: PipelineQcData) -> Self {
        self.output.qc = Some(QcOutput::from(qc));
        self
    }

    /// Set QC data from a fully-formed QcOutput
    pub fn with_qc_output(mut self, qc: QcOutput) -> Self {
        self.output.qc = Some(qc);
        self
    }

    /// Add reads_aligned count to QC data
    ///
    /// If QC data doesn't exist yet, creates a minimal QcOutput with just reads_aligned
    pub fn with_reads_aligned(mut self, count: u64) -> Self {
        if let Some(ref mut qc) = self.output.qc {
            qc.reads_aligned = Some(count);
        } else {
            self.output.qc = Some(QcOutput {
                reads_aligned: Some(count),
                ..Default::default()
            });
        }
        self
    }

    /// Add mean coverage to QC data
    pub fn with_mean_coverage(mut self, coverage: f64) -> Self {
        if let Some(ref mut qc) = self.output.qc {
            qc.mean_coverage = Some(coverage);
        } else {
            self.output.qc = Some(QcOutput {
                mean_coverage: Some(coverage),
                ..Default::default()
            });
        }
        self
    }

    /// Set per-target average coverage on QC data
    pub fn with_target_coverage(mut self, coverage: HashMap<String, f64>) -> Self {
        if let Some(ref mut qc) = self.output.qc {
            qc.target_coverage = Some(coverage);
        } else {
            self.output.qc = Some(QcOutput {
                target_coverage: Some(coverage),
                ..Default::default()
            });
        }
        self
    }

    /// Set fusion calling results
    pub fn with_fusions(mut self, fusions: FusionsOutput) -> Self {
        self.output.fusions = Some(fusions);
        self
    }

    /// Set karyotype analysis results
    pub fn with_karyotype(mut self, karyotype: KaryotypeOutput) -> Self {
        self.output.karyotype = Some(karyotype);
        self
    }

    /// Set CNV calling results
    pub fn with_cnv(mut self, cnv: CnvOutput) -> Self {
        self.output.cnv = Some(cnv);
        self
    }

    /// Set SNV calling results
    pub fn with_snv(mut self, snv: SnvOutput) -> Self {
        self.output.snv = Some(snv);
        self
    }

    /// Set ITD calling results
    pub fn with_itd(mut self, itd: ItdOutput) -> Self {
        self.output.itd = Some(itd);
        self
    }

    /// Set fusion breakpoint consensus results
    pub fn with_breakpoint_consensus(mut self, bp: BreakpointConsensusOutput) -> Self {
        self.output.breakpoint_consensus = Some(bp);
        self
    }

    /// Build and return the final unified output
    pub fn build(self) -> UnifiedOutput {
        self.output
    }

    /// Get a reference to the current output (for inspection)
    pub fn output(&self) -> &UnifiedOutput {
        &self.output
    }

    /// Write unified JSON to the specified path
    ///
    /// The path should be the full filename (e.g., "sample.result.json")
    pub fn write_json(&self, path: &str) -> std::io::Result<()> {
        if super::schema::should_validate() {
            let value = serde_json::to_value(&self.output)
                .map_err(std::io::Error::other)?;
            if let Err(msg) = super::schema::validate(&value) {
                log::warn!("Schema validation failed for {}: {}", path, msg);
                if cfg!(debug_assertions) {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, msg));
                }
            }
        }
        let file = File::create(path)?;
        serde_json::to_writer_pretty(file, &self.output)
            .map_err(std::io::Error::other)
    }

    /// Write unified JSON using the output prefix
    ///
    /// Creates "{prefix}.result.json"
    pub fn write_to_prefix(&self, prefix: &str) -> std::io::Result<()> {
        let path = format!("{}.result.json", prefix);
        self.write_json(&path)
    }
}

impl Default for OutputCollector {
    fn default() -> Self {
        Self::new()
    }
}

impl UnifiedOutput {
    /// Write this output to a JSON file
    pub fn write_json(&self, path: &str) -> std::io::Result<()> {
        if super::schema::should_validate() {
            let value = serde_json::to_value(self)
                .map_err(std::io::Error::other)?;
            if let Err(msg) = super::schema::validate(&value) {
                log::warn!("Schema validation failed for {}: {}", path, msg);
                if cfg!(debug_assertions) {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, msg));
                }
            }
        }
        let file = File::create(path)?;
        serde_json::to_writer_pretty(file, self)
            .map_err(std::io::Error::other)
    }

    /// Load unified output from a JSON file
    pub fn load_json(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let output: Self = serde_json::from_reader(file)?;
        Ok(output)
    }

    /// Check if any results are present
    pub fn has_results(&self) -> bool {
        self.fusions.is_some()
            || self.karyotype.is_some()
            || self.cnv.is_some()
            || self.snv.is_some()
            || self.itd.is_some()
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collector_builder() {
        let collector = OutputCollector::new()
            .with_reads_aligned(1000)
            .with_mean_coverage(45.5);

        let output = collector.build();

        assert!(output.qc.is_some());
        let qc = output.qc.unwrap();
        assert_eq!(qc.reads_aligned, Some(1000));
        assert_eq!(qc.mean_coverage, Some(45.5));
    }
}
