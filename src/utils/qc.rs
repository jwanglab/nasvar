//! QC data structures for pipeline metrics

use serde::{Deserialize, Serialize};

/// QC metrics accumulated during pipeline execution
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PipelineQcData {
    pub nt_on_target: f64,
    pub reads_on_target: f64,
    pub target_regions_nt: f64,
}
