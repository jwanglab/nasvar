//! Unified output module for pipeline results
//!
//! This module provides:
//! - `UnifiedOutput`: A single structure containing all pipeline results
//! - `OutputCollector`: A builder pattern for collecting results from modules
//!
//! # Example
//!
//! ```ignore
//! use nasvar::output::{OutputCollector, FusionsOutput};
//!
//! let collector = OutputCollector::new()
//!     .with_metadata(metadata)
//!     .with_qc(qc_data)
//!     .with_fusions(fusions_output);
//!
//! // Write unified output
//! collector.write_to_prefix("sample")?;
//! ```

pub mod collector;
pub mod schema;
pub mod types;

// Re-export main types for convenience
pub use collector::OutputCollector;
pub use types::{
    BreakpointConsensusOutput, CnvGeneResult, CnvOutput, FusionBreakpoint,
    FusionBreakpointConsensus, FusionEvent, FusionsOutput, GeneInfo, ItdEvent, ItdOutput,
    KaryotypeOutput, QcOutput, SnvGeneResult, SnvOutput, StructuralVariant, UnifiedOutput,
};
