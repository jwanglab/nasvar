//! BAM/SAM header metadata extraction utilities

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Sequencing run metadata extracted from BAM headers
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonSchema)]
pub struct SequencingMetaData {
    #[serde(default)]
    pub run_start_time: Option<String>,
    #[serde(default)]
    pub run_id: Option<String>,
    #[serde(default)]
    pub basecall_model: Option<String>,
    #[serde(default)]
    pub library_id: Option<String>,
    #[serde(default)]
    pub sequencer_id: Option<String>,
    #[serde(default)]
    pub flow_cell_id: Option<String>,
}

/// Extract sequencing metadata from BAM/SAM header text
///
/// Parses @RG (read group) and @CO (comment) lines to extract:
/// - run_start_time: from DT field
/// - run_id: from DS field (runid=...) or ID field
/// - basecall_model: from DS field (basecall_model=...) or @CO lines
/// - library_id: from LB field
/// - sequencer_id: from PM field (if not "NONE")
/// - flow_cell_id: from PU field
pub fn extract_from_header(header_text: &str) -> Option<SequencingMetaData> {
    let mut run_start_time = None;
    let mut run_id = None;
    let mut basecall_model = None;
    let mut library_id = None;
    let mut sequencer_id = None;
    let mut flow_cell_id = None;

    for line in header_text.lines() {
        if line.starts_with("@RG") {
            for field in line.split('\t') {
                if let Some(val) = field.strip_prefix("DS:") {
                    // Parse space-separated key=value pairs from DS field
                    // e.g. "runid=abc123 basecall_model=dna_r10.4.1_e8.2_400bps_fast@v4.3.0"
                    for kv in val.split(' ') {
                        if let Some(idx) = kv.find('=') {
                            let key = &kv[..idx];
                            let value = &kv[idx + 1..];
                            match key {
                                "runid" => run_id = Some(value.to_string()),
                                "basecall_model" => basecall_model = Some(value.to_string()),
                                _ => {}
                            }
                        }
                    }
                } else if let Some(val) = field.strip_prefix("ID:") {
                    if run_id.is_none() {
                        run_id = Some(val.to_string());
                    }
                } else if let Some(val) = field.strip_prefix("LB:") {
                    library_id = Some(val.to_string());
                } else if let Some(val) = field.strip_prefix("DT:") {
                    run_start_time = Some(val.to_string());
                } else if let Some(val) = field.strip_prefix("PM:") {
                    // Skip "NONE" placeholder
                    if val != "NONE" {
                        sequencer_id = Some(val.to_string());
                    }
                } else if let Some(val) = field.strip_prefix("PU:") {
                    flow_cell_id = Some(val.to_string());
                }
            }
        } else if line.starts_with("@CO") {
            let content = line.trim_start_matches("@CO").trim();
            if let Some(val) = content.strip_prefix("basecalling_model=") {
                basecall_model = Some(val.to_string());
            } else if let Some(val) = content.strip_prefix("basecaller:model=") {
                basecall_model = Some(val.to_string());
            } else if let Some(val) = content.strip_prefix("run_id=") {
                run_id = Some(val.to_string());
            } else if let Some(val) = content.strip_prefix("flow_cell_id=") {
                flow_cell_id = Some(val.to_string());
            } else if let Some(val) = content.strip_prefix("sample_id=")
                && library_id.is_none() {
                    library_id = Some(val.to_string());
                }
        }
    }

    if run_start_time.is_some()
        || run_id.is_some()
        || basecall_model.is_some()
        || library_id.is_some()
        || sequencer_id.is_some()
        || flow_cell_id.is_some()
    {
        Some(SequencingMetaData {
            run_start_time,
            run_id,
            basecall_model,
            library_id,
            sequencer_id,
            flow_cell_id,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_from_rg_line() {
        let header = "@RG\tID:run123\tLB:sample1\tDT:2025-05-13T16:00:00\tPU:flowcell1\tPM:GridION\tDS:basecall_model=sup@v5.0 runid=abc123";
        let meta = extract_from_header(header).unwrap();

        assert_eq!(meta.run_id, Some("abc123".to_string()));
        assert_eq!(meta.basecall_model, Some("sup@v5.0".to_string()));
        assert_eq!(meta.library_id, Some("sample1".to_string()));
        assert_eq!(meta.run_start_time, Some("2025-05-13T16:00:00".to_string()));
        assert_eq!(meta.flow_cell_id, Some("flowcell1".to_string()));
        assert_eq!(meta.sequencer_id, Some("GridION".to_string()));
    }

    #[test]
    fn test_skip_none_sequencer() {
        let header = "@RG\tID:run123\tPM:NONE";
        let meta = extract_from_header(header).unwrap();

        assert_eq!(meta.sequencer_id, None);
    }
}
