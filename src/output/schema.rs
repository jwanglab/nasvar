//! JSON Schema generation and runtime validation for unified output

use std::sync::LazyLock;

use schemars::schema_for;
use serde_json::Value;

use super::types::UnifiedOutput;

/// Cached JSON Schema for UnifiedOutput.
static SCHEMA: LazyLock<schemars::Schema> = LazyLock::new(|| schema_for!(UnifiedOutput));

/// Returns the JSON Schema as a pretty-printed JSON string.
pub fn schema_json_pretty() -> String {
    serde_json::to_string_pretty(&*SCHEMA).expect("schema serialization should not fail")
}

/// Validate a JSON value against the UnifiedOutput schema.
///
/// Returns `Ok(())` if valid, or `Err` with a description of all validation errors.
pub fn validate(value: &Value) -> Result<(), String> {
    let schema_val = serde_json::to_value(&*SCHEMA).expect("schema serialization should not fail");
    let validator = jsonschema::validator_for(&schema_val)
        .map_err(|e| format!("Failed to compile schema: {}", e))?;

    let errors: Vec<String> = validator
        .iter_errors(value)
        .map(|e| format!("  - {}: {}", e.instance_path, e))
        .collect();

    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "Output JSON failed schema validation ({} errors):\n{}",
            errors.len(),
            errors.join("\n")
        ))
    }
}

/// Returns true if runtime validation should be performed.
///
/// Always true in debug builds. In release builds, true only if `NASVAR_VALIDATE_OUTPUT=1`.
pub fn should_validate() -> bool {
    if cfg!(debug_assertions) {
        true
    } else {
        std::env::var("NASVAR_VALIDATE_OUTPUT")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::OutputCollector;

    #[test]
    fn test_schema_generation() {
        let schema = schema_json_pretty();
        let parsed: Value = serde_json::from_str(&schema).unwrap();
        assert_eq!(
            parsed.get("type").and_then(|v| v.as_str()),
            Some("object")
        );
    }

    #[test]
    fn test_validate_minimal_output() {
        let collector = OutputCollector::new();
        let output = collector.build();
        let value = serde_json::to_value(&output).unwrap();
        assert!(validate(&value).is_ok());
    }

    #[test]
    fn test_validate_invalid_output() {
        let bad_json: Value = serde_json::json!({"bogus": true});
        let result = validate(&bad_json);
        assert!(result.is_err());
    }
}
