//! Pack pipeline outputs into a `<prefix>.nasvar.zip` result package.
//!
//! The package is a zip containing:
//!   - `manifest.json` — small header describing the run (see [`Manifest`]).
//!   - every sibling file matching `<basename>.*` in `out_prefix`'s directory
//!     (the unified output JSON, coverage TSV, MAF, CNV TSV, slice BAM/BAI,
//!     etc.). The previously-written zip itself is excluded so re-runs don't
//!     nest packages.
//!
//! `.bam` entries are stored uncompressed because they're already bgzf; every
//! other entry is deflated. Intended for portable sharing of a run between
//! machines (the browser frontend imports the same format directly).

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{anyhow, Result};
use log::info;
use serde::Serialize;
use zip::CompressionMethod;
use zip::write::{SimpleFileOptions, ZipWriter};

const FORMAT: &str = "nasvar-result-package";
const VERSION: u32 = 1;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Manifest {
    format: &'static str,
    version: u32,
    /// Unix milliseconds the package was written. Browser side renders this
    /// with `new Date(savedAtMs).toISOString()` if it wants an ISO string.
    saved_at_ms: u64,
    /// Basename of the input BAM (no directory), or null if unknown.
    bam_name: Option<String>,
    /// Basename of the run's `out_prefix` (no directory).
    out_prefix: String,
    /// Sum of bundled file sizes (uncompressed). Manifest itself excluded.
    total_bytes: u64,
    /// Number of bundled output files. Manifest itself excluded.
    file_count: u32,
}

/// Write `<out_prefix>.nasvar.zip` containing the manifest and every sibling
/// file matching `<basename>.*`. Returns the zip path on success.
pub fn pack_results(out_prefix: &str, bam_path: Option<&str>) -> Result<String> {
    let zip_path = format!("{}.nasvar.zip", out_prefix);

    let prefix_path = Path::new(out_prefix);
    let dir = match prefix_path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let prefix_base = prefix_path
        .file_name()
        .ok_or_else(|| anyhow!("out_prefix '{}' has no basename", out_prefix))?
        .to_string_lossy()
        .to_string();

    let dot_prefix = format!("{}.", prefix_base);
    let mut files: Vec<(String, PathBuf)> = Vec::new();
    for entry in std::fs::read_dir(&dir)
        .map_err(|e| anyhow!("read_dir {}: {}", dir.display(), e))?
    {
        let entry = entry?;
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.starts_with(&dot_prefix) {
            continue;
        }
        // Skip the zip we're about to write (and any prior one with the same name).
        if name.ends_with(".nasvar.zip") {
            continue;
        }
        if !entry.file_type()?.is_file() {
            continue;
        }
        files.push((name, entry.path()));
    }
    files.sort_by(|a, b| a.0.cmp(&b.0));

    let total_bytes: u64 = files
        .iter()
        .map(|(_, p)| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0))
        .sum();
    let saved_at_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis() as u64)
        .unwrap_or(0);

    let manifest = Manifest {
        format: FORMAT,
        version: VERSION,
        saved_at_ms,
        bam_name: bam_path.and_then(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
        }),
        out_prefix: prefix_base.clone(),
        total_bytes,
        file_count: files.len() as u32,
    };

    let zip_file = std::fs::File::create(&zip_path)
        .map_err(|e| anyhow!("create {}: {}", zip_path, e))?;
    let mut zw = ZipWriter::new(zip_file);

    zw.start_file(
        "manifest.json",
        SimpleFileOptions::default().compression_method(CompressionMethod::Deflated),
    )?;
    serde_json::to_writer_pretty(&mut zw, &manifest)
        .map_err(|e| anyhow!("write manifest: {}", e))?;

    for (name, path) in &files {
        // Already-compressed payloads (bgzf BAM, gzipped TSVs) get stored to
        // avoid redundant deflate work.
        let method = if name.ends_with(".bam") || name.ends_with(".gz") {
            CompressionMethod::Stored
        } else {
            CompressionMethod::Deflated
        };
        zw.start_file(
            name,
            SimpleFileOptions::default().compression_method(method),
        )?;
        let mut f = std::fs::File::open(path)
            .map_err(|e| anyhow!("open {}: {}", path.display(), e))?;
        std::io::copy(&mut f, &mut zw)
            .map_err(|e| anyhow!("copy {}: {}", path.display(), e))?;
    }
    zw.finish().map_err(|e| anyhow!("finish zip: {}", e))?;

    info!(
        "[pack] wrote {} ({} files, {} bytes uncompressed)",
        zip_path,
        files.len(),
        total_bytes
    );
    Ok(zip_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use zip::ZipArchive;

    #[test]
    fn pack_round_trips_outputs_and_skips_prior_zip() {
        let tmp = std::env::temp_dir().join(format!("nasvar_pack_{}", std::process::id()));
        std::fs::create_dir_all(&tmp).unwrap();
        let prefix = tmp.join("sample1");
        let prefix_str = prefix.to_str().unwrap();

        // Emit a TSV, a JSON, a fake bgzf-stored BAM, and a previous zip.
        std::fs::write(format!("{}.coverage.tsv", prefix_str), b"chrom\tstart\tend\n").unwrap();
        std::fs::write(format!("{}.json", prefix_str), b"{\"k\":1}").unwrap();
        std::fs::write(format!("{}.slice.bam", prefix_str), b"\x1f\x8bfake-bgzf").unwrap();
        std::fs::write(format!("{}.nasvar.zip", prefix_str), b"stale").unwrap();
        // Sibling that doesn't belong (different basename) must be ignored.
        std::fs::write(tmp.join("other.tsv"), b"nope").unwrap();

        let zip_path = pack_results(prefix_str, Some("/in/sample1.bam")).unwrap();
        assert_eq!(zip_path, format!("{}.nasvar.zip", prefix_str));

        let f = std::fs::File::open(&zip_path).unwrap();
        let mut zip = ZipArchive::new(f).unwrap();
        let names: Vec<String> = (0..zip.len())
            .map(|i| zip.by_index(i).unwrap().name().to_string())
            .collect();
        assert!(names.contains(&"manifest.json".to_string()));
        assert!(names.contains(&"sample1.coverage.tsv".to_string()));
        assert!(names.contains(&"sample1.json".to_string()));
        assert!(names.contains(&"sample1.slice.bam".to_string()));
        // The stale prior zip and the unrelated sibling must NOT be packed.
        assert!(!names.contains(&"sample1.nasvar.zip".to_string()));
        assert!(!names.contains(&"other.tsv".to_string()));

        let mut buf = String::new();
        zip.by_name("manifest.json").unwrap().read_to_string(&mut buf).unwrap();
        let v: serde_json::Value = serde_json::from_str(&buf).unwrap();
        assert_eq!(v["format"], "nasvar-result-package");
        assert_eq!(v["version"], 1);
        assert_eq!(v["bamName"], "sample1.bam");
        assert_eq!(v["outPrefix"], "sample1");
        assert_eq!(v["fileCount"], 3);

        std::fs::remove_dir_all(&tmp).ok();
    }
}

