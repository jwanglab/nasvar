use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;
use log::warn;

use crate::utils::contig::ContigMapper;

#[derive(Debug, Clone)]
pub struct BedRegion {
    pub segment: String,
    pub start: u32,
    pub end: u32,
    pub name: String,
}

/// Drop BED regions whose contig is absent from the alignment header (checked
/// via the `ContigMapper`, which understands chr<->accession translation),
/// logging a strong warning naming the skipped genes.
///
/// This guards against a BED built for a different assembly (patch/alt/unplaced
/// contigs, or a version mismatch): such regions would otherwise be silently
/// dropped by some stages and hard-abort others (a region query on a contig not
/// in the header errors). Returns the regions that ARE present in the alignment.
pub fn filter_regions_present_in_bam(
    regions: Vec<BedRegion>,
    mapper: &ContigMapper,
    label: &str,
) -> Vec<BedRegion> {
    let (kept, dropped): (Vec<BedRegion>, Vec<BedRegion>) =
        regions.into_iter().partition(|r| mapper.exists_in_bam(&r.segment));
    if !dropped.is_empty() {
        let names: Vec<String> = dropped
            .iter()
            .map(|r| format!("{} [{}]", r.name, r.segment))
            .collect();
        warn!(
            "{}: {} region(s) reference a contig ABSENT from the alignment and will be SKIPPED (not assessed) \
             — verify the BED matches the reference assembly the BAM was aligned to: {}",
            label, dropped.len(), names.join(", ")
        );
    }
    kept
}

pub fn read_bed(bed_path: &str) -> Result<Vec<BedRegion>, Box<dyn std::error::Error>> {
    let mut file = File::open(bed_path)
        .map_err(|e| std::io::Error::other(format!("Error opening BED file {}: {}", bed_path, e)))?;

    // Sniff the first two bytes for the gzip magic (0x1f 0x8b). This catches
    // both plain gzip (produced by `gzip`) and BGZF (blocked gzip produced by
    // `bgzip`), since BGZF is a concatenation of gzip members that
    // MultiGzDecoder handles transparently.
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file = File::open(bed_path)?; // re-open to reset position

    let reader: Box<dyn BufRead> = if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    read_bed_from_reader(reader)
}

pub fn read_bed_from_reader<R: BufRead>(reader: R) -> Result<Vec<BedRegion>, Box<dyn std::error::Error>> {
    let mut targets: Vec<BedRegion> = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
             return Err(format!("Malformed BED line {} (expected at least 3 columns, got {}): {}", i+1, parts.len(), line).into());
        }

        let chrom = parts[0].to_string();
        let start: u32 = parts[1].parse().map_err(|e| format!("Invalid start at line {}: {}", i+1, e))?;
        let end: u32 = parts[2].parse().map_err(|e| format!("Invalid end at line {}: {}", i+1, e))?;
        let name = if parts.len() > 3 { parts[3].to_string() } else { format!("{}:{}-{}", chrom, start, end) };

        targets.push(BedRegion {
            segment: chrom,
            start,
            end,
            name
        });
    }
    Ok(targets)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_bed_parsing() {
        let data = "chr1\t100\t200\tgene1\nchr2\t500\t600\n#Comment\n";
        let cursor = Cursor::new(data);
        let regions = read_bed_from_reader(cursor).unwrap();

        assert_eq!(regions.len(), 2);

        // 4 column line
        assert_eq!(regions[0].segment, "chr1");
        assert_eq!(regions[0].start, 100);
        assert_eq!(regions[0].end, 200);
        assert_eq!(regions[0].name, "gene1");

        // 3 column line (name auto-generated)
        assert_eq!(regions[1].segment, "chr2");
        assert_eq!(regions[1].start, 500);
        assert_eq!(regions[1].end, 600);
        assert_eq!(regions[1].name, "chr2:500-600");
    }

    #[test]
    fn test_bed_gzipped() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write;

        let data = "chr1\t100\t200\tgene1\nchr2\t500\t600\n";
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data.as_bytes()).unwrap();
        let gz_bytes = enc.finish().unwrap();

        let path = std::env::temp_dir().join("nasvar_bed_gzip_test.bed.gz");
        std::fs::write(&path, &gz_bytes).unwrap();

        let regions = read_bed(path.to_str().unwrap()).unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].segment, "chr1");
        assert_eq!(regions[0].start, 100);
        assert_eq!(regions[0].end, 200);
        assert_eq!(regions[1].segment, "chr2");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_bed_malformed() {
        // Line 2 is malformed (only 2 cols) — should be a fatal error
        let data = "chr1\t100\t200\nchr2\t500\nchr3\t1000\t2000";
        let cursor = Cursor::new(data);
        let result = read_bed_from_reader(cursor);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Malformed BED line 2"));
    }
}
