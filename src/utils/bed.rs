use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct BedRegion {
    pub segment: String,
    pub start: u32,
    pub end: u32,
    pub name: String,
}

pub fn read_bed(bed_path: &str) -> Result<Vec<BedRegion>, Box<dyn std::error::Error>> {
    let file = File::open(bed_path)
        .map_err(|e| std::io::Error::other(format!("Error opening BED file {}: {}", bed_path, e)))?;
    let reader = BufReader::new(file);
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
    fn test_bed_malformed() {
        // Line 2 is malformed (only 2 cols) — should be a fatal error
        let data = "chr1\t100\t200\nchr2\t500\nchr3\t1000\t2000";
        let cursor = Cursor::new(data);
        let result = read_bed_from_reader(cursor);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Malformed BED line 2"));
    }
}
