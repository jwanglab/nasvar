use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    None,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => Strand::None,
        }
    }
}

/// Parsed gene annotation from a GFF file.
pub struct GeneAnnotation {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub mrna_id: String,
    pub mrna_id_bytes: Option<Vec<u8>>,
    pub chrom_bytes: Option<Vec<u8>>,
    pub strand: Option<Strand>,
    pub cds_list: Vec<(u32, u32)>,
}

// Minimal GFF parser for gene annotation
pub fn get_gene_annotation(gene_name: &str, gff_path: &str, transcript_id: Option<&str>) -> Result<GeneAnnotation, Box<dyn std::error::Error>> {
    let file = File::open(gff_path).map_err(|e| format!("Error opening GFF {}: {}", gff_path, e))?;
    let reader = BufReader::new(file);
    get_gene_annotation_from_reader(gene_name, reader, transcript_id)
}

pub fn get_gene_annotation_from_reader<R: BufRead>(gene_name: &str, reader: R, target_transcript_id: Option<&str>) -> Result<GeneAnnotation, Box<dyn std::error::Error>> {
    let mut gene_id = String::new();

    let mut chrom: Option<String> = None;
    let mut strand: Option<Strand> = None;

    // CDS list
    let mut cds_list: Vec<(u32, u32)> = Vec::new();

    let mut relevant_lines = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        // We cannot filter by gene_name here because children (mRNA, CDS) might not have it.
        // We must buffer all lines or perform multi-pass.
        // For now, buffer all (memory intensive but correct).
        relevant_lines.push(line);
    }

    // Parse gene line
    for line in &relevant_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 { continue; }

        let feat_type = parts[2];
        let info = parts[8];

        if feat_type == "gene" {
            // Check Name=...
            if info.contains(&format!("Name={}", gene_name)) {
                // Found it
                chrom = Some(parts[0].to_string());
                strand = Some(Strand::from_char(parts[6].chars().next().unwrap_or('.')));

                // Get ID
                for tag in info.split(';') {
                    let tag = tag.trim();
                    if let Some(stripped) = tag.strip_prefix("ID=") {
                         gene_id = stripped.to_string();
                    }
                }
            }
        }
    }

    if gene_id.is_empty() {
        return Err(format!("Gene {} not found in GFF", gene_name).into());
    }

    // Find ALL mRNA children of gene_id
    let mut candidate_transcripts: Vec<String> = Vec::new();

    for line in &relevant_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 { continue; }
        let feat_type = parts[2];
        let info = parts[8];

        if (feat_type == "mRNA" || feat_type == "transcript")
             && info.contains(&format!("Parent={}", gene_id)) {
                 // Get ID
                 for tag in info.split(';') {
                    let tag = tag.trim();
                    if let Some(stripped) = tag.strip_prefix("ID=") {
                         candidate_transcripts.push(stripped.to_string());
                         break;
                    }
                 }
             }
    }

    if candidate_transcripts.is_empty() {
         return Err(format!("No mRNA found for gene ID {}", gene_id).into());
    }

    // Collect CDS for all candidates
    let mut transcript_cds: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for t in &candidate_transcripts {
         transcript_cds.insert(t.clone(), Vec::new());
    }

    for line in &relevant_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 { continue; }
        let feat_type = parts[2];
        let info = parts[8];

        if feat_type == "CDS" {
             // Check parent
             let mut parent = "";
             for tag in info.split(';') {
                 let tag = tag.trim();
                 if let Some(stripped) = tag.strip_prefix("Parent=") {
                     parent = stripped;
                     break;
                 }
             }

             // Handle comma-separated parents in some GFFs? Assuming single parent for now or primary.
             // If parent in candidates
             if let Some(list) = transcript_cds.get_mut(parent) {
                  let start: u32 = parts[3].parse().unwrap_or(0);
                  let end: u32 = parts[4].parse().unwrap_or(0);
                  if start > 0 && end > 0 {
                      list.push((start, end));
                  }
             }
        }
    }

    // Select best transcript (longest CDS) OR specific target
    let mut best_mrna = String::new();
    let mut best_len: u32 = 0;

    // If target transcript specified, verify it exists
    if let Some(target) = target_transcript_id {
         if transcript_cds.contains_key(target) {
              best_mrna = target.to_string();
              cds_list = transcript_cds.get(target).unwrap().clone();
         } else {
             // Fallback or Error?
             // Let's error if specific ID requested but not found (under this gene).
             // However, maybe it wasn't a child of this gene ID?
             // The loop above only finding children of 'gene_id'.
             // If user provides a transcript ID, it MUST be a child of the gene name provided.
             return Err(format!("Transcript {} not found for gene {}", target, gene_name).into());
         }
    } else {
        // Default to first if none have CDS
        if !candidate_transcripts.is_empty() {
             best_mrna = candidate_transcripts[0].clone();
             if let Some(bg) = transcript_cds.get(&best_mrna) {
                 cds_list = bg.clone();
                 best_len = cds_list.iter().map(|(s,e)| e - s + 1).sum();
             }
        }

        for (tid, cds) in &transcript_cds {
             let len: u32 = cds.iter().map(|(s,e)| e - s + 1).sum();
             if len > best_len {
                  best_len = len;
                  best_mrna = tid.clone();
                  cds_list = cds.clone();
             }
        }
    }

    let mrna_id = best_mrna;

    Ok(GeneAnnotation {
        gene_id,
        gene_name: Some(gene_name.to_string()),
        mrna_id: mrna_id.clone(),
        mrna_id_bytes: Some(mrna_id.as_bytes().to_vec()),
        chrom_bytes: chrom.map(|s| s.as_bytes().to_vec()),
        strand,
        cds_list,
    })
}

// ============================================================================
// Partner Gene Index for One-Sided Fusion Filtering
// ============================================================================

/// Gene boundary information (chromosome, start, end)
#[derive(Debug, Clone)]
pub struct GeneBounds {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

/// Index of partner gene coordinates for one-sided fusion filtering
///
/// This struct caches gene coordinates loaded from a GFF file,
/// allowing efficient checking of whether a fusion breakpoint
/// lands within margin of an allowed partner gene.
#[derive(Debug, Clone, Default)]
pub struct PartnerGeneIndex {
    /// gene_name -> (chromosome, start, end)
    genes: HashMap<String, GeneBounds>,
}

impl PartnerGeneIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self { genes: HashMap::new() }
    }

    /// Load partner gene coordinates from a GFF file
    ///
    /// Only extracts coordinates for genes in the provided set.
    /// Returns an index with all found genes; genes not found in GFF are silently skipped.
    pub fn load_from_gff(gff_path: &str, partner_genes: &HashSet<String>) -> Result<Self, Box<dyn std::error::Error>> {
        if partner_genes.is_empty() {
            return Ok(Self::new());
        }

        let file = File::open(gff_path)?;
        let reader = BufReader::new(file);
        Self::load_from_reader(reader, partner_genes)
    }

    /// Load from a reader (for testing)
    pub fn load_from_reader<R: BufRead>(reader: R, partner_genes: &HashSet<String>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut genes = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') { continue; }

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 9 { continue; }

            let feat_type = parts[2];
            if feat_type != "gene" { continue; }

            let info = parts[8];

            // Extract gene name from Name= attribute
            let mut gene_name: Option<&str> = None;
            for tag in info.split(';') {
                let tag = tag.trim();
                if let Some(stripped) = tag.strip_prefix("Name=") {
                    gene_name = Some(stripped);
                    break;
                }
            }

            if let Some(name) = gene_name
                && partner_genes.contains(name) {
                    let chrom = parts[0].to_string();
                    let start: u32 = parts[3].parse().unwrap_or(0);
                    let end: u32 = parts[4].parse().unwrap_or(0);

                    if start > 0 && end > 0 {
                        genes.insert(name.to_string(), GeneBounds { chrom, start, end });
                    }
                }
        }

        Ok(Self { genes })
    }

    /// Check if a breakpoint falls within margin of a specific gene.
    /// Normalizes chromosome names to handle NC_* vs chr* mismatches.
    pub fn is_near_gene(&self, gene_name: &str, chrom: &str, pos: u32, margin: u32) -> bool {
        if let Some(bounds) = self.genes.get(gene_name) {
            let mapper = crate::bam::ContigMapper::new();
            let norm_chrom = mapper.to_chr_name(chrom);
            let norm_bounds_chrom = mapper.to_chr_name(&bounds.chrom);
            if norm_chrom == norm_bounds_chrom {
                let g_st = bounds.start.saturating_sub(margin);
                let g_en = bounds.end + margin;
                return pos >= g_st && pos <= g_en;
            }
        }
        false
    }

    /// Check if a breakpoint falls within margin of ANY of the provided partner genes
    pub fn is_near_any_partner(&self, partners: &[String], chrom: &str, pos: u32, margin: u32) -> bool {
        partners.iter().any(|p| self.is_near_gene(p, chrom, pos, margin))
    }

    /// Get gene coordinates if present
    pub fn get_gene(&self, gene_name: &str) -> Option<&GeneBounds> {
        self.genes.get(gene_name)
    }

    /// Check if any genes were loaded
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    /// Get number of genes in the index
    pub fn len(&self) -> usize {
        self.genes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_gff_parsing_simple() {
        let gff_data = "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=TargetGene\n\
                        chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna1;Parent=gene1\n\
                        chr1\t.\tCDS\t1100\t1200\t.\t+\t0\tID=cds1;Parent=mrna1\n\
                        chr1\t.\tCDS\t1300\t1400\t.\t+\t1\tID=cds2;Parent=mrna1\n";

        let cursor = Cursor::new(gff_data);
        let result = get_gene_annotation_from_reader("TargetGene", cursor, None).unwrap();

        assert_eq!(result.gene_id, "gene1");
        assert_eq!(result.gene_name.unwrap(), "TargetGene");
        assert_eq!(result.mrna_id, "mrna1");
        assert_eq!(result.chrom_bytes.unwrap(), b"chr1");
        assert_eq!(result.strand.unwrap(), Strand::Forward);
        assert_eq!(result.cds_list.len(), 2);
        assert_eq!(result.cds_list[0], (1100, 1200));
        assert_eq!(result.cds_list[1], (1300, 1400));
    }

    #[test]
    fn test_gff_missing_gene() {
        let gff_data = "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=OtherGene\n";
        let cursor = Cursor::new(gff_data);
        let result = get_gene_annotation_from_reader("TargetGene", cursor, None);
        assert!(result.is_err());
    }
    #[test]
    fn test_gff_multiple_isoforms() {
        // Gene with 2 transcripts.
        // mrna1: Non-coding (no CDS) or weird. Appears first.
        // mrna2: Coding (has CDS).
        // Current logic picks mrna1 and returns 0 CDS.
        // Desired logic picks mrna2.
        let gff_data = "chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=MultiIsoformGene\n\
                        chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mrna1;Parent=gene1\n\
                        chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mrna2;Parent=gene1\n\
                        chr1\t.\tCDS\t1100\t1200\t.\t+\t0\tID=cds1;Parent=mrna2\n\
                        chr1\t.\tCDS\t1300\t1400\t.\t+\t1\tID=cds2;Parent=mrna2\n";

        // Auto select best
        let cursor = Cursor::new(gff_data);
        let result = get_gene_annotation_from_reader("MultiIsoformGene", cursor, None).unwrap();
        println!("Picked mRNA: {}", result.mrna_id);
        assert!(!result.cds_list.is_empty(), "Should pick the isoform with CDS");
        assert_eq!(result.mrna_id, "mrna2");

         // Explicit select mrna1 (bad one)
        let cursor2 = Cursor::new(gff_data);
        let result2 = get_gene_annotation_from_reader("MultiIsoformGene", cursor2, Some("mrna1")).unwrap();
        assert_eq!(result2.mrna_id, "mrna1");
        assert!(result2.cds_list.is_empty());

    }

    #[test]
    fn test_partner_gene_index() {
        let gff_data = "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=GATA2\n\
                        chr1\t.\tgene\t5000\t6000\t.\t+\t.\tID=gene2;Name=ETV6\n\
                        chr2\t.\tgene\t10000\t20000\t.\t-\t.\tID=gene3;Name=RUNX1\n";

        let partner_genes: HashSet<String> = ["GATA2", "ETV6", "RUNX1", "UNKNOWN"]
            .iter().map(|s| s.to_string()).collect();

        let cursor = Cursor::new(gff_data);
        let index = PartnerGeneIndex::load_from_reader(cursor, &partner_genes).unwrap();

        // Should have loaded 3 genes (UNKNOWN not in GFF)
        assert_eq!(index.len(), 3);

        // Test is_near_gene
        assert!(index.is_near_gene("GATA2", "chr1", 1500, 100));  // Inside gene
        assert!(index.is_near_gene("GATA2", "chr1", 900, 200));   // Within margin before
        assert!(index.is_near_gene("GATA2", "chr1", 2100, 200));  // Within margin after
        assert!(!index.is_near_gene("GATA2", "chr1", 3000, 100)); // Outside margin
        assert!(!index.is_near_gene("GATA2", "chr2", 1500, 100)); // Wrong chromosome

        // Test is_near_any_partner
        let partners = vec!["GATA2".to_string(), "ETV6".to_string()];
        assert!(index.is_near_any_partner(&partners, "chr1", 1500, 100));  // Near GATA2
        assert!(index.is_near_any_partner(&partners, "chr1", 5500, 100));  // Near ETV6
        assert!(!index.is_near_any_partner(&partners, "chr2", 15000, 100)); // Near RUNX1 but not in partners list

        // Test get_gene
        let gata2 = index.get_gene("GATA2").unwrap();
        assert_eq!(gata2.chrom, "chr1");
        assert_eq!(gata2.start, 1000);
        assert_eq!(gata2.end, 2000);
    }
}
