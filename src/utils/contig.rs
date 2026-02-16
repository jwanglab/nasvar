//! Chromosome/contig name mapping utilities.
//!
//! This module provides `ContigMapper` for translating between different naming
//! conventions used in BAM files (e.g., chr1 vs NC_060925.1).

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use log::info;

/// T2T-CHM13v2.0 accession IDs in order: chr1-22, chrX, chrY
const T2T_ACCESSIONS: [&str; 24] = [
    "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1",
    "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1",
    "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1",
    "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"
];

/// Chromosome names in order: chr1-22, chrX, chrY
const CHR_NAMES: [&str; 24] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
];

/// The naming convention used in a BAM file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NamingConvention {
    /// Chromosome names like chr1, chr2, ..., chrX, chrY
    ChrNames,
    /// NCBI accession IDs like NC_060925.1 (T2T-CHM13v2.0)
    Accession,
    /// Mixed or unrecognized convention
    Unknown,
}

/// Maps between chromosome names and accession IDs.
///
/// Provides bidirectional translation between naming conventions
/// and auto-detects the convention used in a BAM header.
#[derive(Debug, Clone)]
pub struct ContigMapper {
    /// The naming convention detected in the BAM header
    pub detected_convention: NamingConvention,

    /// chr name -> accession (e.g., "chr1" -> "NC_060925.1")
    chr_to_acc: HashMap<String, String>,

    /// accession -> chr name (e.g., "NC_060925.1" -> "chr1")
    acc_to_chr: HashMap<String, String>,

    /// chr name -> index (0-23 for chr1-22, X, Y)
    chr_to_idx: HashMap<String, usize>,

    /// accession -> index (0-23)
    acc_to_idx: HashMap<String, usize>,

    /// Reference names actually present in BAM header
    bam_refs: HashSet<String>,
}

impl Default for ContigMapper {
    fn default() -> Self {
        Self::new()
    }
}

impl ContigMapper {
    /// Create a new ContigMapper with T2T-CHM13v2.0 mappings.
    pub fn new() -> Self {
        let mut chr_to_acc = HashMap::new();
        let mut acc_to_chr = HashMap::new();
        let mut chr_to_idx = HashMap::new();
        let mut acc_to_idx = HashMap::new();

        for (i, (chr, acc)) in CHR_NAMES.iter().zip(T2T_ACCESSIONS.iter()).enumerate() {
            chr_to_acc.insert(chr.to_string(), acc.to_string());
            acc_to_chr.insert(acc.to_string(), chr.to_string());
            chr_to_idx.insert(chr.to_string(), i);
            acc_to_idx.insert(acc.to_string(), i);
        }

        Self {
            detected_convention: NamingConvention::Unknown,
            chr_to_acc,
            acc_to_chr,
            chr_to_idx,
            acc_to_idx,
            bam_refs: HashSet::new(),
        }
    }

    /// Create a ContigMapper from BAM header reference names, auto-detecting convention.
    pub fn from_refs(refs: &[String]) -> Self {
        let mut mapper = Self::new();
        mapper.bam_refs = refs.iter().cloned().collect();
        mapper.detected_convention = mapper.detect_convention();
        mapper
    }

    /// Create a ContigMapper from a FASTA index (.fai) file, auto-detecting convention.
    ///
    /// Reads sequence names from the first column of the .fai file and uses them
    /// to detect the FASTA's naming convention.
    pub fn from_fai(fai_path: &str) -> std::io::Result<Self> {
        let file = std::fs::File::open(fai_path)?;
        let reader = BufReader::new(file);
        let mut refs = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if let Some(name) = line.split('\t').next()
                && !name.is_empty() {
                    refs.push(name.to_string());
                }
        }
        let mapper = Self::from_refs(&refs);
        info!("FASTA index: detected {} sequences, convention: {:?}",
            refs.len(), mapper.detected_convention);
        Ok(mapper)
    }

    /// Detect the naming convention based on BAM header references.
    fn detect_convention(&self) -> NamingConvention {
        let mut chr_count = 0;
        let mut acc_count = 0;

        for r in &self.bam_refs {
            if self.chr_to_acc.contains_key(r) {
                chr_count += 1;
            }
            if self.acc_to_chr.contains_key(r) {
                acc_count += 1;
            }
        }

        // If we find at least half the canonical chromosomes in one convention
        if chr_count >= 12 && chr_count > acc_count {
            NamingConvention::ChrNames
        } else if acc_count >= 12 && acc_count > chr_count {
            NamingConvention::Accession
        } else if chr_count > 0 && acc_count == 0 {
            NamingConvention::ChrNames
        } else if acc_count > 0 && chr_count == 0 {
            NamingConvention::Accession
        } else {
            NamingConvention::Unknown
        }
    }

    /// Translate any chromosome name to the BAM's convention (for queries).
    ///
    /// If the input is already in the BAM's convention or the name is unrecognized,
    /// returns the input unchanged.
    pub fn to_bam_name(&self, name: &str) -> String {
        // If it already exists in BAM refs, use as-is
        if self.bam_refs.contains(name) {
            return name.to_string();
        }

        match self.detected_convention {
            NamingConvention::ChrNames => {
                // BAM uses chr names; convert accession -> chr
                if let Some(chr) = self.acc_to_chr.get(name)
                    && self.bam_refs.contains(chr) {
                        return chr.clone();
                    }
            }
            NamingConvention::Accession => {
                // BAM uses accessions; convert chr -> accession
                if let Some(acc) = self.chr_to_acc.get(name)
                    && self.bam_refs.contains(acc) {
                        return acc.clone();
                    }
            }
            NamingConvention::Unknown => {
                // Try both directions
                if let Some(chr) = self.acc_to_chr.get(name)
                    && self.bam_refs.contains(chr) {
                        return chr.clone();
                    }
                if let Some(acc) = self.chr_to_acc.get(name)
                    && self.bam_refs.contains(acc) {
                        return acc.clone();
                    }
            }
        }

        // Return as-is if no translation found
        name.to_string()
    }

    /// Translate any name to chr convention (for output/display).
    ///
    /// Converts accession IDs to chr names. If already in chr format
    /// or unrecognized, returns the input unchanged.
    pub fn to_chr_name(&self, name: &str) -> String {
        // If it's an accession, convert to chr
        if let Some(chr) = self.acc_to_chr.get(name) {
            return chr.clone();
        }
        // Already chr or unknown
        name.to_string()
    }

    /// Normalize a name to canonical chr format.
    ///
    /// This is an alias for `to_chr_name` for clarity in contexts
    /// where normalization is the intent.
    pub fn normalize(&self, name: &str) -> String {
        self.to_chr_name(name)
    }

    /// Check if a name (in any convention) exists in the BAM header.
    pub fn exists_in_bam(&self, name: &str) -> bool {
        if self.bam_refs.contains(name) {
            return true;
        }
        // Try translation
        let translated = self.to_bam_name(name);
        self.bam_refs.contains(&translated)
    }

    /// Get the chromosome index (0-23) for a name in any convention.
    ///
    /// Returns None for unrecognized names.
    pub fn get_chr_index(&self, name: &str) -> Option<usize> {
        if let Some(&idx) = self.chr_to_idx.get(name) {
            return Some(idx);
        }
        if let Some(&idx) = self.acc_to_idx.get(name) {
            return Some(idx);
        }
        None
    }

    /// Get the chr name for a given index (0-23).
    pub fn chr_name_from_index(idx: usize) -> Option<&'static str> {
        CHR_NAMES.get(idx).copied()
    }

    /// Get the accession for a given index (0-23).
    pub fn accession_from_index(idx: usize) -> Option<&'static str> {
        T2T_ACCESSIONS.get(idx).copied()
    }

    /// Get all T2T accession IDs.
    pub fn t2t_accessions() -> &'static [&'static str; 24] {
        &T2T_ACCESSIONS
    }

    /// Get all chromosome names.
    pub fn chr_names() -> &'static [&'static str; 24] {
        &CHR_NAMES
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detection_chr_names() {
        let refs: Vec<String> = CHR_NAMES.iter().map(|s| s.to_string()).collect();
        let mapper = ContigMapper::from_refs(&refs);
        assert_eq!(mapper.detected_convention, NamingConvention::ChrNames);
    }

    #[test]
    fn test_detection_accession() {
        let refs: Vec<String> = T2T_ACCESSIONS.iter().map(|s| s.to_string()).collect();
        let mapper = ContigMapper::from_refs(&refs);
        assert_eq!(mapper.detected_convention, NamingConvention::Accession);
    }

    #[test]
    fn test_to_bam_name_chr_convention() {
        let refs: Vec<String> = CHR_NAMES.iter().map(|s| s.to_string()).collect();
        let mapper = ContigMapper::from_refs(&refs);

        // Accession -> chr (for BAM query)
        assert_eq!(mapper.to_bam_name("NC_060925.1"), "chr1");
        assert_eq!(mapper.to_bam_name("NC_060947.1"), "chrX");

        // chr -> chr (already correct)
        assert_eq!(mapper.to_bam_name("chr1"), "chr1");
    }

    #[test]
    fn test_to_bam_name_accession_convention() {
        let refs: Vec<String> = T2T_ACCESSIONS.iter().map(|s| s.to_string()).collect();
        let mapper = ContigMapper::from_refs(&refs);

        // chr -> accession (for BAM query)
        assert_eq!(mapper.to_bam_name("chr1"), "NC_060925.1");
        assert_eq!(mapper.to_bam_name("chrX"), "NC_060947.1");

        // accession -> accession (already correct)
        assert_eq!(mapper.to_bam_name("NC_060925.1"), "NC_060925.1");
    }

    #[test]
    fn test_to_chr_name() {
        let mapper = ContigMapper::new();

        assert_eq!(mapper.to_chr_name("NC_060925.1"), "chr1");
        assert_eq!(mapper.to_chr_name("NC_060948.1"), "chrY");
        assert_eq!(mapper.to_chr_name("chr1"), "chr1");
        assert_eq!(mapper.to_chr_name("unknown"), "unknown");
    }

    #[test]
    fn test_get_chr_index() {
        let mapper = ContigMapper::new();

        assert_eq!(mapper.get_chr_index("chr1"), Some(0));
        assert_eq!(mapper.get_chr_index("chr22"), Some(21));
        assert_eq!(mapper.get_chr_index("chrX"), Some(22));
        assert_eq!(mapper.get_chr_index("chrY"), Some(23));
        assert_eq!(mapper.get_chr_index("NC_060925.1"), Some(0));
        assert_eq!(mapper.get_chr_index("NC_060948.1"), Some(23));
        assert_eq!(mapper.get_chr_index("unknown"), None);
    }

    #[test]
    fn test_exists_in_bam() {
        let refs: Vec<String> = T2T_ACCESSIONS.iter().map(|s| s.to_string()).collect();
        let mapper = ContigMapper::from_refs(&refs);

        // Direct match
        assert!(mapper.exists_in_bam("NC_060925.1"));
        // Via translation
        assert!(mapper.exists_in_bam("chr1"));
        // Non-existent
        assert!(!mapper.exists_in_bam("chrM"));
    }
}
