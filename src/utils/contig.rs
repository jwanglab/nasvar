//! Chromosome/contig name mapping utilities.
//!
//! This module provides `ContigMapper` for translating between different naming
//! conventions used in BAM/FASTA files (e.g., chr1 vs NC_000001.11).
//!
//! The mapping is **not** hardcoded to any assembly: it is built from the ordered
//! `contigs` list in the reference config (`ReferenceConfig::contigs`). The position
//! of each contig in that list defines the canonical chromosome index (0-based) used
//! throughout the pipeline (karyotype, MAF, coverage).

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use log::info;

use crate::config::Contig;

/// The naming convention used in a BAM/FASTA file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NamingConvention {
    /// Chromosome names like chr1, chr2, ..., chrX, chrY
    ChrNames,
    /// NCBI accession IDs like NC_000001.11
    Accession,
    /// Mixed or unrecognized convention
    Unknown,
}

/// Maps between chromosome names and accession IDs for a specific assembly.
///
/// Built from the reference config's ordered contig list. Provides bidirectional
/// translation between naming conventions and auto-detects the convention used in a
/// BAM header or FASTA index.
#[derive(Debug, Clone)]
pub struct ContigMapper {
    /// The naming convention detected in the BAM/FASTA references
    pub detected_convention: NamingConvention,

    /// chr name -> accession (e.g., "chr1" -> "NC_000001.11")
    chr_to_acc: HashMap<String, String>,

    /// accession -> chr name (e.g., "NC_000001.11" -> "chr1")
    acc_to_chr: HashMap<String, String>,

    /// chr name -> canonical index (position in the config contig list)
    chr_to_idx: HashMap<String, usize>,

    /// accession -> canonical index
    acc_to_idx: HashMap<String, usize>,

    /// short bare name (e.g., "1", "X", "MT") -> canonical index. Populated
    /// from each contig's chr name by stripping the `chr` prefix. `MT` is
    /// registered as an alias for `chrM`. Lets the MAF writer emit
    /// Ensembl-style short names to save disk space, and lets any reader
    /// interpret them via the same `to_chr_name` / `get_chr_index` funnel.
    short_to_idx: HashMap<String, usize>,

    /// canonical index -> chr name
    idx_to_chr: Vec<String>,

    /// canonical index -> accession
    idx_to_acc: Vec<String>,

    /// canonical index -> short bare name (e.g., "1", "X").
    idx_to_short: Vec<String>,

    /// Reference names actually present in the BAM/FASTA header
    bam_refs: HashSet<String>,
}

impl Default for ContigMapper {
    /// An empty mapper: all translations are identity, `num_chromosomes()` is 0.
    /// Useful for tools that only need pass-through name handling.
    fn default() -> Self {
        Self::from_contigs(&[])
    }
}

impl ContigMapper {
    /// Build a mapper from an ordered canonical contig list (no BAM refs yet).
    ///
    /// The order of `contigs` defines the canonical chromosome index.
    pub fn from_contigs(contigs: &[Contig]) -> Self {
        let mut chr_to_acc = HashMap::new();
        let mut acc_to_chr = HashMap::new();
        let mut chr_to_idx = HashMap::new();
        let mut acc_to_idx = HashMap::new();
        let mut short_to_idx = HashMap::new();
        let mut idx_to_chr = Vec::with_capacity(contigs.len());
        let mut idx_to_acc = Vec::with_capacity(contigs.len());
        let mut idx_to_short = Vec::with_capacity(contigs.len());

        for (i, c) in contigs.iter().enumerate() {
            chr_to_acc.insert(c.name.clone(), c.accession.clone());
            acc_to_chr.insert(c.accession.clone(), c.name.clone());
            chr_to_idx.insert(c.name.clone(), i);
            acc_to_idx.insert(c.accession.clone(), i);
            idx_to_chr.push(c.name.clone());
            idx_to_acc.push(c.accession.clone());

            // Short bare name: strip `chr` prefix. For mitochondrion we also
            // register `MT` as an alias since it's the Ensembl convention.
            let short = c.name.strip_prefix("chr").unwrap_or(&c.name).to_string();
            short_to_idx.insert(short.clone(), i);
            if short == "M" { short_to_idx.insert("MT".to_string(), i); }
            idx_to_short.push(short);
        }

        Self {
            detected_convention: NamingConvention::Unknown,
            chr_to_acc,
            acc_to_chr,
            chr_to_idx,
            acc_to_idx,
            short_to_idx,
            idx_to_chr,
            idx_to_acc,
            idx_to_short,
            bam_refs: HashSet::new(),
        }
    }

    /// Build from a contig list and detect the convention from BAM/FASTA reference names.
    pub fn from_contigs_and_refs(contigs: &[Contig], refs: &[String]) -> Self {
        let mut mapper = Self::from_contigs(contigs);
        mapper.bam_refs = refs.iter().cloned().collect();
        mapper.detected_convention = mapper.detect_convention();
        mapper
    }

    /// Build from a contig list, detecting the convention from a FASTA index (.fai) file.
    ///
    /// Reads sequence names from the first column of the .fai file.
    pub fn from_fai(fai_path: &str, contigs: &[Contig]) -> std::io::Result<Self> {
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
        let mapper = Self::from_contigs_and_refs(contigs, &refs);
        info!("FASTA index: detected {} sequences, convention: {:?}",
            refs.len(), mapper.detected_convention);
        Ok(mapper)
    }

    /// Number of canonical chromosomes defined by the reference config.
    pub fn num_chromosomes(&self) -> usize {
        self.idx_to_chr.len()
    }

    /// Detect the naming convention based on the BAM/FASTA references.
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
        let half = self.num_chromosomes() / 2;
        if chr_count >= half && chr_count > acc_count {
            NamingConvention::ChrNames
        } else if acc_count >= half && acc_count > chr_count {
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
    /// Accepts accession IDs (`NC_060925.1`), short bare names (`1`, `X`,
    /// `MT`), or chr-style names (`chr1`); returns the canonical
    /// chr-style form. Unrecognized names pass through unchanged.
    pub fn to_chr_name(&self, name: &str) -> String {
        if let Some(chr) = self.acc_to_chr.get(name) {
            return chr.clone();
        }
        if let Some(&idx) = self.short_to_idx.get(name) {
            return self.idx_to_chr[idx].clone();
        }
        name.to_string()
    }

    /// Translate any name to its short bare form (e.g., `1`, `X`, `M`).
    /// Used by the MAF writer to keep per-row chrom bytes minimal.
    /// Unrecognized names pass through unchanged.
    pub fn to_short_name(&self, name: &str) -> String {
        if let Some(&idx) = self.chr_to_idx.get(name) {
            return self.idx_to_short[idx].clone();
        }
        if let Some(&idx) = self.acc_to_idx.get(name) {
            return self.idx_to_short[idx].clone();
        }
        if let Some(&idx) = self.short_to_idx.get(name) {
            return self.idx_to_short[idx].clone();
        }
        name.to_string()
    }

    /// Normalize a name to canonical chr format. Alias for `to_chr_name`.
    pub fn normalize(&self, name: &str) -> String {
        self.to_chr_name(name)
    }

    /// Check if a name (in any convention) exists in the BAM header.
    pub fn exists_in_bam(&self, name: &str) -> bool {
        if self.bam_refs.contains(name) {
            return true;
        }
        let translated = self.to_bam_name(name);
        self.bam_refs.contains(&translated)
    }

    /// Get the canonical chromosome index for a name in any convention.
    ///
    /// Accepts chr-style, accession, or short bare (`1`, `X`, `MT`) forms.
    /// Returns None for unrecognized names.
    pub fn get_chr_index(&self, name: &str) -> Option<usize> {
        if let Some(&idx) = self.chr_to_idx.get(name) {
            return Some(idx);
        }
        if let Some(&idx) = self.acc_to_idx.get(name) {
            return Some(idx);
        }
        if let Some(&idx) = self.short_to_idx.get(name) {
            return Some(idx);
        }
        None
    }

    /// Get the chr name for a given canonical index.
    pub fn chr_name_from_index(&self, idx: usize) -> Option<&str> {
        self.idx_to_chr.get(idx).map(|s| s.as_str())
    }

    /// Get the accession for a given canonical index.
    pub fn accession_from_index(&self, idx: usize) -> Option<&str> {
        self.idx_to_acc.get(idx).map(|s| s.as_str())
    }

    /// Map a chromosome string to its canonical index.
    ///
    /// Accepts the configured names/accessions directly, and also tolerates bare
    /// forms like "1"-"22", "X", "Y" (matched against "chr"-prefixed config names).
    /// Returns `None` for unrecognized names.
    pub fn parse_chr_index(&self, name: &str) -> Option<usize> {
        if let Some(&idx) = self.chr_to_idx.get(name) {
            return Some(idx);
        }
        if let Some(&idx) = self.acc_to_idx.get(name) {
            return Some(idx);
        }
        // Tolerate bare chromosome tokens ("1", "X", ...) by trying a "chr" prefix.
        let bare = name.strip_prefix("chr").unwrap_or(name);
        let candidate = format!("chr{}", bare);
        self.chr_to_idx.get(&candidate).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A T2T-CHM13v2.0 style contig list (chr1..chr22, chrX, chrY).
    fn t2t_contigs() -> Vec<Contig> {
        let accs = [
            "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1",
            "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1",
            "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1",
            "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1",
            "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1",
        ];
        contig_list(&accs)
    }

    /// A GRCh38 style contig list (chr1..chr22, chrX, chrY).
    fn grch38_contigs() -> Vec<Contig> {
        let accs = [
            "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10",
            "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11",
            "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10",
            "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
            "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10",
        ];
        contig_list(&accs)
    }

    fn contig_list(accs: &[&str]) -> Vec<Contig> {
        let names = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX", "chrY",
        ];
        names.iter().zip(accs.iter())
            .map(|(n, a)| Contig { name: n.to_string(), accession: a.to_string() })
            .collect()
    }

    fn refs_of<F: Fn(&Contig) -> String>(contigs: &[Contig], f: F) -> Vec<String> {
        contigs.iter().map(f).collect()
    }

    #[test]
    fn test_detection_chr_names() {
        let c = grch38_contigs();
        let refs = refs_of(&c, |x| x.name.clone());
        let mapper = ContigMapper::from_contigs_and_refs(&c, &refs);
        assert_eq!(mapper.detected_convention, NamingConvention::ChrNames);
    }

    #[test]
    fn test_detection_accession() {
        let c = grch38_contigs();
        let refs = refs_of(&c, |x| x.accession.clone());
        let mapper = ContigMapper::from_contigs_and_refs(&c, &refs);
        assert_eq!(mapper.detected_convention, NamingConvention::Accession);
    }

    #[test]
    fn test_to_bam_name_grch38_accession_convention() {
        let c = grch38_contigs();
        let refs = refs_of(&c, |x| x.accession.clone());
        let mapper = ContigMapper::from_contigs_and_refs(&c, &refs);

        // chr -> GRCh38 accession (for BAM query)
        assert_eq!(mapper.to_bam_name("chr1"), "NC_000001.11");
        assert_eq!(mapper.to_bam_name("chr21"), "NC_000021.9");
        assert_eq!(mapper.to_bam_name("chrY"), "NC_000024.10");
        // already correct
        assert_eq!(mapper.to_bam_name("NC_000001.11"), "NC_000001.11");
    }

    #[test]
    fn test_to_bam_name_t2t_chr_convention() {
        let c = t2t_contigs();
        let refs = refs_of(&c, |x| x.name.clone());
        let mapper = ContigMapper::from_contigs_and_refs(&c, &refs);

        assert_eq!(mapper.to_bam_name("NC_060925.1"), "chr1");
        assert_eq!(mapper.to_bam_name("NC_060947.1"), "chrX");
        assert_eq!(mapper.to_bam_name("chr1"), "chr1");
    }

    #[test]
    fn test_to_chr_name_and_index() {
        let c = grch38_contigs();
        let mapper = ContigMapper::from_contigs(&c);

        assert_eq!(mapper.to_chr_name("NC_000001.11"), "chr1");
        assert_eq!(mapper.to_chr_name("NC_000024.10"), "chrY");
        assert_eq!(mapper.to_chr_name("chr1"), "chr1");
        assert_eq!(mapper.to_chr_name("unknown"), "unknown");

        assert_eq!(mapper.num_chromosomes(), 24);
        assert_eq!(mapper.get_chr_index("chr1"), Some(0));
        assert_eq!(mapper.get_chr_index("chr22"), Some(21));
        assert_eq!(mapper.get_chr_index("chrX"), Some(22));
        assert_eq!(mapper.get_chr_index("chrY"), Some(23));
        assert_eq!(mapper.get_chr_index("NC_000021.9"), Some(20));
        assert_eq!(mapper.get_chr_index("unknown"), None);

        assert_eq!(mapper.chr_name_from_index(0), Some("chr1"));
        assert_eq!(mapper.accession_from_index(20), Some("NC_000021.9"));
        assert_eq!(mapper.chr_name_from_index(99), None);
    }

    #[test]
    fn test_parse_chr_index() {
        let c = grch38_contigs();
        let mapper = ContigMapper::from_contigs(&c);

        // Bare tokens
        assert_eq!(mapper.parse_chr_index("1"), Some(0));
        assert_eq!(mapper.parse_chr_index("22"), Some(21));
        assert_eq!(mapper.parse_chr_index("X"), Some(22));
        assert_eq!(mapper.parse_chr_index("Y"), Some(23));

        // With chr prefix
        assert_eq!(mapper.parse_chr_index("chr1"), Some(0));
        assert_eq!(mapper.parse_chr_index("chrX"), Some(22));

        // Accession
        assert_eq!(mapper.parse_chr_index("NC_000001.11"), Some(0));

        // Invalid
        assert_eq!(mapper.parse_chr_index("chrM"), None);
        assert_eq!(mapper.parse_chr_index("unknown"), None);
    }

    #[test]
    fn test_exists_in_bam() {
        let c = grch38_contigs();
        let refs = refs_of(&c, |x| x.accession.clone());
        let mapper = ContigMapper::from_contigs_and_refs(&c, &refs);

        assert!(mapper.exists_in_bam("NC_000001.11")); // direct
        assert!(mapper.exists_in_bam("chr1"));         // via translation
        assert!(!mapper.exists_in_bam("chrM"));
    }

    #[test]
    fn test_empty_mapper_is_identity() {
        let mapper = ContigMapper::default();
        assert_eq!(mapper.num_chromosomes(), 0);
        assert_eq!(mapper.to_chr_name("NC_000001.11"), "NC_000001.11");
        assert_eq!(mapper.to_bam_name("chr1"), "chr1");
        assert_eq!(mapper.get_chr_index("chr1"), None);
    }
}
