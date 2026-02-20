use std::collections::HashMap;
use std::fs;
use std::path::Path;

use clap::Parser;
use log::warn;
use nasvar::config::AggregateConfig;
use nasvar::output::{StructuralVariant, UnifiedOutput};
use nasvar::utils::contig::ContigMapper;

#[derive(Parser)]
#[command(name = "aggregate")]
#[command(about = "Make aggregate report (TSV) from analysis output directories")]
struct Cli {
    /// Analysis version (string)
    version: String,
    /// Path to aggregate configuration file (JSON). If not provided, uses defaults.
    #[arg(long)]
    config: Option<String>,
    /// Analysis output directories
    #[arg(required = true)]
    dirs: Vec<String>,
}

fn extract_sample_name(run_dir: &str) -> String {
    let path = Path::new(run_dir);
    let filename = path.file_name().and_then(|f| f.to_str()).unwrap_or(run_dir);

    // Extract name up to first '.'
    if let Some(idx) = filename.find('.') {
        filename[..idx].to_string()
    } else {
        filename.to_string()
    }
}

/// Search a directory for a file ending with any of the given suffixes
fn find_file_with_suffix(dir: &str, suffixes: &[&str]) -> Option<String> {
    let dir_path = Path::new(dir);
    if let Ok(entries) = fs::read_dir(dir_path) {
        for entry in entries.flatten() {
            if let Some(name) = entry.file_name().to_str() {
                for suffix in suffixes {
                    if name.ends_with(suffix) {
                        return Some(entry.path().to_string_lossy().into_owned());
                    }
                }
            }
        }
    }
    None
}

fn format_gene_text(name: &str, chr: &str, pos: u32, mapper: &ContigMapper) -> String {
    let chr_display = mapper.to_chr_name(chr);
    let pos_mbp = pos as f64 / 1_000_000.0;
    format!("{} ({} ~{:.2}Mbp)", name, chr_display, pos_mbp)
}

fn format_sv_list(variants: &[StructuralVariant]) -> String {
    variants
        .iter()
        .map(|d| {
            let chrom = d.chrom.as_deref().unwrap_or("?");
            format!(
                "{}: {} - {} ({} nt): {} supporting reads",
                chrom,
                d.start,
                d.end,
                d.end - d.start,
                d.reads
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn process_directory(
    run_dir: &str,
    config: &AggregateConfig,
) -> Option<HashMap<String, String>> {
    let name = extract_sample_name(run_dir);
    let mut result: HashMap<String, String> = HashMap::new();
    result.insert("name".to_string(), name.clone());
    result.insert("fusions".to_string(), String::new());
    result.insert("one_sided_fusions".to_string(), String::new());

    // Create mapper for chromosome name conversion
    let mapper = ContigMapper::new();

    // Find and load result.json
    let unified_path = match find_file_with_suffix(run_dir, &["result.json"]) {
        Some(p) => p,
        None => {
            warn!(
                "No *result.json file found in {} -- skipping this sample",
                run_dir
            );
            return None;
        }
    };

    let unified = match UnifiedOutput::load_json(&unified_path) {
        Ok(u) => u,
        Err(e) => {
            warn!(
                "Error loading {}: {} -- skipping this sample",
                unified_path, e
            );
            return None;
        }
    };

    // ----------------------- Fusions ----------------------------
    if let Some(ref fusions_data) = unified.fusions {
        let min_reads = config.min_supporting_reads as usize;
        let min_bp_reads = config.min_breakpoint_reads as usize;

        // Filter fusions
        let mut filtered_fusions = Vec::new();
        for fus in &fusions_data.fusions {
            let g0_name = &fus.gene1.name;
            let g1_name = &fus.gene2.name;

            // Skip blacklisted pairs
            if !g0_name.is_empty()
                && !g1_name.is_empty()
                && config.is_blacklisted_pair(g0_name, g1_name)
            {
                continue;
            }

            // Skip if supporting reads below threshold
            if fus.supporting_reads < min_reads {
                continue;
            }

            // Check breakpoint filter (special genes bypass this)
            let bypass = config.bypass_breakpoint_filter(g0_name)
                || config.bypass_breakpoint_filter(g1_name);
            let has_valid_breakpoint = fus.breakpoints.iter().any(|b| b.n_reads >= min_bp_reads);

            if has_valid_breakpoint || bypass {
                filtered_fusions.push(fus);
            }
        }

        // Separate two-sided and one-sided fusions
        let two_sided: Vec<_> = filtered_fusions
            .iter()
            .filter(|f| !f.gene1.name.is_empty() && !f.gene2.name.is_empty())
            .collect();

        let one_sided: Vec<_> = filtered_fusions
            .iter()
            .filter(|f| f.gene1.name.is_empty() || f.gene2.name.is_empty())
            .collect();

        // Format two-sided fusions
        if !two_sided.is_empty() {
            let mut fusions_text = String::new();
            for fus in &two_sided {
                let g0_txt =
                    format_gene_text(&fus.gene1.name, &fus.gene1.chr, fus.gene1.pos, &mapper);
                let g1_txt =
                    format_gene_text(&fus.gene2.name, &fus.gene2.chr, fus.gene2.pos, &mapper);
                let n_reads = fus.supporting_reads;
                fusions_text
                    .push_str(&format!("{} :: {} ({} reads); ", g0_txt, g1_txt, n_reads));
            }
            result.insert("fusions".to_string(), fusions_text);
        }

        // Format spike-in controls
        if let Some(spike_ins) = &fusions_data.spike_in
            && !spike_ins.is_empty() {
                let mut spike_text = String::new();
                for fus in spike_ins {
                    let g0_txt = format_gene_text(
                        &fus.gene1.name,
                        &fus.gene1.chr,
                        fus.gene1.pos,
                        &mapper,
                    );
                    let g1_txt = format_gene_text(
                        &fus.gene2.name,
                        &fus.gene2.chr,
                        fus.gene2.pos,
                        &mapper,
                    );
                    let n_reads = fus.supporting_reads;
                    spike_text
                        .push_str(&format!("{} :: {} ({} reads); ", g0_txt, g1_txt, n_reads));
                }
                result.insert("spike_in".to_string(), spike_text);
            }
        if !result.contains_key("spike_in") {
            result.insert("spike_in".to_string(), String::new());
        }

        // Format one-sided fusions
        if !one_sided.is_empty() {
            let mut fusions_text = String::new();
            for fus in &one_sided {
                let g0_txt =
                    format_gene_text(&fus.gene1.name, &fus.gene1.chr, fus.gene1.pos, &mapper);
                let g1_txt =
                    format_gene_text(&fus.gene2.name, &fus.gene2.chr, fus.gene2.pos, &mapper);
                let n_reads = fus.supporting_reads;
                fusions_text
                    .push_str(&format!("{} :: {} ({} reads); ", g0_txt, g1_txt, n_reads));
            }
            result.insert("one_sided_fusions".to_string(), fusions_text);
        }
    }

    // ----------------------- Karyotype ----------------------------
    if let Some(ref k) = unified.karyotype {
        let warnings = k
            .warnings
            .as_deref()
            .unwrap_or("")
            .split('\n')
            .collect::<Vec<_>>()
            .join(";");
        result.insert("warnings".to_string(), warnings);
        result.insert("karyotype".to_string(), k.karyotype_string.clone());
        result.insert("ISCN".to_string(), k.iscn_string.clone());

        // reads_aligned from karyotype
        if let Some(ra) = k.reads_aligned {
            result.insert("reads_aligned".to_string(), format!("{}", ra));
        }
    }

    // ----------------------- CNVs ----------------------------
    if let Some(ref cnv_data) = unified.cnv {
        let cnvs = &cnv_data.genes;

        // Local CN genes
        for gene in &config.local_cn_genes {
            if let Some(g) = cnvs.get(gene.as_str())
                && let Some(local) = g.local {
                    result.insert(
                        gene.clone(),
                        format!("{}x ({:.1}x)", local.round() as i64, local),
                    );
                }
        }

        // Focal CN genes
        for gene in &config.cnv_genes {
            if let Some(g) = cnvs.get(gene.as_str())
                && let Some(focal) = g.focal {
                    result.insert(
                        gene.clone(),
                        format!("{}x ({:.1}x)", focal.round() as i64, focal),
                    );
                }
        }

        // Deletion genes
        for gene in &config.deletion_genes {
            let del_key = format!("{}_deletions", gene);
            if let Some(g) = cnvs.get(gene.as_str()) {
                if let Some(deletions) = &g.deletions {
                    if !deletions.is_empty() {
                        result.insert(del_key.clone(), format_sv_list(deletions));
                    } else if let Some(focal) = g.focal {
                        // Fallback: show focal CN if no deletions detected
                        result.insert(
                            del_key.clone(),
                            format!("{}x  ({:.1}x)", focal.round() as i64, focal),
                        );
                    }
                } else if let Some(focal) = g.focal {
                    result.insert(
                        del_key.clone(),
                        format!("{}x  ({:.1}x)", focal.round() as i64, focal),
                    );
                }
            }
            result.entry(del_key).or_default();
        }

        // Duplication genes (CNV-detected intragenic duplications)
        for gene in &config.duplication_genes {
            let dup_key = format!("{}_ITD", gene);
            if let Some(g) = cnvs.get(gene.as_str())
                && let Some(dups) = &g.duplications
                    && !dups.is_empty() {
                        result.insert(dup_key.clone(), format_sv_list(dups));
                    }
            result.entry(dup_key).or_default();
        }
    }

    // ----------------------- ITD / insertions ----------------------------
    if let Some(ref itd_data) = unified.itd {
        for (gene, label) in &config.itd_genes {
            let key = format!("{}_{}", gene, label);
            if let Some(itds) = itd_data.genes.get(gene.as_str()) {
                if !itds.is_empty() {
                    let text = itds
                        .iter()
                        .map(|ins| {
                            let ar = if ins.coverage > ins.merged {
                                ins.merged as f64 / (ins.coverage - ins.merged) as f64
                            } else {
                                0.0
                            };
                            format!(
                                "{} nt insertion at position {}: {} reads (of {}); {:.2} AR",
                                ins.length, ins.position, ins.merged, ins.coverage, ar
                            )
                        })
                        .collect::<Vec<_>>()
                        .join("; ");
                    result.insert(key, text);
                } else {
                    result.insert(key, String::new());
                }
            } else {
                result.insert(key, String::new());
            }
        }
    }

    // ----------------------- SNVs / genotypes ----------------------------
    if let Some(ref snv_data) = unified.snv {
        for gene in &config.snv_genes {
            if let Some(g) = snv_data.genes.get(gene.as_str()) {
                if let Some(cov) = g.coverage {
                    result.insert(format!("{}_depth", gene), format!("{:.2}", cov));
                }

                let mutations_text: String = g
                    .mutations
                    .iter()
                    .map(|(k, v)| format!("{}: {} AF", k, v))
                    .collect::<Vec<_>>()
                    .join("; ");
                result.insert(format!("{}_mutations", gene), mutations_text);
                result.insert(format!("{}_genotype", gene), g.genotype.clone());
            }
        }
    }

    // QC metrics from unified output
    if let Some(ref qc) = unified.qc {
        // Override reads_aligned from QC if available (more authoritative than karyotype)
        if let Some(ra) = qc.reads_aligned {
            result.insert("reads_aligned".to_string(), format!("{}", ra));
        }
        if qc.reads_on_target > 0.0 {
            let read_length = (qc.nt_on_target / qc.reads_on_target) as i64;
            result.insert(
                "on_target_read_length".to_string(),
                format!("{}", read_length),
            );
        }
        if qc.target_regions_nt > 0.0 {
            let depth = qc.nt_on_target / qc.target_regions_nt;
            result.insert("on_target_depth".to_string(), format!("{:.2}", depth));
        }
    }

    // Metadata from unified output
    if let Some(ref meta) = unified.metadata {
        if let Some(v) = &meta.run_start_time {
            result.insert("seq_date".to_string(), v.clone());
        }
        if let Some(v) = &meta.run_id {
            result.insert("seq_ID".to_string(), v.clone());
        }
        if let Some(v) = &meta.basecall_model {
            result.insert("basecalling_model".to_string(), v.clone());
        }
        if let Some(v) = &meta.library_id {
            result.insert("library_ID".to_string(), v.clone());
        }
        if let Some(v) = &meta.sequencer_id {
            result.insert("sequencer_ID".to_string(), v.clone());
        }
        if let Some(v) = &meta.flow_cell_id {
            result.insert("flow_cell_ID".to_string(), v.clone());
        }
    }

    Some(result)
}

/// Build output columns from config gene lists
fn build_columns(config: &AggregateConfig) -> Vec<String> {
    if !config.columns.is_empty() {
        return config.columns.clone();
    }

    let mut cols: Vec<String> = vec![
        "name".into(),
        "karyotype".into(),
        "ISCN".into(),
        "warnings".into(),
        "fusions".into(),
        "one_sided_fusions".into(),
        "spike_in".into(),
        "on_target_depth".into(),
        "reads_aligned".into(),
        "on_target_read_length".into(),
    ];

    // SNV genes: {gene}_mutations, {gene}_genotype, {gene}_depth
    for gene in &config.snv_genes {
        cols.push(format!("{}_mutations", gene));
        cols.push(format!("{}_genotype", gene));
        cols.push(format!("{}_depth", gene));
    }

    // ITD / insertion genes (from ITD caller): {gene}_{label}
    let mut sorted_itd: Vec<_> = config.itd_genes.iter().collect();
    sorted_itd.sort_by_key(|(name, _)| *name);
    for (gene, label) in sorted_itd {
        cols.push(format!("{}_{}", gene, label));
    }

    // Duplication genes (from CNV caller): {gene}_ITD (skip if already present)
    for gene in &config.duplication_genes {
        let col = format!("{}_ITD", gene);
        if !cols.contains(&col) {
            cols.push(col);
        }
    }

    // Local CN genes: {gene}
    for gene in &config.local_cn_genes {
        cols.push(gene.clone());
    }

    // Focal CN genes: {gene}
    for gene in &config.cnv_genes {
        cols.push(gene.clone());
    }

    // Deletion genes: {gene}_deletions
    for gene in &config.deletion_genes {
        cols.push(format!("{}_deletions", gene));
    }

    // Metadata columns
    cols.extend(vec![
        "seq_date".into(),
        "library_ID".into(),
        "basecalling_model".into(),
        "sequencer_ID".into(),
        "flow_cell_ID".into(),
        "seq_ID".into(),
    ]);

    cols
}

fn main() {
    let cli = Cli::parse();

    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .format_module_path(false)
        .init();

    let agg_config = match &cli.config {
        Some(path) => match AggregateConfig::load(path) {
            Ok(c) => c,
            Err(e) => {
                warn!(
                    "Could not load config {}: {} -- using defaults",
                    path, e
                );
                AggregateConfig::default()
            }
        },
        None => AggregateConfig::default(),
    };

    let mut aggregate_results: Vec<HashMap<String, String>> = Vec::new();

    for run_dir in &cli.dirs {
        if let Some(result) = process_directory(run_dir, &agg_config) {
            aggregate_results.push(result);
        }
    }

    // Build output columns from config
    let cols = build_columns(&agg_config);

    // Print header
    println!("{}", cols.join("\t"));

    // Print rows
    for r in &aggregate_results {
        let row: Vec<&str> = cols
            .iter()
            .map(|c| r.get(c.as_str()).map(|s| s.as_str()).unwrap_or("MISSING"))
            .collect();
        println!("{}", row.join("\t"));
    }
}
