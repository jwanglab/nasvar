use log::{info, warn};
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::Path;

use crate::config::PipelineConfig;
use crate::output::{
    CnvGeneResult, CnvOutput, FusionBreakpoint, FusionEvent, FusionsOutput, GeneInfo,
    ItdOutput, KaryotypeOutput, QcOutput, SnvOutput, UnifiedOutput,
};
use crate::utils::metadata::SequencingMetaData;

// ---------------------------------------------------------------------------
// JSON loading helpers
// ---------------------------------------------------------------------------

fn load_json<T: serde::de::DeserializeOwned>(path: &str) -> Option<T> {
    if !Path::new(path).exists() {
        return None;
    }
    match File::open(path) {
        Ok(file) => {
            let reader = BufReader::new(file);
            match serde_json::from_reader(reader) {
                Ok(v) => Some(v),
                Err(e) => {
                    warn!("Failed to parse {}: {}", path, e);
                    None
                }
            }
        }
        Err(e) => {
            warn!("Failed to open {}: {}", path, e);
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn generate_report(out_prefix: &str, config: Option<&PipelineConfig>, targets_bed: Option<&str>) -> Result<(), Box<dyn std::error::Error>> {
    // Load unified output (sole source for all data)
    let mut unified: UnifiedOutput = load_json(&format!("{}.result.json", out_prefix))
        .unwrap_or_default();

    // Always use the current binary version for the report header
    unified.version = env!("CARGO_PKG_VERSION").to_string();

    // Load analysis targets BED file if provided
    let target_genes: Vec<String> = if let Some(bed_path) = targets_bed {
        match crate::utils::bed::read_bed(bed_path) {
            Ok(regions) => {
                let mut genes: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
                for region in regions {
                    genes.insert(region.name.clone());
                }
                genes.into_iter().collect()
            }
            Err(e) => {
                log::warn!("Could not load targets BED file {}: {}", bed_path, e);
                Vec::new()
            }
        }
    } else {
        Vec::new()
    };

    // Build markdown
    let md = build_markdown(out_prefix, &unified, config, &target_genes);

    // Write markdown
    let md_path = format!("{}.report.md", out_prefix);
    {
        let mut f = File::create(&md_path)?;
        f.write_all(md.as_bytes())?;
    }
    info!("Markdown report written to: {}", md_path);

    // Write HTML
    let html_path = format!("{}.report.html", out_prefix);
    {
        let html = md_to_html(&md);
        let mut f = File::create(&html_path)?;
        f.write_all(html.as_bytes())?;
    }
    info!("HTML report written to: {}", html_path);

    // Stdout summary
    let karyo = &unified.karyotype;
    let default_snvs = SnvOutput::default();
    let default_cnvs = CnvOutput::default();
    let default_itds = ItdOutput::default();
    let default_fus = FusionsOutput::default();
    let snvs = unified.snv.as_ref().unwrap_or(&default_snvs);
    let cnvs = unified.cnv.as_ref().unwrap_or(&default_cnvs);
    let itds = unified.itd.as_ref().unwrap_or(&default_itds);
    let fusions = unified.fusions.as_ref().unwrap_or(&default_fus);
    let qc = unified.qc.as_ref();
    print_stdout_summary(out_prefix, karyo, snvs, cnvs, itds, fusions, qc);

    Ok(())
}

// ---------------------------------------------------------------------------
// Markdown builder
// ---------------------------------------------------------------------------

fn build_markdown(
    out_prefix: &str,
    unified: &UnifiedOutput,
    config: Option<&PipelineConfig>,
    enriched_genes: &[String],
) -> String {
    let default_snvs = SnvOutput::default();
    let default_cnvs = CnvOutput::default();
    let default_itds = ItdOutput::default();
    let default_fus = FusionsOutput::default();

    let karyo = &unified.karyotype;
    let snvs = unified.snv.as_ref().unwrap_or(&default_snvs);
    let cnvs = unified.cnv.as_ref().unwrap_or(&default_cnvs);
    let itds = unified.itd.as_ref().unwrap_or(&default_itds);
    let fus = unified.fusions.as_ref().unwrap_or(&default_fus);
    let qc = unified.qc.as_ref();
    let metadata = unified.metadata.as_ref();
    let version = &unified.version;

    let mut md = String::with_capacity(8192);

    // Header
    md.push_str("# Adaptive Whole Genome Sequencing Report\n\n");
    md.push_str("**FOR RESEARCH USE ONLY**\n\n");

    // Run info
    let run_id = out_prefix
        .rsplit('/')
        .next()
        .unwrap_or(out_prefix);
    let now = crate::utils::time::utc_now_iso8601();
    md.push_str(&format!("**Run ID:** {}\n\n", run_id));
    md.push_str(&format!("**Date/time of report:** {}\n\n", now));
    md.push_str(&format!("**Analysis version:** {}\n\n", version));

    // Specimen details
    md.push_str("## Specimen details\n\n");
    md.push_str("**Specimen ID:**\n\n");
    md.push_str("**Date of specimen collection:**\n\n");
    md.push_str("**Specimen source:**\n\n");

    // ======================= Results =======================
    md.push_str("# Results\n\n");

    // --- Gene Fusions ---
    build_fusions_section(&mut md, fus, config);

    // --- Karyotype ---
    build_karyotype_section(&mut md, karyo, out_prefix);

    // --- CNVs ---
    build_cnv_section(&mut md, cnvs, config);

    // --- ITDs ---
    build_itd_section(&mut md, itds, config);

    // --- SNVs ---
    build_snv_section(&mut md, snvs, config);

    // --- QC ---
    build_qc_section(&mut md, karyo, qc);

    // --- Sequencing details ---
    build_seq_details_section(&mut md, metadata);

    // --- Methodology ---
    build_methodology_section(&mut md, config);

    // --- Gene enrichment set ---
    build_enrichment_table(&mut md, enriched_genes);

    md
}

// --- Section builders ---

fn build_fusions_section(md: &mut String, fus: &FusionsOutput, config: Option<&PipelineConfig>) {
    md.push_str("## Gene Fusions\n\n");

    let min_reads = config.map_or(3, |c| c.thresholds.fusions.min_supporting_reads) as usize;
    let min_bp_reads = config.map_or(3, |c| c.thresholds.fusions.min_breakpoint_reads) as usize;

    // Filter fusions: blacklist + min reads + breakpoint filter
    let filtered: Vec<&FusionEvent> = fus.fusions.iter().filter(|f| {
        if config.is_some_and(|c| c.is_blacklisted_pair(&f.gene1.name, &f.gene2.name)) {
            return false;
        }
        if f.supporting_reads < min_reads {
            return false;
        }
        // Breakpoint filter: genes in skip_self_fusions bypass this check
        let bypass = config.is_some_and(|c| {
            c.skip_self_fusion(&f.gene1.name) || c.skip_self_fusion(&f.gene2.name)
        });
        if bypass {
            return true;
        }
        if f.breakpoints.is_empty() {
            true // no breakpoint data available, pass through
        } else {
            f.breakpoints.iter().any(|bp| bp.n_reads >= min_bp_reads)
        }
    }).collect();

    let two_sided: Vec<&&FusionEvent> = filtered.iter()
        .filter(|f| !f.gene1.name.is_empty() && !f.gene2.name.is_empty())
        .collect();
    let one_sided: Vec<&&FusionEvent> = filtered.iter()
        .filter(|f| f.gene1.name.is_empty() || f.gene2.name.is_empty())
        .collect();

    if !two_sided.is_empty() {
        md.push_str("| Gene/region 1 | Gene/region 2 | Supporting reads | Breakpoint(s) |\n");
        md.push_str("|---|---|---|---|\n");
        for f in &two_sided {
            let g1_info = format_gene_info(&f.gene1);
            let g2_info = format_gene_info(&f.gene2);
            let bp_str = format_breakpoints(&f.breakpoints);
            let reads_str = if f.supporting_reads >= 3 {
                format!("**{}**", f.supporting_reads)
            } else {
                format!("{}", f.supporting_reads)
            };
            md.push_str(&format!("| {} | {} | {} | {} |\n", g1_info, g2_info, reads_str, bp_str));
        }
        md.push('\n');
    } else {
        md.push_str("**No fusions detected**\n\n");
    }

    md.push_str("### Rearrangements with only one side in the panel\n\n");
    if !one_sided.is_empty() {
        md.push_str("| Gene/region 1 | Gene/region 2 | Supporting reads | Breakpoint(s) |\n");
        md.push_str("|---|---|---|---|\n");
        for f in &one_sided {
            let g1_info = format_gene_info(&f.gene1);
            let g2_info = format_gene_info(&f.gene2);
            let bp_str = format_breakpoints(&f.breakpoints);
            let reads_str = if f.supporting_reads >= 3 {
                format!("**{}**", f.supporting_reads)
            } else {
                format!("{}", f.supporting_reads)
            };
            md.push_str(&format!("| {} | {} | {} | {} |\n", g1_info, g2_info, reads_str, bp_str));
        }
        md.push('\n');
    } else {
        md.push_str("**None detected**\n\n");
    }
}

fn format_gene_info(g: &GeneInfo) -> String {
    if !g.chr.is_empty() && g.pos != 0 {
        format!("{} ({} ~{:.2}Mbp)", g.name, g.chr, g.pos as f64 / 1_000_000.0)
    } else {
        g.name.clone()
    }
}

fn format_breakpoints(bps: &[FusionBreakpoint]) -> String {
    if bps.is_empty() {
        return String::new();
    }
    bps.iter()
        .map(|bp| {
            let d0 = bp.gene0_dir.arrow();
            let d1 = bp.gene1_dir.arrow();
            let overlap_info = match bp.overlap_gap {
                Some(v) if v > 0 => format!(" [{}bp hom]", v),
                Some(v) if v < 0 => format!(" [{}bp ins]", -v),
                _ => String::new(),
            };
            format!(
                "{}:{}{} -- {}:{}{}  ({} reads{})",
                bp.gene0_chr, bp.gene0_pos, d0,
                bp.gene1_chr, bp.gene1_pos, d1,
                bp.n_reads, overlap_info,
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn build_karyotype_section(md: &mut String, karyo: &Option<KaryotypeOutput>, out_prefix: &str) {
    md.push_str("## Inferred Full Chromosome or Arm Level Gain or Loss\n\n");

    if let Some(k) = karyo {
        if let Some(w) = &k.warnings {
            for line in w.split('\n') {
                md.push_str(&format!("> **WARNING:** {}\n", line));
            }
            md.push('\n');
            md.push_str(&format!("Digital karyotype: **{}**\n\n", k.karyotype_string));
            md.push_str(&format!("ISCN nomenclature: **{}**\n\n", k.iscn_string));
        } else {
            md.push_str(&format!("Digital karyotype: **{}**\n\n", k.karyotype_string));
            md.push_str(&format!("ISCN nomenclature: **{}**\n\n", k.iscn_string));
        }
        if let Some(br) = k.blast_ratio {
            md.push_str(&format!("Blast ratio (tumor fraction): **{:.2}%**\n\n", br * 100.0));
        }

        // Karyotype + BAF plot (prefer GC-corrected with BAF)
        let baf_plot = format!("{}.karyotype_baf.gc_corrected.svg", out_prefix);
        let karyo_gc_plot = format!("{}.karyotype.gc_corrected.svg", out_prefix);
        if Path::new(&baf_plot).exists() {
            md.push_str(&format!("![Karyotype + BAF]({})\n\n", baf_plot));
        } else if Path::new(&karyo_gc_plot).exists() {
            md.push_str(&format!("![Karyotype]({})\n\n", karyo_gc_plot));
        }
    } else {
        md.push_str("Not run\n\n");
    }
}

fn build_cnv_section(md: &mut String, cnvs: &CnvOutput, config: Option<&PipelineConfig>) {
    md.push_str("## Focal copy-number variation\n\n");

    if let Some(cfg) = config {
        let del_set: std::collections::HashSet<&str> = cfg.genes.cnv.deletion_genes.iter().map(|s| s.as_str()).collect();

        // 1. Local CN genes
        for gene in &cfg.genes.cnv.local_cn_genes {
            md.push_str(&format!("### {}\n\n", gene));
            if let Some(c) = cnvs.genes.get(gene.as_str()) {
                if let Some(local) = c.local {
                    md.push_str(&format!(
                        "Local read-depth based copy number estimate: **{}x ({:.1}x)**\n\n",
                        local.round() as i64, local
                    ));
                }
            } else {
                md.push_str("Not tested\n\n");
            }
        }

        // 2. Focal CN genes (skip local_cn_genes; for genes also in deletion_genes, show focal CN only if no deletions)
        for gene in &cfg.genes.cnv.focal_genes {
            if cfg.genes.cnv.local_cn_genes.contains(gene) {
                continue;
            }
            if del_set.contains(gene.as_str()) {
                // Dual-role gene: show focal CN only if no deletions detected
                if let Some(c) = cnvs.genes.get(gene.as_str()) {
                    let has_dels = c.deletions.as_ref().map(|d| !d.is_empty()).unwrap_or(false);
                    if !has_dels {
                        md.push_str(&format!("### {}\n\n", gene));
                        if let Some(focal) = c.focal {
                            md.push_str(&format!(
                                "Focal sequencing-depth based copy number estimate: **{}x ({:.1}x)**\n\n",
                                focal.round() as i64, focal
                            ));
                        }
                    }
                }
            } else {
                md.push_str(&format!("### {}\n\n", gene));
                if let Some(c) = cnvs.genes.get(gene.as_str()) {
                    if let Some(focal) = c.focal {
                        md.push_str(&format!(
                            "Focal sequencing-depth based copy number estimate: **{}x ({:.1}x)**\n\n",
                            focal.round() as i64, focal
                        ));
                    }
                } else {
                    md.push_str("Not tested\n\n");
                }
            }
        }

        // 3. Deletion genes
        for gene in &cfg.genes.cnv.deletion_genes {
            md.push_str(&format!("### {} intragenic deletions\n\n", gene));
            format_structural_variants(md, cnvs.genes.get(gene.as_str()), true);
        }

        // 4. Duplication genes
        for gene in &cfg.genes.cnv.duplication_genes {
            md.push_str(&format!("### {} intragenic tandem duplications\n\n", gene));
            format_structural_variants(md, cnvs.genes.get(gene.as_str()), false);
        }
    } else {
        // No config: iterate all genes in the result data
        let mut genes: Vec<&String> = cnvs.genes.keys().collect();
        genes.sort();
        for gene in genes {
            let c = &cnvs.genes[gene];
            md.push_str(&format!("### {}\n\n", gene));
            if let Some(local) = c.local {
                md.push_str(&format!(
                    "Local read-depth based copy number estimate: **{}x ({:.1}x)**\n\n",
                    local.round() as i64, local
                ));
            }
            if let Some(focal) = c.focal {
                md.push_str(&format!(
                    "Focal sequencing-depth based copy number estimate: **{}x ({:.1}x)**\n\n",
                    focal.round() as i64, focal
                ));
            }
            let has_dels = c.deletions.as_ref().map(|d| !d.is_empty()).unwrap_or(false);
            if has_dels {
                md.push_str("**Intragenic deletions:**\n\n");
                format_structural_variants(md, Some(c), true);
            }
            let has_dups = c.duplications.as_ref().map(|d| !d.is_empty()).unwrap_or(false);
            if has_dups {
                md.push_str("**Intragenic duplications:**\n\n");
                format_structural_variants(md, Some(c), false);
            }
        }
    }
}

fn format_structural_variants(md: &mut String, gene: Option<&CnvGeneResult>, deletions: bool) {
    match gene {
        Some(c) => {
            let variants = if deletions { &c.deletions } else { &c.duplications };
            match variants {
                Some(vs) if !vs.is_empty() => {
                    for d in vs {
                        let chrom = d.chrom.as_deref().unwrap_or("?");
                        md.push_str(&format!(
                            "{}: {} - {} ({} nt): **{} supporting reads**\n\n",
                            chrom, d.start, d.end, d.end - d.start, d.reads
                        ));
                    }
                }
                _ => md.push_str("None detected\n\n"),
            }
        }
        None => md.push_str("Not tested\n\n"),
    }
}

fn build_itd_section(md: &mut String, itds: &ItdOutput, config: Option<&PipelineConfig>) {
    md.push_str("## Internal/partial tandem duplications\n\n");

    let genes: Vec<String> = if let Some(cfg) = config {
        let mut g: Vec<String> = cfg.genes.itd.keys().cloned().collect();
        g.sort();
        g
    } else {
        let mut g: Vec<String> = itds.genes.keys().cloned().collect();
        g.sort();
        g
    };

    for gene in &genes {
        md.push_str(&format!("### {} internal tandem duplication (ITD)\n\n", gene));
        match itds.genes.get(gene.as_str()) {
            Some(entries) if !entries.is_empty() => {
                for e in entries {
                    let ar = if e.coverage > e.merged {
                        e.merged as f64 / (e.coverage - e.merged) as f64
                    } else {
                        e.merged as f64 / e.coverage as f64
                    };
                    md.push_str(&format!(
                        "{} nt insertion at position {}: {} reads (of {}); {:.2} AR\n\n",
                        e.length, e.position, e.merged, e.coverage, ar
                    ));
                }
            }
            _ => md.push_str("None detected\n\n"),
        }
    }
}

fn build_snv_section(md: &mut String, snvs: &SnvOutput, config: Option<&PipelineConfig>) {
    md.push_str("## SNV Genotypes\n\n");

    if let Some(cfg) = config {
        let pathogenic_genes: Vec<&str> = cfg.genes.snv.pathogenic.iter().map(|s| s.as_str()).collect();
        let pgx_genes: Vec<&str> = cfg.genes.snv.pharmacogenomics.iter().map(|s| s.as_str()).collect();

        // Disease/pathogenic variants
        if !pathogenic_genes.is_empty() {
            md.push_str("### Disease variants\n\n");
            for gene_name in &pathogenic_genes {
                format_snv_gene(md, snvs, gene_name);
            }
        }

        // Pharmacogenomic variants
        if !pgx_genes.is_empty() {
            md.push_str("### Pharmacogenomic variants\n\n");
            for gene_name in &pgx_genes {
                format_snv_gene(md, snvs, gene_name);
            }
        }
    } else {
        // No config: show all genes from result data
        let mut all_genes: Vec<&String> = snvs.genes.keys().collect();
        all_genes.sort();
        for gene_name in &all_genes {
            format_snv_gene(md, snvs, gene_name);
        }
    }
}

fn format_snv_gene(md: &mut String, snvs: &SnvOutput, gene_name: &str) {
    md.push_str(&format!("**{}**\n\n", gene_name));
    if let Some(g) = snvs.genes.get(gene_name) {
        if let Some(cov) = g.coverage {
            md.push_str(&format!("Sequencing depth: {:.2}x\n\n", cov));
        }
        md.push_str("Mutations:\n\n");
        if g.mutations.is_empty() {
            md.push_str("(none detected)\n\n");
        } else {
            for (m, af) in &g.mutations {
                md.push_str(&format!("- {}: {} AF\n", m, af));
            }
            md.push('\n');
        }
        if !g.aa_changes.is_empty() {
            md.push_str("Amino acid changes:\n\n");
            for aa in &g.aa_changes {
                md.push_str(&format!("- {}\n", aa));
            }
            md.push('\n');
        }
        md.push_str(&format!("Genotype: **{}**\n\n", g.genotype));
    } else {
        md.push_str("Not tested\n\n");
    }
}

fn format_snv_methodology_gene(md: &mut String, gene: &str, transcript_cfg: Option<&crate::config::SnvTranscriptConfig>) {
    if let Some(cfg) = transcript_cfg {
        let mut variants: Vec<&String> = cfg.variants.keys().collect();
        variants.sort();
        let variants_str = variants.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ");
        if let Some(transcript) = &cfg.transcript {
            md.push_str(&format!("- **{}** ({}): {}\n", gene, transcript, variants_str));
        } else {
            md.push_str(&format!("- **{}**: {}\n", gene, variants_str));
        }
    } else {
        md.push_str(&format!("- **{}**\n", gene));
    }
}

fn build_qc_section(md: &mut String, karyo: &Option<KaryotypeOutput>, qc: Option<&QcOutput>) {
    md.push_str("## Quality control metrics\n\n");

    md.push_str("- Blast percentage reported by hematopathology: \n");

    // Try to get reads_aligned from QC first, fall back to karyotype
    let reads_aligned = qc.and_then(|q| q.reads_aligned)
        .or_else(|| karyo.as_ref().and_then(|k| k.reads_aligned));
    if let Some(ra) = reads_aligned {
        md.push_str(&format!("- Genome-wide reads aligned: **{:}**\n", format_number(ra)));
    }

    if let Some(q) = qc {
        if q.reads_on_target > 0.0 {
            let mean_len = q.nt_on_target / q.reads_on_target;
            md.push_str(&format!("- Mean on-target read length: **{} nt**\n", mean_len as u64));
        }
        if q.target_regions_nt > 0.0 {
            let mean_cov = q.nt_on_target / q.target_regions_nt;
            md.push_str(&format!("- Mean coverage over target regions: **{:.2}x**\n", mean_cov));
        }
    }
    md.push('\n');
}

fn build_seq_details_section(md: &mut String, metadata: Option<&SequencingMetaData>) {
    if let Some(meta) = metadata {
        // Only show if we have at least one field
        let has_any = meta.run_start_time.is_some()
            || meta.run_id.is_some()
            || meta.basecall_model.is_some()
            || meta.library_id.is_some()
            || meta.sequencer_id.is_some()
            || meta.flow_cell_id.is_some();
        if !has_any {
            return;
        }
        md.push_str("## Sequencing details\n\n");
        if let Some(v) = &meta.run_start_time {
            md.push_str(&format!("**Date/time of sequencing start:** {}\n\n", v));
        }
        if let Some(v) = &meta.run_id {
            md.push_str(&format!("**Sequencing run ID:** {}\n\n", v));
        }
        if let Some(v) = &meta.basecall_model {
            md.push_str(&format!("**Basecalling model:** {}\n\n", v));
        }
        if let Some(v) = &meta.library_id {
            md.push_str(&format!("**Library ID:** {}\n\n", v));
        }
        if let Some(v) = &meta.sequencer_id {
            md.push_str(&format!("**Sequencer ID:** {}\n\n", v));
        }
        if let Some(v) = &meta.flow_cell_id {
            md.push_str(&format!("**Flow cell ID:** {}\n\n", v));
        }
    }
}

fn build_methodology_section(md: &mut String, config: Option<&PipelineConfig>) {
    md.push_str("---\n\n");
    md.push_str("## Methodology\n\n");

    // Digital Karyotype (static - no config needed)
    md.push_str("**Digital Karyotype and Copy Number Variants:** ");
    md.push_str("Copy number across the genome at the chromosome level ('digital karyotype') is inferred based on relative sequencing depth and minor allele frequency (MAF). To avoid the potentially confounding effect of adaptive sampling on relative sequencing depth assessment, coarse-scale depth is constructed as a function of reads per million base pairs (Mbp), where each read contributes a count of one to the bin in which the center of the read aligns. Depth is corrected for GC bias using a linear correction model.\n\n");

    // ITD - dynamic from config.genes.itd
    md.push_str("**Internal Tandem Duplications**\n\n");
    if let Some(cfg) = config {
        let mut genes: Vec<_> = cfg.genes.itd.iter().collect();
        genes.sort_by_key(|(name, _)| *name);
        for (gene, region) in genes {
            md.push_str(&format!(
                "**{}:** Insertions are called among reads aligning to {}:{}-{}",
                gene, region.chrom, region.start, region.end
            ));
            if region.min_length > 0 {
                md.push_str(&format!(" (min length: {}bp", region.min_length));
                if region.min_frequency > 0.0 {
                    md.push_str(&format!(", min frequency: {:.0}%)", region.min_frequency * 100.0));
                } else {
                    md.push(')');
                }
            }
            md.push_str(". Called insertions are clustered together if their respective lengths are within 10% of one another and their position in the reference is within 110% of the length of the insertion from one another. Allelic ratio (AR) is calculated as the ratio of the number of reads supporting the insertion cluster divided by the number of reads aligning to the insertion site that are not in the insertion cluster.\n\n");
        }
    } else {
        md.push_str("*No ITD configuration available.*\n\n");
    }

    // SNV - dynamic from config.genes.snv
    md.push_str("**Single Nucleotide Variants**\n\n");
    if let Some(cfg) = config {
        // Collect genes by category
        let mut pgx_genes: Vec<&String> = cfg.genes.snv.pharmacogenomics.iter().collect();
        pgx_genes.sort();

        // Pharmacogenomics genes with their variants
        if !pgx_genes.is_empty() {
            md.push_str("Pharmacogenomics:\n\n");
            for gene in &pgx_genes {
                format_snv_methodology_gene(md, gene, cfg.genes.snv.transcripts.get(*gene));
            }
            md.push('\n');
        }

        // Pathogenic variants
        let mut pathogenic_genes: Vec<&String> = cfg.genes.snv.pathogenic.iter().collect();
        pathogenic_genes.sort();

        if !pathogenic_genes.is_empty() {
            md.push_str("Pathogenic variants:\n\n");
            for gene in &pathogenic_genes {
                format_snv_methodology_gene(md, gene, cfg.genes.snv.transcripts.get(*gene));
            }
            md.push('\n');
        }

        // Disease/subtype defining variants (not in pharmacogenomics or pathogenic)
        let mut disease_genes: Vec<(&String, &crate::config::SnvTranscriptConfig)> = cfg.genes.snv.transcripts
            .iter()
            .filter(|(gene, _)| {
                !cfg.genes.snv.pharmacogenomics.contains(*gene) &&
                !cfg.genes.snv.pathogenic.contains(*gene)
            })
            .collect();
        disease_genes.sort_by_key(|(gene, _)| *gene);

        if !disease_genes.is_empty() {
            md.push_str("Disease/subtype defining variants:\n\n");
            for (gene, transcript_cfg) in disease_genes {
                format_snv_methodology_gene(md, gene, Some(transcript_cfg));
            }
            md.push('\n');
        }
    } else {
        md.push_str("*No SNV configuration available.*\n\n");
    }

    // CNV - dynamic from config.genes.cnv
    md.push_str("**Focal Copy Number Analysis**\n\n");
    if let Some(cfg) = config {
        if !cfg.genes.cnv.focal_genes.is_empty() {
            let mut focal = cfg.genes.cnv.focal_genes.clone();
            focal.sort();
            md.push_str(&format!("Focal genes: {}\n\n", focal.join(", ")));
        }
        if !cfg.genes.cnv.deletion_genes.is_empty() {
            let mut del = cfg.genes.cnv.deletion_genes.clone();
            del.sort();
            md.push_str(&format!("Deletion screening: {}\n\n", del.join(", ")));
        }
        if !cfg.genes.cnv.duplication_genes.is_empty() {
            let mut dup = cfg.genes.cnv.duplication_genes.clone();
            dup.sort();
            md.push_str(&format!("Duplication screening: {}\n\n", dup.join(", ")));
        }
    } else {
        md.push_str("*No CNV configuration available.*\n\n");
    }

    // Fusions - dynamic one-sided genes
    md.push_str("**Gene Fusions**\n\n");
    if let Some(cfg) = config
        && !cfg.genes.fusions.one_sided_genes.is_empty() {
            let mut one_sided = cfg.genes.fusions.one_sided_genes.clone();
            one_sided.sort();
            md.push_str("One-sided fusion detection enabled for:\n\n");
            for gene in &one_sided {
                if let Some(partners) = cfg.genes.fusions.one_sided_partners.get(gene) {
                    let mut partners_sorted = partners.clone();
                    partners_sorted.sort();
                    md.push_str(&format!("- {} (partners: {})\n", gene, partners_sorted.join(", ")));
                } else {
                    md.push_str(&format!("- {} (any partner)\n", gene));
                }
            }
            md.push('\n');
        }
}

fn build_enrichment_table(md: &mut String, enriched_genes: &[String]) {
    if !enriched_genes.is_empty() {
        md.push_str("**Gene / Regions Enrichment set:**\n\n");
        let cols = 6;
        // Header row (empty cells)
        md.push_str(&format!("|{}|\n", vec![" "; cols].join("|")));
        // Separator row
        md.push_str(&format!("|{}|\n", vec!["---"; cols].join("|")));
        // Data rows
        for chunk in enriched_genes.chunks(cols) {
            let mut cells: Vec<String> = chunk.iter().map(|g| format!(" {} ", g)).collect();
            while cells.len() < cols {
                cells.push(String::from(" "));
            }
            md.push_str(&format!("|{}|\n", cells.join("|")));
        }
        md.push('\n');
    }
}

// ---------------------------------------------------------------------------
// Markdown → HTML converter
// ---------------------------------------------------------------------------

fn md_to_html(md: &str) -> String {
    let mut html = String::with_capacity(md.len() * 2);

    // Boilerplate
    html.push_str("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"UTF-8\">\n");
    html.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
    html.push_str("<title>Adaptive Whole Genome Sequencing Report</title>\n");
    html.push_str("<style>\n");
    html.push_str("body { font-family: Calibri, Arial, sans-serif; margin: 0.5in; line-height: 1.4; }\n");
    html.push_str("h1 { font-size: 1.6em; }\n");
    html.push_str("h2 { font-size: 1.3em; border-bottom: 1px solid #ccc; padding-bottom: 4px; }\n");
    html.push_str("h3 { font-size: 1.1em; }\n");
    html.push_str("table { border-collapse: collapse; margin: 8px 0 16px 0; }\n");
    html.push_str("th, td { border: 1px solid #999; padding: 4px 8px; text-align: left; }\n");
    html.push_str("th { font-weight: bold; background: #f0f0f0; }\n");
    html.push_str(".warning { color: #cc0000; font-weight: bold; }\n");
    html.push_str(".research-only { color: #ff0000; font-weight: bold; font-size: 1.1em; }\n");
    html.push_str(".karyotype-plot { width: 100%; max-width: 1600px; margin: 12px 0; }\n");
    html.push_str("hr { border: none; border-bottom: 1px solid #999; margin: 24px 0; }\n");
    html.push_str("</style>\n</head><body>\n");

    let lines: Vec<&str> = md.lines().collect();
    let mut i = 0;
    let mut in_table = false;
    let mut in_ul = false;
    let mut in_blockquote = false;

    while i < lines.len() {
        let line = lines[i];

        // Close open blockquote if line doesn't start with >
        if in_blockquote && !line.starts_with('>') {
            html.push_str("</div>\n");
            in_blockquote = false;
        }

        // Close open list if line doesn't start with -
        if in_ul && !line.starts_with("- ") {
            html.push_str("</ul>\n");
            in_ul = false;
        }

        // Close open table if line doesn't start with |
        if in_table && !line.starts_with('|') {
            html.push_str("</table>\n");
            in_table = false;
        }

        if line.is_empty() {
            // Skip blank lines (they separate paragraphs, handled by block elements)
            i += 1;
            continue;
        }

        // Horizontal rule
        if line == "---" {
            html.push_str("<hr>\n");
            i += 1;
            continue;
        }

        // Headings
        if let Some(stripped) = line.strip_prefix("### ") {
            html.push_str(&format!("<h3>{}</h3>\n", inline_format(&html_escape(stripped))));
            i += 1;
            continue;
        }
        if let Some(stripped) = line.strip_prefix("## ") {
            html.push_str(&format!("<h2>{}</h2>\n", inline_format(&html_escape(stripped))));
            i += 1;
            continue;
        }
        if let Some(stripped) = line.strip_prefix("# ") {
            html.push_str(&format!("<h1>{}</h1>\n", inline_format(&html_escape(stripped))));
            i += 1;
            continue;
        }

        // Blockquotes (warnings)
        if let Some(content) = line.strip_prefix("> ") {
            if !in_blockquote {
                html.push_str("<div class=\"warning\">\n");
                in_blockquote = true;
            }
            html.push_str(&format!("<p>{}</p>\n", inline_format(&html_escape(content))));
            i += 1;
            continue;
        }

        // Images → inline SVG or <img>
        if line.starts_with("![") {
            if let Some(path) = extract_image_path(line) {
                if path.ends_with(".svg") && Path::new(&path).exists() {
                    if let Ok(svg) = std::fs::read_to_string(&path) {
                        let svg = strip_svg_dimensions(&svg);
                        html.push_str(&format!("<div class=\"karyotype-plot\">{}</div>\n", svg));
                    }
                } else {
                    html.push_str(&format!("<img src=\"{}\" style=\"max-width:100%\">\n", html_escape(&path)));
                }
            }
            i += 1;
            continue;
        }

        // Tables
        if line.starts_with('|') {
            if !in_table {
                html.push_str("<table>\n");
                in_table = true;
                // First row is header
                let cells = parse_table_row(line);
                html.push_str("<tr>");
                for c in &cells {
                    html.push_str(&format!("<th>{}</th>", inline_format(&html_escape(c))));
                }
                html.push_str("</tr>\n");
                i += 1;
                continue;
            }
            // Skip separator row
            if line.contains("---|") || line.contains("--- |") {
                i += 1;
                continue;
            }
            let cells = parse_table_row(line);
            html.push_str("<tr>");
            for c in &cells {
                html.push_str(&format!("<td>{}</td>", inline_format(&html_escape(c))));
            }
            html.push_str("</tr>\n");
            i += 1;
            continue;
        }

        // Bullet lists
        if let Some(stripped) = line.strip_prefix("- ") {
            if !in_ul {
                html.push_str("<ul>\n");
                in_ul = true;
            }
            html.push_str(&format!("<li>{}</li>\n", inline_format(&html_escape(stripped))));
            i += 1;
            continue;
        }

        // Regular paragraph
        html.push_str(&format!("<p>{}</p>\n", inline_format(&html_escape(line))));
        i += 1;
    }

    // Close any open blocks
    if in_blockquote { html.push_str("</div>\n"); }
    if in_ul { html.push_str("</ul>\n"); }
    if in_table { html.push_str("</table>\n"); }

    html.push_str("</body></html>\n");
    html
}

/// Apply inline formatting: **bold**, *italic*
fn inline_format(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if i + 1 < chars.len() && chars[i] == '*' && chars[i + 1] == '*' {
            // Find closing **
            if let Some(end) = find_closing_double_star(&chars, i + 2) {
                result.push_str("<b>");
                let inner: String = chars[i + 2..end].iter().collect();
                result.push_str(&inner);
                result.push_str("</b>");
                i = end + 2;
                continue;
            }
        }
        if chars[i] == '*' && (i + 1 < chars.len()) && chars[i + 1] != '*' {
            // Find closing *
            if let Some(end) = find_closing_star(&chars, i + 1) {
                result.push_str("<i>");
                let inner: String = chars[i + 1..end].iter().collect();
                result.push_str(&inner);
                result.push_str("</i>");
                i = end + 1;
                continue;
            }
        }
        result.push(chars[i]);
        i += 1;
    }
    result
}

fn find_closing_double_star(chars: &[char], start: usize) -> Option<usize> {
    let mut i = start;
    while i + 1 < chars.len() {
        if chars[i] == '*' && chars[i + 1] == '*' {
            return Some(i);
        }
        i += 1;
    }
    None
}

fn find_closing_star(chars: &[char], start: usize) -> Option<usize> {
    let mut i = start;
    while i < chars.len() {
        if chars[i] == '*' {
            return Some(i);
        }
        i += 1;
    }
    None
}

fn parse_table_row(line: &str) -> Vec<String> {
    let trimmed = line.trim().trim_matches('|');
    trimmed.split('|').map(|s| s.trim().to_string()).collect()
}

fn extract_image_path(line: &str) -> Option<String> {
    // ![alt](path)
    let start = line.find("](")?;
    let end = line[start + 2..].find(')')?;
    Some(line[start + 2..start + 2 + end].to_string())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn strip_svg_dimensions(svg: &str) -> String {
    let mut result = svg.to_string();
    if let Some(start) = result.find("<svg")
        && let Some(end) = result[start..].find('>') {
            let tag = result[start..start + end + 1].to_string();
            let mut new_tag = tag.clone();
            // Remove width="..."
            while let Some(w) = new_tag.find(" width=\"") {
                if let Some(close) = new_tag[w + 7..].find('"') {
                    new_tag = format!("{}{}", &new_tag[..w], &new_tag[w + 8 + close..]);
                } else {
                    break;
                }
            }
            // Remove height="..."
            while let Some(w) = new_tag.find(" height=\"") {
                if let Some(close) = new_tag[w + 9..].find('"') {
                    new_tag = format!("{}{}", &new_tag[..w], &new_tag[w + 10 + close..]);
                } else {
                    break;
                }
            }
            result = format!("{}{}{}", &result[..start], new_tag, &result[start + end + 1..]);
        }
    result
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

// ---------------------------------------------------------------------------
// Stdout summary (abbreviated)
// ---------------------------------------------------------------------------

fn print_stdout_summary(
    out_prefix: &str,
    karyo: &Option<KaryotypeOutput>,
    snvs: &SnvOutput,
    cnvs: &CnvOutput,
    itds: &ItdOutput,
    fus: &FusionsOutput,
    qc: Option<&QcOutput>,
) {
    println!("--------------------------------------------------------------------------------");
    println!("                           NASVAR AGGREGATE REPORT                              ");
    println!("--------------------------------------------------------------------------------");
    println!("Sample Prefix: {}", out_prefix);
    println!();

    if let Some(k) = karyo {
        println!("=== Karyotype ===");
        if let Some(w) = &k.warnings {
            println!("WARNING: {}", w);
        }
        println!("Karyotype: {}", k.karyotype_string);
        println!("ISCN:      {}", k.iscn_string);
        if let Some(br) = k.blast_ratio {
            println!("Blast Ratio: {:.2}%", br * 100.0);
        }
        println!();
    }

    if !fus.fusions.is_empty() {
        println!("=== Gene Fusions ===");
        for f in &fus.fusions {
            println!("  {}-{}: {} reads", f.gene1.name, f.gene2.name, f.supporting_reads);
        }
        println!();
    }

    if !snvs.genes.is_empty() {
        println!("=== SNVs ===");
        let mut keys: Vec<&String> = snvs.genes.keys().collect();
        keys.sort();
        for gene in keys {
            let d = &snvs.genes[gene];
            println!("  {}: {}", gene, d.genotype);
        }
        println!();
    }

    if !cnvs.genes.is_empty() {
        println!("=== CNVs ===");
        let mut keys: Vec<&String> = cnvs.genes.keys().collect();
        keys.sort();
        for gene in keys {
            let d = &cnvs.genes[gene];
            let f_str = d.focal.map(|x| format!("{:.1}x", x)).unwrap_or("-".into());
            let l_str = d.local.map(|x| format!("{:.1}x", x)).unwrap_or("-".into());
            println!("  {}: focal={} local={}", gene, f_str, l_str);
        }
        println!();
    }

    let has_itd = itds.genes.values().any(|v| !v.is_empty());
    if has_itd {
        println!("=== ITDs ===");
        for (gene, entries) in &itds.genes {
            for e in entries {
                let ar = if e.coverage > e.merged {
                    e.merged as f64 / (e.coverage - e.merged) as f64
                } else {
                    e.merged as f64 / e.coverage as f64
                };
                println!("  {}: {}nt at pos {} ({} reads, {:.2} AR)", gene, e.length, e.position, e.merged, ar);
            }
        }
        println!();
    }

    if let Some(q) = qc
        && let Some(ra) = q.reads_aligned {
            println!("=== QC ===");
            println!("  Reads aligned: {}", ra);
            println!();
        }

    println!("--------------------------------------------------------------------------------");
}
