use indexmap::IndexMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::bam::{ContigMapper, AlignmentCoords, fix_sam_coords};
use crate::input::AlignmentInput;
use crate::config::{Contig, PipelineConfig, ReferenceConfig};
use crate::output::{CnvOutput, CnvGeneResult, FocalNbFit, KaryotypeOutput, StructuralVariant};
use crate::var::focal_nb::{fit_nb, test_gene, GeneReadCount, MIN_NB_BASELINE};
use crate::utils::bed::BedRegion;
use log::{info, warn, debug};

// get_segment using reference config for centromere/PAR coordinates
fn get_segment(chrom: &str, start: u32, end: u32, ref_config: &ReferenceConfig) -> Option<String> {
    // Extract chromosome number/letter (remove "chr" prefix)
    let chr_num = chrom.strip_prefix("chr").unwrap_or(chrom);

    if let Some(&(cen_start, cen_end)) = ref_config.centromeres.get(chrom) {
        if end < cen_start {
            // p-arm (before centromere)
            Some(format!("{}p", chr_num))
        } else if start > cen_end {
            // q-arm (after centromere)
            Some(format!("{}q", chr_num))
        } else {
            // Overlaps centromere - skip
            None
        }
    } else if let Some(&(par_start, par_end)) = ref_config.par1.get(chrom) {
        if start < par_end && end > par_start {
            // Overlaps PAR1
            if start >= par_start && end <= par_end {
                // Completely contained in PAR1
                Some("PAR1".to_string())
            } else {
                None
            }
        } else if let Some(&(par2_start, par2_end)) = ref_config.par2.get(chrom) {
            if end > par2_start && start < par2_end {
                // Overlaps PAR2 - skip
                None
            } else {
                Some(chr_num.to_string())
            }
        } else {
            Some(chr_num.to_string())
        }
    } else {
        // Chromosomes without centromeres in the dict (13, 14, 15, 21, 22)
        Some(chr_num.to_string())
    }
}

struct DelDupResult {
    deletions: HashMap<(u32, u32), usize>,
    duplications: HashMap<(u32, u32), usize>,
}

fn internal_del_dup(
    bam: &mut AlignmentInput,
    _gene: &str,
    target: &BedRegion,
    config: &PipelineConfig,
) -> Result<DelDupResult, Box<dyn std::error::Error>>
{
    let gene_margin = config.thresholds.cnv.gene_margin;
    let breakpoint_margin = config.thresholds.cnv.breakpoint_margin;
    let min_anchor = config.thresholds.cnv.min_anchor;

    // use IndexMap to preserve insertion order
    let mut hits: IndexMap<String, Vec<AlignmentCoords>> = IndexMap::new();

    let s = target.start.saturating_sub(gene_margin);
    let e = target.end + gene_margin;

    // resolve chromosome name
    let chrom = bam.contig_mapper.to_bam_name(&target.segment);
    let region_str = format!("{}:{}-{}", chrom, s, e);

    let query = bam.query(&region_str)?;

    for result in query {
        let record = result?;

        let name = record.name().unwrap().to_string();
        if name.contains(';') {
            continue;
        } // filter out duplex reads

        let start = match record.alignment_start() {
            Some(p) => p as i32,
            _ => continue,
        };
        let end = start + record.alignment_span() as i32;

        // fix SAM coords to match canonical query
        let sa = fix_sam_coords(start as u32, end as u32, &record.cigar_raw(), record.flags())?;
        hits.entry(name).or_default().push(sa);
    }

    let mut dels: HashMap<(u32, u32), usize> = HashMap::new();
    let mut dups: HashMap<(u32, u32), usize> = HashMap::new();

    for (_name, alns) in hits {
        if alns.len() < 2 {
            continue;
        }

        for (i, h0) in alns.iter().enumerate() {
            for h1 in alns.iter().skip(i + 1) {

                // Cast to i64 to avoid underflow when alignments don't overlap
                let q_overlap =
                    std::cmp::min(h0.qe, h1.qe) as i64 - std::cmp::max(h0.qs, h1.qs) as i64;
                if q_overlap > 100 {
                    continue;
                }
                if h0.is_reverse != h1.is_reverse {
                    continue;
                }

                // min_anchor check
                if (h0.te - h0.ts) > min_anchor && (h1.te - h1.ts) > min_anchor {
                    /* tandem duplications:
                             -------->
                            | aln 1  |
                      -------===--------
                        -------->
                       | aln 0  |
                      -------===--------
                    */
                    if (!h1.is_reverse && h0.qe < h1.qe && h0.te.saturating_sub(h1.ts) > 1000)
                        || (h0.is_reverse && h1.qe < h0.qe && h0.te.saturating_sub(h1.ts) > 1000)
                    {
                        let cand_s = h1.ts;
                        let cand_e = h0.te;

                        // Find the closest existing breakpoint within margin (for deterministic behavior)
                        let mut best_match: Option<(u32, u32)> = None;
                        let mut best_dist = u32::MAX;
                        for (k, _v) in dups.iter() {
                            let dist_s = k.0.abs_diff(cand_s);
                            let dist_e = k.1.abs_diff(cand_e);
                            if dist_s < breakpoint_margin && dist_e < breakpoint_margin {
                                let total_dist = dist_s + dist_e;
                                if total_dist < best_dist {
                                    best_dist = total_dist;
                                    best_match = Some(*k);
                                }
                            }
                        }

                        if let Some(k) = best_match {
                            *dups.get_mut(&k).unwrap() += 1;
                        } else {
                            dups.insert((cand_s, cand_e), 1);
                        }
                    } else if (h1.is_reverse && h0.qe < h1.qe && h1.te.saturating_sub(h0.ts) > 1000)
                        || (!h0.is_reverse && h1.qe < h0.qe && h1.te.saturating_sub(h0.ts) > 1000)
                    {
                        let cand_s = h0.ts;
                        let cand_e = h1.te;

                        // Find the closest existing breakpoint within margin
                        let mut best_match: Option<(u32, u32)> = None;
                        let mut best_dist = u32::MAX;
                        for (k, _v) in dups.iter() {
                            let dist_s = k.0.abs_diff(cand_s);
                            let dist_e = k.1.abs_diff(cand_e);
                            if dist_s < breakpoint_margin && dist_e < breakpoint_margin {
                                let total_dist = dist_s + dist_e;
                                if total_dist < best_dist {
                                    best_dist = total_dist;
                                    best_match = Some(*k);
                                }
                            }
                        }

                        if let Some(k) = best_match {
                            *dups.get_mut(&k).unwrap() += 1;
                        } else {
                            dups.insert((cand_s, cand_e), 1);
                        }
                    }
                    // Dels
                    else if ((!h0.is_reverse && h0.qe < h1.qe)
                        || (h0.is_reverse && h0.qe > h1.qe))
                        && (h0.te < h1.ts)
                    {
                        let cand_s = h0.te;
                        let cand_e = h1.ts;

                        // Find the closest existing breakpoint within margin
                        let mut best_match: Option<(u32, u32)> = None;
                        let mut best_dist = u32::MAX;
                        for (k, _v) in dels.iter() {
                            let dist_s = k.0.abs_diff(cand_s);
                            let dist_e = k.1.abs_diff(cand_e);
                            if dist_s < breakpoint_margin && dist_e < breakpoint_margin {
                                let total_dist = dist_s + dist_e;
                                if total_dist < best_dist {
                                    best_dist = total_dist;
                                    best_match = Some(*k);
                                }
                            }
                        }

                        if let Some(k) = best_match {
                            *dels.get_mut(&k).unwrap() += 1;
                        } else {
                            dels.insert((cand_s, cand_e), 1);
                        }
                    }
                }
            }
        }
    }
    Ok(DelDupResult { deletions: dels, duplications: dups })
}

fn calculate_focal_depths(
    bam: &mut AlignmentInput,
    targets: &[BedRegion],
) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>>
{
    let mut depths = HashMap::new();
    for (ti, t) in targets.iter().enumerate() {
        // Convert chromosome name to BAM format (e.g., chr9 -> NC_060933.1)
        let chrom = bam.contig_mapper.to_bam_name(&t.segment);
        let region = format!("{}:{}-{}", chrom, t.start, t.end);
        eprint!("\r  Querying {}/{}: {} ({})...", ti + 1, targets.len(), t.name, region);
        std::io::Write::flush(&mut std::io::stderr()).ok();

        // Resilient: a region on a contig absent from the BAM header errors here.
        // Warn and skip that gene instead of aborting the whole CNV stage.
        let query = match bam.query(&region) {
            Ok(q) => q,
            Err(e) => {
                warn!("CNV focal depth: skipping {} ({}) — region query failed: {}", t.name, region, e);
                continue;
            }
        };
        let mut total_bases = 0;

        for result in query {
            let record = result?;
            let a_start = match record.alignment_start() {
                Some(p) => p,
                _ => continue,
            };
            let a_len = record.alignment_span();
            let a_end = a_start + a_len;

            // Overlap
            let o_start = std::cmp::max(a_start, t.start as usize);
            let o_end = std::cmp::min(a_end, t.end as usize);

            if o_end > o_start {
                total_bases += o_end - o_start;
            }
        }

        let len = (t.end - t.start) as usize;
        let depth = if len > 0 {
            total_bases as f64 / len as f64
        } else {
            0.0
        };
        depths.insert(t.name.clone(), depth);
    }
    eprintln!("\r  Focal depths done ({} targets).          ", targets.len());
    Ok(depths)
}

/// Returns (raw_depths, median_coverage) - raw depths are NOT normalized
fn calculate_local_depths_raw(
    covg_path: &str,
    targets: &[BedRegion],
    contigs: &[Contig],
) -> Result<(HashMap<String, f64>, f64), Box<dyn std::error::Error>> {
    // Read coverage TSV: chrom start end value
    let file = File::open(covg_path).map_err(|e| {
        std::io::Error::other(format!("Error opening coverage file {}: {}", covg_path, e))
    })?;
    let reader = BufReader::new(file);

    let mut bins: Vec<(String, usize, usize, f64)> = Vec::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with("chromosome") { continue; }
        let p: Vec<&str> = l.split('\t').collect();
        if p.len() < 4 {
            continue;
        }
        // Unified TSV: column 6 is `n_reads_gc_adj` (GC-corrected coverage).
        // Prefer it when present and finite; fall back to column 4 (raw)
        // otherwise, which also covers legacy 5-column files and the
        // no-karyotype standalone-coverage path.
        let raw: f64 = p[3].parse().unwrap_or(0.0);
        let val: f64 = if p.len() >= 6 {
            p[5].parse::<f64>().ok().filter(|v| v.is_finite()).unwrap_or(raw)
        } else {
            raw
        };
        let start: usize = p[1].parse().unwrap_or(0);
        let end: usize = p[2].parse().unwrap_or(0);
        bins.push((p[0].to_string(), start, end, val));
    }

    let mut local_depths = HashMap::new();

    // Calculate median coverage (fallback baseline)
    let mut all_vals: Vec<f64> = bins.iter().map(|x| x.3).filter(|v| *v > 0.0).collect();
    all_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if !all_vals.is_empty() {
        all_vals[all_vals.len() / 2]
    } else {
        1.0
    };
    debug!(
        "Median local coverage (from {} bins): {:.2}",
        all_vals.len(),
        median
    );

    let chr_mapper = ContigMapper::from_contigs(contigs);
    for t in targets {
        let mut sum_val = 0.0;
        let mut count = 0;
        let t_chr = chr_mapper.to_chr_name(&t.segment);

        for (chr, s, e, v) in &bins {
            if chr == &t_chr && *s < t.end as usize && *e > t.start as usize {
                sum_val += v;
                count += 1;
            }
        }
        let avg_raw = if count > 0 {
            sum_val / count as f64
        } else {
            0.0
        };
        // Return raw depth, not normalized
        local_depths.insert(t.name.clone(), avg_raw);
    }

    Ok((local_depths, median))
}

/// Blast ratios below this make the focal note ask to check other blast
/// estimates, and the report shows the note in red.
pub const LOW_FOCAL_BLAST_RATIO: f64 = 0.5;

/// Bundled optional parameters for CNV calling.
pub struct CnvCallParams<'a> {
    pub coverage_file: Option<&'a str>,
    pub blast_ratio: Option<f64>,
    pub karyotype: Option<&'a KaryotypeOutput>,
    pub ref_config: &'a ReferenceConfig,
    pub precomputed_focal_depths: Option<HashMap<String, f64>>,
    /// Per-gene target read counts (from the pipeline's QC pass or
    /// `{prefix}.target_counts.tsv`); None falls back to coverage focal CN.
    pub target_read_counts: Option<HashMap<String, GeneReadCount>>,
}

pub fn call_cnvs(
    bam: &mut AlignmentInput,
    targets: &[BedRegion],
    config: &PipelineConfig,
    params: CnvCallParams<'_>,
) -> Result<CnvOutput, Box<dyn std::error::Error>>
{
    let CnvCallParams { coverage_file, blast_ratio, karyotype, ref_config, precomputed_focal_depths, target_read_counts } = params;
    let mut result: HashMap<String, CnvGeneResult> = HashMap::new();
    let chr_mapper = ContigMapper::from_contigs(&ref_config.contigs);

    // Map targets by name
    let mut t_map: HashMap<&str, &BedRegion> = HashMap::new();
    for t in targets {
        t_map.insert(t.name.as_str(), t);
    }

    // Extract diploid segments and blast ratio from karyotype if available
    let (diploid_segments, bin_cov_2n, karyo_blast_ratio) = if let Some(karyo) = karyotype {
        let mut diploid_segs: Vec<String> = Vec::new();
        let mut diploid_covs: Vec<f64> = Vec::new();

        for (seg, &cn) in &karyo.karyotype {
            if cn == 2 {
                diploid_segs.push(seg.clone());
                if let Some(&med) = karyo.medians.get(seg) {
                    diploid_covs.push(med);
                }
            }
        }

        let cov_2n = if !diploid_covs.is_empty() {
            let sum: f64 = diploid_covs.iter().sum();
            Some(sum / diploid_covs.len() as f64)
        } else {
            None
        };

        // Cap blast_ratio at 1.0 as sanity check
        let br = karyo.blast_ratio.map(|v| if v > 1.0 { 1.0 } else { v });

        info!(
            "Diploid segments from karyotype: {} (mean cov: {:?}, blast_ratio: {:?})",
            diploid_segs.len(),
            cov_2n,
            br
        );
        (diploid_segs, cov_2n, br)
    } else {
        (Vec::new(), None, None)
    };

    // Use passed blast_ratio, or fall back to karyotype's blast_ratio
    let effective_blast_ratio = blast_ratio.or(karyo_blast_ratio).map(|v| v.min(1.0));
    if let Some(br) = effective_blast_ratio {
        info!("Using blast ratio: {:.4}", br);
    }
    // Focal CN assumes all cells are blasts when no blast ratio is available
    let focal_br = effective_blast_ratio.filter(|&b| b > 1e-6).unwrap_or(1.0);
    // Blast percentage assumed: a note when a blast ratio is available, a
    // warning when 100% is assumed for lack of one
    let mut notes: Vec<String> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();
    if effective_blast_ratio.is_some() {
        let mut note = format!(
            "Focal copy number calculation assumes a blast percentage of {:.0}%.",
            focal_br * 100.0
        );
        if focal_br < LOW_FOCAL_BLAST_RATIO {
            note.push_str(
                " Low blast percentage: verify against an independent blast estimate (e.g. flow cytometry) \
                 and rerun with --blast-ratio if needed.",
            );
        }
        notes.push(note);
    } else {
        warnings.push(
            "Focal copy number calculation assumes a blast percentage of 100% (no blast ratio estimated)."
                .to_string(),
        );
    }

    // Determine which genes are in diploid segments (for focal depth filtering)
    let mut diploid_genes: Vec<String> = Vec::new();
    for t in targets {
        let chr_name = chr_mapper.to_chr_name(&t.segment);
        if let Some(seg) = get_segment(&chr_name, t.start, t.end, ref_config)
            && diploid_segments.contains(&seg) {
                diploid_genes.push(t.name.clone());
            }
    }
    info!("Genes in diploid segments: {}", diploid_genes.len());

    // 1. Calculate Focal (Gene) Depths from BAM (or use precomputed)
    let focal_depths: HashMap<String, f64> = if let Some(fd) = precomputed_focal_depths {
        info!("Using precomputed focal depths ({} targets).", fd.len());
        fd
    } else {
        info!("Calculating focal depths from BAM...");
        calculate_focal_depths(bam, targets)?
    };

    // Calculate median of focal depths from DIPLOID genes only
    // Use exclude_from_baseline from config
    let exclude_list = &config.genes.cnv.exclude_from_baseline;

    // Baseline genes: not excluded, and in a diploid segment (or, without
    // karyotype info, not on a sex chromosome)
    let is_baseline = |name: &str| -> bool {
        if exclude_list.iter().any(|e| name.contains(e.as_str())) {
            return false;
        }
        if !diploid_segments.is_empty() {
            diploid_genes.iter().any(|g| g == name)
        } else if let Some(t) = t_map.get(name) {
            let chr_name = chr_mapper.to_chr_name(&t.segment);
            chr_name != "chrX" && chr_name != "chrY"
        } else {
            true
        }
    };

    let mut filtered_depths: Vec<f64> = Vec::new();
    for (name, depth) in &focal_depths {
        if is_baseline(name) {
            filtered_depths.push(*depth);
        }
    }

    filtered_depths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_focal = if !filtered_depths.is_empty() {
        filtered_depths[filtered_depths.len() / 2]
    } else {
        // Fallback to all if filtration removes everything
        warn!("No diploid genes found for focal baseline, using all genes");
        let mut all: Vec<f64> = focal_depths.values().cloned().collect();
        all.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if !all.is_empty() {
            all[all.len() / 2]
        } else {
            1.0
        }
    };
    info!(
        "Median focal depth (baseline from {} genes): {:.2}",
        filtered_depths.len(),
        median_focal
    );

    // 2. Calculate Local (Bin) Depths - get raw depths
    let (local_depths_raw, fallback_median) = if let Some(cov_path) = coverage_file {
        info!("Reading local depths from {}...", cov_path);
        match calculate_local_depths_raw(cov_path, targets, &ref_config.contigs) {
            Ok((d, m)) => (d, m),
            Err(e) => {
                warn!("Could not read coverage file: {}", e);
                (HashMap::new(), 1.0)
            }
        }
    } else {
        (HashMap::new(), 1.0)
    };

    // Use bin_cov_2n from karyotype if available, otherwise fallback to median
    let local_baseline = bin_cov_2n.unwrap_or(fallback_median);
    info!(
        "Local depth baseline (2n): {:.2} (source: {})",
        local_baseline,
        if bin_cov_2n.is_some() {
            "karyotype diploid segments"
        } else {
            "coverage median"
        }
    );

    // Focal CN from read counts: NB fit of the baseline genes' counts
    let cnv_genes = &config.genes.cnv.focal_genes;
    let mut nb_baseline: Vec<(&str, GeneReadCount)> = target_read_counts
        .as_ref()
        .map(|c| c.iter().filter(|(g, _)| is_baseline(g)).map(|(g, &n)| (g.as_str(), n)).collect())
        .unwrap_or_default();
    // Fixed order, so the fit does not depend on HashMap iteration order
    nb_baseline.sort_by(|a, b| a.0.cmp(b.0));
    let nb_fit = match &target_read_counts {
        None => {
            warn!("No target read counts; focal CN from coverage.");
            warnings.push(
                "Negative binomial test not run (no target read counts). Fall back to gene coverage CN estimates.".to_string(),
            );
            None
        }
        Some(_) if nb_baseline.len() < MIN_NB_BASELINE => {
            warn!("Only {} diploid baseline genes for the NB fit; focal CN from coverage.", nb_baseline.len());
            warnings.push(format!(
                "Negative binomial test not run ({} diploid baseline genes). Fall back to gene coverage CN estimates.",
                nb_baseline.len()
            ));
            None
        }
        Some(_) => {
            let counts: Vec<GeneReadCount> = nb_baseline.iter().map(|(_, c)| *c).collect();
            let fit = fit_nb(&counts);
            match fit {
                Some(f) => info!(
                    "Focal NB fit on {} baseline genes: {:.2} reads/10kb, offset {:.0} bp, alpha {:.4}",
                    counts.len(), f.r * 1e4, f.l, f.alpha
                ),
                None => {
                    warn!("Focal NB fit failed; focal CN from coverage.");
                    warnings.push(format!(
                        "Negative binomial test not run (fit failed on {} diploid baseline genes). Fall back to gene coverage CN estimates.",
                        counts.len()
                    ));
                }
            }
            fit
        }
    };

    for g in cnv_genes {
        let mut gene_result = CnvGeneResult::default();

        // Coverage focal CN: cn = (depth - med) / (med * 0.5 * br) + 2
        let cov = focal_depths
            .get(g)
            .map(|&depth| (2.0 + 2.0 * (depth / median_focal - 1.0) / focal_br).max(0.0));
        if let Some(cn) = cov {
            gene_result.focal_coverage = Some(cn);
        }

        let count = target_read_counts.as_ref().and_then(|c| c.get(g).copied());
        match (nb_fit, count) {
            (Some(fit), Some(count)) => {
                // A focal gene in the baseline (e.g. INO80) is tested against a
                // fit without it
                let fit = if nb_baseline.iter().any(|(b, _)| *b == g.as_str()) {
                    let others: Vec<GeneReadCount> =
                        nb_baseline.iter().filter(|(b, _)| *b != g.as_str()).map(|(_, c)| *c).collect();
                    fit_nb(&others).unwrap_or(fit)
                } else {
                    fit
                };
                let t = test_gene(&fit, count);
                // cn = 2 + 2 (reads / mu - 1) / br
                let cn = (2.0 + 2.0 * (t.reads as f64 / t.mu - 1.0) / focal_br).max(0.0);
                gene_result.focal = Some(cn);
                gene_result.focal_p_2n = Some(t.p_2n);
                gene_result.focal_reads = Some(t.reads);
                gene_result.focal_expected_reads = Some(t.mu);
                if cn.round() as i64 != 2 && t.p_2n >= 0.05 {
                    gene_result.focal_flag = Some(
                        "Low confidence in this copy number change.".to_string(),
                    );
                }
            }
            (Some(_), None) => {
                warn!("No read count for focal gene {}", g);
            }
            (None, _) => {
                if let Some(cn) = cov {
                    gene_result.focal = Some(cn);
                }
            }
        }

        // Local CN - using raw bin depth and karyotype-based baseline
        if let Some(bin_depth) = local_depths_raw.get(g) {
            let cn = if let Some(br) = effective_blast_ratio {
                // Blast ratio formula:
                // expected_3n = br*3*covg_2n/2 + (1-br)*covg_2n
                // delta = expected_3n - covg_2n
                // cn = (bin_depth - covg_2n) / delta + 2
                let expected_3n = br * 3.0 * local_baseline / 2.0 + (1.0 - br) * local_baseline;
                let delta = expected_3n - local_baseline;
                if delta.abs() > 1e-6 {
                    (bin_depth - local_baseline) / delta + 2.0
                } else {
                    bin_depth / local_baseline * 2.0
                }
            } else {
                // Simple 2n normalization
                bin_depth / local_baseline * 2.0
            };
            // Clamp negative CN values to 0
            let cn = if cn < 0.0 { 0.0 } else { cn };
            gene_result.local = Some(cn);
        }

        result.insert(g.clone(), gene_result);
    }

    // Build list of genes to check for internal deletions/duplications
    let mut del_dup_genes: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for g in &config.genes.cnv.deletion_genes {
        del_dup_genes.insert(g.as_str());
    }
    for g in &config.genes.cnv.duplication_genes {
        del_dup_genes.insert(g.as_str());
    }
    let chr_mapper = ContigMapper::from_contigs(&ref_config.contigs);

    for g in del_dup_genes {
        if let Some(t) = t_map.get(g) {
            info!("Calling internal deletions/duplications for {}...", g);
            // Convert NC_ accession to chr name
            let chr_name = chr_mapper.to_chr_name(&t.segment);
            match internal_del_dup(bam, g, t, config) {
                Ok(DelDupResult { deletions: dels, duplications: dups }) => {
                    let min_reads = config.thresholds.cnv.min_supporting_reads;
                    let mut del_list: Vec<StructuralVariant> = Vec::new();
                    for ((s, e), count) in dels {
                        if count >= min_reads {
                            // Convert from 1-based (internal) to 0-based (output)
                            let s0 = s - 1;
                            let e0 = e - 1;
                            debug!("{} deletion {}-{}: {} reads", g, s0, e0, count);
                            del_list.push(StructuralVariant {
                                chrom: Some(chr_name.clone()),
                                start: s0 as i64,
                                end: e0 as i64,
                                reads: count,
                            });
                        }
                    }

                    let mut dup_list: Vec<StructuralVariant> = Vec::new();
                    for ((s, e), count) in dups {
                        if count >= min_reads {
                            // Convert from 1-based (internal) to 0-based (output)
                            let s0 = s - 1;
                            let e0 = e - 1;
                            debug!("{} duplication {}-{}: {} reads", g, s0, e0, count);
                            dup_list.push(StructuralVariant {
                                chrom: Some(chr_name.clone()),
                                start: s0 as i64,
                                end: e0 as i64,
                                reads: count,
                            });
                        }
                    }

                    // Get or create gene entry
                    let gene_entry = result.entry(g.to_string()).or_default();
                    if !del_list.is_empty() {
                        gene_entry.deletions = Some(del_list);
                    }
                    if !dup_list.is_empty() {
                        gene_entry.duplications = Some(dup_list);
                    }
                }
                Err(e) => {
                    warn!("Error calling del/dup for {}: {}", g, e);
                }
            }
        }
    }

    Ok(CnvOutput {
        genes: result,
        focal_method: Some(if nb_fit.is_some() { "read_count_nb" } else { "coverage" }.to_string()),
        focal_blast_ratio: Some(focal_br),
        focal_baseline_genes: Some(if nb_fit.is_some() { nb_baseline.len() } else { filtered_depths.len() }),
        focal_nb_fit: nb_fit.map(|f| FocalNbFit { reads_per_base: f.r, length_offset: f.l, alpha: f.alpha }),
        notes: if notes.is_empty() { None } else { Some(notes) },
        warnings: if warnings.is_empty() { None } else { Some(warnings) },
    })
}
