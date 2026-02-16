use std::collections::{HashMap, HashSet};
use std::io::Write;

use crate::bam::{ContigMapper, AlignmentCoords, fix_sam_coords};
use crate::input::{AlignmentInput, AlignmentHeader, AlignmentRecord, CigarKind};
use crate::config::PipelineConfig;
use crate::utils::bed::BedRegion;
use crate::utils::annotation::PartnerGeneIndex;
use log::info;

// Re-export output types for use by this module
pub use crate::output::{FusionBreakpoint, FusionEvent, FusionsOutput, GeneInfo};
use crate::output::types::BreakpointDirection;

// ==================== Shared FusionAccumulator ====================
// Used by both standalone call_fusions() and pipeline mode

/// Shared fusion candidate accumulator - used by both CLI and pipeline.
/// Processes records one at a time, accumulating hit reads and gene depths.
pub struct FusionAccumulator {
    target_map: HashMap<usize, Vec<BedRegion>>,  // ref_id -> targets (owned)
    pub hit_reads: HashSet<String>,
    gene_depth_sums: HashMap<String, u64>,       // gene -> total overlap bases
    min_anchor: u32,
    /// Per-gene margins from config, keyed by gene name
    margins: HashMap<String, u32>,
    /// Default margin for genes not in config
    default_margin: u32,
}

impl FusionAccumulator {
    /// Create new accumulator with targets indexed by BAM ref_id.
    /// Uses config for per-gene fusion margins.
    pub fn new(
        header: &AlignmentHeader,
        targets: &[BedRegion],
        config: &PipelineConfig,
    ) -> Self {
        let mut target_map: HashMap<usize, Vec<BedRegion>> = HashMap::new();
        let mut margins: HashMap<String, u32> = HashMap::new();

        let mapper = ContigMapper::from_refs(&header.refs);

        for t in targets {
            let bam_chrom = mapper.to_bam_name(&t.segment);
            if let Some(id) = header.refs.iter().position(|r| r == &bam_chrom) {
                target_map.entry(id).or_default().push(t.clone());
            }
            // Cache per-gene margin
            margins.insert(t.name.clone(), config.get_fusion_margin(&t.name));
        }

        Self {
            target_map,
            hit_reads: HashSet::new(),
            gene_depth_sums: HashMap::new(),
            min_anchor: config.thresholds.fusions.min_anchor,
            margins,
            default_margin: config.thresholds.fusions.default_margin,
        }
    }

    /// Process a single record - called per-record during BAM scan.
    /// Accumulates both hit reads and gene depth information.
    pub fn process(&mut self, record: &AlignmentRecord) {
        if record.ref_id < 0 {
            return;
        }
        let ref_id = record.ref_id as usize;

        let start = record.alignment_start().unwrap_or(0) as u32;
        let mut end = start;
        for &(op, len) in record.cigar_ops() {
            match op {
                CigarKind::Match
                | CigarKind::SequenceMatch
                | CigarKind::SequenceMismatch
                | CigarKind::Deletion
                | CigarKind::Skip => {
                    end += len as u32;
                }
                _ => {}
            }
        }

        if let Some(tgt_list) = self.target_map.get(&ref_id) {
            for tgt in tgt_list {
                // Accumulate depth: calculate actual overlap with gene (no margin)
                let overlap_start = std::cmp::max(tgt.start, start);
                let overlap_end = std::cmp::min(tgt.end, end);
                if overlap_end > overlap_start {
                    *self.gene_depth_sums.entry(tgt.name.clone()).or_insert(0)
                        += (overlap_end - overlap_start) as u64;
                }

                // Hit detection: check overlap with per-gene margin
                let margin = self.margins.get(&tgt.name).copied().unwrap_or(self.default_margin);
                let g_st = tgt.start.saturating_sub(margin);
                let g_en = tgt.end + margin;

                if std::cmp::min(g_en, end).saturating_sub(std::cmp::max(g_st, start)) > self.min_anchor {
                    if let Some(name) = record.name() {
                        self.hit_reads.insert(name.to_string());
                    }
                    break;
                }
            }
        }
    }

    /// Finalize and return results after all records processed.
    /// Returns (hit_reads, gene_depths) where gene_depths is average depth per gene.
    pub fn finalize(self, targets: &[BedRegion]) -> (HashSet<String>, HashMap<String, f64>) {
        // Convert accumulated depths to average depth per gene
        let gene_depths: HashMap<String, f64> = targets.iter()
            .map(|t| {
                let sum = self.gene_depth_sums.get(&t.name).copied().unwrap_or(0);
                let gene_len = (t.end - t.start) as f64;
                (t.name.clone(), if gene_len > 0.0 { sum as f64 / gene_len } else { 0.0 })
            })
            .collect();

        (self.hit_reads, gene_depths)
    }
}

// ==================== End FusionAccumulator ====================

/// Bundled optional parameters for fusion calling.
pub struct FusionCallParams<'a> {
    pub one_sided: Option<&'a HashSet<String>>,
    pub partner_index: Option<&'a PartnerGeneIndex>,
    pub gene_depths: Option<&'a HashMap<String, f64>>,
}

pub fn call_fusions(
    bam: &mut AlignmentInput,
    bam_path: &str,
    targets: &Vec<BedRegion>,
    config: &PipelineConfig,
    params: FusionCallParams<'_>,
) -> std::result::Result<FusionsOutput, Box<dyn std::error::Error>>
{
    info!("Calling fusions...");

    // Get file size for progress (use block count for consistent units with virtual position)
    let file_size = std::fs::metadata(bam_path)?.len();
    let file_blocks = (file_size >> 16) + 1; // Convert to BGZF block units

    // Create accumulator (uses shared code for both CLI and pipeline)
    let mut accumulator = FusionAccumulator::new(&bam.header, targets, config);

    // Pass 1: Scan BAM with progress reporting
    bam.seek(bam.start_pos)?;
    let mut last_percentage = 0;

    eprint!("Scanning for candidates... 0%");
    std::io::stderr().flush()?;

    loop {
        // Progress check (use block units for both)
        let curr = bam.tell() >> 16;
        let percentage = curr * 100 / file_blocks;
        if percentage > last_percentage {
            eprint!("\rScanning for candidates... {}%", percentage);
            std::io::stderr().flush()?;
            last_percentage = percentage;
        }

        match bam.read_record()? {
            Some(record) => accumulator.process(&record),
            None => break,
        }
    }
    eprintln!("\rScanning for candidates... 100%");

    // Finalize Pass 1
    let (hit_reads, gene_depths) = accumulator.finalize(targets);

    info!(
        "Found {} reads overlapping targets. Collecting paired alignments...",
        hit_reads.len()
    );

    // Pass 2: Breakpoint detection
    call_fusions_from_hits(
        bam, bam_path, targets, hit_reads, config,
        FusionCallParams { one_sided: params.one_sided, partner_index: params.partner_index, gene_depths: Some(&gene_depths) },
    )
}

pub fn call_fusions_from_hits(
    bam: &mut AlignmentInput,
    bam_path: &str,
    targets: &Vec<BedRegion>,
    hit_reads: HashSet<String>,
    config: &PipelineConfig,
    params: FusionCallParams<'_>,
) -> std::result::Result<FusionsOutput, Box<dyn std::error::Error>>
{
    let FusionCallParams { one_sided, partner_index, gene_depths } = params;
    let min_anchor = config.thresholds.fusions.min_anchor;
    let min_supporting_reads = config.thresholds.fusions.min_supporting_reads;
    let min_breakpoint_reads = config.thresholds.fusions.min_breakpoint_reads;
    let start_time = std::time::Instant::now();
    let file_size = std::fs::metadata(bam_path)?.len();
    let file_blocks = (file_size >> 16) + 1; // Convert to BGZF block units

    // Create a ContigMapper to translate BED chromosome names to BAM chromosome names
    let mapper = ContigMapper::from_refs(&bam.header.refs);

    let mut target_ref_map: HashMap<String, Vec<&BedRegion>> = HashMap::new();
    for target in targets {
        // Translate the BED file's chromosome name to the BAM's naming convention
        let bam_chrom = mapper.to_bam_name(&target.segment);
        target_ref_map
            .entry(bam_chrom)
            .or_default()
            .push(target);
    }

    // Pass 2: Collect all alignments for Hit Reads
    bam.seek(bam.start_pos)?;
    let mut hits: HashMap<String, Vec<(AlignmentCoords, String)>> = HashMap::new(); // SA, ref_name

    // Reset progress
    let mut last_percentage = 0;
    eprint!("Collecting alignments... 0%");

    let mut record;

    loop {
        // Progress check (use block units for both)
        let curr = bam.tell() >> 16;
        let percentage = curr * 100 / file_blocks;
        if percentage > last_percentage {
            eprint!("\rCollecting alignments... {}%", percentage);
            std::io::stderr().flush()?;
            last_percentage = percentage;
        }

        match bam.read_record()? {
            Some(rec) => record = rec,
            None => break,
        }

        let name = record.name().unwrap().to_string();

        if hit_reads.contains(&name) {
            if name.contains(';') {
                continue;
            } // Skip duplex

            // let ref_id = match record.reference_sequence_id() { Some(Ok(i)) => i, _ => continue };
            let ref_id = if record.ref_id >= 0 {
                record.ref_id as usize
            } else {
                continue;
            };

            // let (name_b, _) = bam.header.reference_sequences().get_full(&BString::from(target.segment.as_str()))?;
            let name_b: &String = bam.header.refs.get(ref_id).ok_or("Invalid ref id")?;

            // Check if chromosome is target
            // name_b is &String (or &BString? Custom header has String).
            let target_refs: HashSet<&String> = target_ref_map.keys().collect();
            if target_refs.contains(name_b) {
                // Mate checking logic would go here
            }

            let ref_name_b: &[u8] = name_b.as_bytes();
            let ref_name = String::from_utf8(ref_name_b.to_vec())?;

            let start = record.alignment_start().unwrap_or(0) as i32;
            let mut end = start;
            for &(op, len) in record.cigar_ops() {
                match op {
                    CigarKind::Match
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
                    | CigarKind::Deletion
                    | CigarKind::Skip => {
                        end += len as i32;
                    }
                    _ => {}
                }
            }
            let sa = fix_sam_coords(start as u32, end as u32, &record.cigar_raw(), record.flags())?;
            hits.entry(name).or_default().push((sa, ref_name));
        }
    }
    info!("Collecting alignments... 100%");

    // Fusion Processing
    #[derive(Hash, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
    struct Locus {
        chrom: String,
        bin: u32,
        gene: String,
    }

    struct FusionReadPair {
        gene0: (AlignmentCoords, String),
        gene1: (AlignmentCoords, String),
    }

    struct FusionCluster {
        read_count: usize,
        pairs: Vec<FusionReadPair>,
    }

    // Key: (g0, g1) sorted. Val: cluster of supporting reads.
    let mut fusions: HashMap<(Locus, Locus), FusionCluster> = HashMap::new();

    // Gene range lookup for validation
    let mut gene_ranges: HashMap<String, (u32, u32)> = HashMap::new();
    for t in targets {
        gene_ranges.insert(t.name.clone(), (t.start, t.end));
    }

    for (_read, alns) in hits {
        if alns.len() < 2 {
            continue;
        }

        let mut found_fusions: HashSet<(Locus, Locus)> = HashSet::new();

        let check_ambiguity = alns.len() > 2;

        for i in 0..alns.len() {
            let (h0, ref0) = &alns[i];

            // Skip if this alignment is ambiguously mapped (same query region
            // aligns to a different genomic locus elsewhere in this read)
            if check_ambiguity && is_query_ambiguous(h0, ref0, &alns, i, min_anchor, &target_ref_map, config) {
                continue;
            }

            // Identify Gene 0 Candidates
            let mut end0s: Vec<Locus> = Vec::new();
            if let Some(tgt_list) = target_ref_map.get::<String>(ref0) {
                for g in tgt_list {
                    let margin = config.get_fusion_margin(&g.name);
                    let g_st = g.start.saturating_sub(margin);
                    let g_en = g.end + margin;
                    if std::cmp::min(g_en, h0.te).saturating_sub(std::cmp::max(g_st, h0.ts)) > min_anchor {
                        end0s.push(Locus {
                            chrom: ref0.clone(),
                            bin: (g.start + (g.end - g.start) / 2) / 10000,
                            gene: g.name.clone(),
                        });
                        break; // Assign to first matching gene
                    }
                }
            }
            if end0s.is_empty() {
                let s_bin = h0.ts / 10000;
                let e_bin = h0.te / 10000;
                for b in s_bin..=e_bin + 1 {
                    end0s.push(Locus {
                        chrom: ref0.clone(),
                        bin: b,
                        gene: "".to_string(),
                    });
                }
            }

            for j in (i + 1)..alns.len() {
                let (h1, ref1) = &alns[j];

                // Skip if this alignment is ambiguously mapped
                if check_ambiguity && is_query_ambiguous(h1, ref1, &alns, j, min_anchor, &target_ref_map, config) {
                    continue;
                }

                // Identify Gene 1 Candidates
                let mut end1s: Vec<Locus> = Vec::new();
                if let Some(tgt_list) = target_ref_map.get::<String>(ref1) {
                    for g in tgt_list {
                        let margin = config.get_fusion_margin(&g.name);
                        let g_st = g.start.saturating_sub(margin);
                        let g_en = g.end + margin;
                        if std::cmp::min(g_en, h1.te).saturating_sub(std::cmp::max(g_st, h1.ts)) > min_anchor {
                            end1s.push(Locus {
                                chrom: ref1.clone(),
                                bin: (g.start + (g.end - g.start) / 2) / 10000,
                                gene: g.name.clone(),
                            });
                            break;
                        }
                    }
                }
                if end1s.is_empty() {
                    // If one end is not a gene, the other MUST be.
                    if !end0s.is_empty() && end0s[0].gene.is_empty() {
                        continue;
                    }
                    let s_bin = h1.ts / 10000;
                    let e_bin = h1.te / 10000;
                    for b in s_bin..=e_bin + 1 {
                        end1s.push(Locus {
                            chrom: ref1.clone(),
                            bin: b,
                            gene: "".to_string(),
                        });
                    }
                }

                for g0 in &end0s {
                    for g1 in &end1s {
                        if g0 == g1 {
                            continue;
                        }
                        // Skip self-fusions for configured genes (e.g., DUX4)
                        if !g0.gene.is_empty()
                            && config.skip_self_fusion(&g0.gene)
                            && !g1.gene.is_empty()
                            && config.skip_self_fusion(&g1.gene)
                        {
                            continue;
                        }

                        let g0_is_gene = !g0.gene.is_empty();
                        let g1_is_gene = !g1.gene.is_empty();

                        let allowed = if g0_is_gene && g1_is_gene {
                            true
                        } else {
                            // One sided fusion
                            if let Some(list) = one_sided {
                                // Identify the known gene and the unknown side
                                let (known_gene, unknown_chrom, unknown_pos) = if g0_is_gene {
                                    (&g0.gene, ref1, h1.ts)
                                } else {
                                    (&g1.gene, ref0, h0.ts)
                                };

                                if list.contains(known_gene) {
                                    // Check if this gene has configured partners
                                    if let Some(partners) = config.get_one_sided_partners(known_gene) {
                                        // Must land near an allowed partner
                                        if let Some(index) = partner_index {
                                            let margin = config.thresholds.fusions.default_margin;
                                            index.is_near_any_partner(partners, unknown_chrom, unknown_pos, margin)
                                        } else {
                                            // No partner index loaded - cannot validate partners
                                            false
                                        }
                                    } else {
                                        // No partners configured for this gene - allow any
                                        true
                                    }
                                } else {
                                    false
                                }
                            } else {
                                false
                            }
                        };

                        if !allowed {
                            continue;
                        }

                        let ((k0, k1), pair) = if g0 < g1 {
                            (
                                (g0.clone(), g1.clone()),
                                FusionReadPair { gene0: (h0.clone(), ref0.clone()), gene1: (h1.clone(), ref1.clone()) },
                            )
                        } else {
                            (
                                (g1.clone(), g0.clone()),
                                FusionReadPair { gene0: (h1.clone(), ref1.clone()), gene1: (h0.clone(), ref0.clone()) },
                            )
                        };

                        if found_fusions.contains(&(k0.clone(), k1.clone())) {
                            continue;
                        }

                        let entry =
                            fusions
                                .entry((k0.clone(), k1.clone()))
                                .or_insert(FusionCluster { read_count: 0, pairs: Vec::new() });
                        entry.read_count += 1;
                        entry.pairs.push(pair);
                        found_fusions.insert((k0, k1));
                    }
                }
            }
        }
    }

    // Process Fusions & Breakpoints
    let mut final_fusions: Vec<FusionEvent> = Vec::new();

    for ((g0, g1), cluster) in fusions {
        if (cluster.read_count as i64) < min_supporting_reads {
            continue;
        } // Min supporting reads threshold

        // Find Breakpoints — range-based clustering using overlap/gap intervals
        let mut clusters: Vec<BreakpointCluster> = Vec::new();

        for FusionReadPair { gene0: (p0, r0), gene1: (p1, r1) } in &cluster.pairs {
            // p0 corresponds to g0, p1 to g1

            let (g0_bpt, g0_dir, g1_bpt, g1_dir) = if p0.qe > p1.qe {
                let gb0 = if p0.is_reverse { p0.te } else { p0.ts };
                let gb1 = if p1.is_reverse { p1.ts } else { p1.te };
                let d0 = if !p0.is_reverse { BreakpointDirection::Right } else { BreakpointDirection::Left };
                let d1 = if p1.is_reverse { BreakpointDirection::Right } else { BreakpointDirection::Left };
                (gb0, d0, gb1, d1)
            } else {
                let gb1 = if p1.is_reverse { p1.te } else { p1.ts };
                let gb0 = if p0.is_reverse { p0.ts } else { p0.te };
                let d0 = if p0.is_reverse { BreakpointDirection::Right } else { BreakpointDirection::Left };
                let d1 = if !p1.is_reverse { BreakpointDirection::Right } else { BreakpointDirection::Left };
                (gb0, d0, gb1, d1)
            };

            // Compute overlap/gap in query coordinates.
            // Positive = microhomology (overextended alignment).
            // Negative = gap/insertion (alignment stopped short of breakpoint).
            let overlap_gap = std::cmp::min(p0.qe, p1.qe) as i32
                            - std::cmp::max(p0.qs, p1.qs) as i32;

            // Require each alignment's non-overlapping query portion to exceed min_anchor.
            // This filters out short or heavily overlapping alignments that don't
            // contribute meaningful breakpoint evidence.
            let shared = std::cmp::max(0, overlap_gap) as u32;
            let unique0 = (p0.qe - p0.qs).saturating_sub(shared);
            let unique1 = (p1.qe - p1.qs).saturating_sub(shared);
            if unique0 < min_anchor || unique1 < min_anchor {
                continue;
            }

            // Compute per-read breakpoint intervals.
            // The interval extends from the observed position toward the junction
            // (by gap amount) or toward the gene body (by overlap amount).
            let gap = std::cmp::max(0, -overlap_gap) as u32;
            let overlap = std::cmp::max(0, overlap_gap) as u32;

            let g0_interval = match g0_dir {
                BreakpointDirection::Left  => (g0_bpt.saturating_sub(overlap), g0_bpt + gap),
                BreakpointDirection::Right => (g0_bpt.saturating_sub(gap), g0_bpt + overlap),
            };
            let g1_interval = match g1_dir {
                BreakpointDirection::Left  => (g1_bpt.saturating_sub(overlap), g1_bpt + gap),
                BreakpointDirection::Right => (g1_bpt.saturating_sub(gap), g1_bpt + overlap),
            };

            // Cluster via interval overlap
            let mut found = false;
            for c in clusters.iter_mut() {
                if intervals_overlap(c.interval0, g0_interval, INTERVAL_MARGIN)
                    && intervals_overlap(c.interval1, g1_interval, INTERVAL_MARGIN)
                    && c.dir0 == g0_dir
                    && c.dir1 == g1_dir
                {
                    c.count += 1;
                    // Expand cluster intervals (union)
                    c.interval0.0 = std::cmp::min(c.interval0.0, g0_interval.0);
                    c.interval0.1 = std::cmp::max(c.interval0.1, g0_interval.1);
                    c.interval1.0 = std::cmp::min(c.interval1.0, g1_interval.0);
                    c.interval1.1 = std::cmp::max(c.interval1.1, g1_interval.1);
                    c.raw_positions0.push(g0_bpt);
                    c.raw_positions1.push(g1_bpt);
                    c.overlaps.push(overlap_gap);
                    found = true;
                    break;
                }
            }
            if !found {
                clusters.push(BreakpointCluster {
                    ref0: r0.clone(),
                    interval0: g0_interval,
                    dir0: g0_dir,
                    ref1: r1.clone(),
                    interval1: g1_interval,
                    dir1: g1_dir,
                    count: 1,
                    raw_positions0: vec![g0_bpt],
                    raw_positions1: vec![g1_bpt],
                    overlaps: vec![overlap_gap],
                });
            }
        }

        // Filter Clusters
        // 1. Min breakpoint reads (use one-sided threshold for one-sided fusions)
        //    Genes with bypass_breakpoint_filter (e.g., DUX4) skip this step
        //    because reads align to multiple DUX4 cassettes, splitting read support
        //    across breakpoint clusters. These rely on min_supporting_reads instead.
        // 2. Breakpoint Consistency with Locus (Bin/Gene) - applied on median positions
        let is_one_sided = g0.gene.is_empty() || g1.gene.is_empty();
        let bypass_bp_filter = config.bypass_breakpoint_filter(&g0.gene)
            || config.bypass_breakpoint_filter(&g1.gene);

        if !bypass_bp_filter {
            let breakpoint_threshold = if is_one_sided {
                config.get_one_sided_min_breakpoint_reads()
            } else {
                min_breakpoint_reads
            };
            clusters.retain(|c| (c.count as i64) >= breakpoint_threshold);
        }

        // Build FusionBreakpoints using median positions, then apply gene-boundary filter
        let mut bpts_out: Vec<FusionBreakpoint> = clusters.into_iter().map(|mut c| {
            FusionBreakpoint {
                gene0_name: g0.gene.clone(),
                gene0_chr: c.ref0,
                gene0_pos: median_u32(&mut c.raw_positions0),
                gene0_dir: c.dir0,
                gene1_name: g1.gene.clone(),
                gene1_chr: c.ref1,
                gene1_pos: median_u32(&mut c.raw_positions1),
                gene1_dir: c.dir1,
                n_reads: c.count,
                overlap_gap: Some(median_i32(&mut c.overlaps)),
            }
        }).collect();

        // Gene-boundary filter (always applied, using median positions)
        bpts_out.retain(|bp| {
            let g0_ok = if g0.gene.is_empty() {
                bp.gene0_pos / 10000 == g0.bin
            } else {
                let margin = config.get_fusion_margin(&g0.gene);
                if let Some((st, en)) = gene_ranges.get(&g0.gene) {
                    bp.gene0_pos > st.saturating_sub(margin) && bp.gene0_pos < en + margin
                } else {
                    false
                }
            };

            let g1_ok = if g1.gene.is_empty() {
                bp.gene1_pos / 10000 == g1.bin
            } else {
                let margin = config.get_fusion_margin(&g1.gene);
                if let Some((st, en)) = gene_ranges.get(&g1.gene) {
                    bp.gene1_pos > st.saturating_sub(margin) && bp.gene1_pos < en + margin
                } else {
                    false
                }
            };

            g0_ok && g1_ok
        });

        if !bpts_out.is_empty() {
            // Apply min_fraction filter for one-sided fusions
            if is_one_sided
                && let Some(min_frac) = config.get_one_sided_min_fraction() {
                    // Get the known gene (the one that's not empty)
                    let known_gene = if g0.gene.is_empty() { &g1.gene } else { &g0.gene };
                    if let Some(depths) = gene_depths
                        && let Some(&avg_depth) = depths.get(known_gene)
                            && avg_depth > 0.0 {
                                let fraction = cluster.read_count as f64 / avg_depth;
                                if fraction < min_frac {
                                    info!(
                                        "Filtering one-sided fusion {}-{}: fraction {:.3} < min_fraction {:.3}",
                                        g0.gene, g1.gene, fraction, min_frac
                                    );
                                    continue; // Skip this fusion
                                }
                            }
                }

            final_fusions.push(FusionEvent {
                gene1: GeneInfo {
                    name: g0.gene.clone(),
                    chr: bpts_out[0].gene0_chr.clone(),
                    pos: bpts_out[0].gene0_pos,
                },
                gene2: GeneInfo {
                    name: g1.gene.clone(),
                    chr: bpts_out[0].gene1_chr.clone(),
                    pos: bpts_out[0].gene1_pos,
                },
                supporting_reads: cluster.read_count,
                repetitive_reads: 0,
                breakpoints: bpts_out,
            });
        }
    }

    // Separate spike-in control fusions from real fusions
    let (spike_in_fusions, real_fusions): (Vec<_>, Vec<_>) = final_fusions
        .into_iter()
        .partition(|f| config.is_spike_in_fusion(&f.gene1.name, &f.gene2.name));

    // Build output structure
    let fusions_output = FusionsOutput {
        fusions: real_fusions.clone(),
        spike_in: if spike_in_fusions.is_empty() {
            None
        } else {
            Some(spike_in_fusions.clone())
        },
    };

    let duration = start_time.elapsed();
    info!("Done.");
    info!("--------------------------------------------------");
    info!("Summary:");
    info!("  Time elapsed: {:.2?}", duration);
    info!("  Total fusions found: {}", real_fusions.len());
    if !spike_in_fusions.is_empty() {
        info!("  Spike-in controls detected: {}", spike_in_fusions.len());
    }
    info!("--------------------------------------------------");

    Ok(fusions_output)
}

// ==================== Breakpoint Clustering Helpers ====================

/// Small error margin for interval overlap comparisons (in bp).
const INTERVAL_MARGIN: u32 = 10;

/// A cluster of breakpoints defined by interval overlap rather than fixed position tolerance.
struct BreakpointCluster {
    ref0: String,
    interval0: (u32, u32),
    dir0: BreakpointDirection,
    ref1: String,
    interval1: (u32, u32),
    dir1: BreakpointDirection,
    count: usize,
    raw_positions0: Vec<u32>,
    raw_positions1: Vec<u32>,
    overlaps: Vec<i32>,
}

/// Check if two intervals overlap, with a small margin of error.
fn intervals_overlap(a: (u32, u32), b: (u32, u32), margin: u32) -> bool {
    a.0 <= b.1 + margin && b.0 <= a.1 + margin
}

fn median_u32(vals: &mut [u32]) -> u32 {
    vals.sort_unstable();
    vals[vals.len() / 2]
}

fn median_i32(vals: &mut [i32]) -> i32 {
    vals.sort_unstable();
    vals[vals.len() / 2]
}

/// Check if an alignment is ambiguously mapped — i.e., another alignment from the
/// same read covers the same query region but maps to a different genomic locus.
/// This detects repetitive elements data-drivenly without needing an external BED file.
/// Alignments within the same gene target region (e.g., DUX4 cassettes) are not
/// considered competing even if they don't directly overlap in reference coordinates.
fn is_query_ambiguous(
    target: &AlignmentCoords,
    target_ref: &str,
    all_alns: &[(AlignmentCoords, String)],
    exclude_idx: usize,
    min_anchor: u32,
    target_ref_map: &HashMap<String, Vec<&BedRegion>>,
    config: &PipelineConfig,
) -> bool {
    let target_qlen = target.qe - target.qs;
    for (j, (other, other_ref)) in all_alns.iter().enumerate() {
        if j == exclude_idx {
            continue;
        }
        // Skip if same locus (same chrom AND overlapping reference coords) — not competing
        let same_locus = other_ref == target_ref
            && std::cmp::min(target.te, other.te) > std::cmp::max(target.ts, other.ts);
        if same_locus {
            continue;
        }
        // Skip if both alignments fall within the same gene target region.
        // This handles genes like DUX4 with multiple cassettes spread across
        // a region — multi-mapping within one gene is expected, not ambiguous.
        if other_ref == target_ref
            && let Some(genes) = target_ref_map.get(target_ref) {
                let same_gene = genes.iter().any(|g| {
                    let margin = config.get_fusion_margin(&g.name);
                    let g_st = g.start.saturating_sub(margin);
                    let g_en = g.end + margin;
                    target.ts < g_en && target.te > g_st
                        && other.ts < g_en && other.te > g_st
                });
                if same_gene {
                    continue;
                }
            }
        // Compute query-coordinate overlap
        let q_overlap = std::cmp::min(target.qe, other.qe)
            .saturating_sub(std::cmp::max(target.qs, other.qs));
        if q_overlap == 0 {
            continue;
        }
        let other_qlen = other.qe - other.qs;
        let shorter = std::cmp::min(target_qlen, other_qlen);
        if q_overlap > shorter / 2 && q_overlap > min_anchor {
            return true;
        }
    }
    false
}
