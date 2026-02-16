//! Fusion breakpoint consensus sequence reconstruction.
//!
//! Given detected fusion breakpoints, re-scans the BAM using indexed queries
//! to find supporting split reads, then builds a majority-vote consensus
//! spanning each breakpoint junction.

use std::collections::HashMap;
use std::io::Write;

use log::info;

use crate::bam::{fix_sam_coords, AlignmentCoords};
use crate::input::AlignmentInput;
use crate::output::types::{
    BreakpointConsensusOutput, BreakpointDirection, FusionBreakpoint, FusionBreakpointConsensus, FusionsOutput,
};
use crate::utils::align::{needleman_wunsch, AlignOp, AlignScoring};

// ============================================================================
// Configuration
// ============================================================================

/// Configuration for breakpoint consensus calling.
pub struct BreakpointConsensusConfig {
    /// Number of flanking bases on each side of the breakpoint (default: 200)
    pub flank_size: usize,
    /// Minimum number of reads to attempt consensus (default: 3)
    pub min_reads: usize,
    /// Minimum per-position coverage to call a consensus base (default: 3)
    pub min_coverage: u32,
    /// Margin around breakpoint position for BAM queries in bp (default: 500)
    pub query_margin: u32,
    /// Tolerance for matching read alignments to breakpoint position (default: 30)
    pub cluster_tolerance: u32,
    /// Minimum fraction of reads with an insertion for it to be called (default: 0.5)
    pub insertion_fraction: f64,
}

impl Default for BreakpointConsensusConfig {
    fn default() -> Self {
        Self {
            flank_size: 200,
            min_reads: 3,
            min_coverage: 3,
            query_margin: 500,
            cluster_tolerance: 30,
            insertion_fraction: 0.5,
        }
    }
}

// ============================================================================
// Internal data structures
// ============================================================================

/// A single read's extracted breakpoint-spanning data.
struct BreakpointRead {
    /// Query span for the "left" alignment (earlier in query coords)
    left_span: (u32, u32),
    /// Query span for the "right" alignment (later in query coords)
    right_span: (u32, u32),
    /// Whether the left alignment corresponds to gene0
    left_is_gene0: bool,
    /// Full read sequence (original strand, ASCII A/C/G/T)
    sequence: Vec<u8>,
    /// Junction sequence between left_span.1 and right_span.0
    junction_seq: Vec<u8>,
}

/// Result of consensus building for one breakpoint.
struct ConsensusResult {
    /// Full consensus sequence
    consensus_seq: Vec<u8>,
    /// Index where the left gene flank ends
    left_end: usize,
    /// Index where the right gene flank begins
    right_start: usize,
    /// Inserted bases at junction (may be empty)
    inserted_bases: Vec<u8>,
    /// Whether the left side is gene0
    left_is_gene0: bool,
    /// Number of reads used
    n_reads: usize,
    /// Per-position coverage
    coverage: Vec<u32>,
}

// ============================================================================
// Public API
// ============================================================================

/// Build breakpoint consensus for all breakpoints in a FusionsOutput.
///
/// Returns structured output and optionally writes a FASTA file.
pub fn call_breakpoint_consensus(
    bam: &mut AlignmentInput,
    fusions_output: &FusionsOutput,
    config: &BreakpointConsensusConfig,
    out_prefix: Option<&str>,
) -> Result<BreakpointConsensusOutput, Box<dyn std::error::Error>> {
    let mut results: Vec<FusionBreakpointConsensus> = Vec::new();
    let mut fasta_entries: Vec<(String, String)> = Vec::new();

    for fusion in &fusions_output.fusions {
        for bp in &fusion.breakpoints {
            info!(
                "Building consensus for {} {}:{} -- {} {}:{}",
                bp.gene0_name, bp.gene0_chr, bp.gene0_pos, bp.gene1_name, bp.gene1_chr,
                bp.gene1_pos
            );

            match build_consensus_for_breakpoint(bam, bp, config) {
                Ok(Some(cons)) => {
                    let seq_str =
                        String::from_utf8(cons.consensus_seq.clone()).unwrap_or_default();
                    let (bracketed, inserted_str) = format_bracketed(&cons, bp);

                    let mean_cov = if cons.coverage.is_empty() {
                        0.0
                    } else {
                        cons.coverage.iter().sum::<u32>() as f64 / cons.coverage.len() as f64
                    };

                    let fasta_header = format!(
                        "{}-{}_{}_{}_{}_{}_n{}",
                        bp.gene0_name,
                        bp.gene1_name,
                        bp.gene0_chr,
                        bp.gene0_pos,
                        bp.gene1_chr,
                        bp.gene1_pos,
                        cons.n_reads,
                    );

                    fasta_entries.push((fasta_header, seq_str.clone()));

                    results.push(FusionBreakpointConsensus {
                        gene0_name: bp.gene0_name.clone(),
                        gene0_chr: bp.gene0_chr.clone(),
                        gene0_pos: bp.gene0_pos,
                        gene1_name: bp.gene1_name.clone(),
                        gene1_chr: bp.gene1_chr.clone(),
                        gene1_pos: bp.gene1_pos,
                        consensus_sequence: seq_str,
                        bracketed_sequence: bracketed,
                        inserted_bases: inserted_str,
                        n_reads: cons.n_reads,
                        mean_coverage: mean_cov,
                    });
                    info!(
                        "  Consensus: {} bp, {} reads, {:.1}x mean coverage",
                        cons.consensus_seq.len(),
                        cons.n_reads,
                        mean_cov,
                    );
                }
                Ok(None) => {
                    info!(
                        "  Insufficient reads for consensus (< {})",
                        config.min_reads
                    );
                }
                Err(e) => {
                    log::warn!("  Error building consensus: {}", e);
                }
            }
        }
    }

    // Write FASTA if prefix provided
    if let Some(prefix) = out_prefix
        && !fasta_entries.is_empty() {
            let fasta_path = format!("{}.breakpoints.fa", prefix);
            write_fasta(&fasta_path, &fasta_entries)?;
            info!("Breakpoint FASTA written to {}", fasta_path);
        }

    Ok(BreakpointConsensusOutput {
        breakpoints: results,
    })
}

// ============================================================================
// Core algorithm
// ============================================================================

// ============================================================================
// Template-based polishing helpers
// ============================================================================

/// Extract the full breakpoint region from a read: left_flank + junction + right_flank.
/// Returns (region_sequence, left_end_offset, right_start_offset).
fn extract_breakpoint_region(
    read: &BreakpointRead,
    flank_size: usize,
) -> Option<(Vec<u8>, usize, usize)> {
    let left_start = read.left_span.0 as usize;
    let left_end = read.left_span.1 as usize;
    let right_start = read.right_span.0 as usize;
    let right_end = read.right_span.1 as usize;

    if left_end <= left_start || right_end <= right_start
        || left_end > read.sequence.len() || right_end > read.sequence.len()
    {
        return None;
    }

    // Take up to flank_size bases from left flank (rightmost portion, near junction)
    let left_seq_full = &read.sequence[left_start..left_end];
    let left_take_start = left_seq_full.len().saturating_sub(flank_size);
    let left_flank = &left_seq_full[left_take_start..];

    // Take up to flank_size bases from right flank (leftmost portion, near junction)
    let right_seq_full = &read.sequence[right_start..right_end];
    let right_take = std::cmp::min(right_seq_full.len(), flank_size);
    let right_flank = &right_seq_full[..right_take];

    let left_end_offset = left_flank.len();
    let right_start_offset = left_end_offset + read.junction_seq.len();

    let mut region = Vec::with_capacity(left_flank.len() + read.junction_seq.len() + right_flank.len());
    region.extend_from_slice(left_flank);
    region.extend_from_slice(&read.junction_seq);
    region.extend_from_slice(right_flank);

    Some((region, left_end_offset, right_start_offset))
}


/// Majority-vote base call from counts [A, C, G, T]. Ties broken by ACGT order.
fn majority_base(counts: &[u32; 4]) -> u8 {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut best_idx = 0;
    for i in 1..4 {
        if counts[i] > counts[best_idx] {
            best_idx = i;
        }
    }
    bases[best_idx]
}

fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Build consensus by aligning all reads to a template and doing column-based
/// majority vote. Handles indels naturally through the alignment.
///
/// Returns (consensus_seq, coverage, new_left_end, new_right_start).
fn build_pileup_consensus(
    template: &[u8],
    read_regions: &[Vec<u8>],
    template_left_end: usize,
    template_right_start: usize,
    scoring: &AlignScoring,
    min_coverage: u32,
) -> (Vec<u8>, Vec<u32>, usize, usize) {
    let n_reads = read_regions.len();
    if n_reads == 0 || template.is_empty() {
        return (template.to_vec(), vec![0; template.len()], template_left_end, template_right_start);
    }

    let t_len = template.len();

    // Phase 1: Align all reads and collect per-template-position data.
    // For each read, we get: aligned base (or None=gap) at each template position,
    // plus insertions between template positions.
    // insertions_after[read_idx][template_pos] = inserted bases after that position.
    let mut aligned_bases: Vec<Vec<Option<u8>>> = Vec::with_capacity(n_reads);
    // insertions_after[template_pos] = Vec of (read_idx, inserted_bases)
    let mut insertions_after: Vec<Vec<(usize, Vec<u8>)>> = vec![Vec::new(); t_len + 1];

    for (r_idx, read_seq) in read_regions.iter().enumerate() {
        let aln = needleman_wunsch(template, read_seq, scoring);

        let mut bases_at_pos = vec![None; t_len];
        let mut t_pos: usize = 0;
        let mut q_pos: usize = 0;

        for &(op, len) in &aln.ops {
            match op {
                AlignOp::Match => {
                    for _ in 0..len {
                        if t_pos < t_len {
                            bases_at_pos[t_pos] = Some(read_seq[q_pos]);
                        }
                        t_pos += 1;
                        q_pos += 1;
                    }
                }
                AlignOp::InsRef => {
                    // Gap in query — template consumed, no query base
                    for _ in 0..len {
                        // bases_at_pos[t_pos] remains None (gap)
                        t_pos += 1;
                    }
                }
                AlignOp::InsQry => {
                    // Gap in template — query bases inserted
                    let ins_start = q_pos;
                    q_pos += len;
                    let ins_bases = read_seq[ins_start..q_pos].to_vec();
                    insertions_after[t_pos].push((r_idx, ins_bases));
                }
            }
        }
        aligned_bases.push(bases_at_pos);
    }

    // Phase 2: Decide which insertions to include.
    // For each template position, check if >50% of reads have an insertion there.
    let half = n_reads.div_ceil(2);
    let mut accepted_ins_len: Vec<usize> = vec![0; t_len + 1]; // 0 = no insertion accepted

    for pos in 0..=t_len {
        let n_with_ins = insertions_after[pos].len();
        if n_with_ins >= half {
            // Use median insertion length
            let mut lens: Vec<usize> = insertions_after[pos].iter().map(|(_, b)| b.len()).collect();
            lens.sort_unstable();
            accepted_ins_len[pos] = lens[lens.len() / 2];
        }
    }

    // Phase 3: Build consensus with junction tracking.
    let mut consensus: Vec<u8> = Vec::new();
    let mut coverage: Vec<u32> = Vec::new();
    let mut new_left_end = template_left_end;
    let mut new_right_start = template_right_start;

    // Helper: emit insertion columns after a given template position
    let emit_insertion = |pos: usize, consensus: &mut Vec<u8>, coverage: &mut Vec<u32>,
                               new_left_end: &mut usize, new_right_start: &mut usize| {
        let ins_len = accepted_ins_len[pos];
        if ins_len == 0 {
            return;
        }
        // Build per-position consensus for this insertion
        for k in 0..ins_len {
            let mut counts = [0u32; 4];
            let mut total = 0u32;
            for (_r_idx, ins_bases) in &insertions_after[pos] {
                if k < ins_bases.len()
                    && let Some(idx) = base_to_idx(ins_bases[k]) {
                        counts[idx] += 1;
                        total += 1;
                    }
            }
            if total >= min_coverage {
                consensus.push(majority_base(&counts));
                coverage.push(total);
                // Insertion before tracked positions shifts them
                if pos < template_left_end {
                    *new_left_end += 1;
                }
                if pos < template_right_start {
                    *new_right_start += 1;
                }
            }
        }
    };

    // Leading insertions (before first template position)
    emit_insertion(0, &mut consensus, &mut coverage, &mut new_left_end, &mut new_right_start);

    // Template columns + trailing insertions
    for t_pos in 0..t_len {
        // Template column: majority vote among aligned bases
        let mut counts = [0u32; 4];
        let mut n_gap = 0u32;
        for read_bases in &aligned_bases {
            match read_bases[t_pos] {
                Some(b) => {
                    if let Some(idx) = base_to_idx(b) {
                        counts[idx] += 1;
                    }
                }
                None => n_gap += 1,
            }
        }
        let n_base: u32 = counts.iter().sum();

        // Skip this column if majority of reads have a gap (deletion from template)
        if n_gap > n_base {
            // Adjust junction tracking if a column before the boundary is deleted
            if t_pos < template_left_end && new_left_end > 0 {
                new_left_end -= 1;
            }
            if t_pos < template_right_start && new_right_start > 0 {
                new_right_start -= 1;
            }
        } else if n_base >= min_coverage {
            consensus.push(majority_base(&counts));
            coverage.push(n_base);
        }

        // Insertions after this template position
        emit_insertion(t_pos + 1, &mut consensus, &mut coverage, &mut new_left_end, &mut new_right_start);
    }

    // Clamp junction boundaries
    new_left_end = new_left_end.min(consensus.len());
    new_right_start = new_right_start.min(consensus.len()).max(new_left_end);

    (consensus, coverage, new_left_end, new_right_start)
}

fn build_consensus_for_breakpoint(
    bam: &mut AlignmentInput,
    bp: &FusionBreakpoint,
    config: &BreakpointConsensusConfig,
) -> Result<Option<ConsensusResult>, Box<dyn std::error::Error>> {
    let split_reads = collect_split_reads(bam, bp, config)?;

    if split_reads.len() < config.min_reads {
        return Ok(None);
    }

    log::debug!(
        "  {} supporting reads, junction sizes: {:?}",
        split_reads.len(),
        split_reads.iter().take(5).map(|r| r.junction_seq.len()).collect::<Vec<_>>(),
    );

    // Verify orientation consistency — majority vote on left_is_gene0
    let n_left_g0 = split_reads.iter().filter(|r| r.left_is_gene0).count();
    let left_is_gene0 = n_left_g0 * 2 >= split_reads.len();

    let flank = config.flank_size;
    let total_reads = split_reads.len();

    // Extract full breakpoint regions from all reads (with junction boundaries)
    let extracted: Vec<(Vec<u8>, usize, usize)> = split_reads
        .iter()
        .filter_map(|r| extract_breakpoint_region(r, flank))
        .collect();

    let n_failed = total_reads - extracted.len();
    if n_failed > 0 {
        log::debug!("  {} of {} reads failed region extraction", n_failed, total_reads);
    }

    if extracted.len() < config.min_reads {
        return Ok(None);
    }

    // Select template: longest extracted region
    let template_idx = extracted
        .iter()
        .enumerate()
        .max_by_key(|(_, (seq, _, _))| seq.len())
        .map(|(i, _)| i)
        .unwrap(); // safe: extracted is non-empty

    let read_regions: Vec<Vec<u8>> = extracted.iter().map(|(seq, _, _)| seq.clone()).collect();
    let (mut current_template, mut current_left_end, mut current_right_start) =
        extracted[template_idx].clone();

    let scoring = AlignScoring::default();
    let max_rounds = 5;
    let mut prev_templates: Vec<Vec<u8>> = Vec::new();
    let mut final_coverage = vec![0u32; current_template.len()];

    for round in 0..max_rounds {
        let (consensus, coverage, new_left_end, new_right_start) = build_pileup_consensus(
            &current_template,
            &read_regions,
            current_left_end,
            current_right_start,
            &scoring,
            config.min_coverage,
        );

        // Check convergence
        if consensus == current_template {
            log::debug!("  Polishing converged after {} round(s)", round + 1);
            final_coverage = coverage;
            current_template = consensus;
            current_left_end = new_left_end;
            current_right_start = new_right_start;
            break;
        }

        // Check for oscillation
        if prev_templates.iter().any(|t| t == &consensus) {
            log::debug!("  Polishing oscillating, stopping at round {}", round + 1);
            final_coverage = coverage;
            current_template = consensus;
            current_left_end = new_left_end;
            current_right_start = new_right_start;
            break;
        }

        prev_templates.push(current_template);
        current_template = consensus;
        current_left_end = new_left_end;
        current_right_start = new_right_start;
        final_coverage = coverage;
    }

    // Extract junction insertion
    let inserted_bases = if current_right_start > current_left_end {
        current_template[current_left_end..current_right_start].to_vec()
    } else {
        Vec::new()
    };

    Ok(Some(ConsensusResult {
        consensus_seq: current_template,
        left_end: current_left_end,
        right_start: current_right_start,
        inserted_bases,
        left_is_gene0,
        n_reads: total_reads,
        coverage: final_coverage,
    }))
}

/// Collect split reads supporting a breakpoint via indexed BAM queries.
fn collect_split_reads(
    bam: &mut AlignmentInput,
    bp: &FusionBreakpoint,
    config: &BreakpointConsensusConfig,
) -> Result<Vec<BreakpointRead>, Box<dyn std::error::Error>> {
    let margin = config.query_margin;

    // Translate chromosome names to BAM naming convention
    let g0_chr_bam = bam.contig_mapper.to_bam_name(&bp.gene0_chr);
    let g1_chr_bam = bam.contig_mapper.to_bam_name(&bp.gene1_chr);

    let region_a = format!(
        "{}:{}-{}",
        g0_chr_bam,
        bp.gene0_pos.saturating_sub(margin),
        bp.gene0_pos + margin
    );
    let region_b = format!(
        "{}:{}-{}",
        g1_chr_bam,
        bp.gene1_pos.saturating_sub(margin),
        bp.gene1_pos + margin
    );

    // Collect all records from region A, keyed by read name
    let mut reads_a: HashMap<String, Vec<RecordData>> = HashMap::new();
    for result in bam.query(&region_a)? {
        let rec = result?;
        if rec.flag & 4 != 0 {
            continue;
        } // skip unmapped
        if let Some(name) = &rec.name {
            if name.contains(';') {
                continue;
            } // skip duplex
            reads_a.entry(name.clone()).or_default().push(RecordData {
                ref_id: rec.ref_id,
                pos: rec.pos,
                flag: rec.flag,
                seq: rec.seq.clone(),
                cigar_raw: rec.cigar_raw(),
                alignment_span: rec.alignment_span(),
            });
        }
    }

    // Query region B, find reads present in both regions
    let mut split_reads: Vec<BreakpointRead> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();

    for result in bam.query(&region_b)? {
        let rec = result?;
        if rec.flag & 4 != 0 {
            continue;
        }
        if let Some(name) = &rec.name {
            if seen.contains(name.as_str()) || name.contains(';') {
                continue;
            }
            if let Some(a_records) = reads_a.get(name.as_str()) {
                seen.insert(name.clone());

                let b_data = RecordData {
                    ref_id: rec.ref_id,
                    pos: rec.pos,
                    flag: rec.flag,
                    seq: rec.seq.clone(),
                    cigar_raw: rec.cigar_raw(),
                    alignment_span: rec.alignment_span(),
                };

                if let Some(br) = extract_breakpoint_read(a_records, &b_data, bp, config) {
                    split_reads.push(br);
                }
            }
        }
    }

    // Also check: same-chromosome fusions may have both segments in region A
    // (reads_a may contain reads with multiple alignments to the same region)
    if bp.gene0_chr == bp.gene1_chr {
        for (name, records) in &reads_a {
            if seen.contains(name.as_str()) {
                continue;
            }
            if records.len() >= 2 {
                // Try all pairs
                for i in 0..records.len() {
                    for j in (i + 1)..records.len() {
                        if let Some(br) =
                            extract_breakpoint_read(&[records[i].clone()], &records[j], bp, config)
                        {
                            seen.insert(name.clone());
                            split_reads.push(br);
                        }
                    }
                }
            }
        }
    }

    Ok(split_reads)
}

/// Intermediate record data for matching across regions.
#[derive(Clone)]
struct RecordData {
    #[allow(dead_code)]
    ref_id: i32,
    pos: i32,
    flag: u16,
    seq: Vec<u8>,
    cigar_raw: Vec<u32>,
    alignment_span: usize,
}

/// Extract a BreakpointRead from a pair of alignment records matching a breakpoint.
///
/// `a_records` are from the gene0 region query, `b_record` from gene1.
/// Returns None if no valid pair is found.
fn extract_breakpoint_read(
    a_records: &[RecordData],
    b_record: &RecordData,
    bp: &FusionBreakpoint,
    config: &BreakpointConsensusConfig,
) -> Option<BreakpointRead> {
    let tol = config.cluster_tolerance;

    for a_rec in a_records {
        // Compute AlignmentCoords for both
        let a_start = if a_rec.pos >= 0 {
            a_rec.pos as u32
        } else {
            continue;
        };
        let a_end = a_start + a_rec.alignment_span as u32;
        let a_coords = fix_sam_coords(a_start, a_end, &a_rec.cigar_raw, a_rec.flag).ok()?;

        let b_start = if b_record.pos >= 0 {
            b_record.pos as u32
        } else {
            continue;
        };
        let b_end = b_start + b_record.alignment_span as u32;
        let b_coords = fix_sam_coords(b_start, b_end, &b_record.cigar_raw, b_record.flag).ok()?;

        // Check: a_coords should be near gene0_pos, b_coords near gene1_pos
        // (or vice versa — we check both assignments)
        // Use direction-aware matching to check only the correct alignment end.
        let (g0_coords, g0_rec, g1_coords, g1_rec, swapped) =
            if is_near_breakpoint_directed(&a_coords, bp.gene0_pos, bp.gene0_dir, tol)
                && is_near_breakpoint_directed(&b_coords, bp.gene1_pos, bp.gene1_dir, tol)
            {
                (&a_coords, a_rec, &b_coords, b_record, false)
            } else if is_near_breakpoint_directed(&a_coords, bp.gene1_pos, bp.gene1_dir, tol)
                && is_near_breakpoint_directed(&b_coords, bp.gene0_pos, bp.gene0_dir, tol)
            {
                (&b_coords, b_record, &a_coords, a_rec, true)
            } else {
                continue;
            };

        // Determine which record is primary (no 0x800 flag) vs supplementary
        let g0_is_supp = g0_rec.flag & 0x800 != 0;
        let g1_is_supp = g1_rec.flag & 0x800 != 0;

        // Get full sequence from the PRIMARY alignment (no hard clips)
        // Supplementary alignments have hard clips so their seq is truncated
        let (seq_rec, is_rev) = if !g0_is_supp && g1_is_supp {
            (g0_rec, g0_coords.is_reverse)
        } else if !g1_is_supp && g0_is_supp {
            (g1_rec, g1_coords.is_reverse)
        } else {
            // Both primary or both supplementary — pick the one with more sequence
            if g0_rec.seq.len() >= g1_rec.seq.len() {
                (g0_rec, g0_coords.is_reverse)
            } else {
                (g1_rec, g1_coords.is_reverse)
            }
        };

        let sequence = if is_rev {
            reverse_complement(&seq_rec.seq)
        } else {
            seq_rec.seq.clone()
        };

        if sequence.is_empty() {
            continue;
        }

        // Compute adjusted junction points in query coords based on expected
        // breakpoint positions. This corrects for alignment variability where
        // the aligner extends a segment past the actual breakpoint.
        let g0_junction_q = compute_junction_q(g0_coords, bp.gene0_pos, bp.gene0_dir);
        let g1_junction_q = compute_junction_q(g1_coords, bp.gene1_pos, bp.gene1_dir);

        // Sanity check: junction points should be within the sequence
        let seq_len = sequence.len() as i64;
        if g0_junction_q < 0 || g0_junction_q > seq_len
            || g1_junction_q < 0 || g1_junction_q > seq_len
        {
            continue;
        }
        let g0_jq = g0_junction_q as u32;
        let g1_jq = g1_junction_q as u32;

        // Determine gene content spans based on direction + strand.
        // For Left (breakpoint at te): gene extends LEFTWARD from breakpoint.
        //   Forward: content = qs..junction_q (junction on right)
        //   Reverse: content = junction_q..qe (junction on left)
        // For Right (breakpoint at ts): gene extends RIGHTWARD from breakpoint.
        //   Forward: content = junction_q..qe (junction on left)
        //   Reverse: content = qs..junction_q (junction on right)
        let (g0_content_start, g0_content_end) = gene_content_span(
            g0_coords, g0_jq, bp.gene0_dir,
        );
        let (g1_content_start, g1_content_end) = gene_content_span(
            g1_coords, g1_jq, bp.gene1_dir,
        );

        // Determine which gene is left (earlier in query) vs right
        let (left_start, left_end, right_start, right_end, left_is_gene0) =
            if g0_content_start <= g1_content_start {
                (g0_content_start, g0_jq, g1_jq, g1_content_end, true)
            } else {
                (g1_content_start, g1_jq, g0_jq, g0_content_end, false)
            };

        // Validate coordinates
        let left_end_usize = left_end as usize;
        let right_start_usize = right_start as usize;
        if left_end_usize > sequence.len() || right_start_usize > sequence.len() {
            continue;
        }

        let mut junction_seq = if right_start > left_end {
            sequence[left_end_usize..right_start_usize].to_vec()
        } else {
            Vec::new()
        };

        let _ = swapped;

        // Normalize orientation: ensure all reads have the same gene on the
        // same side. The expected orientation is gene with "<-|" direction on
        // the left (it extends leftward from breakpoint) and "|->'" on right.
        // For reverse-strand reads, fix_sam_coords puts genes in the opposite
        // order (original-strand = sequencer direction, not reference direction).
        // We normalize by RC'ing the sequence and swapping spans when needed.
        let expected_gene0_left = bp.gene0_dir == BreakpointDirection::Left;
        let mut final_sequence = sequence;
        let mut final_left_start = left_start;
        let mut final_left_end = left_end;
        let mut final_right_start = right_start;
        let mut final_right_end = right_end;
        let mut final_left_is_gene0 = left_is_gene0;

        if left_is_gene0 != expected_gene0_left {
            // Orientation is flipped — RC sequence and swap spans
            let slen = final_sequence.len() as u32;
            final_sequence = reverse_complement(&final_sequence);
            junction_seq = reverse_complement(&junction_seq);

            // Span (a,b) in original → (slen-b, slen-a) in RC'd sequence
            let new_left_start = slen - right_end;
            let new_left_end = slen - right_start;
            let new_right_start = slen - left_end;
            let new_right_end = slen - left_start;

            final_left_start = new_left_start;
            final_left_end = new_left_end;
            final_right_start = new_right_start;
            final_right_end = new_right_end;
            final_left_is_gene0 = !left_is_gene0;
        }

        return Some(BreakpointRead {
            left_span: (final_left_start, final_left_end),
            right_span: (final_right_start, final_right_end),
            left_is_gene0: final_left_is_gene0,
            sequence: final_sequence,
            junction_seq,
        });
    }

    None
}

/// Compute the junction point in query coordinates for a gene's alignment,
/// adjusted for the offset between the actual alignment endpoint and the
/// expected breakpoint position.
///
/// For direction `<-|` (breakpoint at te):
///   Forward: te maps to qe. Junction at qe adjusted by (gene_pos - te).
///   Reverse: te maps to qs. Junction at qs adjusted by -(gene_pos - te).
/// For direction `|->` (breakpoint at ts):
///   Forward: ts maps to qs. Junction at qs adjusted by (gene_pos - ts).
///   Reverse: ts maps to qe. Junction at qe adjusted by -(gene_pos - ts).
fn compute_junction_q(coords: &AlignmentCoords, bp_pos: u32, direction: BreakpointDirection) -> i64 {
    match direction {
        BreakpointDirection::Left => {
            // Breakpoint at te
            let delta = bp_pos as i64 - coords.te as i64;
            if coords.is_reverse {
                // Reverse: te maps to qs
                coords.qs as i64 - delta
            } else {
                // Forward: te maps to qe
                coords.qe as i64 + delta
            }
        }
        BreakpointDirection::Right => {
            // Breakpoint at ts
            let delta = bp_pos as i64 - coords.ts as i64;
            if coords.is_reverse {
                // Reverse: ts maps to qe
                coords.qe as i64 - delta
            } else {
                // Forward: ts maps to qs
                coords.qs as i64 + delta
            }
        }
    }
}

/// Determine the gene content span in query coordinates.
///
/// Returns (content_start, content_end) where the gene's actual sequence is.
/// The junction point is at the boundary — either content_end (junction on right)
/// or content_start (junction on left).
fn gene_content_span(coords: &AlignmentCoords, junction_q: u32, direction: BreakpointDirection) -> (u32, u32) {
    match direction {
        BreakpointDirection::Left => {
            // Gene extends leftward from breakpoint
            if coords.is_reverse {
                // Reverse: junction at qs side → content is junction_q..qe
                (junction_q, coords.qe)
            } else {
                // Forward: junction at qe side → content is qs..junction_q
                (coords.qs, junction_q)
            }
        }
        BreakpointDirection::Right => {
            // Gene extends rightward from breakpoint
            if coords.is_reverse {
                // Reverse: junction at qe side → content is qs..junction_q
                (coords.qs, junction_q)
            } else {
                // Forward: junction at qs side → content is junction_q..qe
                (junction_q, coords.qe)
            }
        }
    }
}

/// Check if an alignment's template coordinates match a breakpoint position,
/// using the breakpoint direction to check only the relevant alignment end.
///
/// Check if an alignment's template coordinates match a breakpoint position,
/// using the breakpoint direction to check only the relevant alignment end.
///
/// `Right` → breakpoint at template start (ts)
/// `Left`  → breakpoint at template end (te)
fn is_near_breakpoint_directed(
    coords: &AlignmentCoords,
    bp_pos: u32,
    direction: BreakpointDirection,
    tolerance: u32,
) -> bool {
    match direction {
        BreakpointDirection::Right => coords.ts.abs_diff(bp_pos) <= tolerance,
        BreakpointDirection::Left => coords.te.abs_diff(bp_pos) <= tolerance,
    }
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

// ============================================================================
// Output formatting
// ============================================================================

/// Format the bracketed primer-design-ready string.
///
/// Returns (bracketed_string, inserted_bases_string).
fn format_bracketed(cons: &ConsensusResult, bp: &FusionBreakpoint) -> (String, String) {
    let left_seq = String::from_utf8_lossy(&cons.consensus_seq[..cons.left_end]);
    let right_seq = String::from_utf8_lossy(&cons.consensus_seq[cons.right_start..]);
    let ins_str = String::from_utf8_lossy(&cons.inserted_bases).to_string();

    // Label with gene names based on orientation
    let (left_gene, right_gene) = if cons.left_is_gene0 {
        (&bp.gene0_name, &bp.gene1_name)
    } else {
        (&bp.gene1_name, &bp.gene0_name)
    };

    let bracketed = if cons.inserted_bases.is_empty() {
        format!("{}({})|({}){}",left_seq, left_gene, right_gene, right_seq)
    } else {
        format!(
            "{}({})[{}]({}){}",
            left_seq, left_gene, ins_str, right_gene, right_seq
        )
    };

    (bracketed, ins_str)
}

/// Write consensus sequences as FASTA.
fn write_fasta(
    path: &str,
    entries: &[(String, String)],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut f = std::fs::File::create(path)?;
    for (header, seq) in entries {
        writeln!(f, ">{}", header)?;
        for chunk in seq.as_bytes().chunks(80) {
            f.write_all(chunk)?;
            f.write_all(b"\n")?;
        }
    }
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
        assert_eq!(reverse_complement(b""), b"" as &[u8]);
    }

    #[test]
    fn test_extract_breakpoint_region_basic() {
        let read = BreakpointRead {
            left_span: (0, 10),
            right_span: (12, 22),
            left_is_gene0: true,
            sequence: b"AAAAACCCCCXXGGGGGTTTTT".to_vec(),
            junction_seq: b"XX".to_vec(),
        };
        let (region, left_end, right_start) = extract_breakpoint_region(&read, 200).unwrap();
        assert_eq!(region, b"AAAAACCCCCXXGGGGGTTTTT");
        assert_eq!(left_end, 10);
        assert_eq!(right_start, 12);
    }

    #[test]
    fn test_extract_breakpoint_region_short_flanks() {
        let read = BreakpointRead {
            left_span: (0, 5),
            right_span: (5, 10),
            left_is_gene0: true,
            sequence: b"AAAAACCCCC".to_vec(),
            junction_seq: Vec::new(),
        };
        let (region, left_end, right_start) = extract_breakpoint_region(&read, 3).unwrap();
        // Takes last 3 of left (AAA), no junction, first 3 of right (CCC)
        assert_eq!(region, b"AAACCC");
        assert_eq!(left_end, 3);
        assert_eq!(right_start, 3);
    }

    #[test]
    fn test_pileup_consensus_identical_reads() {
        let template = b"ACGTACGT".to_vec();
        let reads = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
        ];
        let scoring = AlignScoring::default();
        let (cons, cov, le, rs) = build_pileup_consensus(&template, &reads, 4, 4, &scoring, 1);
        assert_eq!(cons, b"ACGTACGT");
        assert_eq!(le, 4);
        assert_eq!(rs, 4);
        assert!(cov.iter().all(|&c| c == 3));
    }

    #[test]
    fn test_pileup_consensus_homopolymer_variation() {
        // Template: AAACCC (3A + 3C)
        // 2 reads with AAACCC (same), 1 read with AAAACCC (4A + 3C)
        // Should produce AAACCC (majority) — the extra A is an insertion in 1/3 reads
        let template = b"AAACCC".to_vec();
        let reads = vec![
            b"AAACCC".to_vec(),
            b"AAACCC".to_vec(),
            b"AAAACCC".to_vec(),
        ];
        let scoring = AlignScoring::default();
        let (cons, _cov, _le, _rs) = build_pileup_consensus(&template, &reads, 3, 3, &scoring, 1);
        assert_eq!(cons, b"AAACCC");
    }

    #[test]
    fn test_pileup_consensus_majority_insertion() {
        // Template: AAACCC
        // 2 reads with AAAACCC (extra A), 1 read with AAACCC
        // Majority has insertion → should produce AAAACCC
        let template = b"AAACCC".to_vec();
        let reads = vec![
            b"AAAACCC".to_vec(),
            b"AAAACCC".to_vec(),
            b"AAACCC".to_vec(),
        ];
        let scoring = AlignScoring::default();
        let (cons, _cov, _le, _rs) = build_pileup_consensus(&template, &reads, 3, 3, &scoring, 1);
        assert_eq!(cons, b"AAAACCC");
    }

    #[test]
    fn test_format_bracketed_no_insertion() {
        let bp = FusionBreakpoint {
            gene0_name: "BCR".to_string(),
            gene0_chr: "chr22".to_string(),
            gene0_pos: 23632000,
            gene0_dir: BreakpointDirection::Right,
            gene1_name: "ABL1".to_string(),
            gene1_chr: "chr9".to_string(),
            gene1_pos: 133729000,
            gene1_dir: BreakpointDirection::Left,
            n_reads: 10,
            overlap_gap: None,
        };
        let cons = ConsensusResult {
            consensus_seq: b"AAAAACCCCC".to_vec(),
            left_end: 5,
            right_start: 5,
            inserted_bases: Vec::new(),
            left_is_gene0: true,
            n_reads: 10,
            coverage: vec![10; 10],
        };
        let (bracketed, ins) = format_bracketed(&cons, &bp);
        assert_eq!(bracketed, "AAAAA(BCR)|(ABL1)CCCCC");
        assert_eq!(ins, "");
    }

    #[test]
    fn test_format_bracketed_with_insertion() {
        let bp = FusionBreakpoint {
            gene0_name: "EBF1".to_string(),
            gene0_chr: "chr5".to_string(),
            gene0_pos: 150662290,
            gene0_dir: BreakpointDirection::Left,
            gene1_name: "PDGFRB".to_string(),
            gene1_chr: "chr5".to_string(),
            gene1_pos: 159246527,
            gene1_dir: BreakpointDirection::Right,
            n_reads: 15,
            overlap_gap: None,
        };
        let cons = ConsensusResult {
            consensus_seq: b"AAAAAGGGCCCCC".to_vec(),
            left_end: 5,
            right_start: 8,
            inserted_bases: b"GGG".to_vec(),
            left_is_gene0: false, // gene1 (PDGFRB) is on left
            n_reads: 15,
            coverage: vec![15; 13],
        };
        let (bracketed, ins) = format_bracketed(&cons, &bp);
        assert_eq!(bracketed, "AAAAA(PDGFRB)[GGG](EBF1)CCCCC");
        assert_eq!(ins, "GGG");
    }

    #[test]
    fn test_config_defaults() {
        let config = BreakpointConsensusConfig::default();
        assert_eq!(config.flank_size, 200);
        assert_eq!(config.min_reads, 3);
        assert_eq!(config.min_coverage, 3);
        assert_eq!(config.query_margin, 500);
        assert_eq!(config.cluster_tolerance, 30);
        assert!((config.insertion_fraction - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_is_near_breakpoint_directed() {
        let coords = AlignmentCoords {
            qs: 0,
            qe: 1000,
            ts: 5000,
            te: 6000,
            is_reverse: false,
        };
        // Right: breakpoint at ts=5000
        assert!(is_near_breakpoint_directed(&coords, 5020, BreakpointDirection::Right, 30));
        assert!(!is_near_breakpoint_directed(&coords, 6020, BreakpointDirection::Right, 30));
        // Left: breakpoint at te=6000
        assert!(is_near_breakpoint_directed(&coords, 6020, BreakpointDirection::Left, 30));
        assert!(!is_near_breakpoint_directed(&coords, 5020, BreakpointDirection::Left, 30));
    }
}
