use log::info;
use crate::bam::ContigMapper;
use crate::utils::contig::NUM_CHROMOSOMES;
use crate::input::{AlignmentInput, AlignmentHeader, AlignmentRecord, CigarKind};
use crate::config::PipelineConfig;
use crate::output::{OutputCollector, FusionsOutput, UnifiedOutput};
use crate::utils::bed::BedRegion;
use crate::utils::qc::PipelineQcData;
use crate::var::maf::{Site, filter_enriched_sites};
use crate::var::fusions::FusionAccumulator;  // Shared fusion accumulator
use crate::var::coverage::CoverageAccumulator;  // Shared coverage accumulator
use crate::utils::annotation::PartnerGeneIndex;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

/// Result from PipelineRunner::run().
pub struct PipelineResult {
    pub output: UnifiedOutput,
    pub reads_aligned: Option<u64>,
    pub focal_depths: Option<HashMap<String, f64>>,
}

pub struct PipelineRunner<'a> {
    bam_path: String,
    out_prefix: String,
    ref_path: Option<String>,
    // MAF data
    maf_sites: Option<Vec<Vec<Site>>>,
    maf_regions: Option<Vec<BedRegion>>,
    // Coverage data
    coverage_repeats: Option<Vec<BedRegion>>,
    // Fusion data
    fusion_targets: Option<Vec<BedRegion>>,
    fusion_one_sided: Option<HashSet<String>>,
    fusion_partner_index: Option<PartnerGeneIndex>,
    // Pipeline config
    config: Option<&'a PipelineConfig>,
    // Adaptive sampling alignments — replaces main BAM for coverage
    as_alignments: Option<String>,
    // GFF path used for variant-region BAM slicing (gene bounds for fusion
    // partners + chrom lookup for SNV/ITD genes). Required for slicing.
    gff_path: Option<String>,
    // When true (default), emit a <prefix>.slice.bam + .bai with reads
    // around called variants for inclusion in result-package exports.
    export_slice_bam: bool,
}

impl<'a> PipelineRunner<'a> {
    pub fn new(bam_path: &str, out_prefix: &str) -> Self {
        Self {
            bam_path: bam_path.to_string(),
            out_prefix: out_prefix.to_string(),
            ref_path: None,
            maf_sites: None,
            maf_regions: None,
            coverage_repeats: None,
            fusion_targets: None,
            fusion_one_sided: None,
            fusion_partner_index: None,
            config: None,
            as_alignments: None,
            gff_path: None,
            export_slice_bam: true,
        }
    }

    pub fn with_reference(mut self, ref_path: Option<&str>) -> Self {
        self.ref_path = ref_path.map(|s| s.to_string());
        self
    }

    pub fn with_maf(mut self, sites: Vec<Vec<Site>>, regions: Vec<BedRegion>) -> Self {
        self.maf_sites = Some(sites);
        self.maf_regions = Some(regions);
        self
    }

    pub fn with_coverage(mut self, repeats: Vec<BedRegion>) -> Self {
        self.coverage_repeats = Some(repeats);
        self
    }

    pub fn with_fusions(mut self, targets: Vec<BedRegion>) -> Self {
        self.fusion_targets = Some(targets);
        self
    }

    pub fn with_one_sided(mut self, genes: Option<HashSet<String>>) -> Self {
        self.fusion_one_sided = genes;
        self
    }

    pub fn with_config(mut self, config: &'a PipelineConfig) -> Self {
        self.config = Some(config);
        self
    }

    pub fn with_partner_index(mut self, index: Option<PartnerGeneIndex>) -> Self {
        self.fusion_partner_index = index;
        self
    }

    pub fn with_as_alignments(mut self, path: Option<String>) -> Self {
        self.as_alignments = path;
        self
    }

    pub fn with_gff_path(mut self, gff_path: Option<String>) -> Self {
        self.gff_path = gff_path;
        self
    }

    pub fn with_export_slice_bam(mut self, on: bool) -> Self {
        self.export_slice_bam = on;
        self
    }

    pub fn run(self, mut bam: &mut AlignmentInput) -> Result<PipelineResult, Box<dyn std::error::Error>> {
        info!("Starting Pipeline...");

        // Snapshot fusion_targets here — fusion calling later in this
        // function consumes self.fusion_targets, but the variant-region
        // slice at the end needs them for the partner-gene bounds lookup.
        let fusion_targets_for_slice: Option<Vec<BedRegion>> = self.fusion_targets.clone();

        // If no on-disk BAI was supplied, build one inline from pass 1's
        // sequential scan so the SNV / ITD / fusion-consensus / CNV stages
        // that need `query()` can run afterward.
        if !bam.has_index() {
            info!(
                "No index for {}; will build one on the fly during pass 1 (BAM must be coordinate-sorted).",
                self.bam_path
            );
            bam.enable_inline_index();
        }

        let header = bam.header.clone();

        // Initialize output collector with metadata from BAM header
        let mut collector = OutputCollector::new();
        if let Some(meta) = extract_metadata_from_header(&header.text) {
            collector = collector.with_metadata(meta);
        }

        // 1. Initialize Accumulators

        // Coverage Accumulator (skip if AS alignments will replace it)
        let use_as_for_coverage = self.as_alignments.is_some();
        let mut cov_acc = if !use_as_for_coverage {
            if let Some(reps) = &self.coverage_repeats {
                Some(CoverageAccumulator::new(&header, reps, self.ref_path.as_deref()))
            } else {
                None
            }
        } else {
            None
        };

        // MAF accumulator
        let mut maf_acc = self.maf_sites.as_ref().map(|sites| {
            let mut filtered = sites.clone();
            if let Some(regions) = &self.maf_regions {
                filter_enriched_sites(&mut filtered, regions);
            }
            MafAccumulator::new(&header, filtered)
        });

        // Fusion accumulator requires config for per-gene margins
        let mut fusion_acc = if let Some(targets) = &self.fusion_targets {
            let config = self.config.expect("PipelineConfig required for fusion calling - use .with_config()");
            Some(FusionAccumulator::new(&header, targets, config))
        } else {
            None
        };

        // QC Accumulator (uses same targets as fusions)
        let mut qc_acc = self
            .fusion_targets
            .as_ref()
            .map(|targets| QcAccumulator::new(&header, targets));


        // Main Loop
        info!("Pass 1: Main Scan...");
        let _start_time = std::time::Instant::now();
        let mut record;
        let mut i = 0;

        loop {
            if let Some(r) = bam.read_record()? {
                record = r;
            } else {
                break;
            }
            i += 1;
            if i % 10000 == 0 {
                eprint!("\rProcessed {} reads...", i);
                std::io::Write::flush(&mut std::io::stderr())?;
            }

            // Feed to accumulators
            if let Some(c) = &mut cov_acc {
                c.process(&record);
            }
            if let Some(m) = &mut maf_acc {
                m.process(&record);
            }
            if let Some(f) = &mut fusion_acc {
                f.process(&record);
            }
            if let Some(q) = &mut qc_acc {
                q.process(&record);
            }
        }
        info!("Processed {} reads. Done.", i);

        // Pass-1 stream is done. If we were building an inline BAI, finalize
        // it now so the subsequent random-access stages can query.
        bam.finalize_inline_index()?;

        // Adaptive sampling coverage scan (replaces main BAM for coverage and reads_aligned)
        let as_reads_aligned: Option<u64> = if let Some(ref as_path) = self.as_alignments
            && let Some(reps) = &self.coverage_repeats
        {
            use crate::var::coverage::scan_as_alignments;
            info!("Scanning adaptive sampling alignments for coverage...");
            let as_acc = scan_as_alignments(as_path, &header, reps, self.ref_path.as_deref())?;
            as_acc.write_output(&format!("{}.coverage.tsv", self.out_prefix), true)?;
            Some(as_acc.reads_aligned())
        } else {
            None
        };

        // Write coverage from main BAM (when not using AS alignments)
        if let Some(c) = cov_acc {
            c.write_output(&format!("{}.coverage.tsv", self.out_prefix), true)?;
        }

        // Write MAF
        if let Some(m) = maf_acc {
            m.write_output(&format!("{}.maf", self.out_prefix), &header)?;
        }

        // Finalize QC and get reads_aligned.
        // reads_aligned comes from the AS BAM when --as-alignments is provided,
        // otherwise from the main-BAM QcAccumulator. Target QC stats
        // (nt_on_target, reads_on_target, focal_depths) always come from the
        // main BAM regardless of AS.
        let (reads_aligned, focal_depths) = if let Some(q) = qc_acc {
            let count = as_reads_aligned.unwrap_or_else(|| q.reads_aligned());
            if as_reads_aligned.is_some() {
                info!("Total aligned reads (primary, from AS): {}", count);
            } else {
                info!("Total aligned reads (primary): {}", count);
            }
            let fd = q.focal_depths();
            collector = collector.with_qc(q.to_qc_data());
            collector = collector.with_reads_aligned(count);
            collector = collector.with_target_coverage(fd.clone());
            (Some(count), Some(fd))
        } else if let Some(count) = as_reads_aligned {
            info!("Total aligned reads (primary, from AS): {}", count);
            collector = collector.with_reads_aligned(count);
            (Some(count), None)
        } else {
            (None, None)
        };

        // Fusion Pass 2
        if let Some(f) = fusion_acc {
            // Finalize the accumulator to get hit_reads and gene_depths
            let targets = self.fusion_targets.as_ref().unwrap(); // safe because accumulator existed
            let (hit_reads, gene_depths) = f.finalize(targets);

            if !hit_reads.is_empty() {
                info!(
                    "Pass 2: Calling Fusions ({} candidates)...",
                    hit_reads.len()
                );

                use crate::var::fusions::{call_fusions_from_hits, FusionCallParams};
                let targets = self.fusion_targets.unwrap(); // move ownership

                // Get pipeline config (required for fusion calling)
                let config = self.config.expect("PipelineConfig required for fusion calling - use .with_config()");

                let fusions_output = call_fusions_from_hits(
                    &mut bam,
                    self.bam_path.as_str(),
                    &targets,
                    hit_reads,
                    config,
                    FusionCallParams {
                        one_sided: self.fusion_one_sided.as_ref(),
                        partner_index: self.fusion_partner_index.as_ref(),
                        gene_depths: Some(&gene_depths),
                    },
                )?;

                // Breakpoint consensus (requires indexed BAM)
                if !fusions_output.fusions.is_empty() && bam.has_index() {
                    use crate::var::fusion_consensus;
                    let bp_config = fusion_consensus::BreakpointConsensusConfig::default();
                    match fusion_consensus::call_breakpoint_consensus(
                        &mut bam,
                        &fusions_output,
                        &bp_config,
                        Some(&self.out_prefix),
                    ) {
                        Ok(bp_output) => {
                            info!(
                                "Breakpoint consensus: {} breakpoints resolved.",
                                bp_output.breakpoints.len()
                            );
                            collector = collector.with_breakpoint_consensus(bp_output);
                        }
                        Err(e) => {
                            log::error!("Breakpoint consensus error (non-fatal): {}", e);
                        }
                    }
                }

                collector = collector.with_fusions(fusions_output);
            } else {
                // No fusions found - add empty fusions output
                let empty_fusions = FusionsOutput {
                    fusions: vec![],
                    spike_in: None,
                };
                collector = collector.with_fusions(empty_fusions);
            }
        }

        // Write unified output
        collector.write_to_prefix(&self.out_prefix)?;

        let unified = collector.build();

        // Variant-region slice BAM. Best-effort: log and continue on failure
        // since this is a packaging convenience, not a pipeline result.
        if self.export_slice_bam {
            match self.gff_path.as_deref() {
                Some(gff) => {
                    if let Err(e) = crate::var::slice::dump_variant_slice(
                        &unified, &self.bam_path, gff,
                        fusion_targets_for_slice.as_deref(),
                        &self.out_prefix,
                    ) {
                        log::warn!("[slice] dump failed: {}", e);
                    }
                }
                None => {
                    log::info!("[slice] skipping — no GFF path supplied (use with_gff_path)");
                }
            }
        }

        Ok(PipelineResult {
            output: unified,
            reads_aligned,
            focal_depths,
        })
    }
}

// Use shared metadata extraction from utils
use crate::utils::metadata::extract_from_header as extract_metadata_from_header;

// ------------------ MAF Accumulator ------------------

pub struct MafAccumulator {
    sites: Vec<Vec<Site>>,
    counts: Vec<Vec<[u32; 5]>>,
}

impl MafAccumulator {
    pub fn new(header: &AlignmentHeader, sites_input: Vec<Vec<Site>>) -> Self {
        let mut sites_by_ref = vec![Vec::new(); header.refs.len()];
        let mut counts = vec![Vec::new(); header.refs.len()];

        // Map sites (indexed by canonical chromosome) to BAM ref_ids.
        // Support both NC_* accession and chr* naming conventions.
        for i in 0..NUM_CHROMOSOMES {
            let acc = ContigMapper::accession_from_index(i).unwrap();
            let chr = ContigMapper::chr_name_from_index(i).unwrap();
            let ref_id = header.refs.iter().position(|r| r == acc)
                .or_else(|| header.refs.iter().position(|r| r == chr));
            if let Some(ref_id) = ref_id
                && i < sites_input.len() {
                    sites_by_ref[ref_id] = sites_input[i].clone();
                    counts[ref_id] = vec![[0; 5]; sites_input[i].len()];
                }
        }

        Self {
            sites: sites_by_ref,
            counts,
        }
    }

    pub fn process(&mut self, record: &AlignmentRecord) {
        if record.ref_id < 0 {
            return;
        }
        let id = record.ref_id as usize;

        if id >= self.sites.len() || self.sites[id].is_empty() {
            return;
        }

        let ref_sites = &self.sites[id];
        let ref_counts = &mut self.counts[id];

        let start = record.alignment_start().unwrap_or(0);
        let seq = record.sequence();

        // Find relevant sites
        let mut site_idx = ref_sites.partition_point(|s| s.pos < start);

        // Iterate pileup
        let mut t = start;
        let mut q = 0; // query pos

        for &(op, len) in record.cigar_ops() {
            match op {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    let t_end = t + len;

                    while site_idx < ref_sites.len() && ref_sites[site_idx].pos < t_end {
                        let s = &ref_sites[site_idx];
                        if s.pos >= t {
                            let offset = s.pos - t;
                            if q + offset < seq.len() {
                                let base = seq.get(q + offset).copied().unwrap_or(b'N');
                                let al = match base {
                                    b'A' => 0,
                                    b'C' => 1,
                                    b'G' => 2,
                                    b'T' => 3,
                                    _ => 4,
                                };
                                ref_counts[site_idx][al] += 1;
                            }
                        }
                        site_idx += 1;
                    }

                    t += len;
                    q += len;
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    q += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    let t_end = t + len;
                    site_idx = ref_sites[site_idx..].partition_point(|s| s.pos < t_end) + site_idx;
                    t += len;
                }
                _ => {}
            }
        }
    }

    pub fn write_output(&self, path: &str, header: &AlignmentHeader) -> std::io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "chrom\tpos\tref\talt")?;

        for (id, sites) in self.sites.iter().enumerate() {
            if sites.is_empty() {
                continue;
            }
            let name = &header.refs[id];

            for (i, site) in sites.iter().enumerate() {
                let cnt = self.counts[id][i];
                let ref_ct = cnt[site.al0 as usize];
                let alt_ct = cnt[site.al1 as usize];
                writeln!(file, "{}\t{}\t{}\t{}", name, site.pos, ref_ct, alt_ct)?;
            }
        }
        Ok(())
    }
}

// ------------------ QC Accumulator ------------------

struct PerTarget {
    name: String,
    start: i32,
    end: i32,
    total_bases: u64,
}

pub struct QcAccumulator {
    /// Per-target tracking (flat list), indexed by `by_ref`
    targets: Vec<PerTarget>,
    /// ref_id -> indices into `targets`
    by_ref: HashMap<usize, Vec<usize>>,
    margin: i32,
    pub nt_on_target: u64,
    pub reads_on_target: u64,
    pub target_regions_nt: u64,
    pub reads_aligned: u64,
}

impl QcAccumulator {
    pub fn new(header: &AlignmentHeader, targets: &[BedRegion]) -> Self {
        let mapper = ContigMapper::from_refs(&header.refs);
        let mut per_targets = Vec::with_capacity(targets.len());
        let mut by_ref: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut target_regions_nt: u64 = 0;

        for t in targets {
            let bam_chrom = mapper.to_bam_name(&t.segment);
            if let Some(ref_id) = header.refs.iter().position(|r| r == &bam_chrom) {
                let idx = per_targets.len();
                per_targets.push(PerTarget {
                    name: t.name.clone(),
                    start: t.start as i32,
                    end: t.end as i32,
                    total_bases: 0,
                });
                by_ref.entry(ref_id).or_default().push(idx);
            }
            target_regions_nt += (t.end - t.start) as u64;
        }
        Self {
            targets: per_targets,
            by_ref,
            margin: 5000,
            nt_on_target: 0,
            reads_on_target: 0,
            target_regions_nt,
            reads_aligned: 0,
        }
    }

    pub fn process(&mut self, record: &AlignmentRecord) {
        // Count primary aligned reads (not unmapped, not secondary, not supplementary)
        let flag = record.flags();
        if (flag & 0x904) == 0 {
            self.reads_aligned += 1;
        }

        if record.ref_id < 0 {
            return;
        }
        let ref_id = record.ref_id as usize;

        let indices = match self.by_ref.get(&ref_id) {
            Some(v) => v,
            None => return,
        };

        let start = record.alignment_start().unwrap_or(0) as i32;
        let end = start + record.alignment_span() as i32;
        let mut counted_read = false;

        for &idx in indices {
            let t = &mut self.targets[idx];

            // Use margin for read detection, exact boundaries for base counting
            let g_st = t.start - self.margin;
            let g_en = t.end + self.margin;
            if start > g_en || end < g_st {
                continue;
            }

            // Walk CIGAR, count only M/=/X bases clipped to the target interval
            let mut ref_pos = start;
            let mut on_target_bases: u64 = 0;
            for &(op, len) in record.cigar_ops() {
                let len_i = len as i32;
                match op {
                    CigarKind::Match
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch => {
                        let ov_start = ref_pos.max(t.start);
                        let ov_end = (ref_pos + len_i).min(t.end);
                        if ov_start < ov_end {
                            on_target_bases += (ov_end - ov_start) as u64;
                        }
                        ref_pos += len_i;
                    }
                    CigarKind::Deletion | CigarKind::Skip => {
                        ref_pos += len_i;
                    }
                    _ => {}
                }
            }

            t.total_bases += on_target_bases;
            self.nt_on_target += on_target_bases;
            if !counted_read && on_target_bases > 0 {
                self.reads_on_target += 1;
                counted_read = true;
            }
        }
    }

    pub fn to_qc_data(&self) -> PipelineQcData {
        PipelineQcData {
            nt_on_target: self.nt_on_target as f64,
            reads_on_target: self.reads_on_target as f64,
            target_regions_nt: self.target_regions_nt as f64,
        }
    }

    pub fn reads_aligned(&self) -> u64 {
        self.reads_aligned
    }

    /// Per-gene average depth: total aligned bases / total region length,
    /// aggregated across all BED entries sharing the same gene name.
    pub fn focal_depths(&self) -> HashMap<String, f64> {
        let mut gene_bases: HashMap<String, u64> = HashMap::new();
        let mut gene_lengths: HashMap<String, u64> = HashMap::new();
        for t in &self.targets {
            *gene_bases.entry(t.name.clone()).or_default() += t.total_bases;
            *gene_lengths.entry(t.name.clone()).or_default() += (t.end - t.start) as u64;
        }
        gene_bases
            .iter()
            .map(|(name, &bases)| {
                let len = gene_lengths[name];
                let depth = if len > 0 { bases as f64 / len as f64 } else { 0.0 };
                (name.clone(), depth)
            })
            .collect()
    }
}
