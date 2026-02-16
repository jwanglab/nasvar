use log::info;
use crate::bam::ContigMapper;
use crate::input::{AlignmentInput, AlignmentHeader, AlignmentRecord, CigarKind};
use crate::config::PipelineConfig;
use crate::output::{OutputCollector, FusionsOutput, UnifiedOutput};
use crate::utils::bed::BedRegion;
use crate::utils::qc::PipelineQcData;
use crate::var::maf::Site;
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
    // Focal depth targets (CNV genes)
    focal_targets: Option<Vec<BedRegion>>,
    // Pipeline config
    config: Option<&'a PipelineConfig>,
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
            focal_targets: None,
            config: None,
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

    pub fn with_focal_targets(mut self, targets: Vec<BedRegion>) -> Self {
        self.focal_targets = Some(targets);
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

    pub fn run(self) -> Result<PipelineResult, Box<dyn std::error::Error>> {
        info!("Starting Pipeline...");

        let mut bam = AlignmentInput::open(&self.bam_path, self.ref_path.as_deref())?;
        let header = bam.header.clone();

        // Initialize output collector with metadata from BAM header
        let mut collector = OutputCollector::new();
        if let Some(meta) = extract_metadata_from_header(&header.text) {
            collector = collector.with_metadata(meta);
        }

        // 1. Initialize Accumulators

        // Coverage Accumulator
        let mut cov_acc = if let Some(reps) = &self.coverage_repeats {
            let acc = CoverageAccumulator::new(&header, reps, self.ref_path.as_deref());
            Some(acc)
        } else {
            None
        };

        // MAF accumulator
        let mut maf_acc = self
            .maf_sites
            .as_ref()
            .map(|sites| MafAccumulator::new(&header, sites.clone()));

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

        // focal depth accumulator (to be used for CNV calls later)
        let mut focal_acc = self
            .focal_targets
            .as_ref()
            .map(|targets| FocalDepthAccumulator::new(&header, targets));

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
            if let Some(fd) = &mut focal_acc {
                fd.process(&record);
            }
        }
        info!("Processed {} reads. Done.", i);
        // Finalize
        // Write coverage (with GC content for pipeline mode)
        if let Some(c) = cov_acc {
            c.write_output(&format!("{}.coverage.tsv", self.out_prefix), true)?;
        }

        // Write MAF
        if let Some(m) = maf_acc {
            m.write_output(&format!("{}.maf", self.out_prefix), &header)?;
        }

        // Finalize QC and get reads_aligned
        let reads_aligned = if let Some(q) = qc_acc {
            let count = q.reads_aligned();
            info!("Total aligned reads (primary): {}", count);
            collector = collector.with_qc(q.to_qc_data());
            collector = collector.with_reads_aligned(count);
            Some(count)
        } else {
            None
        };

        // Finalize focal depths
        let focal_depths = focal_acc.map(|fd| fd.finalize());

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

        Ok(PipelineResult {
            output: collector.build(),
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
    min_depth: u32,
}

impl MafAccumulator {
    pub fn new(header: &AlignmentHeader, sites_input: Vec<Vec<Site>>) -> Self {
        let mut sites_by_ref = vec![Vec::new(); header.refs.len()];
        let mut counts = vec![Vec::new(); header.refs.len()];

        // Map sites (indexed 0-23 for chr1-22,X,Y) to BAM ref_ids.
        // Support both NC_* accession and chr* naming conventions.
        for i in 0..24 {
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
            min_depth: 20,
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
                if ref_ct + alt_ct >= self.min_depth {
                    writeln!(file, "{}\t{}\t{}\t{}", name, site.pos, ref_ct, alt_ct)?;
                }
            }
        }
        Ok(())
    }
}

// ------------------ QC Accumulator ------------------

pub struct QcAccumulator {
    target_map: HashMap<usize, Vec<BedRegion>>,
    margin: i32,
    pub nt_on_target: u64,
    pub reads_on_target: u64,
    pub target_regions_nt: u64,
    pub reads_aligned: u64,
}

impl QcAccumulator {
    pub fn new(header: &AlignmentHeader, targets: &[BedRegion]) -> Self {
        let mut target_map = HashMap::new();
        let mut target_regions_nt: u64 = 0;

        let mapper = ContigMapper::from_refs(&header.refs);

        for t in targets {
            let bam_chrom = mapper.to_bam_name(&t.segment);
            if let Some(id) = header.refs.iter().position(|r| r == &bam_chrom) {
                target_map.entry(id).or_insert(Vec::new()).push(t.clone());
            }
            target_regions_nt += (t.end - t.start) as u64;
        }
        Self {
            target_map,
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
        let id = record.ref_id as usize;

        if let Some(tgt_list) = self.target_map.get(&id) {
            let start = record.alignment_start().unwrap_or(0) as i32;
            let mut end = start;
            let mut read_len: u64 = 0;
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
                match op {
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                        read_len += len as u64;
                    }
                    _ => {}
                }
            }

            for tgt in tgt_list {
                let g_st = tgt.start as i32 - self.margin;
                let g_en = tgt.end as i32 + self.margin;

                if start <= g_en && end >= g_st {
                    self.reads_on_target += 1;
                    self.nt_on_target += read_len;
                    break;
                }
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
}

// ------------------ Focal Depth Accumulator ------------------

struct FocalTarget {
    name: String,
    start: usize,
    end: usize,
    total_bases: u64,
}

pub struct FocalDepthAccumulator {
    /// Targets sorted by (ref_id, start) for efficient lookup
    targets: Vec<FocalTarget>,
    /// Index into `targets` grouped by ref_id for fast filtering
    by_ref: HashMap<usize, Vec<usize>>,
}

impl FocalDepthAccumulator {
    pub fn new(header: &AlignmentHeader, targets: &[BedRegion]) -> Self {
        let mapper = ContigMapper::from_refs(&header.refs);
        let mut focal_targets = Vec::with_capacity(targets.len());
        let mut by_ref: HashMap<usize, Vec<usize>> = HashMap::new();

        for t in targets {
            let bam_name = mapper.to_bam_name(&t.segment);
            if let Some(ref_id) = header.refs.iter().position(|r| r == &bam_name) {
                let idx = focal_targets.len();
                focal_targets.push(FocalTarget {
                    name: t.name.clone(),
                    start: t.start as usize,
                    end: t.end as usize,
                    total_bases: 0,
                });
                by_ref.entry(ref_id).or_default().push(idx);
            }
        }

        Self {
            targets: focal_targets,
            by_ref,
        }
    }

    pub fn process(&mut self, record: &AlignmentRecord) {
        if record.ref_id < 0 {
            return;
        }
        let ref_id = record.ref_id as usize;

        let indices = match self.by_ref.get(&ref_id) {
            Some(v) => v,
            None => return,
        };

        let a_start = match record.alignment_start() {
            Some(p) => p,
            None => return,
        };
        let a_end = a_start + record.alignment_span();

        for &idx in indices {
            let t = &mut self.targets[idx];
            let o_start = std::cmp::max(a_start, t.start);
            let o_end = std::cmp::min(a_end, t.end);
            if o_end > o_start {
                t.total_bases += (o_end - o_start) as u64;
            }
        }
    }

    pub fn finalize(self) -> HashMap<String, f64> {
        let mut depths = HashMap::new();
        for t in &self.targets {
            let len = t.end - t.start;
            let depth = if len > 0 {
                t.total_bases as f64 / len as f64
            } else {
                0.0
            };
            depths.insert(t.name.clone(), depth);
        }
        depths
    }
}
