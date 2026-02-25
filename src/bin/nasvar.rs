use std::path::Path;
use clap::{Parser, Subcommand, ValueEnum};
use log::{info, warn, error};

#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use nasvar::config::{PipelineConfig, ReferenceConfig};
use nasvar::input::AlignmentInput;
use nasvar::output::{OutputCollector, UnifiedOutput};
use nasvar::utils::bed::read_bed;
use nasvar::var::cnv::{call_cnvs, CnvCallParams};
use nasvar::var::coverage::read_depth;
use nasvar::var::fusions::{call_fusions, FusionCallParams};
use nasvar::var::itd;
use nasvar::var::maf::{calc_maf, read_sites};
use nasvar::var::snv::call_snvs;
use nasvar::utils::annotation::PartnerGeneIndex;

#[derive(Parser)]
#[command(name = "nasvar")]
#[command(
    about = "Nanopore Adaptive Sample VARiant caller",
    long_about = "A targeted variant caller for long-read sequencing with a focus on human tumor-only data, supporting SNVs, CNVs, fusions, and more."
)]
struct Cli {
    /// Log verbosity level
    #[arg(long, global = true, default_value = "info")]
    log_level: LogLevel,
    /// Write log output to a file instead of stderr
    #[arg(long, global = true)]
    log_file: Option<String>,
    /// Append to log file instead of truncating
    #[arg(long, global = true)]
    append_log: bool,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, ValueEnum)]
enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl LogLevel {
    fn to_level_filter(&self) -> log::LevelFilter {
        match self {
            LogLevel::Error => log::LevelFilter::Error,
            LogLevel::Warn => log::LevelFilter::Warn,
            LogLevel::Info => log::LevelFilter::Info,
            LogLevel::Debug => log::LevelFilter::Debug,
            LogLevel::Trace => log::LevelFilter::Trace,
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    /// Call Single Nucleotide Variants
    Snv {
        /// Sorted and indexed BAM file containing aligned long reads. Must have an associated .bai index file.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA file. Must be indexed (i.e., a corresponding .fai file must exist).
        #[arg(long, required = true)]
        fasta: String,
        /// GFF3 file containing gene annotations. Ensure it contains 'gene' features with 'Name' or 'ID' attributes.
        #[arg(long, required = true)]
        gff: String,
        /// Optional list of specific genes to restrict calling to (comma-separated, e.g., "BRCA1,TP53").
        #[arg(long, value_delimiter = ',')]
        genes: Option<Vec<String>>,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named `<prefix>.result.json`.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Path to pipeline configuration JSON file (genes, thresholds, SNV transcripts).
        #[arg(long, required = true)]
        config: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
    },
    /// Call Copy Number Variants
    Cnv {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// BED file containing target regions (e.g., exons, capture targets).
        #[arg(long, required = true)]
        targets: String,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named <prefix>.result.json.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
        /// Path to coverage.tsv file (from 'coverage' command) for local CN calculation.
        #[arg(long)]
        coverage: Option<String>,
        /// BLAST ratio (tumor fraction) for CN adjustment (0.0 - 1.0).
        #[arg(long)]
        blast_ratio: Option<f64>,
        /// Path to pipeline configuration JSON file (genes, thresholds).
        #[arg(long, required = true)]
        config: String,
        /// Path to reference genome configuration JSON file (centromeres, PAR regions).
        #[arg(long, required = true)]
        reference: String,
    },
    /// Call Karyotype (LRDK)
    Karyotype {
        /// Path to coverage.tsv file.
        #[arg(long, required = true)]
        coverage: String,
        /// Path to MAF file.
        #[arg(long)]
        maf: Option<String>,
        /// Prefix for output files.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite.
        #[arg(short, long)]
        force: bool,
        /// Path to pipeline configuration JSON file (thresholds).
        #[arg(long, required = true)]
        config: String,
        /// Path to reference genome configuration JSON file (centromeres, PAR regions).
        #[arg(long, required = true)]
        reference: String,
        /// GC bias correction method for karyotype inference.
        #[arg(long, value_enum, default_value_t = nasvar::karyotype::GcCorrectionMethod::Linear)]
        gc_correction: nasvar::karyotype::GcCorrectionMethod,
    },
    /// Call Gene Fusions
    Fusions {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// BED file containing target regions (e.g., genes involved in fusions).
        #[arg(long, required = true)]
        targets: String,
        /// BED file containing repetitive elements to mask (e.g., generated by RepeatMasker).
        #[arg(long, required = true)]
        repeats: String,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named <prefix>.result.json.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
        /// Path to a text file containing genes to look for one-sided fusions (optional).
        #[arg(long)]
        one_sided: Option<String>,
        /// GFF3 file for partner gene coordinates (required for one-sided partner filtering).
        #[arg(long)]
        gff: Option<String>,
        /// Path to pipeline configuration JSON file (genes, thresholds).
        #[arg(long, required = true)]
        config: String,
    },
    /// Calculate Minor Allele Frequencies
    Maf {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// BED file containing enriched regions (targets) to prioritize.
        #[arg(long, required = true)]
        enriched: String,
        /// Text file containing sites of interest for MAF calculation (Format: CHR POS REF ALT).
        #[arg(long, required = true)]
        sites: String,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named <prefix>.maf.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
    },
    /// Calculate Read Depth / Coverage
    Coverage {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// BED file containing repetitive elements to mask.
        #[arg(long, required = true)]
        repeats: String,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named <prefix>.coverage.tsv.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
    },
    /// Run Full Pipeline
    Pipeline {
        #[arg(required = true)]
        bam: String,
        #[arg(required = true)]
        repeats: String,
        #[arg(required = true)]
        enriched: String,
        #[arg(required = true)]
        sites: String,
        #[arg(required = true)]
        targets: String,
        #[arg(required = true)]
        fasta: String,
        #[arg(required = true)]
        gff: String,
        #[arg(required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
        /// BLAST ratio (tumor fraction) for CN adjustment (0.0 - 1.0).
        #[arg(long)]
        blast_ratio: Option<f64>,
        /// Path to pipeline configuration JSON file (genes, thresholds, SNV transcripts, one-sided genes).
        #[arg(long, required = true)]
        config: String,
        /// Path to reference genome configuration JSON file (centromeres, PAR regions).
        #[arg(long, required = true)]
        reference: String,
        /// GC bias correction method for karyotype inference.
        #[arg(long, value_enum, default_value_t = nasvar::karyotype::GcCorrectionMethod::Linear)]
        gc_correction: nasvar::karyotype::GcCorrectionMethod,
    },
    /// Call Internal Tandem Duplications (ITDs)
    Itd {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// Path to ITD configuration JSON file specifying genes and regions to search.
        #[arg(long, required = true)]
        itd_config: String,
        /// Prefix for output files (e.g., "results/sample1"). Output files will be named <prefix>.result.json.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
    },
    /// Build consensus sequences for fusion breakpoints
    Breakpoints {
        /// Sorted and indexed BAM/CRAM file containing aligned long reads.
        #[arg(long, required = true)]
        bam: String,
        /// Reference genome FASTA (required for CRAM input).
        #[arg(long)]
        ref_fasta: Option<String>,
        /// Prefix for input/output files. Reads fusions from <prefix>.result.json.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Force overwrite of existing output files.
        #[arg(short, long)]
        force: bool,
        /// Number of flanking bases on each side of the breakpoint.
        #[arg(long, default_value = "200")]
        flank_size: usize,
        /// Minimum number of reads required for consensus.
        #[arg(long, default_value = "3")]
        min_reads: usize,
        /// Minimum per-position coverage to call a consensus base.
        #[arg(long, default_value = "3")]
        min_coverage: u32,
    },
    /// Print JSON Schema for unified output format
    Schema {
        /// Write schema to file instead of stdout
        #[arg(long)]
        output: Option<String>,
    },
    /// Generate Aggregate Report
    Report {
        /// Prefix for input files (e.g., "results/sample1"). Looks for <prefix>.result.json.
        #[arg(long, required = true)]
        out_prefix: String,
        /// Path to pipeline configuration file (for methodology section details)
        #[arg(long)]
        config: Option<String>,
        /// BED file containing analysis target regions (for enrichment table in report)
        #[arg(long)]
        targets: Option<String>,
    },
}

// Helper to check output paths and create directories
fn check_output_paths(
    prefix: &str,
    suffixes: &[&str],
    force: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(prefix);
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty() && !parent.exists() {
            info!("Creating output directory: {:?}", parent);
            std::fs::create_dir_all(parent)?;
        }

    if !force {
        for suffix in suffixes {
            let p = format!("{}{}", prefix, suffix);
            if Path::new(&p).exists() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::AlreadyExists,
                    format!(
                        "Output file {} already exists. Use --force to overwrite.",
                        p
                    ),
                )));
            }
        }
    }

    Ok(())
}

// Helper to open reader correctly typed
// noodles::bam::io::reader::Builder returns Reader<File> (which is Reader<std::io::Reader<File>>)

fn get_alignment_reader(path: &str, ref_path: Option<&str>) -> Result<AlignmentInput, Box<dyn std::error::Error>> {
    let reader = AlignmentInput::open(path, ref_path).map_err(|e| {
        std::io::Error::other(format!("Error opening alignment file {}: {}", path, e))
    })?;
    Ok(reader)
}

struct StepTimer {
    total_start: std::time::Instant,
    step_start: std::time::Instant,
}

impl StepTimer {
    fn new() -> Self {
        let now = std::time::Instant::now();
        Self {
            total_start: now,
            step_start: now,
        }
    }
    fn start(&mut self, name: &str) {
        info!("===== [STAGE] {} =====", name);
        self.step_start = std::time::Instant::now();
    }
    fn end(&self) {
        let now = std::time::Instant::now();
        info!("----- Stage Time: {:.2?} -----", now.duration_since(self.step_start));
        info!("----- Total Time: {:.2?} -----", now.duration_since(self.total_start));
    }
}

fn main() {
    // -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let mut log_builder = env_logger::Builder::from_default_env();
    log_builder
        .filter_level(cli.log_level.to_level_filter())
        .format_module_path(false);
    if let Some(ref path) = cli.log_file {
        let file = if cli.append_log {
            std::fs::File::options().create(true).append(true).open(path)
        } else {
            std::fs::File::create(path)
        }
        .unwrap_or_else(|e| panic!("Could not open log file '{}': {}", path, e));
        log_builder.target(env_logger::Target::Pipe(Box::new(file)));
    }
    log_builder.init();

    match &cli.command {
        Commands::Snv {
            bam,
            fasta,
            gff,
            genes,
            out_prefix,
            force,
            config,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".result.json"], *force) {
                error!("{}", e);
                return;
            }

            let pipeline_config = match PipelineConfig::load(config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading pipeline config {}: {}", config, e);
                    return;
                }
            };

            let mut br = match get_alignment_reader(bam, Some(fasta)) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }
            let g_list = genes.clone().unwrap_or(Vec::new());
            match call_snvs(&mut br, fasta, gff, &g_list, &pipeline_config) {
                Ok(snv_output) => {
                    let collector = OutputCollector::new().with_snv(snv_output);
                    if let Err(e) = collector.write_to_prefix(out_prefix) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => error!("Error calling SNVs: {}", e),
            }
        }
        Commands::Cnv {
            bam,
            ref_fasta,
            targets,
            out_prefix,
            force,
            coverage,
            blast_ratio,
            config,
            reference,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".result.json"], *force) {
                error!("{}", e);
                return;
            }

            // Load configs
            let pipeline_config = match PipelineConfig::load(config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading pipeline config {}: {}", config, e);
                    return;
                }
            };
            let ref_config = match ReferenceConfig::load(reference) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading reference config {}: {}", reference, e);
                    return;
                }
            };

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }
            let t_vec = match read_bed(targets) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading targets: {}", e);
                    return;
                }
            };
            // Try to load karyotype from output
            let unified_path = format!("{}.result.json", out_prefix);
            let karyo_data = if std::path::Path::new(&unified_path).exists() {
                UnifiedOutput::load_json(&unified_path)
                    .ok()
                    .and_then(|u| u.karyotype)
            } else {
                None
            };
            match call_cnvs(
                &mut br,
                &t_vec,
                &pipeline_config,
                CnvCallParams {
                    coverage_file: coverage.as_deref(),
                    blast_ratio: *blast_ratio,
                    karyotype: karyo_data.as_ref(),
                    ref_config: &ref_config,
                    precomputed_focal_depths: None,
                },
            ) {
                Ok(cnv_output) => {
                    let collector = OutputCollector::new().with_cnv(cnv_output);
                    if let Err(e) = collector.write_to_prefix(out_prefix) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => error!("Error calling CNVs: {}", e),
            }
        }
        Commands::Karyotype {
            coverage,
            maf,
            out_prefix,
            force,
            config,
            reference,
            gc_correction,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".result.json"], *force) {
                error!("{}", e);
                return;
            }

            // Load configs
            let pipeline_config = match PipelineConfig::load(config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading pipeline config {}: {}", config, e);
                    return;
                }
            };
            let ref_config = match ReferenceConfig::load(reference) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading reference config {}: {}", reference, e);
                    return;
                }
            };

            match nasvar::karyotype::call_karyotype_gc_corrected(coverage, maf.as_deref(), out_prefix, None, &ref_config, &pipeline_config.thresholds.karyotype, *gc_correction, None) {
                Ok(karyo_output) => {
                    let collector = OutputCollector::new().with_karyotype(karyo_output);
                    if let Err(e) = collector.write_to_prefix(out_prefix) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => error!("Error calling karyotype: {}", e),
            }
        }
        Commands::Fusions {
            bam,
            ref_fasta,
            targets,
            repeats,
            out_prefix,
            force,
            one_sided,
            gff,
            config,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".result.json"], *force) {
                error!("{}", e);
                return;
            }

            // Load pipeline config
            let pipeline_config = match PipelineConfig::load(config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading pipeline config {}: {}", config, e);
                    return;
                }
            };

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            let t_vec = match read_bed(targets) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading targets: {}", e);
                    return;
                }
            };
            // Note: --repeats is accepted but not used by fusion calling (kept for CLI compat).
            // Ambiguous alignments are now detected via query-coordinate overlap instead.
            let _r_vec = match read_bed(repeats) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading repeats: {}", e);
                    return;
                }
            };

            // Get one-sided genes: from file if provided, otherwise from config
            let os_set = if let Some(p) = one_sided {
                match std::fs::read_to_string(p) {
                    Ok(c) => Some(
                        c.lines()
                            .map(|s| s.trim().to_string())
                            .collect::<std::collections::HashSet<_>>(),
                    ),
                    Err(e) => {
                        error!("Error reading one_sided file: {}", e);
                        return;
                    }
                }
            } else {
                // Fall back to config's one-sided genes
                let genes = pipeline_config.get_one_sided_genes();
                if genes.is_empty() {
                    None
                } else {
                    Some(genes)
                }
            };

            // Load partner gene index if GFF provided
            let partner_index = if let Some(gff_path) = gff {
                let partner_genes = pipeline_config.get_all_partner_genes();
                if !partner_genes.is_empty() {
                    match PartnerGeneIndex::load_from_gff(gff_path, &partner_genes) {
                        Ok(idx) => {
                            info!("Loaded {} partner gene coordinates from GFF", idx.len());
                            Some(idx)
                        }
                        Err(e) => {
                            warn!("Could not load partner genes from GFF: {}", e);
                            None
                        }
                    }
                } else {
                    None
                }
            } else {
                None
            };

            match call_fusions(
                &mut br,
                bam,
                &t_vec,
                &pipeline_config,
                FusionCallParams {
                    one_sided: os_set.as_ref(),
                    partner_index: partner_index.as_ref(),
                    gene_depths: None,
                },
            ) {
                Ok(fusions_output) => {
                    info!("Fusion calling complete. Found {} fusions.", fusions_output.fusions.len());
                    let collector = OutputCollector::new().with_fusions(fusions_output);
                    if let Err(e) = collector.write_to_prefix(out_prefix) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => {
                    error!("Error calling Fusions: {}", e);
                }
            }
        }
        Commands::Maf {
            bam,
            ref_fasta,
            enriched,
            sites,
            out_prefix,
            force,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".maf"], *force) {
                error!("{}", e);
                return;
            }

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            let e_vec = match read_bed(enriched) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading enriched: {}", e);
                    return;
                }
            };
            let s_vec = match read_sites(sites) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading sites: {}", e);
                    return;
                }
            };
            if let Err(e) = calc_maf(&mut br, &e_vec, &s_vec, out_prefix) {
                error!("Error calculating MAF: {}", e);
            }
        }
        Commands::Coverage {
            bam,
            ref_fasta,
            repeats,
            out_prefix,
            force,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".coverage.tsv"], *force) {
                error!("{}", e);
                return;
            }

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }
            let r_vec = match read_bed(repeats) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading repeats: {}", e);
                    return;
                }
            };
            match read_depth(&mut br, &r_vec, out_prefix) {
                Ok(reads_aligned) => info!("Aligned reads counted: {}", reads_aligned),
                Err(e) => error!("Error calculating coverage: {}", e),
            }
        }
        Commands::Pipeline {
            bam,
            repeats,
            enriched,
            sites,
            targets,
            fasta,
            gff,
            out_prefix,
            force,
            blast_ratio,
            config,
            reference,
            gc_correction,
        } => {
            let mut timer = StepTimer::new();

            timer.start("Initialization");
            if let Err(e) = check_output_paths(
                out_prefix,
                &[
                    ".result.json",
                    ".maf",
                    ".coverage.tsv",
                ],
                *force,
            ) {
                error!("{}", e);
                return;
            }

            let mut br = match get_alignment_reader(bam, Some(fasta.as_str())) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }
            let r_vec = match read_bed(repeats) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading repeats: {}", e);
                    return;
                }
            };
            let e_vec = match read_bed(enriched) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading enriched: {}", e);
                    return;
                }
            };
            let s_vec = match read_sites(sites) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading sites: {}", e);
                    return;
                }
            };
            let t_vec = match read_bed(targets) {
                Ok(v) => v,
                Err(e) => {
                    error!("Error reading targets: {}", e);
                    return;
                }
            };

            // Load pipeline config
            let pipeline_config = match PipelineConfig::load(config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading pipeline config {}: {}", config, e);
                    return;
                }
            };

            // Load reference config
            let ref_config = match ReferenceConfig::load(reference) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading reference config {}: {}", reference, e);
                    return;
                }
            };

            // Load ITD config from pipeline config's itd section
            let itd_cfg = pipeline_config.genes.itd.clone();

            // Get one-sided fusion genes from config
            let os_set = {
                let genes = pipeline_config.get_one_sided_genes();
                if genes.is_empty() {
                    None
                } else {
                    Some(genes)
                }
            };

            // Load partner gene index for one-sided fusion filtering
            let partner_index = {
                let partner_genes = pipeline_config.get_all_partner_genes();
                if !partner_genes.is_empty() {
                    match PartnerGeneIndex::load_from_gff(gff, &partner_genes) {
                        Ok(idx) => {
                            info!("Loaded {} partner gene coordinates from GFF", idx.len());
                            Some(idx)
                        }
                        Err(e) => {
                            warn!("Could not load partner genes from GFF: {}", e);
                            None
                        }
                    }
                } else {
                    None
                }
            };

            timer.end();

            timer.start("Coverage, MAF, Fusions");
            let runner = nasvar::pipeline::PipelineRunner::new(bam, out_prefix)
                .with_reference(Some(fasta.as_str()))
                .with_coverage(r_vec.clone())
                .with_maf(s_vec.clone(), e_vec.clone())
                .with_fusions(t_vec.clone())
                .with_one_sided(os_set.clone())
                .with_partner_index(partner_index.clone())
                .with_focal_targets(t_vec.clone())
                .with_config(&pipeline_config);

            let (pipeline_output, reads_aligned, focal_depths) = match runner.run() {
                Ok(r) => (r.output, r.reads_aligned, r.focal_depths),
                Err(e) => {
                    error!("Error in PipelineRunner: {}", e);
                    return;
                }
            };
            timer.end();

            // Start building output from pipeline results (includes fusions, QC, metadata)
            let mut collector = OutputCollector::new();
            if let Some(ref fusions) = pipeline_output.fusions {
                collector = collector.with_fusions(fusions.clone());
            }
            if let Some(qc) = pipeline_output.qc.clone() {
                collector = collector.with_qc_output(qc);
            }
            if let Some(ref meta) = pipeline_output.metadata {
                collector = collector.with_metadata(meta.clone());
            }

            let cov_file = format!("{}.coverage.tsv", out_prefix);
            let maf_file = format!("{}.maf", out_prefix);

            // Karyotype (two-pass with GC bias correction)
            timer.start("Karyotype Analysis");
            let seg_bases = nasvar::karyotype::compute_seg_bases(&s_vec, &e_vec, &ref_config);
            let (karyo_result, est_blast_ratio) = match nasvar::karyotype::call_karyotype_gc_corrected(
                &cov_file,
                Some(&maf_file),
                out_prefix,
                reads_aligned,
                &ref_config,
                &pipeline_config.thresholds.karyotype,
                *gc_correction,
                Some(&seg_bases),
            ) {
                Ok(karyo_output) => {
                    let br = karyo_output.blast_ratio;
                    (Some(karyo_output), br)
                }
                Err(e) => {
                    error!("Error in call_karyotype: {}", e);
                    (None, None)
                }
            };
            let effective_blast_ratio = blast_ratio.or(est_blast_ratio);
            if let Some(br) = effective_blast_ratio {
                info!("Using effective blast ratio: {:.4}", br);
            }
            timer.end();

            // CNV (use GC-adjusted coverage if available)
            timer.start("CNV Calling");
            let gc_adjusted_cov_file = format!("{}.coverage.gc_adjusted.tsv", out_prefix);
            let cnv_cov_file = if std::path::Path::new(&gc_adjusted_cov_file).exists() {
                &gc_adjusted_cov_file
            } else {
                &cov_file
            };
            let cnv_result = match call_cnvs(
                &mut br,
                &t_vec,
                &pipeline_config,
                CnvCallParams {
                    coverage_file: Some(cnv_cov_file),
                    blast_ratio: effective_blast_ratio,
                    karyotype: karyo_result.as_ref(),
                    ref_config: &ref_config,
                    precomputed_focal_depths: focal_depths,
                },
            ) {
                Ok(cnv_output) => {
                    Some(cnv_output)
                }
                Err(e) => {
                    error!("Error in call_cnvs: {}", e);
                    None
                }
            };
            timer.end();

            // SNV
            timer.start("SNV Calling");
            let snv_result = match call_snvs(&mut br, fasta, gff, &Vec::new(), &pipeline_config) {
                Ok(snv_output) => {
                    Some(snv_output)
                }
                Err(e) => {
                    error!("Error in call_snvs: {}", e);
                    None
                }
            };
            timer.end();

            // ITD
            timer.start("ITD Detection");
            let itd_result = match itd::call_itds(&mut br, &itd_cfg) {
                Ok(itd_output) => {
                    Some(itd_output)
                }
                Err(e) => {
                    error!("Error in call_itds: {}", e);
                    None
                }
            };
            timer.end();

            // Build and write output
            if let Some(karyo) = karyo_result {
                collector = collector.with_karyotype(karyo);
            }
            if let Some(cnv) = cnv_result {
                collector = collector.with_cnv(cnv);
            }
            if let Some(snv) = snv_result {
                collector = collector.with_snv(snv);
            }
            if let Some(itd) = itd_result {
                collector = collector.with_itd(itd);
            }

            if let Err(e) = collector.write_to_prefix(out_prefix) {
                error!("Error writing output: {}", e);
            } else {
                info!("Output written to {}.result.json", out_prefix);
            }

            // Report
            timer.start("Report Generation");
            if let Err(e) = nasvar::report::generate_report(out_prefix, Some(&pipeline_config), Some(targets)) {
                error!("Error generating report: {}", e);
            }
            timer.end();
        }

        Commands::Itd {
            bam,
            ref_fasta,
            itd_config,
            out_prefix,
            force,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".result.json"], *force) {
                error!("{}", e);
                return;
            }

            let itd_cfg = match itd::load_itd_config(itd_config) {
                Ok(c) => c,
                Err(e) => {
                    error!("Error loading ITD config {}: {}", itd_config, e);
                    return;
                }
            };

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }
            match itd::call_itds(&mut br, &itd_cfg) {
                Ok(itd_output) => {
                    let collector = OutputCollector::new().with_itd(itd_output);
                    if let Err(e) = collector.write_to_prefix(out_prefix) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => error!("Error calling ITDs: {}", e),
            }
        }

        Commands::Breakpoints {
            bam,
            ref_fasta,
            out_prefix,
            force,
            flank_size,
            min_reads,
            min_coverage,
        } => {
            if let Err(e) = check_output_paths(out_prefix, &[".breakpoints.fa"], *force) {
                error!("{}", e);
                return;
            }

            // Load existing unified output to get fusion results
            let unified_path = format!("{}.result.json", out_prefix);
            let unified = match UnifiedOutput::load_json(&unified_path) {
                Ok(u) => u,
                Err(e) => {
                    error!(
                        "Error loading {}: {}. Run fusion calling first.",
                        unified_path, e
                    );
                    return;
                }
            };

            let fusions = match &unified.fusions {
                Some(f) => f.clone(),
                None => {
                    error!(
                        "No fusion results found in {}. Run fusion calling first.",
                        unified_path
                    );
                    return;
                }
            };

            if fusions.fusions.is_empty() {
                info!("No fusions detected -- nothing to do.");
                return;
            }

            let mut br = match get_alignment_reader(bam, ref_fasta.as_deref()) {
                Ok(b) => b,
                Err(e) => {
                    error!("Error opening alignment file: {}", e);
                    return;
                }
            };
            if let Err(e) = br.require_index(bam) {
                error!("{}", e);
                return;
            }

            let bp_config = nasvar::var::fusion_consensus::BreakpointConsensusConfig {
                flank_size: *flank_size,
                min_reads: *min_reads,
                min_coverage: *min_coverage,
                ..Default::default()
            };

            match nasvar::var::fusion_consensus::call_breakpoint_consensus(
                &mut br,
                &fusions,
                &bp_config,
                Some(out_prefix),
            ) {
                Ok(bp_output) => {
                    info!(
                        "Breakpoint consensus complete. {} breakpoints resolved.",
                        bp_output.breakpoints.len()
                    );
                    // Update existing unified output with consensus data
                    let mut updated = unified;
                    updated.breakpoint_consensus = Some(bp_output);
                    if let Err(e) = updated.write_json(&unified_path) {
                        error!("Error writing output: {}", e);
                    }
                }
                Err(e) => error!("Error building breakpoint consensus: {}", e),
            }
        }

        Commands::Schema { output } => {
            let schema = nasvar::output::schema::schema_json_pretty();
            if let Some(path) = output {
                if let Err(e) = std::fs::write(path, &schema) {
                    error!("Error writing schema: {}", e);
                } else {
                    info!("Schema written to {}", path);
                }
            } else {
                println!("{}", schema);
            }
        }
        Commands::Report {
            out_prefix,
            config,
            targets,
        } => {
            let pipeline_config = config.as_ref().and_then(|p| {
                match PipelineConfig::load(p) {
                    Ok(cfg) => Some(cfg),
                    Err(e) => {
                        warn!("Could not load config file {}: {}", p, e);
                        None
                    }
                }
            });
            if let Err(e) = nasvar::report::generate_report(out_prefix, pipeline_config.as_ref(), targets.as_deref()) {
                error!("Error generating report: {}", e);
            }
        }
    }
}
