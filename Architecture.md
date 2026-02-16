# nasvar Architecture & Technical Reference

## Overview

The package produces two binaries:

- **`nasvar`** -- the main CLI with subcommands for each analysis type and a combined `pipeline` mode
- **`aggregate`** -- a multi-sample reporting tool that combines per-sample JSON outputs into a TSV

## Module Organization

```
src/
├── lib.rs                  # Library root; re-exports and declares submodules
├── bam.rs                  # BAM utility helpers
├── config.rs               # JSON configuration loading (PipelineConfig, ReferenceConfig, AggregateConfig)
├── input.rs                # Unified BAM/CRAM reader abstraction (AlignmentInput)
├── bin/
│   ├── nasvar.rs           # Main CLI binary (clap subcommands)
│   └── aggregate.rs        # Multi-sample TSV report generator
├── var/
│   ├── mod.rs              # Variant module declarations
│   ├── snv.rs              # SNV calling with pharmacogenomics annotation
│   ├── cnv.rs              # CNV calling (local/focal CN + split-read deletions/duplications)
│   ├── fusions.rs          # Gene fusion detection (two-stage) with ambiguity filtering
│   ├── fusion_consensus.rs # Breakpoint consensus sequence reconstruction
│   ├── itd.rs              # Internal tandem duplication calling
│   ├── coverage.rs         # Genome-wide binned read depth
│   └── maf.rs              # Minor allele frequency at specified sites
├── pipeline/
│   └── mod.rs              # PipelineRunner: single-pass accumulator orchestration
├── karyotype/
│   ├── mod.rs              # Karyotype inference with GC bias correction
│   └── plot.rs             # SVG karyotype visualization
├── output/
│   ├── mod.rs              # Output module declarations and re-exports
│   ├── types.rs            # All output structs (UnifiedOutput, per-module outputs)
│   ├── collector.rs        # OutputCollector builder pattern
│   └── schema.rs           # JSON Schema generation (via schemars)
├── report/
│   └── mod.rs              # Report generation (Markdown + HTML) from unified JSON
├── plotting/               # Custom SVG plotting library
│   ├── mod.rs              # Module root and re-exports
│   ├── axes.rs             # Axis rendering
│   ├── data.rs             # Data structures
│   ├── error.rs            # Error types
│   ├── figure.rs           # Figure container
│   ├── backend/            # SVG rendering backend
│   ├── element/            # Plot elements (axis, grid, legend, text)
│   ├── plot/               # Plot types (line, scatter)
│   ├── scale/              # Scaling (linear)
│   └── style/              # Colors, markers, themes, line/fill styles
└── utils/
    ├── mod.rs              # Utility module declarations
    ├── align.rs            # Affine-gap Needleman-Wunsch sequence aligner
    ├── annotation.rs       # Partner gene indexing from GFF (PartnerGeneIndex)
    ├── bed.rs              # BED file parsing
    ├── contig.rs           # Contig/chromosome name mapping
    ├── metadata.rs         # Sequencing metadata extraction from BAM headers
    ├── qc.rs               # QC data structures
    └── time.rs             # ISO 8601 timestamp utilities
```

## Core Design Patterns

### 1. Accumulator Pattern (Single-Pass Scanning)

The `PipelineRunner` (`src/pipeline/mod.rs`) processes the BAM file once and feeds each record to multiple accumulators simultaneously:

- **CoverageAccumulator** -- bins reads into 1 Mb windows, masks repeats via `BitVec`, tracks GC content
- **MafAccumulator** -- counts ref/alt alleles at pre-specified SNP sites using binary search
- **FusionScanner** -- identifies reads overlapping fusion target regions, collects read names for Pass 2
- **FocalDepthAccumulator** -- computes per-gene depth for CNV focal genes

This avoids redundant full-genome scans, reducing I/O from 7+ passes to 2.

**Accumulator internals** (`PipelineRunner::run()`):

```
for each record:
    coverage_accumulator.process(&record)
    maf_accumulator.process(&record)
    fusion_scanner.process(&record)
    focal_depth_accumulator.process(&record)
```

**CoverageAccumulator**:
- Builds a `BitVec` mask from repeat BED regions (1 bit per genome position)
- Bins reads into configurable windows (default 1 Mb)
- Tracks per-bin GC content queried from the reference FASTA at initialization
- Only counts primary alignments (skips supplementary/secondary)

**MafAccumulator**:
- Pre-indexes sites by chromosome and position
- Uses `partition_point` (binary search) for efficient site lookup within each read span
- Counts `[A, C, G, T, N]` at each site position via CIGAR walk
- Filters output by `min_depth` threshold (default: 20)

**FusionScanner**:
- Builds a `HashMap` of target regions by reference ID
- Checks read overlap with targets (including configurable margin per gene)
- Collects read names in a `HashSet` for Pass 2 re-query
- Tracks QC metrics: `nt_on_target`, `reads_on_target`, `target_regions_nt`

**FocalDepthAccumulator**:
- Computes exact depth over each focal gene region
- Returns `HashMap<gene_name, depth>` consumed by the CNV caller

### 2. Builder Pattern (PipelineRunner)

`PipelineRunner` uses a fluent builder API:

```rust
PipelineRunner::new(bam_path, out_prefix)
    .with_reference(Some(fasta))
    .with_coverage(repeats)
    .with_maf(sites, enriched)
    .with_fusions(targets)
    .with_focal_targets(targets)
    .with_one_sided(one_sided_genes)
    .with_partner_index(partner_index)
    .with_config(&config)
    .run()
```

Components are optional -- the runner skips accumulators that aren't configured.

### 3. Unified Output (OutputCollector)

All pipeline results flow through `OutputCollector` (`src/output/collector.rs`) into a single `UnifiedOutput` struct (`src/output/types.rs`), which is serialized to `{prefix}.result.json`. The collector uses a builder pattern:

```rust
OutputCollector::new()
    .with_metadata(meta)
    .with_qc(qc_data)
    .with_fusions(fusions_output)
    .with_karyotype(karyo_output)
    .with_cnv(cnv_output)
    .with_snv(snv_output)
    .with_itd(itd_output)
    .with_breakpoint_consensus(bp_output)
    .write_to_prefix(&out_prefix)
```

TSV outputs (`.coverage.tsv`, `.maf`), SVG plots, and reports are written separately by their respective modules, not through the unified JSON.

The `UnifiedOutput` struct derives `JsonSchema` (via `schemars`) for automatic schema generation. The `nasvar schema` subcommand prints the schema to stdout or file. `UnifiedOutput` can be loaded back from disk via `UnifiedOutput::load_json()` for downstream tools (report, breakpoints, aggregate).

### 4. Unified Input Abstraction

`AlignmentInput` (`src/input.rs`) wraps both BAM and CRAM readers behind a single interface:

- **Format detection**: Reads the first 4 bytes of the file. If they match `b"CRAM"`, opens as CRAM; otherwise opens as BAM.
- **Index handling**: Automatically looks for `.bai` (BAM) or `.crai` (CRAM) index files. The `require_index()` method validates that an index exists before region queries.
- **Record decoding**: Noodles SAM records are converted to owned `AlignmentRecord` structs containing name, ref\_id, position (0-based), flags, MAPQ, sequence, quality, and pre-parsed CIGAR operations.
- **CRAM buffering**: CRAM files are read container-by-container into a `VecDeque<AlignmentRecord>` buffer. Only `seek(0)` is supported (re-opens the file).
- **CRAM reference**: Requires a FASTA reference for decoding. A `noodles::fasta::Repository` is constructed from an indexed FASTA reader.

## Data Flow

### Full Pipeline Execution (`nasvar pipeline`)

```
Initialization
  ├── Load PipelineConfig + ReferenceConfig from JSON
  ├── Open BAM/CRAM via AlignmentInput
  ├── Parse BED files (targets, repeats, enriched) and MAF sites
  └── Load partner gene index from GFF (for one-sided fusions)
       │
       ▼
Pass 1: PipelineRunner.run()
  ├── CoverageAccumulator  ──→  {prefix}.coverage.tsv
  ├── MafAccumulator       ──→  {prefix}.maf
  ├── FusionScanner        ──→  candidate read names (HashSet)
  └── FocalDepthAccumulator ──→  HashMap<gene, depth>
       │
       ▼
Pass 2: Fusion Refinement
  └── Re-query BAM for candidate reads → FusionsOutput
       │
       ▼
Breakpoint Consensus (optional, requires indexed BAM)
  └── Re-query BAM per breakpoint → iterative polishing → BreakpointConsensusOutput
      Also writes {prefix}.breakpoints.fa
       │
       ▼
Karyotype Analysis
  ├── Two-pass GC bias correction (linear or LOESS, configurable via --gc-correction)
  ├── MAF-based ploidy estimation
  ├── SVG plots ({prefix}.karyotype*.svg, {prefix}.gc_vs_coverage*.svg)
  └── KaryotypeOutput + {prefix}.coverage.gc_adjusted.tsv
       │
       ▼
CNV Calling
  ├── Uses GC-adjusted coverage for local CN
  ├── Uses focal depths from Pass 1 for focal CN
  ├── Split-read deletion/duplication detection
  ├── Blast ratio adjustment (if provided or estimated)
  └── CnvOutput (accepts KaryotypeOutput directly for chromosomal context)
       │
       ▼
SNV Calling → SnvOutput
ITD Detection → ItdOutput
       │
       ▼
OutputCollector
  └── Collects all *Output structs → {prefix}.result.json
       │
       ▼
Report Generation (if run)
  └── Reads {prefix}.result.json → {prefix}.report.md + {prefix}.report.html
```

### Standalone Subcommand Execution

Each analysis can be run independently via its subcommand (`nasvar snv`, `nasvar cnv`, etc.). These read the BAM directly and produce their respective output via `OutputCollector` → `.result.json`.

## Algorithms

### Karyotype Inference (LRDK)

**File**: `src/karyotype/mod.rs`

The entry point is `call_karyotype_gc_corrected`, which runs the full two-pass pipeline:

#### Two-Pass Architecture

- **Pass 1** (`call_karyotype` on raw coverage): Produces an initial karyotype estimate. This is needed to identify which chromosome arms are at majority ploidy (typically diploid), so the GC correction can be trained on segments with consistent copy number.
- **GC bias correction** (`gc_correct_coverage`): Fits a regression model on majority-ploidy autosomal bins, then adjusts all bins. Outputs `{prefix}.coverage.gc_adjusted.tsv`. The correction method is controlled by `--gc-correction` (default: `linear`).
- **Pass 2** (`call_karyotype` on corrected coverage): Re-runs the full karyotype inference on GC-adjusted coverage for the final result.

When `--gc-correction none` is specified, GC correction is skipped entirely and Pass 1 results are returned directly.

#### GC Bias Correction

The correction method is configurable via the `--gc-correction` CLI flag (`GcCorrectionMethod` enum):

| Value | Method | Description |
|-------|--------|-------------|
| `linear` | Linear regression (default) | Fits `predicted = m * gc + b` on majority-ploidy bins |
| `loess` | LOESS local regression | Non-parametric smooth fit; better for non-linear GC bias |
| `none` | Disabled | Skips GC correction entirely |

**Common steps (both methods):**

1. Identify the majority autosomal ploidy from Pass 1 (most common CN among autosomal arms)
2. Collect `(gc_content, coverage)` pairs from bins belonging to majority-ploidy segments only (excludes sex chromosomes and non-majority arms to avoid confounding)
3. Fit the selected regression model on those pairs
4. Correct each bin: `corrected = observed * reference / predicted`, where `reference` is the model's predicted coverage at GC = 0.41 and `predicted` is the model's prediction at the bin's actual GC content. This normalizes all bins to the expected coverage at 41% GC.
5. Clamp negative values to 0

**Linear regression**: Fits a single slope and intercept via ordinary least squares. The fit curve is a straight line. This is adequate when GC bias is approximately linear across the 0.3-0.7 GC range.

**LOESS (Locally Estimated Scatterplot Smoothing)**: Fits a smooth, non-parametric curve that can capture non-linear GC bias patterns. Implementation details:

- **Tricube kernel**: Neighbor weights use `w(u) = (1 - |u|^3)^3` where `u = distance / max_distance`
- **Bandwidth**: 0.3 (30% of data points used as neighbors for each prediction)
- **Local model**: Weighted linear regression at each query point using the k nearest neighbors
- **Curve generation**: 50 evenly-spaced points from GC 0.3 to 0.7 are evaluated for plotting
- **Reference point**: Same GC = 0.41 normalization as linear, evaluated via `loess_predict_single`

The LOESS approach is recommended when GC bias appears curved rather than linear in the `gc_vs_coverage` diagnostic plots.

#### Level Finding (`find_levels`)

Histogram-based peak detection matching scipy `find_peaks` behavior:

1. Compute adaptive histogram bin width from the coverage data
2. Build histogram of all coverage values (autosomes + sex chromosomes), trimmed at 95th percentile per segment
3. Find local maxima with plateau handling: scan for runs of equal values, check that both sides are strictly lower, take the midpoint as peak position
4. Distance post-filter: sort peaks by height (descending), greedily keep peaks that are at least 3 bins apart
5. Refine each peak: take the median of all data points within +/- 1.5 bin widths of the peak center, then recount support within +/- 1 bin width of the refined position
6. Height filter: discard peaks with < 5% of the tallest peak's count
7. Autosomal assignment filter: assign each autosomal segment to its closest level (by ratio, within 0.3 tolerance); keep only levels that have at least one autosomal segment assigned
8. Duplicate filter: merge peaks within 9.5% of each other (keep the taller one)

#### CN State Resolution (`resolve_cn_states`)

Determines which coverage level corresponds to CN=1, CN=2, CN=3. Uses MAF data when available to identify the diploid level:

- **MAF peak > 0.4** at a level indicates that level is likely diploid (heterozygous SNPs cluster near 0.5 in diploid regions)
- **Coverage ratio check**: if `v2/v1 >= 1.6`, the two levels are likely 1n/2n; if < 1.6, likely 2n/3n
- When both levels have MAF > 0.4, the algorithm compares MAF peak heights and coverage ratios to disambiguate
- **Heuristic fallback** (no MAF data): if `v2 < v1 * 1.51`, assume 2n/3n; otherwise assume 1n/2n

#### CN Assignment

Each chromosome arm's median coverage is assigned to the closest CN state (CN=1, 2, or 3) within half the inter-level distance. Arms with coverage above CN=3 are extrapolated (`((med - cn2) / delta + 2).round()`). Arms below `delta/2` are assigned CN=0.

#### Post-Assignment Adjustments

After initial CN assignment, three corrections handle edge cases:

1. **X/Y sex chromosome adjustment**: Y chromosome coverage is compared to the expected germline 1n level (`cn2_avg / 2`). If Y falls near germline 1n, it's set to CN=1. X is then adjusted using a preliminary blast ratio estimate -- if X coverage matches the expected tumor 2n level accounting for blast fraction, it's set to CN=2.

2. **Low-penetrance hypodiploid detection**: When Y is the only segment at CN=1, and the "CN=3" average is approximately 2x the "CN=1" average, the algorithm detects a misidentified level structure. The supposed diploid level is actually monosomic, and the supposed triploid is actually diploid. All CN assignments are shifted down (2→1, 3→2). This correction is suppressed when MAF data confirms the lower level is genuinely diploid (peak > 0.4).

3. **High trisomy count**: When >25 segments are at CN=3 with none at CN=1, all levels are shifted down by 1 (indicating the level finder identified the wrong base ploidy).

#### Blast Ratio Estimation

Estimated from the spacing between CN depth levels (autosomal only):
- From 1n and 2n: `blast_ratio = (1 - CN1/CN2) * 2`
- From 2n and 3n: `blast_ratio = (CN3/CN2 * 2) - 2`
- When all three are available, uses the formula with more supporting segments

#### Karyotype String

Sum autosomal copies across all arms, determine sex chromosomes, format ISCN modifications (gains, losses, deletions). Special handling for acrocentric chromosomes (13, 14, 15, 21, 22) which lack a p arm.

#### QC Warnings

- Within-segment spread > `spread_warning` (0.08): indicates noisy or wavy coverage data
- Insufficient MAF sites (< 20 segments with > 1% MAF density): warns that MAF-based ploidy estimation may be unreliable

#### Visualization

Generate SVG plots for pre/post GC-correction karyotype views, combined karyotype + BAF plot, and GC bias diagnostics (`src/karyotype/plot.rs` using `src/plotting/`). The GC vs coverage plots (`plot_gc_vs_coverage`) accept a fitted curve as a series of `(gc, predicted)` points -- 2 points for linear (rendered as a straight line) or ~50 points for LOESS (rendered as a smooth curve). The plot label automatically indicates the fit type.

### Fusion Detection (Two-Stage with Ambiguity Filtering)

**File**: `src/var/fusions.rs`

**Stage 1 (FusionScanner in PipelineRunner)**:
- For each primary alignment, check if it overlaps any fusion target region (expanded by `default_margin` or gene-specific `special_margins`)
- Collect names of overlapping reads into a `HashSet`

**Stage 2 (call_fusions)**:
- Re-query BAM for collected reads by name
- Parse SA (supplementary alignment) tags to find split reads
- Identify fusion candidates: reads with alignments in two different target gene regions
- Cluster breakpoints and count supporting reads

**Query-coordinate ambiguity detection** (`is_query_ambiguous`):
- Replaces the previous BED-based repetitive element filtering with data-driven detection
- Checks if the same query region maps to multiple genomic loci
- Criteria: overlap > 50% of shorter alignment AND overlap > `min_anchor`
- Same-locus exemption: if both alignments map to overlapping reference regions, not ambiguous
- Same-gene exemption: multi-mapping within one gene target (e.g., DUX4 cassettes) is expected
- The `--repeats` CLI flag is still accepted for compatibility but does not affect fusion calling

**Range-based breakpoint clustering** (`BreakpointCluster`):
- Each read's breakpoint is represented as an interval `[observed_bpt, observed_bpt+gap]` (or `[bpt-overlap, bpt]` for microhomology)
- `overlap_gap = min(qe_A, qe_B) - max(qs_A, qs_B)`: negative = gap (inserted bases), positive = microhomology (shared bases)
- Gap direction extends toward the junction (Left: rightward, Right: leftward)
- Clusters merge when both gene intervals overlap (with 10bp `INTERVAL_MARGIN`) AND directions match
- Final positions: median of raw observed positions across reads
- Output: `FusionBreakpoint.overlap_gap: Option<i32>` -- median overlap/gap, shown in reports as `[Nbp hom]` or `[Nbp ins]`

**Filtering pipeline**:
1. Query ambiguity check -- multi-mapping reads rejected
2. Blacklist pairs (e.g., BCR-IGL, CYP2C19-CYP2C9) -- rejected
3. Self-fusions -- rejected unless gene is in `skip_self_fusions` list
4. Min supporting reads >= `min_supporting_reads` (3)
5. Min breakpoint reads >= `min_breakpoint_reads` (3), unless gene is in `bypass_breakpoint_filter`
6. One-sided fusions: if gene is in `one_sided_genes`, validate that the partner lands near a configured partner gene coordinate (from GFF). Uses separate `one_sided.min_breakpoint_reads` and optional `one_sided.min_fraction` thresholds.
7. Spike-in detection: fusions matching `spike_in.gene1`/`spike_in.gene2` are separated into a `spike_in` array in the output JSON

### Breakpoint Consensus (Iterative Polishing)

**File**: `src/var/fusion_consensus.rs`

Reconstructs exact sequences over fusion breakpoints using an iterative template-based polishing approach (Racon-like).

**Configuration** (`BreakpointConsensusConfig`):
- `flank_size`: 200 -- flanking bases on each side of breakpoint
- `min_reads`: 3 -- minimum reads required
- `min_coverage`: 3 -- minimum per-position coverage
- `query_margin`: 500 -- BAM query margin around breakpoint
- `cluster_tolerance`: 30 -- tolerance for matching read alignments to expected breakpoint
- `insertion_fraction`: 0.5 -- minimum prevalence for accepting an insertion column

**Algorithm**:
1. For each detected fusion breakpoint, query the indexed BAM for split reads in the breakpoint region (flank_size on each side)
2. Use direction-aware matching (`is_near_breakpoint_directed`) to find reads that span the specific breakpoint -- checking only the relevant end (ts or te) based on `BreakpointDirection` to avoid cross-contamination between nearby breakpoints
3. Extract left-flank and right-flank sequences from each read, normalize orientation (RC + swap spans when reverse-strand read has opposite gene order)
4. Select the longest read as the initial template
5. Iterative polishing (max 5 rounds):
   a. Align all reads to the current template via Needleman-Wunsch
   b. Build pileup consensus: majority-vote per column, delete columns with >50% gaps, accept insert columns if >50% of reads have insertion
   c. Track junction boundary (`left_end`, `right_start`) through column insertions/deletions
   d. Use consensus as the next template; stop if consensus matches previous round
6. Adjust junction position via `compute_junction_q` based on expected breakpoint position (the aligner can extend alignment past the actual breakpoint)
7. Output: bracketed primer format `GENE0_SEQ(gene0)[INSERTED](gene1)GENE1_SEQ`

**Integration**:
- Pipeline mode: auto-runs after fusion pass 2 if BAM is indexed (non-fatal on failure)
- Standalone: `nasvar breakpoints --bam <BAM> --out-prefix <PFX>` reads from `.result.json`
- Outputs: JSON in unified output (`breakpoint_consensus` field) + `.breakpoints.fa` FASTA

### Needleman-Wunsch Aligner

**File**: `src/utils/align.rs`

Affine-gap Needleman-Wunsch global sequence aligner used by the breakpoint consensus module.

**Scoring** (`AlignScoring`):
- `match_score`: +2
- `mismatch_score`: -2
- `gap_open`: -3
- `gap_extend`: -1

**Output**: `PairwiseAlignment` with run-length encoded operations (`AlignOp` enum: Match, InsRef, InsQry). Case-insensitive matching.

### CNV Calling

**File**: `src/var/cnv.rs`

**Copy number estimation**:
- Local CN: `gene_depth / genome_median_depth` (from coverage file)
- Focal CN: `focal_depth / genome_median_depth` (from Pass 1 accumulator)

**Blast ratio adjustment** (when tumor fraction is known):
```
adjusted_cn = (raw_cn * (1.0 - blast_ratio)) + (2.0 * blast_ratio)
```

**Structural events**:
- Query region = gene coordinates +/- `gene_margin` (default: 50 kb)
- Parse CIGAR for split reads; cluster breakpoints within `breakpoint_margin` (50 bp)
- Report deletions/duplications with >= `min_supporting_reads` (3) and >= `min_anchor` (500 bp)

**Integration with karyotype**: The CNV caller accepts `Option<&KaryotypeOutput>` directly (in-memory, not from file) to provide chromosomal context (per-arm copy numbers). Uses centromere coordinates from `ReferenceConfig` for p/q arm assignment.

### SNV Calling

**File**: `src/var/snv.rs`

1. Parse GFF3 to extract gene coordinates and CDS exons
2. For each configured gene (or all genes if `--genes` not specified):
   - Query BAM for the gene region
   - Pileup at each CDS position, counting base frequencies
   - Call variant if VAF >= `min_allele_freq` (0.3) and depth >= `min_coverage` (5.0)
3. Translate CDS variants to amino acid changes using a hardcoded codon table
4. Match against configured variant definitions (e.g., `p.R38H` at CDS position 226)
5. Determine zygosity: VAF >= `homozygosity_vaf` (0.8) = homozygous
6. Phasing: check if multiple variants co-occur on the same reads (`phasing_ratio` threshold)

Pharmacogenomics genes (configured in `snv.pharmacogenomics`) report all variants found, not just pre-configured ones. Pathogenic genes (configured in `snv.pathogenic`) are reported separately in reports.

### ITD Calling

**File**: `src/var/itd.rs`

1. Load ITD config (gene-specific regions with chrom, start, end, min\_length, min\_frequency)
2. For each region:
   - Query BAM for overlapping reads
   - Extract insertions from CIGAR (insertion operations >= `min_length`)
   - Identify tandem patterns in inserted sequences
   - Merge proximal events
3. Calculate allelic ratio: ITD-supporting reads / total coverage at the locus
4. Filter by `min_frequency` threshold

### Coverage Calculation

**File**: `src/var/coverage.rs`

- Iterates through all records in the BAM (sequential scan)
- Bins positions into 1 Mb windows (configurable via `coverage.bin_size`)
- Masks repeat regions using a `BitVec` (bits set for repeat positions are excluded from counts)
- Calculates per-bin GC content from the reference FASTA
- Counts only primary aligned reads
- Returns total aligned read count (used for karyotype normalization)

### MAF Calculation

**File**: `src/var/maf.rs`

- Reads sites file (tab-separated: `CHR POS REF ALT`)
- Groups sites by enriched regions (only counts alleles within enriched BED regions)
- For each overlapping read: walks CIGAR to determine the base at each site position
- Outputs ref/alt counts per site, filtered by `min_depth` (20)

## Configuration System

Three JSON configuration files drive the pipeline.

All configuration structs use `#[derive(Deserialize)]` with serde defaults. Each threshold type implements `Default`, so the entire `thresholds` section can be omitted from JSON and sensible defaults apply. The `PipelineConfig`, `ReferenceConfig`, and `AggregateConfig` structs each have a `load(path)` method that reads and deserializes from a JSON file.

### PipelineConfig (`--config`)

Top-level structure:

| Field | Type | Description |
|-------|------|-------------|
| `genes.cnv` | `CnvGeneConfig` | Focal genes, exclusions, deletion/duplication/local\_cn gene lists |
| `genes.itd` | `HashMap<String, ItdRegion>` | ITD regions: chrom, start, end, min\_length, min\_frequency |
| `genes.snv` | `SnvGeneConfig` | Pharmacogenomics genes, pathogenic genes, transcript/variant definitions |
| `genes.fusions` | `FusionGeneConfig` | One-sided genes/partners, blacklist pairs, bypass\_breakpoint\_filter, spike-in, special margins |
| `thresholds.*` | Various | Per-module thresholds (all have sensible defaults) |

### ReferenceConfig (`--reference`)

| Field | Type | Description |
|-------|------|-------------|
| `name` | `String` | Reference genome name (e.g., "T2T-CHM13v2.0") |
| `centromeres` | `HashMap<String, (u32, u32)>` | Centromere coordinates per chromosome |
| `par1` | `HashMap<String, (u32, u32)>` | PAR1 regions (chrX, chrY) |
| `par2` | `HashMap<String, (u32, u32)>` | PAR2 regions (chrX, chrY) |

### AggregateConfig (`--config` on aggregate binary)

| Field | Type | Description |
|-------|------|-------------|
| `itd_genes` | `Vec<String>` | ITD genes to include in report |
| `cnv_genes` | `Vec<String>` | CNV focal genes |
| `local_cn_genes` | `Vec<String>` | Genes showing local CN instead of focal |
| `deletion_genes` | `Vec<String>` | Genes to report deletions for |
| `duplication_genes` | `Vec<String>` | Genes to report duplications for |
| `snv_genes` | `Vec<String>` | SNV genes to include |
| `fusion_blacklist` | `Vec<(String, String)>` | Gene pairs to exclude from fusion report |
| `fusion_special_genes` | `HashMap` | Genes with special reporting rules (skip\_self\_fusion, bypass\_breakpoint\_filter) |
| `min_supporting_reads` | `i64` | Minimum reads for fusion reporting |
| `min_breakpoint_reads` | `i64` | Minimum breakpoint reads for fusion reporting |
| `columns` | `Vec<String>` | Explicit column order (auto-generated from gene lists when empty) |

### ITD Config (embedded in PipelineConfig `genes.itd`)

Each entry maps a gene name to `ItdRegion`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chrom` | `String` | -- | Chromosome name |
| `start` | `usize` | -- | Region start |
| `end` | `usize` | -- | Region end |
| `min_length` | `isize` | 3 | Minimum ITD length in bp |
| `min_frequency` | `f64` | 0.0 | Minimum allelic frequency |

## Threshold Defaults

| Module | Parameter | Default |
|--------|-----------|---------|
| Fusions | `default_margin` | 5000 bp |
| Fusions | `min_anchor` | 500 bp |
| Fusions | `min_supporting_reads` | 3 |
| Fusions | `min_breakpoint_reads` | 3 |
| SNV | `min_coverage` | 5.0x |
| SNV | `min_allele_freq` | 0.3 |
| SNV | `homozygosity_vaf` | 0.8 |
| SNV | `phasing_ratio` | 10 |
| CNV | `gene_margin` | 50000 bp |
| CNV | `min_supporting_reads` | 3 |
| CNV | `breakpoint_margin` | 50 bp |
| CNV | `min_anchor` | 500 bp |
| MAF | `min_depth` | 20 |
| Coverage | `bin_size` | 1000000 (1 Mb) |
| Karyotype | `spread_warning` | 0.08 |
| Karyotype | `level_tolerance` | 0.3 |
| Karyotype | `maf_peak_diploid` | 0.4 |
| Karyotype | `min_maf_sites` | 50 |
| Breakpoint Consensus | `flank_size` | 200 bp |
| Breakpoint Consensus | `min_reads` | 3 |
| Breakpoint Consensus | `min_coverage` | 3 |
| Breakpoint Consensus | `cluster_tolerance` | 30 bp |

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| `noodles` | 0.104 | BAM/CRAM/SAM/FASTA I/O |
| `clap` | 4.4 | CLI parsing with derive macros |
| `serde` / `serde_json` | 1.0 | JSON config and output serialization |
| `schemars` | 1 | JSON Schema generation from Rust types |
| `jsonschema` | 0.28 | JSON Schema validation |
| `bitvec` | 1 | Memory-efficient repeat masking |
| `tikv-jemallocator` | 0.6 | Performance allocator |
| `anyhow` / `thiserror` | 1.0 | Error handling |
| `log` / `env_logger` | 0.4 / 0.11 | Logging |
| `indexmap` | 2.0 | Ordered HashMap (preserves insertion order) |

## Output Formats

### Unified JSON (`{prefix}.result.json`)

Single JSON file containing all analysis results via `UnifiedOutput`. Fields are omitted when their analysis was not run (e.g., standalone `nasvar snv` only populates the `snv` field).

| Field | Type | Content |
|-------|------|---------|
| `version` | `String` | Pipeline version |
| `timestamp` | `String` | ISO 8601 timestamp |
| `metadata` | `SequencingMetaData` | BAM header metadata (instrument, run, flowcell) |
| `qc` | `QcOutput` | On-target nt, reads, mean coverage, total aligned reads |
| `fusions` | `FusionsOutput` | Fusion events, spike-in controls |
| `karyotype` | `KaryotypeOutput` | Karyotype string, ISCN, per-arm CN, warnings |
| `cnv` | `CnvOutput` | Per-gene CN, deletions, duplications |
| `snv` | `SnvOutput` | Per-gene variants with CDS/AA changes |
| `itd` | `ItdOutput` | ITD events with position, length, frequency |
| `breakpoint_consensus` | `BreakpointConsensusOutput` | Consensus sequences per fusion breakpoint |

### Other Outputs

| File | Format | Content |
|------|--------|---------|
| `{prefix}.breakpoints.fa` | FASTA | Consensus sequences spanning fusion breakpoints |
| `{prefix}.coverage.tsv` | TSV | chromosome, start, end, n\_reads, gc\_content |
| `{prefix}.coverage.gc_adjusted.tsv` | TSV | Same as above, after GC bias correction |
| `{prefix}.maf` | TSV | chromosome, position, ref\_count, alt\_count |
| `{prefix}.report.md` | Markdown | Report with summary tables |
| `{prefix}.report.html` | HTML | Report (embeds karyotype SVG plots) |
| `{prefix}.karyotype.svg` | SVG | Coverage plot before GC correction |
| `{prefix}.karyotype.gc_corrected.svg` | SVG | Coverage plot after GC correction |
| `{prefix}.karyotype_baf.gc_corrected.svg` | SVG | Combined karyotype + B-allele frequency plot |
| `{prefix}.gc_vs_coverage.svg` | SVG | GC content vs. coverage (diagnostic) |
| `{prefix}.gc_vs_coverage.gc_corrected.svg` | SVG | GC content vs. coverage after correction |

## Error Handling

- **`anyhow::Result`** is used in library functions for error propagation with context
- **`thiserror::Error`** is available for custom error types
- The CLI (`src/bin/nasvar.rs`) catches errors at the top level, logs them with `error!()`, and returns (no panic)
- `check_output_paths()` validates output directories exist and files don't already exist (unless `--force`)
- `require_index()` validates BAM/CRAM index availability before region queries

### Common Error Scenarios

| Error | Cause | Resolution |
|-------|-------|------------|
| Index file not found | Missing `.bai` or `.crai` | Run `samtools index` |
| CRAM decoding failure | Missing `--ref-fasta` | Provide reference FASTA |
| Output file exists | Previous run output present | Use `--force` / `-f` |
| JSON parse error | Malformed config file | Validate JSON syntax and field names |

## Testing

Unit tests (`#[test]` functions) cover:

- **Config parsing** (`src/config.rs`): default values, fusion margin lookup, blacklist pair matching, one-sided gene detection
- **Sequence alignment** (`src/utils/align.rs`): identical sequences, mismatches, insertions, deletions, homopolymers, affine gap behavior, empty sequences
- **Karyotype string generation** (`src/karyotype/mod.rs`): normal karyotypes, aneuploidies, ISCN formatting
- **Contig mapping** (`src/utils/contig.rs`): chromosome name normalization
- **Fusion detection** (`src/var/fusions.rs`): breakpoint clustering, ambiguity detection

Run tests with:
```bash
cargo test
```

## Debugging

### Per-Stage Execution

Run individual subcommands instead of `pipeline` to isolate issues:

```bash
nasvar coverage --bam sample.bam --repeats repeats.bed --out-prefix debug
# Inspect debug.coverage.tsv
nasvar karyotype --coverage debug.coverage.tsv --out-prefix debug --config pipeline.json --reference reference.json
# Inspect debug.result.json
```

### Verbose Logging

```bash
RUST_LOG=debug nasvar pipeline ...
```

The pipeline logs stage timing and per-stage progress at `info` level. Set `RUST_LOG=debug` for record-level details.

### Error Handling Style

Library functions return `anyhow::Result`. The CLI catches errors and logs them rather than panicking:

```rust
if let Err(e) = call_snvs(...) {
    error!("Error calling SNVs: {}", e);
}
```

### Configuration Pattern

All configurable parameters use serde `#[serde(default = "...")]` attributes so that JSON configs only need to specify non-default values. Each threshold struct implements `Default`.

### Logging

Uses `env_logger` initialized in `main()` with default level `Info`, no timestamps, no module paths:

```rust
env_logger::Builder::from_default_env()
    .filter_level(log::LevelFilter::Info)
    .format_timestamp(None)
    .format_module_path(false)
    .init();
```

Override with `RUST_LOG` environment variable.
