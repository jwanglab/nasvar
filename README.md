# nasvar = **N**anopore **A**daptive **S**ample **Var**iants

A variant calling pipeline for long-read sequencing with adaptive sampling, with a focus on human cancer genomics (tumor only). nasvar detects SNVs, CNVs, gene fusions, internal tandem duplications (ITDs), and infers karyotypes from aligned BAM or CRAM files.

## Features

- **Gene fusion/rearrangement detection** with breakpoint 
- **Fusion breakpoint consensus** sequences via iterative template-based polishing
- **Karyotype inference** at chromosome/arm level with GC bias correction
- Targeted **SNV calling** in enriched regions
- **CNV detection** in and around enriched regions, with GC bias correction and adjustment for inferred tumor fraction
- **ITD calling** (and other small indels)
- **BAM and CRAM support** with automatic format detection (WARNING: CRAM file processing will be MUCH slower than the corresponding BAM file)
- **Output**: JSON with schema, configurable markdown and HTML reports

## Installation

```bash
git clone https://github.com/jwanglab/nasvar
cd nasvar
cargo build --release
```

Binaries are placed in `target/release/`:
- `nasvar` -- main pipeline
- `aggregate` -- multi-sample reporting

## Quick Start

### Run the full pipeline

```bash
nasvar pipeline \
  sample.bam \
  repeats.bed \
  enriched.bed \
  maf_sites.txt \
  targets.bed \
  reference.fa \
  genes.gff3 \
  output/sample1 \
  --config pipeline.json \
  --reference reference.json
```

This produces:
- `output/sample1.result.json` -- unified JSON with all variant calls, QC, and metadata
- `output/sample1.coverage.tsv` -- binned read depth
- `output/sample1.coverage.gc_adjusted.tsv` -- GC-bias corrected binned coverage
- `output/sample1.maf` -- minor allele frequencies
- `output/sample1.breakpoints.fa` -- fusion breakpoint consensus sequences (if fusions found and BAM is indexed)
- `output/sample1.karyotype.gc_corrected.svg` -- karyotype coverage plot (GC-corrected)
- `output/sample1.karyotype_baf.gc_corrected.svg` -- combined karyotype + BAF plot
- `output/sample1.gc_vs_coverage.svg` -- GC bias diagnostic plot

### Required input files

| File | Description |
|------|-------------|
| BAM/CRAM | Sorted and indexed alignment file (`.bai` or `.crai` required) |
| Reference FASTA | Indexed genome reference (`.fai` required) |
| GFF3 | Gene annotations |
| Pipeline config | JSON file specifying genes and thresholds |
| Reference config | JSON file with centromere and PAR coordinates |
| Enrichment BED | The file used for adaptive sampling |
| Targets BED | Only the enriched genes of interest (no padding) |
| Repeats BED | Repetitive regions for coverage masking |
| Sites file | SNP positions for MAF (format: `CHR POS REF ALT`, tab-separated) |

## Subcommands

### `nasvar pipeline`

Runs all analyses in optimized order (parsing BAM/CRAM as few times as possible). Positional arguments:

```
nasvar pipeline <BAM> <REPEATS> <ENRICHED> <SITES> <TARGETS> <FASTA> <GFF> <OUT_PREFIX> \
  --config <JSON> --reference <JSON> [--blast-ratio <FLOAT>] [-f]
```

### `nasvar snv`

Call SNVs in specified genes.

```
nasvar snv --bam <BAM> --fasta <FASTA> --gff <GFF> --out-prefix <PREFIX> --config <JSON> \
  [--genes GENE1,GENE2] [-f]
```

### `nasvar cnv`

Call copy number variants.

```
nasvar cnv --bam <BAM> --targets <BED> --out-prefix <PREFIX> --config <JSON> --reference <JSON> \
  [--ref-fasta <FASTA>] [--coverage <TSV>] [--blast-ratio <FLOAT>] [-f]
```

### `nasvar fusions`

Detect gene fusions. The `--repeats` flag is accepted for compatibility but does not affect fusion calling; ambiguity filtering is now data-driven.

```
nasvar fusions --bam <BAM> --targets <BED> --repeats <BED> --out-prefix <PREFIX> --config <JSON> \
  [--ref-fasta <FASTA>] [--one-sided <FILE>] [--gff <GFF>] [-f]
```

### `nasvar karyotype`

Infer karyotype from pre-computed coverage.

```
nasvar karyotype --coverage <TSV> --out-prefix <PREFIX> --config <JSON> --reference <JSON> \
  [--maf <FILE>] [-f]
```

### `nasvar coverage`

Calculate binned read depth.

```
nasvar coverage --bam <BAM> --repeats <BED> --out-prefix <PREFIX> \
  [--ref-fasta <FASTA>] [-f]
```

### `nasvar maf`

Calculate minor allele frequencies at specified sites.

```
nasvar maf --bam <BAM> --enriched <BED> --sites <FILE> --out-prefix <PREFIX> \
  [--ref-fasta <FASTA>] [-f]
```

### `nasvar itd`

Call internal tandem duplications.

```
nasvar itd --bam <BAM> --itd-config <JSON> --out-prefix <PREFIX> \
  [--ref-fasta <FASTA>] [-f]
```

### `nasvar breakpoints`

Build consensus sequences for fusion breakpoints. Reads fusions from an existing `{prefix}.result.json` and updates it with breakpoint consensus data.

```
nasvar breakpoints --bam <BAM> --out-prefix <PREFIX> \
  [--ref-fasta <FASTA>] [--flank-size <INT>] [--min-reads <INT>] [--min-coverage <INT>] [-f]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--flank-size` | 200 | Flanking bases on each side of the breakpoint |
| `--min-reads` | 3 | Minimum reads required for consensus |
| `--min-coverage` | 3 | Minimum per-position coverage to call a base |

### `nasvar report`

Generate a report from existing unified output.

```
nasvar report --out-prefix <PREFIX> [--config <JSON>] [--targets <BED>]
```

### `nasvar schema`

Print the JSON Schema for the unified output format.

```
nasvar schema [--output <FILE>]
```

## Aggregate Tool

Combine multiple sample outputs into a single TSV:

```
aggregate <VERSION> [--config <JSON>] <DIR1> [DIR2] [DIR3] ...
```

The tool scans each directory for `{sample}.result.json` and produces a tab-separated report to stdout with columns for karyotype, fusions, CNVs, SNVs, ITDs, and QC metrics. When `--config` is omitted, built-in defaults are used.

## Configuration

### Pipeline Config (`pipeline.json`)

[Example](https://github.com/jwanglab/nasvar/blob/main/peds_leukemia_config.json)

All threshold fields have sensible defaults and can be omitted.

### Reference Config (`reference.json`)

[T2T/CHM13v2.0](https://github.com/jwanglab/nasvar/blob/main/T2T-CHM13v2.0_reference.json)

Provide centromere coordinates for all chromosomes (chr1-chr22, chrX, chrY).

## Common Options

| Option | Description |
|--------|-------------|
| `-f`, `--force` | Overwrite existing output files |
| `--blast-ratio <FLOAT>` | Tumor fraction (0.0-1.0) for CN adjustment. If omitted in pipeline mode, estimated from MAF data. |
| `--genes <LIST>` | Comma-separated gene list (SNV subcommand only) |
| `--ref-fasta <PATH>` | Reference FASTA, required for CRAM input |

### Logging

Control log verbosity via the `RUST_LOG` environment variable:

```bash
RUST_LOG=debug nasvar pipeline ...   # verbose
RUST_LOG=warn nasvar pipeline ...    # warnings only
```

Default level is `info`.

## Output Files

To export full JSON output schema:
```
nasvar schema
```

### `{prefix}.result.json`

Unified JSON output containing all analysis results. Top-level fields (all optional except `version` and `timestamp`):

| Field | Content |
|-------|---------|
| `version` | Pipeline version |
| `timestamp` | Analysis timestamp (ISO 8601) |
| `metadata` | Sequencing run metadata from BAM headers |
| `qc` | On-target metrics and aligned read counts |
| `fusions` | Fusion events with breakpoints, supporting reads, spike-in controls |
| `karyotype` | Karyotype string, ISCN notation, per-arm copy numbers, warnings |
| `cnv` | Per-gene CNV calls with local/focal CN, deletions, duplications |
| `snv` | Per-gene SNV calls with CDS positions, amino acid changes, VAF |
| `itd` | ITD calls with position, length, allelic ratio |
| `breakpoint_consensus` | Consensus sequences over fusion breakpoints (when available) |

### `{prefix}.breakpoints.fa`

FASTA file with consensus sequences spanning each fusion breakpoint. Written when breakpoint consensus succeeds (requires indexed BAM and detected fusions).

### `{prefix}.coverage.tsv`

Tab-separated binned coverage: chromosome, start, end, read count, GC content.

### `{prefix}.coverage.gc_adjusted.tsv`

GC-bias corrected binned coverage, produced during karyotype inference.

### `{prefix}.maf`

Tab-separated allele counts at specified sites: chromosome, position, ref count, alt count.

### `{prefix}.report.md` / `{prefix}.report.html`

Report in Markdown and HTML formats, generated by `nasvar report`. Includes summary tables for fusions, karyotype, CNVs, SNVs, ITDs, and methodology details. The HTML report embeds karyotype SVG plots when available.

### Karyotype plots

Generated during karyotype inference (pipeline or standalone `karyotype` subcommand):

| File | Content |
|------|---------|
| `{prefix}.karyotype.svg` | Coverage plot before GC correction |
| `{prefix}.karyotype.gc_corrected.svg` | Coverage plot after GC correction |
| `{prefix}.karyotype_baf.gc_corrected.svg` | Combined karyotype + B-allele frequency plot |
| `{prefix}.gc_vs_coverage.svg` | GC content vs. coverage (diagnostic) |
| `{prefix}.gc_vs_coverage.gc_corrected.svg` | GC content vs. coverage after correction |

## Troubleshooting

### Index not found

```
Error: Index file not found for 'sample.bam'. Expected 'sample.bam.bai'.
```

Create the index with `samtools index sample.bam`.

### CRAM reference required

CRAM files require a reference FASTA for decoding. Pass `--ref-fasta reference.fa` for any subcommand that reads a CRAM file.

### Output file already exists

Use `--force` / `-f` to overwrite existing output files.
