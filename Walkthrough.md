# nasvar Walkthrough

This guide walks through setting up and running a complete nasvar analysis on a single sample, from obtaining reference data through pipeline execution and interpreting outputs.

## Prerequisites

- **Rust** (edition 2024): https://rustup.rs
- **samtools**: for indexing BAM and FASTA files
- **bcftools** (optional): only needed if building a custom MAF sites file
- **minimap2** (or similar): for aligning reads to the reference genome

## 1. Build and install nasvar

```bash
git clone https://github.com/jwanglab/nasvar
cd nasvar
cargo install --path .
```

## 2. Obtain reference data

nasvar is optimized and tested on the T2T-CHM13v2.0 human reference genome. All coordinates in the included configs use NCBI accession names (e.g. `NC_060925.1` for chr1).

### Reference genome (FASTA)

Download T2T-CHM13v2.0 from NCBI:

```bash
# Download the reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

# Index it
samtools faidx GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
```

### GFF3 gene annotations

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```

### Repeat annotations (BED)

A repeat mask BED is used to exclude repetitive regions from coverage calculations. We use a concatenation of RepeatMasker, centromere/satellite, telomere, and segmental duplication annotations from the [T2T CHM13 project](https://github.com/marbl/CHM13). The file should have at minimum columns: `chrom start end name`.

Download the individual annotation tracks and concatenate them:

```bash
# RepeatMasker
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed

# CenSat (centromere/satellite)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

# Telomere
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed

# Segmental duplications
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed

# Concatenate into a single mask
cat chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed \
    chm13v2.0_censat_v2.1.bed \
    chm13v2.0_telomere.bed \
    chm13v2.0_SD.bed \
    > repeats.bed
```

## 3. Prepare panel-specific files

These files define your assay's target regions and are specific to your adaptive sampling panel design.

### Enrichment BED (`enriched.bed`)

The BED file used for adaptive sampling -- the full set of enriched regions including any flanking/padding (we use 50Kbp on either side of every gene, then merge overlapping regions). This is used to determine which genomic regions should be queried for MAF sites. Format:

```
NC_060925.1	2521244	2993423	PRDM16
NC_060925.1	32937140	33212335	ZNF362
...
```

### Targets BED (`targets.bed`)

A tighter set of analysis regions -- the gene bodies of interest without padding. Used for fusion detection, CNV calling, and on-target QC. Same format as enrichment BED but with smaller intervals.

### MAF sites file (`maf_sites.tsv`)

A tab-separated file of common SNP positions used for minor allele frequency estimation (karyotype MAF-based ploidy inference, BAF plots). Format:

```
chr1	2521312	G	A
chr1	2521441	C	T
chr1	2521476	A	C
```

Columns: chromosome, position, reference allele, alternate allele.

**To generate from 1000 Genomes Project data on T2T:**

The [1000 Genomes Project variant calls on T2T-CHM13v2.0](https://github.com/JosephLalli/phasing_T2T) (3,202 samples, phased with SHAPEIT5) are hosted on the Human Pangenome S3 bucket. Download the whole-genome BCF and extract common biallelic SNPs:

```bash
# Download the whole-genome phased BCF (~65 GB)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz

# Filter to common biallelic SNPs (allele count >= 10):
bcftools view -G -H -O v -AA -a --known -U -v snps -c 10 \
  1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz \
  | awk '(length($4)==1 && length($5)==1){ print $1"\t"$2"\t"$4"\t"$5 }' \
  > maf_sites.tsv
```

Optionally intersect with your enrichment BED to restrict to on-target sites.

## 4. Configuration files

nasvar requires two JSON configuration files. Examples are included in the `config/` directory.

### Reference config (`config/T2T-CHM13v2.0_reference.json`)

Defines centromere coordinates and pseudoautosomal regions (PAR) for the reference genome. The included file covers T2T-CHM13v2.0 and should not need modification unless using a different reference.

```json
{
  "name": "T2T-CHM13v2.0",
  "centromeres": {
    "chr1": [121405145, 144008994],
    "chr2": [87403406, 98345164],
    ...
  },
  "par1": {
    "chrX": [0, 2394410],
    "chrY": [0, 2458320]
  },
  "par2": {
    "chrX": [153925834, 154259566],
    "chrY": [62122809, 62460029]
  }
}
```

Note: Acrocentric chromosomes (13, 14, 15, 21, 22) do not need centromere entries -- they are treated as single-arm chromosomes.

### Pipeline config (`config/peds_leukemia_config.json`)

Defines the genes, thresholds, and parameters for each analysis module. The included pediatric leukemia config is a complete working example. Key sections:

- **`genes.cnv`**: Focal CNV genes (e.g. CDKN2A, IKZF1, PAX5), deletion/duplication gene lists
- **`genes.itd`**: ITD target regions with genomic coordinates (e.g. FLT3, UBTF, NPM1)
- **`genes.snv`**: Pharmacogenomics variants (TPMT, NUDT15) and pathogenic variants, with transcript and CDS position definitions
- **`genes.fusions`**: Fusion partner genes, special margins, blacklists, one-sided/promiscuous partner rules
- **`thresholds`**: Calling parameters for each module (all have sensible defaults and can be omitted)

To adapt for a different assay, copy the example and modify the gene lists and coordinates. All threshold fields are optional.

## 5. Align reads

If you don't already have a sorted, indexed BAM:

```bash
# Align with minimap2
minimap2 -a -x map-ont -y --MD \
  GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  reads.fastq.gz \
  | samtools sort -o sample.T2T.sorted.bam

# Index
samtools index sample.T2T.sorted.bam
```

## 6. Run the pipeline

```bash
nasvar pipeline \
  sample.T2T.sorted.bam \
  repeats.bed \
  enriched.bed \
  maf_sites.tsv \
  targets.bed \
  GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  GCF_009914755.1_T2T-CHM13v2.0_genomic.gff \
  results/sample/sample \
  --config config/peds_leukemia_config.json \
  --reference config/T2T-CHM13v2.0_reference.json \
  --force
```

The output prefix (`results/sample/sample`) determines both the directory and the filename base. Create the output directory first:

```bash
mkdir -p results/sample
```

### Optional flags

| Flag | Description |
|------|-------------|
| `--gc-correction loess` | Use LOESS instead of linear regression for GC bias correction (better for non-linear bias) |
| `--gc-correction none` | Skip GC correction entirely |
| `--blast-ratio 0.5` | Override tumor fraction estimate (0.0-1.0); if omitted, estimated from the data |
| `-f` / `--force` | Overwrite existing output files |

### Verbose logging

```bash
RUST_LOG=debug nasvar pipeline ...
```

## 7. Output files

After a successful run, the output prefix directory will contain:

### Primary results

| File | Description |
|------|-------------|
| `sample.result.json` | Unified JSON with all variant calls, karyotype, QC, and metadata |
| `sample.coverage.tsv` | Binned genome-wide read depth (1 Mb windows) with GC content |
| `sample.coverage.gc_adjusted.tsv` | GC-bias corrected binned coverage |
| `sample.maf` | Per-site allele counts (ref/alt depth) |
| `sample.breakpoints.fa` | Fusion breakpoint consensus sequences (if fusions detected) |

### Plots

| File | Description |
|------|-------------|
| `sample.karyotype.svg` | Genome-wide coverage before GC correction |
| `sample.karyotype.gc_corrected.svg` | Genome-wide coverage after GC correction |
| `sample.karyotype_baf.gc_corrected.svg` | Combined karyotype + B-allele frequency |
| `sample.gc_vs_coverage.svg` | GC content vs coverage scatter (raw) with fitted regression |
| `sample.gc_vs_coverage.gc_corrected.svg` | GC content vs coverage scatter (after correction) |

### Reports

```bash
# Generate markdown and HTML reports from existing results
nasvar report --out-prefix results/sample/sample \
  --config config/peds_leukemia_config.json \
  --targets targets.bed
```

This produces `sample.report.md` and `sample.report.html`. The HTML report embeds karyotype SVG plots.

## 8. Interpreting results

### `result.json` structure

```json
{
  "version": "1.0.0",
  "timestamp": "2026-02-16T10:00:00Z",
  "metadata": { "instrument": "...", "run_id": "...", "flowcell": "..." },
  "qc": { "nt_on_target": 1000000, "reads_on_target": 500, "mean_coverage": 30.0 },
  "fusions": { "fusions": [...], "spike_in": [...] },
  "karyotype": {
    "karyotype_string": "46,XY",
    "per_arm_cn": { "1p": 2, "1q": 2, ... },
    "blast_ratio": 0.85,
    "warnings": []
  },
  "cnv": { "genes": { "CDKN2A": { "focal_cn": 0 }, ... }, "deletions": [...] },
  "snv": { "genes": { "TPMT": { ... } } },
  "itd": { "events": [...] },
  "breakpoint_consensus": { "breakpoints": [...] }
}
```

All sections except `version` and `timestamp` are present only if the corresponding analysis produced results.

### Diagnostic plots

- **`gc_vs_coverage.svg`**: Check the scatter pattern. If the fitted line (linear) doesn't follow the data well, re-run with `--gc-correction loess`. The corrected plot (`gc_vs_coverage.gc_corrected.svg`) should show a flattened relationship.
- **`karyotype.svg` vs `karyotype.gc_corrected.svg`**: Compare to verify GC correction improved the coverage uniformity across chromosome arms.

## 9. Running individual subcommands

The full pipeline runs all modules automatically, but you can also run stages independently:

```bash
# Coverage only
nasvar coverage --bam sample.bam --repeats repeats.bed --out-prefix results/sample/sample

# MAF only
nasvar maf --bam sample.bam --enriched enriched.bed --sites maf_sites.tsv --out-prefix results/sample/sample

# Karyotype from pre-computed coverage/MAF
nasvar karyotype \
  --coverage results/sample/sample.coverage.tsv \
  --maf results/sample/sample.maf \
  --out-prefix results/sample/sample \
  --config config/peds_leukemia_config.json \
  --reference config/T2T-CHM13v2.0_reference.json \
  --gc-correction loess

# Fusions only
nasvar fusions --bam sample.bam --targets targets.bed --repeats repeats.bed \
  --out-prefix results/sample/sample --config config/peds_leukemia_config.json \
  --gff GCF_009914755.1_T2T-CHM13v2.0_genomic.gff

# SNV only
nasvar snv --bam sample.bam --fasta GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --gff GCF_009914755.1_T2T-CHM13v2.0_genomic.gff \
  --out-prefix results/sample/sample --config config/peds_leukemia_config.json

# ITD only
nasvar itd --bam sample.bam --itd-config config/peds_leukemia_config.json \
  --out-prefix results/sample/sample
```

## 10. Multi-sample aggregation

After running multiple samples, aggregate results into a single TSV:

```bash
aggregate v1.0 --config config/peds_leukemia_config.json results/sample1/ results/sample2/ results/sample3/
```

This scans each directory for `*.result.json` files and produces a tab-separated summary to stdout with columns for karyotype, fusions, CNVs, SNVs, ITDs, and QC metrics.
