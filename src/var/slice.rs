//! Variant-region BAM slicing.
//!
//! After the main pipeline has run, this module produces a small BAM file
//! containing only the reads in regions of interest around called variants
//! (SNVs, ITDs, fusions), plus an accompanying .bai. The slice is intended
//! for inclusion in a portable "result package" that lets a downstream
//! viewer show alignments around each call without shipping the full BAM.
//!
//! Region policy (defaults):
//!   - 10 kb around every detected SNV mutation entry.
//!   - 10 kb around every reported ITD position.
//!   - For each fusion call: the entirety of both partner genes (annotated
//!     bounds from the GFF), padded by 50 kb on each side.
//!   - CNV breakpoints are intentionally not included.
//!
//! Regions are sorted and merged before querying so each read is visited
//! once. Reads that span multiple non-adjacent regions are still
//! deduplicated by (qname, ref_id, pos, flag) since indexed queries return
//! a record once per region the read overlaps.

use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::{anyhow, Result};
use log::{info, warn};

use noodles::bam;
use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::csi::binning_index::{
    index::reference_sequence::{bin::Chunk, index::LinearIndex},
    Indexer as CsiIndexer,
};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::header::record::value::{
    map::{program::tag as pg_tag, Program},
    Map,
};

use crate::config::Contig;
use crate::output::{FusionBreakpoint, GeneInfo, UnifiedOutput};
use crate::utils::annotation::{GeneBounds, PartnerGeneIndex};
use crate::utils::bed::BedRegion;
use std::collections::HashMap;

pub const DEFAULT_SNV_PAD: u32 = 10_000;
pub const DEFAULT_ITD_PAD: u32 = 10_000;
pub const DEFAULT_FUSION_PAD: u32 = 50_000;
/// One-sided fusion: pad applied around the *unnamed* (intergenic) side's
/// breakpoint position. The named partner still gets full gene bounds +
/// `DEFAULT_FUSION_PAD` like a two-sided fusion.
pub const FUSION_INTERGENIC_PAD: u32 = 10_000;
/// CNV/SV: pad applied around each breakpoint (start and end separately).
/// The full SV span can be megabases for arm-level events, so we only
/// capture reads flanking the actual junction, not the entire deleted /
/// duplicated interval.
pub const DEFAULT_CNV_BREAKPOINT_PAD: u32 = 10_000;

/// Build a name → bounds map from a fusion-targets BED (the file used by
/// fusion calling). Multiple rows for the same gene name are unioned: the
/// resulting bounds span min(start)..max(end) across all rows. Chrom is
/// taken from the first row seen for each name; rows with a conflicting
/// chrom are skipped with a warning.
fn build_targets_index(targets: &[BedRegion]) -> HashMap<String, GeneBounds> {
    let mut map: HashMap<String, GeneBounds> = HashMap::new();
    for r in targets {
        if r.name.is_empty() { continue; }
        match map.get_mut(&r.name) {
            Some(b) => {
                if b.chrom != r.segment {
                    warn!(
                        "[slice] targets BED has {} on both {} and {} — keeping the first",
                        r.name, b.chrom, r.segment
                    );
                    continue;
                }
                if r.start < b.start { b.start = r.start; }
                if r.end > b.end { b.end = r.end; }
            }
            None => {
                map.insert(r.name.clone(), GeneBounds {
                    chrom: r.segment.clone(),
                    start: r.start,
                    end: r.end,
                });
            }
        }
    }
    map
}

/// Compute non-overlapping (chrom, start, end) windows from a unified
/// pipeline output. Coordinates are 1-based inclusive (suitable for a
/// `noodles::core::Region`). Returns regions sorted and merged.
///
/// `fusion_targets` is the same BED that fusion calling consumed — used as
/// the primary source for fusion partner gene bounds, since
/// `FusionEvent.gene1/2.name` is keyed off this BED's `name` column. We
/// fall back to the GFF only when a partner isn't in the BED.
pub fn compute_slice_regions(
    unified: &UnifiedOutput,
    gff_path: &str,
    fusion_targets: Option<&[BedRegion]>,
    snv_pad: u32,
    itd_pad: u32,
    fusion_pad: u32,
    cnv_breakpoint_pad: u32,
) -> Result<Vec<(String, u32, u32)>> {
    // Collect every gene name we'll need so a single GFF pass picks them up.
    let mut wanted: HashSet<String> = HashSet::new();
    if let Some(snv) = &unified.snv {
        for g in snv.genes.keys() { wanted.insert(g.clone()); }
    }
    if let Some(itd) = &unified.itd {
        for g in itd.genes.keys() { wanted.insert(g.clone()); }
    }
    if let Some(fus) = &unified.fusions {
        for ev in &fus.fusions {
            for n in [&ev.gene1.name, &ev.gene2.name] {
                if !n.is_empty() { wanted.insert(n.clone()); }
            }
        }
    }

    let gene_index = if wanted.is_empty() {
        PartnerGeneIndex::new()
    } else {
        match PartnerGeneIndex::load_from_gff(gff_path, &wanted) {
            Ok(idx) => idx,
            Err(e) => {
                warn!("[slice] couldn't load gene bounds from GFF ({}); fusion windows will fall back to breakpoint pads", e);
                PartnerGeneIndex::new()
            }
        }
    };

    // Targets-BED-derived bounds are the primary source for fusions because
    // FusionEvent.gene1/2.name is keyed off the targets BED's `name` column,
    // not the GFF's `Name=`. The GFF lookup is a backup for genes not in
    // the BED.
    let targets_index = fusion_targets.map(build_targets_index).unwrap_or_default();

    let mut regions: Vec<(String, u32, u32)> = Vec::new();
    let pad = |chrom: &str, pos: u32, p: u32, out: &mut Vec<(String, u32, u32)>| {
        out.push((chrom.to_string(), pos.saturating_sub(p), pos.saturating_add(p)));
    };

    // SNVs: 10 kb around each detected mutation. Mutation keys are
    // formatted "<pos><ref>><alt>" (e.g. "1234567A>G"); the leading numeric
    // prefix is the position. Chrom comes from looking the gene up in the
    // GFF — SnvGeneResult itself doesn't carry it.
    if let Some(snv) = &unified.snv {
        for (gene, result) in &snv.genes {
            let Some(bounds) = gene_index.get_gene(gene) else {
                warn!("[slice] no GFF entry for SNV gene {} — skipping its mutations", gene);
                continue;
            };
            for key in result.mutations.keys() {
                if let Some(pos) = parse_leading_pos(key) {
                    pad(&bounds.chrom, pos, snv_pad, &mut regions);
                }
            }
        }
    }

    // ITDs: 10 kb around each ItdEvent.position.
    if let Some(itd) = &unified.itd {
        for (gene, events) in &itd.genes {
            let Some(bounds) = gene_index.get_gene(gene) else {
                warn!("[slice] no GFF entry for ITD gene {} — skipping its events", gene);
                continue;
            };
            for ev in events {
                if ev.position < 0 { continue; }
                pad(&bounds.chrom, ev.position as u32, itd_pad, &mut regions);
            }
        }
    }

    // CNVs: pad each SV breakpoint (start and end separately) by
    // `cnv_breakpoint_pad`. We deliberately *don't* slice the full SV span
    // -- arm-level deletions can be hundreds of Mb and would explode the
    // slice -- just enough flanking reads at each junction to visualize
    // the breakpoints in the viewer.
    if let Some(cnv) = &unified.cnv {
        for gene_info in cnv.genes.values() {
            let dels = gene_info.deletions.iter().flatten();
            let dups = gene_info.duplications.iter().flatten();
            for sv in dels.chain(dups) {
                let Some(chrom) = sv.chrom.as_deref() else { continue; };
                if sv.start >= 0 {
                    pad(chrom, sv.start as u32, cnv_breakpoint_pad, &mut regions);
                }
                if sv.end >= 0 && sv.end != sv.start {
                    pad(chrom, sv.end as u32, cnv_breakpoint_pad, &mut regions);
                }
            }
        }
    }

    // Fusions: each side independently. If the side has a gene name we can
    // look up in the GFF, take its full bounds + 50 kb. Otherwise (intergenic
    // one-sided side, or a name we couldn't resolve), pad each per-breakpoint
    // position by 10 kb on this side. The event-level GeneInfo.pos is only
    // representative and may not cover all breakpoints, so we always read
    // positions from `ev.breakpoints[]` for the unnamed/fallback path.
    if let Some(fus) = &unified.fusions {
        for ev in &fus.fusions {
            fusion_side_regions(
                &targets_index, &gene_index, &ev.gene1, &ev.breakpoints,
                /* is_left= */ true,
                fusion_pad, FUSION_INTERGENIC_PAD, &mut regions,
            );
            fusion_side_regions(
                &targets_index, &gene_index, &ev.gene2, &ev.breakpoints,
                /* is_left= */ false,
                fusion_pad, FUSION_INTERGENIC_PAD, &mut regions,
            );
        }
    }

    Ok(merge_regions(regions))
}

/// Append regions for one side of a fusion event. `is_left == true` reads
/// gene0_* fields from each breakpoint (corresponding to ev.gene1);
/// `is_left == false` reads gene1_* (corresponding to ev.gene2).
///
/// Lookup order for a named partner: targets BED → GFF → breakpoint pads.
/// The targets BED is the primary source because fusion calling itself
/// derived the gene name from this BED's `name` column, so a name there
/// is guaranteed to match.
fn fusion_side_regions(
    targets_index: &HashMap<String, GeneBounds>,
    gene_index: &PartnerGeneIndex,
    side: &GeneInfo,
    breakpoints: &[FusionBreakpoint],
    is_left: bool,
    full_gene_pad: u32,
    intergenic_pad: u32,
    out: &mut Vec<(String, u32, u32)>,
) {
    if !side.name.is_empty() {
        let from_targets = targets_index.get(&side.name);
        let from_gff = if from_targets.is_none() { gene_index.get_gene(&side.name) } else { None };
        if let Some(b) = from_targets.or(from_gff) {
            out.push((
                b.chrom.clone(),
                b.start.saturating_sub(full_gene_pad),
                b.end.saturating_add(full_gene_pad),
            ));
            return;
        }
        warn!(
            "[slice] no targets-BED or GFF entry for fusion partner {} — falling back to {} kb breakpoint pads",
            side.name, intergenic_pad / 1000,
        );
    }
    // Intergenic side, or a named gene we couldn't resolve: pad each
    // breakpoint on this side by `intergenic_pad`.
    for bp in breakpoints {
        let (chr, pos) = if is_left {
            (bp.gene0_chr.as_str(), bp.gene0_pos)
        } else {
            (bp.gene1_chr.as_str(), bp.gene1_pos)
        };
        if !chr.is_empty() {
            out.push((
                chr.to_string(),
                pos.saturating_sub(intergenic_pad),
                pos.saturating_add(intergenic_pad),
            ));
        }
    }
}

fn parse_leading_pos(key: &str) -> Option<u32> {
    let n: String = key.chars().take_while(|c| c.is_ascii_digit()).collect();
    if n.is_empty() { None } else { n.parse().ok() }
}

fn merge_regions(mut regions: Vec<(String, u32, u32)>) -> Vec<(String, u32, u32)> {
    regions.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<(String, u32, u32)> = Vec::new();
    for r in regions {
        match out.last_mut() {
            Some(last) if last.0 == r.0 && r.1 <= last.2 => {
                if r.2 > last.2 { last.2 = r.2; }
            }
            _ => out.push(r),
        }
    }
    out
}

/// Write a slice BAM containing only reads overlapping `regions`, plus a
/// BAI alongside it. Reads are deduplicated across overlapping queries and
/// across multiple source BAMs by (ref_id, pos, flag, qname).
///
/// `src_indices` is an optional in-memory BAI per source BAM. When the
/// pipeline built the index inline during pass 1, passing it directly
/// here is preferable to the disk round-trip (write to `<path>.bai` then
/// read back). None elements fall back to reading `<path>.bai`
/// from disk.
pub fn write_slice_bam(
    src_bams: &[String],
    src_indices: &[Option<bam::bai::Index>],
    regions: &[(String, u32, u32)],
    out_bam: &str,
    contigs: &[Contig],
) -> Result<usize> {
    if regions.is_empty() {
        info!("[slice] no regions of interest — skipping slice BAM");
        return Ok(0);
    }
    if src_bams.is_empty() {
        return Err(anyhow!("write_slice_bam: no source BAMs supplied"));
    }

    // Sources: one indexed reader per BAM. The first BAM's header is
    // canonical for the output — AlignmentInput::open_multi already verified
    // that every member shares the same reference sequences.
    let mut readers: Vec<bam::io::IndexedReader<bgzf::io::Reader<File>>> =
        Vec::with_capacity(src_bams.len());
    let mut canonical_header: Option<noodles::sam::Header> = None;
    for (i, path) in src_bams.iter().enumerate() {
        let index = match src_indices.get(i).and_then(|o| o.clone()) {
            Some(idx) => idx,
            None => {
                let bai_path = format!("{}.bai", path);
                if !Path::new(&bai_path).exists() {
                    return Err(anyhow!(
                        "slice requires {}; pipeline should have built one inline",
                        bai_path
                    ));
                }
                bam::bai::fs::read(&bai_path)
                    .map_err(|e| anyhow!("read bai {}: {}", bai_path, e))?
            }
        };
        let mut reader = bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(path)
            .map_err(|e| anyhow!("open indexed bam {}: {}", path, e))?;
        let header = reader.read_header()
            .map_err(|e| anyhow!("read header for {}: {}", path, e))?;
        if canonical_header.is_none() {
            canonical_header = Some(header);
        }
        readers.push(reader);
    }
    let mut header = canonical_header.unwrap();

    // Stamp our own @PG record into the slice's header chain so the slice
    // is self-documenting: which BAM(s) were sliced, and by what version of
    // nasvar.
    append_nasvar_pg(&mut header, src_bams, regions.len());

    // Translate region chrom names to the BAM's naming convention before
    // querying. Variant outputs source chroms from a mix of places (targets
    // BED, GFF, BAM directly), so a region might have `chr19` when the BAM
    // uses `NC_060943.1` (or vice versa)
    let bam_refs: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).to_string())
        .collect();
    let contig_mapper = crate::bam::ContigMapper::from_contigs_and_refs(contigs, &bam_refs);

    // Sink: BAM writer over a bgzf writer.
    let out_file = File::create(out_bam)
        .map_err(|e| anyhow!("create {}: {}", out_bam, e))?;
    let bgzf_writer = bgzf::io::Writer::new(out_file);
    let mut writer = bam::io::Writer::from(bgzf_writer);
    writer.write_header(&header)
        .map_err(|e| anyhow!("write header: {}", e))?;

    let mut seen: HashSet<(i32, i32, u16, Vec<u8>)> = HashSet::new();
    let mut written = 0usize;

    // Translate region chroms to BAM convention up front, attach the
    // header ref_id, and resort by (ref_id, start). compute_slice_regions
    // sorts by chrom *string*, which doesn't match BAM ref_id ordering
    // when variants mix naming conventions (e.g. CNV outputs use "chr7"
    // while fusion outputs use "NC_060935.1") -- writing in string order
    // then breaks the BAI build, which needs coordinate-sorted input.
    let mut sorted_regions: Vec<(String, u32, u32, i32)> = regions
        .iter()
        .map(|(c, s, e)| {
            let bam_chrom = contig_mapper.to_bam_name(c);
            let ref_id = bam_refs.iter()
                .position(|r| r == &bam_chrom)
                .map(|i| i as i32)
                .unwrap_or(-1);
            (bam_chrom, *s, *e, ref_id)
        })
        .collect();
    sorted_regions.sort_by(|a, b| a.3.cmp(&b.3).then(a.1.cmp(&b.1)));

    for (bam_chrom, start, end, _ref_id) in &sorted_regions {
        // Region coordinates are 1-based inclusive for Region::new.
        let s_pos = Position::try_from((*start).max(1) as usize)
            .map_err(|e| anyhow!("region start: {}", e))?;
        let e_pos = Position::try_from((*end).max((*start).max(1)) as usize)
            .map_err(|e| anyhow!("region end: {}", e))?;
        let region = Region::new(bam_chrom.as_bytes().to_vec(), s_pos..=e_pos);

        // Pull records from every reader into one bucket, then sort by
        // (ref_id, pos) so the output stays coordinate-sorted (required for
        // the post-pass BAI build). The single-BAM case is essentially free:
        // records come back already in pos order and Rust's sort is adaptive.
        let mut buf: Vec<(i32, i32, u16, Vec<u8>, bam::Record)> = Vec::new();
        for reader in &mut readers {
            let query = match reader.query(&header, &region) {
                Ok(q) => q,
                Err(e) => {
                    warn!(
                        "[slice] query failed for {}:{}-{}: {} — skipping",
                        bam_chrom, start, end, e
                    );
                    continue;
                }
            };
            for result in query.records() {
                let record = result.map_err(|e| anyhow!("read record: {}", e))?;

                // Mirrors decode_bam_record() in src/input.rs for the few
                // fields we need for dedup + sorting. CIGAR/seq/quality pass
                // through verbatim via write_record(&header, &record).
                let ref_id: i32 = match record.reference_sequence_id() {
                    Some(Ok(id)) => id as i32,
                    Some(Err(e)) => return Err(anyhow!("ref id: {}", e)),
                    None => -1,
                };
                let pos: i32 = match record.alignment_start() {
                    Some(Ok(p)) => (p.get() as i32) - 1, // 1-based → 0-based
                    Some(Err(e)) => return Err(anyhow!("pos: {}", e)),
                    None => -1,
                };
                let flag: u16 = record.flags().bits();
                let qname: Vec<u8> = record
                    .name()
                    .map(|n| <_ as AsRef<[u8]>>::as_ref(n).to_vec())
                    .unwrap_or_default();
                buf.push((ref_id, pos, flag, qname, record));
            }
        }
        buf.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        for (ref_id, pos, flag, qname, record) in buf {
            let key = (ref_id, pos, flag, qname);
            if !seen.insert(key) {
                continue;
            }
            writer.write_record(&header, &record)
                .map_err(|e| anyhow!("write record: {}", e))?;
            written += 1;
        }
    }

    // Flush bgzf and EOF-mark the BAM so the second-pass reader doesn't see
    // a truncated file.
    drop(writer);

    // Second pass: build the BAI by reading back the freshly-written BAM and
    // feeding records into the CSI indexer.
    build_bai(out_bam, header.reference_sequences().len())?;

    info!(
        "[slice] wrote {} reads across {} merged regions from {} BAM(s) → {} (+ {}.bai)",
        written, regions.len(), src_bams.len(), out_bam, out_bam
    );
    Ok(written)
}

/// Convenience: compute regions from `unified`, then write the slice BAM
/// (and BAI) to `<out_prefix>.slice.bam`. When `ref_path` is supplied, also
/// emit `<out_prefix>.slice.fa` + `.fai` containing just the bases that
/// cover the slice regions, so a downstream viewer can render mismatch
/// ticks without needing the full reference.
/// Append a `nasvar`-authored `@PG` entry to the slice header so the
/// produced BAM records which input(s) it was sliced from and which version
/// did the slicing. `Programs::add` auto-chains the new entry under the
/// existing leaf via `PP:`, so the source's program history is preserved
/// intact above ours.
fn append_nasvar_pg(
    header: &mut noodles::sam::Header,
    src_bams: &[String],
    region_count: usize,
) {
    let cl = format!(
        "nasvar slice --regions {} --inputs {}",
        region_count,
        src_bams.join(","),
    );
    let map = match Map::<Program>::builder()
        .insert(pg_tag::NAME, "nasvar")
        .insert(pg_tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(pg_tag::COMMAND_LINE, cl)
        .build()
    {
        Ok(m) => m,
        Err(e) => {
            warn!("[slice] could not build @PG map ({}); slice header will omit nasvar provenance", e);
            return;
        }
    };
    if let Err(e) = header.programs_mut().add("nasvar", map) {
        warn!("[slice] could not append @PG ({}); slice header will omit nasvar provenance", e);
    }
}

pub fn dump_variant_slice(
    unified: &UnifiedOutput,
    src_bams: &[String],
    src_indices: &[Option<bam::bai::Index>],
    ref_path: Option<&str>,
    gff_path: &str,
    fusion_targets: Option<&[BedRegion]>,
    out_prefix: &str,
    contigs: &[Contig],
) -> Result<Option<(String, String)>> {
    let regions = compute_slice_regions(
        unified, gff_path, fusion_targets,
        DEFAULT_SNV_PAD, DEFAULT_ITD_PAD, DEFAULT_FUSION_PAD,
        DEFAULT_CNV_BREAKPOINT_PAD,
    )?;
    if regions.is_empty() {
        info!("[slice] no variant regions to slice");
        return Ok(None);
    }
    let total_bp: u64 = regions.iter().map(|r| (r.2 - r.1) as u64).sum();
    info!("[slice] {} merged regions, total span {} bp", regions.len(), total_bp);

    let bam_out = format!("{}.slice.bam", out_prefix);
    write_slice_bam(src_bams, src_indices, &regions, &bam_out, contigs)?;
    let bai_out = format!("{}.bai", bam_out);

    // Best-effort: a slice FASTA is useful but not load-bearing for the BAM.
    if let Some(rp) = ref_path {
        let fa_out = format!("{}.slice.fa.gz", out_prefix);
        if let Err(e) = write_slice_fasta(&regions, rp, &fa_out, contigs) {
            warn!("[slice] fasta dump failed: {}", e);
        }
    }

    Ok(Some((bam_out, bai_out)))
}

/// Write a sparse, gzipped FASTA containing just the slice regions, one
/// record per merged region named `<chrom>:<start>-<end>`. Coordinates
/// match the slice BAM (1-based inclusive). Also writes `<out_fa>.fai`
/// with offsets into the *uncompressed* stream (samtools convention) so
/// a viewer can decode-once-then-slice without needing a .gzi index.
///
/// `out_fa` must end in `.fa.gz` -- the entire output is a single gzip
/// member (not bgzipped). For typical slice sizes (sub-MB) whole-file
/// decompression is fast.
pub fn write_slice_fasta(
    regions: &[(String, u32, u32)],
    ref_path: &str,
    out_fa: &str,
    contigs: &[Contig],
) -> Result<usize> {
    use noodles::fasta;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use crate::bam::ContigMapper;

    if regions.is_empty() {
        return Ok(0);
    }

    // The reference may use accession-style names (NC_*) while the slice
    // regions carry BAM-style names (chr1, chr2). Build a mapper from the
    // reference's .fai so we can translate BAM → FASTA when querying.
    let fai_path = format!("{}.fai", ref_path);
    let fasta_mapper = ContigMapper::from_fai(&fai_path, contigs).ok();

    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(ref_path)
        .map_err(|e| anyhow!("open ref {}: {}", ref_path, e))?;

    let out_file = File::create(out_fa)
        .map_err(|e| anyhow!("create {}: {}", out_fa, e))?;
    let mut out = GzEncoder::new(out_file, Compression::default());
    let mut fai = File::create(format!("{}.fai", out_fa))
        .map_err(|e| anyhow!("create {}.fai: {}", out_fa, e))?;

    let line_bases: usize = 60;
    let line_width: usize = line_bases + 1; // 60 bases + LF
    // Tracks offset into the *uncompressed* stream -- that's what .fai
    // entries reference, regardless of the on-disk gzip framing.
    let mut offset: u64 = 0;
    let mut written = 0usize;

    for (chrom, start, end) in regions {
        let fasta_chrom = fasta_mapper
            .as_ref()
            .map(|m| m.to_bam_name(chrom))
            .unwrap_or_else(|| chrom.clone());

        let region_str = format!("{}:{}-{}", fasta_chrom, start, end);
        let region: Region = match region_str.parse() {
            Ok(r) => r,
            Err(e) => {
                warn!("[slice] region parse {}: {}", region_str, e);
                continue;
            }
        };
        let rec = match reader.query(&region) {
            Ok(r) => r,
            Err(e) => {
                warn!("[slice] ref fetch {}: {}", region_str, e);
                continue;
            }
        };
        let bases = rec.sequence().as_ref();
        let total_bases = bases.len();

        // Header — use the BAM-style chrom so the viewer can look up by it.
        let name = format!("{}:{}-{}", chrom, start, end);
        let header_line = format!(">{}\n", name);
        out.write_all(header_line.as_bytes())?;
        offset += header_line.len() as u64;
        let bases_offset = offset;

        // Bases wrapped at line_bases per line.
        let mut i = 0;
        while i < total_bases {
            let end_i = (i + line_bases).min(total_bases);
            out.write_all(&bases[i..end_i])?;
            out.write_all(b"\n")?;
            i = end_i;
        }
        let full_lines = total_bases / line_bases;
        let last_line_bases = total_bases % line_bases;
        let bytes_written = full_lines * line_width
            + if last_line_bases > 0 { last_line_bases + 1 } else { 0 };
        offset += bytes_written as u64;

        // .fai entry — samtools format.
        let fai_line = format!(
            "{}\t{}\t{}\t{}\t{}\n",
            name, total_bases, bases_offset, line_bases, line_width
        );
        fai.write_all(fai_line.as_bytes())?;

        written += 1;
    }

    out.finish().map_err(|e| anyhow!("gzip finish {}: {}", out_fa, e))?;

    info!(
        "[slice] wrote {} ref regions → {} (+ {}.fai)",
        written, out_fa, out_fa
    );
    Ok(written)
}

/// Read the just-written slice BAM in coordinate order and produce a BAI
/// alongside it at `<bam_path>.bai`. Mirrors the indexer-feeding pattern
/// used by `input.rs::update_inline_index` for the inline-built BAI on the
/// input BAM, so behaviour is consistent between the two paths.
fn build_bai(bam_path: &str, header_n_refs: usize) -> Result<()> {
    use std::io::BufReader;

    let file = File::open(bam_path)
        .map_err(|e| anyhow!("reopen {} for indexing: {}", bam_path, e))?;
    let mut reader = bam::io::Reader::new(BufReader::new(file));
    let _header = reader.read_header()
        .map_err(|e| anyhow!("read header for indexing: {}", e))?;

    let mut indexer: CsiIndexer<LinearIndex> = CsiIndexer::default();
    let mut record = bam::Record::default();
    let mut last_seen: Option<(i32, i32)> = None;
    let mut n_refs_seen = 0usize;

    loop {
        // Capture chunk_start BEFORE reading; the read advances the bgzf
        // reader's virtual position to chunk_end.
        let chunk_start = reader.get_ref().virtual_position();
        let n = reader.read_record(&mut record)
            .map_err(|e| anyhow!("read record for indexing: {}", e))?;
        if n == 0 { break; }
        let chunk_end = reader.get_ref().virtual_position();
        let chunk = Chunk::new(chunk_start, chunk_end);

        let ref_id: i32 = match record.reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            Some(Err(e)) => return Err(anyhow!("ref id during index: {}", e)),
            None => -1,
        };
        let pos: i32 = match record.alignment_start() {
            Some(Ok(p)) => (p.get() as i32) - 1,
            Some(Err(e)) => return Err(anyhow!("pos during index: {}", e)),
            None => -1,
        };
        let flag: u16 = record.flags().bits();

        if ref_id < 0 || pos < 0 {
            indexer.add_record(None, chunk)
                .map_err(|e| anyhow!("indexer add unmapped: {}", e))?;
            continue;
        }

        // Sortedness check — same as input.rs::update_inline_index. If we
        // ever wrote out-of-order records we'd rather know than ship a
        // silently-broken BAI.
        if let Some((pr, pp)) = last_seen
            && (ref_id < pr || (ref_id == pr && pos < pp)) {
            return Err(anyhow!(
                "slice BAM is not coordinate-sorted at {}:{} (after {}:{}) — \
                 cannot build a BAI from it",
                ref_id, pos, pr, pp,
            ));
        }
        last_seen = Some((ref_id, pos));

        let span = cigar_ref_span(&record).max(1);
        let start_1b = ((pos + 1) as usize).max(1);
        let end_1b = start_1b + span - 1;
        let s = Position::try_from(start_1b)
            .map_err(|e| anyhow!("indexer start: {}", e))?;
        let e_ = Position::try_from(end_1b)
            .map_err(|e| anyhow!("indexer end: {}", e))?;
        let is_mapped = (flag & 0x4) == 0;
        let rid = ref_id as usize;
        n_refs_seen = n_refs_seen.max(rid + 1);
        indexer.add_record(Some((rid, s, e_, is_mapped)), chunk)
            .map_err(|e| anyhow!("indexer add: {}", e))?;
    }

    let n_refs = header_n_refs.max(n_refs_seen).max(1);
    let index = indexer.build(n_refs);
    let bai_path = format!("{}.bai", bam_path);
    bam::bai::fs::write(&bai_path, &index)
        .map_err(|e| anyhow!("write bai {}: {}", bai_path, e))?;
    Ok(())
}

/// Sum the reference-consuming op lengths in a noodles bam::Record's CIGAR
/// to get the alignment span (in reference bases). Returns 1 if the CIGAR
/// is empty or unreadable, since the BAI bin computation needs a positive
/// span.
fn cigar_ref_span(record: &bam::Record) -> usize {
    let mut span = 0usize;
    for op_result in record.cigar().iter() {
        let Ok(op) = op_result else { continue };
        match op.kind() {
            CigarKind::Match
            | CigarKind::Deletion
            | CigarKind::Skip
            | CigarKind::SequenceMatch
            | CigarKind::SequenceMismatch => span += op.len(),
            _ => {}
        }
    }
    span
}

/// Suppress unused-import warning when we build without something.
#[allow(dead_code)]
fn _silence_writes(_: &mut dyn Write) {}
