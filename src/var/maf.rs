use crate::config::Contig;
use crate::input::{AlignmentInput, CigarKind};
use crate::utils::bed::BedRegion;
use crate::utils::contig::ContigMapper;
use log::debug;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Clone, Copy, Debug)]
pub struct Site {
    pub pos: usize,
    pub al0: u8,
    pub al1: u8,
}

pub fn read_sites(path: &str, contigs: &[Contig]) -> Result<Vec<Vec<Site>>, Box<dyn std::error::Error>> {
    let mapper = ContigMapper::from_contigs(contigs);
    let mut sites = vec![Vec::new(); mapper.num_chromosomes()];
    let file = File::open(path).map_err(|e| {
        std::io::Error::other(format!("Error opening sites file {}: {}", path, e))
    })?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let l = line?;
        let p: Vec<&str> = l.split('\t').collect();
        if p.len() < 4 {
            continue;
        }

        let idx = match mapper.parse_chr_index(p[0]) {
            Some(i) => i,
            None => continue,
        };

        let pos: usize = p[1].parse()?;

        // Parse alleles A,C,G,T -> 0,1,2,3
        let a0 = match p[2] {
            "A" => 0,
            "C" => 1,
            "G" => 2,
            "T" => 3,
            _ => 4,
        };
        let a1 = match p[3] {
            "A" => 0,
            "C" => 1,
            "G" => 2,
            "T" => 3,
            _ => 4,
        };

        if a0 < 4 && a1 < 4 {
            sites[idx].push(Site {
                pos,
                al0: a0,
                al1: a1,
            });
        }
    }
    Ok(sites)
}

/// Keep only sites that fall within enriched BED regions.
pub fn filter_enriched_sites(sites: &mut [Vec<Site>], enriched: &[BedRegion], contigs: &[Contig]) {
    let mapper = ContigMapper::from_contigs(contigs);
    let n = mapper.num_chromosomes();

    let mut enriched_by_chr: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for r in enriched {
        if let Some(idx) = mapper.get_chr_index(&r.segment) && idx < n {
            enriched_by_chr[idx].push((r.start as usize, r.end as usize));
        }
    }
    for intervals in &mut enriched_by_chr {
        intervals.sort_unstable();
    }

    for (chr_idx, chr_sites) in sites.iter_mut().enumerate() {
        if chr_idx >= n { break; }
        let intervals = &enriched_by_chr[chr_idx];
        if intervals.is_empty() {
            chr_sites.clear();
            continue;
        }
        chr_sites.retain(|site| {
            let i = intervals.partition_point(|&(start, _)| start <= site.pos);
            if i == 0 {
                return false;
            }
            let (_, end) = intervals[i - 1];
            site.pos < end
        });
    }
}

/// Remove sites that fall within repeat/masked BED regions.
pub fn filter_repeat_sites(sites: &mut [Vec<Site>], repeats: &[BedRegion], contigs: &[Contig]) {
    let mapper = ContigMapper::from_contigs(contigs);
    let n = mapper.num_chromosomes();

    // Group repeat intervals by chromosome index, sorted by start
    let mut repeats_by_chr: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for r in repeats {
        if let Some(idx) = mapper.get_chr_index(&r.segment) && idx < n {
            repeats_by_chr[idx].push((r.start as usize, r.end as usize));
        }
    }
    for intervals in &mut repeats_by_chr {
        intervals.sort_unstable();
    }

    for (chr_idx, chr_sites) in sites.iter_mut().enumerate() {
        if chr_idx >= n { break; }
        let intervals = &repeats_by_chr[chr_idx];
        if intervals.is_empty() { continue; }
        chr_sites.retain(|site| {
            // Binary search: find the last interval whose start <= site.pos
            let i = intervals.partition_point(|&(start, _)| start <= site.pos);
            if i == 0 {
                return true; // site is before all intervals
            }
            let (_, end) = intervals[i - 1];
            site.pos >= end // keep if site is at or past the end of the interval
        });
    }
}

pub fn calc_maf(
    bam: &mut AlignmentInput,
    enriched: &[BedRegion],
    sites: &[Vec<Site>],
    out_prefix: &str,
) -> Result<(), Box<dyn std::error::Error>>
{
    use std::io::Write;

    // Output file
    let out_maf = format!("{}.maf", out_prefix);
    let mut file = File::create(out_maf)?;
    writeln!(file, "chrom\tpos\tref\talt")?;

    let mut current_chrom = String::new();

    for target in enriched {
        if target.segment != current_chrom {
            debug!("Processing chromosome {}...", target.segment);
            current_chrom = target.segment.clone();
        }

        let chr_idx = match bam.contig_mapper.get_chr_index(&target.segment) {
            Some(i) => i,
            None => continue,
        };

        // Filter sites in this region
        let target_sites: Vec<&Site> = sites[chr_idx]
            .iter()
            .filter(|s| s.pos >= target.start as usize && s.pos <= target.end as usize)
            .collect();

        if target_sites.is_empty() {
            continue;
        }

        // Setup pileup map: Pos -> [A, C, G, T, N] counts
        let mut counts: HashMap<usize, [u32; 5]> = HashMap::new();
        for s in &target_sites {
            counts.insert(s.pos, [0; 5]);
        }

        // Query BAM (translate chromosome name to BAM convention)
        let bam_chrom = bam.contig_mapper.to_bam_name(&target.segment);
        let region_str = format!("{}:{}-{}", bam_chrom, target.start, target.end);

        let query = bam.query(&region_str)?;

        for result in query {
            let record = result?;
            let start = match record.alignment_start() {
                Some(p) => p,
                None => continue,
            };
            let seq = record.sequence();

            let mut t = start;
            let mut q = 0; // query pos

            for &(op, len) in record.cigar_ops() {
                match op {
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                        for i in 0..len {
                            let pos = t + i;
                            if let Some(cnt) = counts.get_mut(&pos)
                                && q + i < seq.len() {
                                    let base = seq.get(q + i).copied().unwrap_or(b'N');
                                    let al = match base {
                                        b'A' => 0,
                                        b'C' => 1,
                                        b'G' => 2,
                                        b'T' => 3,
                                        _ => 4,
                                    };
                                    cnt[al] += 1;
                                }
                        }
                        t += len;
                        q += len;
                    }
                    CigarKind::Insertion | CigarKind::SoftClip => {
                        q += len;
                    }
                    CigarKind::Deletion | CigarKind::Skip => {
                        t += len;
                    }
                    _ => {}
                }
            }
        }

        // Write results (all sites, including sub-threshold)
        for site in &target_sites {
            if let Some(cnt) = counts.get(&site.pos) {
                let ref_ct = cnt[site.al0 as usize];
                let alt_ct = cnt[site.al1 as usize];
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}",
                    target.segment, site.pos, ref_ct, alt_ct
                )?;
            }
        }
    }

    Ok(())
}
