use crate::input::{AlignmentInput, CigarKind};
use crate::config::ItdRegion;
use crate::output::{ItdOutput, ItdEvent};
use crate::utils::annotation::{get_gene_annotation, Strand};
use crate::bam::ContigMapper;
use log::{info, warn, error};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use noodles::fasta;
use noodles::core::Region;

/// ITD configuration: gene name -> region
pub type ItdConfig = HashMap<String, ItdRegion>;

/// Load ITD configuration from a JSON file
pub fn load_itd_config(path: &str) -> std::io::Result<ItdConfig> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let config: ItdConfig = serde_json::from_reader(reader)?;
    Ok(config)
}

pub fn get_indels(
    bam: &mut AlignmentInput,
    covg: &mut [usize],
    chr: String,
    s: usize,
    e: usize,
) -> Result<HashMap<usize, HashMap<isize, usize>>, Box<dyn std::error::Error>>
{
    let region = format!("{}:{}-{}", chr, s, e);
    let query = bam.query(&region)?;

    let mut indels: HashMap<usize, HashMap<isize, usize>> = HashMap::new();
    for result in query {
        let record = result?;
        let mut t: usize = record.alignment_start().unwrap_or(0);
        for &(op, len) in record.cigar_ops() {
            match op {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    if t <= e && t + len >= s {
                        for i in
                            (if s > t { s } else { t })..(if t + len < e { t + len } else { e })
                        {
                            covg[i - s] += 1;
                        }
                    }
                    t += len;
                }
                CigarKind::Insertion => {
                    if t < e && t >= s {
                        let entry = indels.entry(t).or_default();
                        let count = entry.entry(len as isize).or_insert(0);
                        *count += 1;
                    }
                }
                CigarKind::Deletion | CigarKind::Skip | CigarKind::Pad => {
                    if t <= e && t + len >= s {
                        for i in
                            (if s > t { s } else { t })..(if t + len < e { t + len } else { e })
                        {
                            covg[i - s] += 1;
                        }
                        let entry = indels.entry(t).or_default();
                        let count = entry.entry(-(len as isize)).or_insert(0);
                        *count += 1;
                    }
                    t += len;
                }
                CigarKind::SoftClip | CigarKind::HardClip => {}
            }
        }
    }
    Ok(indels)
}

/// Map a genomic position to a (1-based aa position, 0-based cds offset).
/// `cds` must be sorted ascending. Returns None if pos falls outside all CDS intervals.
fn genomic_pos_to_aa(pos: i64, cds: &[(u32, u32)], strand: Strand) -> Option<(i64, i64)> {
    let mut cds_offset: i64 = 0;
    match strand {
        Strand::Forward => {
            for &(s, e) in cds {
                let (s, e) = (s as i64, e as i64);
                if pos > e {
                    cds_offset += e - s + 1;
                } else if pos >= s {
                    cds_offset += pos - s;
                    return Some((cds_offset / 3 + 1, cds_offset));
                } else {
                    break;
                }
            }
        }
        Strand::Reverse => {
            for &(s, e) in cds.iter().rev() {
                let (s, e) = (s as i64, e as i64);
                if pos < s {
                    cds_offset += e - s + 1;
                } else if pos <= e {
                    cds_offset += e - pos;
                    return Some((cds_offset / 3 + 1, cds_offset));
                } else {
                    break;
                }
            }
        }
        Strand::None => {}
    }
    None
}

/// Fetch the 3 bases of the codon containing `cds_offset` (0-based) from FASTA.
/// Returns the codon as a coding-strand string (5'→3'), ready for translation.
fn fetch_codon(
    cds_offset: i64,
    cds: &[(u32, u32)],
    strand: Strand,
    fasta_path: &str,
    fasta_chrom: &str,
) -> Option<char> {
    let codon_start = (cds_offset / 3 * 3) as usize;

    // Traverse CDS in 5'→3' order to find the 3 genomic positions for this codon
    let ordered: Vec<(u32, u32)> = match strand {
        Strand::Reverse => cds.iter().rev().cloned().collect(),
        _ => cds.to_vec(),
    };

    let mut genomic_positions: Vec<u32> = Vec::with_capacity(3);
    let mut cumulative = 0usize;

    'outer: for &(s, e) in &ordered {
        let interval_len = (e - s + 1) as usize;
        for pos_in_codon in 0..3usize {
            let cds_pos = codon_start + pos_in_codon;
            if genomic_positions.len() == pos_in_codon
                && cds_pos >= cumulative
                && cds_pos < cumulative + interval_len
            {
                let within = (cds_pos - cumulative) as u32;
                let genomic = match strand {
                    Strand::Reverse => e - within,
                    _ => s + within,
                };
                genomic_positions.push(genomic);
                if genomic_positions.len() == 3 {
                    break 'outer;
                }
            }
        }
        cumulative += interval_len;
        if cumulative >= codon_start + 3 {
            break;
        }
    }

    if genomic_positions.len() < 3 {
        return None;
    }

    // Fetch each base from FASTA
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .ok()?;
    let mut codon = String::with_capacity(3);
    for &gpos in &genomic_positions {
        let region_str = format!("{}:{}-{}", fasta_chrom, gpos, gpos);
        let region: Region = region_str.parse().ok()?;
        let rec = reader.query(&region).ok()?;
        let base = std::str::from_utf8(rec.sequence().as_ref())
            .ok()?
            .to_ascii_uppercase()
            .chars()
            .next()
            .unwrap_or('N');
        let coding_base = match strand {
            Strand::Reverse => match base {
                'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', x => x,
            },
            _ => base,
        };
        codon.push(coding_base);
    }

    if codon.len() < 3 {
        return None;
    }

    Some(match codon.as_str() {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "TAA" | "TAG" | "TGA" => '*',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        _ => 'X',
    })
}

pub fn call_itds(
    bam: &mut AlignmentInput,
    config: &ItdConfig,
    gff_path: Option<&str>,
    fasta_path: Option<&str>,
) -> Result<ItdOutput, Box<dyn std::error::Error>>
{
    let mut genes: HashMap<String, Vec<ItdEvent>> = HashMap::new();

    for (name, region) in config {
        info!("Processing ITDs for {}...", name);
        let (chr, s, e) = (&region.chrom, region.start, region.end);
        let bam_chr = bam.contig_mapper.to_bam_name(chr);
        let mut covg = vec![0; e - s + 1];
        let indels = match get_indels(bam, &mut covg, bam_chr, s, e) {
            Ok(x) => x,
            Err(e) => {
                error!("Error processing indels: {}", e);
                return Err(e);
            }
        };

        // Collect indel calls
        let mut events: Vec<ItdEvent> = Vec::new();
        let min_length = region.min_length;
        let min_frequency = region.min_frequency;
        let min_supporting_reads = region.min_supporting_reads;
        let mut sorted_ins_positions: Vec<&usize> = indels.keys().collect();
        sorted_ins_positions.sort_unstable();
        for t in &sorted_ins_positions {
            let cts: Vec<(&isize, &usize)> = indels[t]
                .iter()
                .filter(|(l, count)| **l >= min_length && **count >= min_supporting_reads)
                .collect();
            if cts.is_empty() {
                continue;
            }
            let mut sorted_ins: Vec<&isize> = indels[t].keys().collect();
            sorted_ins.sort_unstable();
            for l in sorted_ins {
                if *l >= min_length && indels[t][l] >= min_supporting_reads {
                    events.push(ItdEvent {
                        chrom: chr.clone(),
                        position: **t as i64,
                        length: *l as i64,
                        merged: indels[t][l] as i64,
                        coverage: covg[**t - s] as i64,
                        aa_position: None,
                        aa_ref: None,
                    });
                }
            }
        }
        // Merge nearby similar-length insertions
        for t in &sorted_ins_positions {
            let mut sorted_ins: Vec<&isize> = indels[t].keys().collect();
            sorted_ins.sort_unstable();
            for l in sorted_ins {
                for event in &mut events {
                    if (event.position - (**t as i64)).abs()
                        < (event.length as f64 * 1.1) as i64
                        && (*l as f64 * 0.9) < (event.length as f64)
                        && (event.length as f64) < (*l as f64 * 1.1)
                        && (event.position != **t as i64 || event.length != *l as i64)
                    {
                        event.merged += indels[t][l] as i64;
                    }
                }
            }
        }
        // Filter by minimum allelic frequency
        if min_frequency > 0.0 {
            events.retain(|e| {
                let merged = e.merged as f64;
                let coverage = e.coverage as f64;
                let denom = coverage - merged;
                denom > 0.0 && merged / denom >= min_frequency
            });
        }

        // Annotate amino acid positions if a transcript is configured
        if let Some(transcript) = region.transcript.as_deref() && let Some(gff) = gff_path {
            match get_gene_annotation(name, gff, Some(transcript)) {
                Ok(ann) => {
                    if let Some(strand) = ann.strand {
                        let mut cds = ann.cds_list.clone();
                        cds.sort_by_key(|&(s, _)| s);
                        let chrom_str = ann.chrom_bytes
                            .as_ref()
                            .and_then(|b| std::str::from_utf8(b).ok())
                            .map(|s| s.to_string());

                        // Derive FASTA chromosome name from FAI if available
                        let fasta_chrom = fasta_path.and_then(|fp| {
                            chrom_str.as_ref().map(|chrom| {
                                let fai_path = format!("{}.fai", fp);
                                let mapper = ContigMapper::from_fai(&fai_path)
                                    .unwrap_or_else(|_| bam.contig_mapper.clone());
                                mapper.to_bam_name(chrom)
                            })
                        });

                        for event in &mut events {
                            if let Some((aa_pos, cds_offset)) =
                                genomic_pos_to_aa(event.position, &cds, strand)
                            {
                                event.aa_position = Some(aa_pos);
                                if let (Some(fp), Some(fchrom)) = (fasta_path, &fasta_chrom) {
                                    event.aa_ref = fetch_codon(cds_offset, &cds, strand, fp, fchrom);
                                }
                            }
                        }
                    }
                }
                Err(e) => warn!("Could not get annotation for {} ({}): {}", name, transcript, e),
            }
        }

        genes.insert(name.clone(), events);
    }

    Ok(ItdOutput { genes })
}
