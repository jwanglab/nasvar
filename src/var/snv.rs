use std::collections::HashMap;

use noodles::core::Region;
use noodles::fasta;

use crate::bam::ContigMapper;
use crate::input::{AlignmentInput, CigarKind};
use crate::config::PipelineConfig;
use crate::output::{SnvOutput, SnvGeneResult};
use crate::utils::annotation::{Strand, get_gene_annotation};
use log::{info, warn};

fn get_codon_table() -> HashMap<String, char> {
    let mut map = HashMap::new();
    let codons = vec![
        ("TCA", "S"),
        ("TCC", "S"),
        ("TCG", "S"),
        ("TCT", "S"),
        ("TTC", "F"),
        ("TTT", "F"),
        ("TTA", "L"),
        ("TTG", "L"),
        ("TAC", "Y"),
        ("TAT", "Y"),
        ("TAA", "*"),
        ("TAG", "*"),
        ("TGC", "C"),
        ("TGT", "C"),
        ("TGA", "*"),
        ("TGG", "W"),
        ("CTA", "L"),
        ("CTC", "L"),
        ("CTG", "L"),
        ("CTT", "L"),
        ("CCA", "P"),
        ("CCC", "P"),
        ("CCG", "P"),
        ("CCT", "P"),
        ("CAC", "H"),
        ("CAT", "H"),
        ("CAA", "Q"),
        ("CAG", "Q"),
        ("CGA", "R"),
        ("CGC", "R"),
        ("CGG", "R"),
        ("CGT", "R"),
        ("ATA", "I"),
        ("ATC", "I"),
        ("ATT", "I"),
        ("ATG", "M"),
        ("ACA", "T"),
        ("ACC", "T"),
        ("ACG", "T"),
        ("ACT", "T"),
        ("AAC", "N"),
        ("AAT", "N"),
        ("AAA", "K"),
        ("AAG", "K"),
        ("AGC", "S"),
        ("AGT", "S"),
        ("AGA", "R"),
        ("AGG", "R"),
        ("GTA", "V"),
        ("GTC", "V"),
        ("GTG", "V"),
        ("GTT", "V"),
        ("GCA", "A"),
        ("GCC", "A"),
        ("GCG", "A"),
        ("GCT", "A"),
        ("GAC", "D"),
        ("GAT", "D"),
        ("GAA", "E"),
        ("GAG", "E"),
        ("GGA", "G"),
        ("GGC", "G"),
        ("GGG", "G"),
        ("GGT", "G"),
    ];
    for (k, v) in codons {
        map.insert(k.to_string(), v.chars().next().unwrap());
    }
    map
}

fn translate(nt: &str, table: &HashMap<String, char>) -> String {
    let mut aa = String::new();
    for i in (0..nt.len()).step_by(3) {
        if i + 3 > nt.len() {
            break;
        }
        let codon = &nt[i..i + 3];
        if codon.contains('N') {
            aa.push(' ');
        } else if codon.contains('-') {
            aa.push('-');
        } else {
            aa.push(*table.get(codon).unwrap_or(&'X'));
        }
    }
    aa
}

fn rev_comp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

pub fn call_snvs(
    bam: &mut AlignmentInput,
    fasta_path: &str,
    gff_path: &str,
    gene_list: &[String],
    pipeline_config: &PipelineConfig,
) -> Result<SnvOutput, Box<dyn std::error::Error>>
{
    // Get thresholds from config
    let min_coverage = pipeline_config.thresholds.snv.min_coverage;
    let min_allele_freq = pipeline_config.thresholds.snv.min_allele_freq;
    let homozygosity_vaf = pipeline_config.thresholds.snv.homozygosity_vaf;
    let phasing_ratio = pipeline_config.thresholds.snv.phasing_ratio;

    // Get SNV transcript configs from pipeline config
    let snv_transcripts = &pipeline_config.genes.snv.transcripts;

    // Open FASTA
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .map_err(|e| {
            std::io::Error::other(format!("Error opening FASTA file {}: {}", fasta_path, e))
        })?;

    // Detect FASTA naming convention (may differ from BAM)
    let fai_path = format!("{}.fai", fasta_path);
    let fasta_mapper = ContigMapper::from_fai(&fai_path).unwrap_or_else(|_| {
        warn!("Could not read FASTA index {}, assuming FASTA matches BAM naming", fai_path);
        bam.contig_mapper.clone()
    });

    let codon_table = get_codon_table();

    // Check if we should use gene list from CLI OR from config
    let search_genes = if !gene_list.is_empty() {
        gene_list.to_vec()
    } else {
        snv_transcripts.keys().cloned().collect()
    };

    let mut genes: HashMap<String, SnvGeneResult> = HashMap::new();

    for gene in search_genes {
        info!("Processing SNVs for {}...", gene);

        let mut target_transcript = None;
        let mut gene_targets: Option<&HashMap<String, Vec<usize>>> = None;

        // Lookup in config
        if let Some(gene_cfg) = snv_transcripts.get(&gene) {
            target_transcript = gene_cfg.transcript.as_deref();
            if !gene_cfg.variants.is_empty() {
                gene_targets = Some(&gene_cfg.variants);
            }
        } else {
            // If gene requested but not in config?
            warn!("Gene {} not found in configuration.", gene);
        }

        let ann = get_gene_annotation(&gene, gff_path, target_transcript)?;
        let mrna_id = ann.mrna_id;
        let chrom = ann.chrom_bytes;
        let gene_strand = ann.strand;
        let mut cds = ann.cds_list;

        if chrom.is_none() || cds.is_empty() {
            warn!("Could not find annotation for {}", gene);
            continue;
        }

        let chrom_str = String::from_utf8(chrom.unwrap().to_vec())?;
        // Translate GFF chromosome name (always NC_*) to BAM and FASTA conventions
        let bam_chrom = bam.contig_mapper.to_bam_name(&chrom_str);
        let fasta_chrom = fasta_mapper.to_bam_name(&chrom_str);
        let strand = gene_strand.unwrap(); // assumed valid

        cds.sort_by(|a, b| a.0.cmp(&b.0));

        // Calculate total CDS length
        let mut cds_len: usize = 0;
        for (s, e) in &cds {
            cds_len += (e - s + 1) as usize;
        }

        info!(
            "Transcript: {}, Exons: {}, CDS Length: {} bp",
            mrna_id,
            cds.len(),
            cds_len
        );

        // Data structs
        let mut covg: Vec<[u32; 7]> = vec![[0; 7]; cds_len]; // Depth, A, C, G, T, N, Del
        let mut reads: HashMap<String, Vec<char>> = HashMap::new(); // Read tracks CDS seq

        let mut offset = 0;
        let mut ref_cds_seq = String::new();

        for (s, e) in &cds {
            let start = *s;
            let end = *e;

            // Get Ref Seq for this exon (FASTA may use different naming than BAM)
            let region_fasta: Region = format!("{}:{}-{}", fasta_chrom, start, end).parse()
                .map_err(|e| std::io::Error::other(format!("Invalid FASTA region: {}", e)))?;
            let rec = fasta_reader.query(&region_fasta)?;
            ref_cds_seq.push_str(
                std::str::from_utf8(rec.sequence().as_ref())?
                    .to_ascii_uppercase()
                    .as_str(),
            );

            // Check BAM (using translated chromosome name)
            let region_str = format!("{}:{}-{}", bam_chrom, s, e);

            let query = bam.query(&region_str)?;

            for result in query {
                let record = result?;
                let name = record.name().unwrap().to_string();

                let rec_start = match record.alignment_start() {
                    Some(p) => p,
                    _ => continue,
                };
                let seq = record.sequence();

                // If not in reads map, init with ' '
                let read_vec = reads.entry(name).or_insert(vec![' '; cds_len]);

                // Pileup logic
                let mut t_pos = rec_start; // 1-based target pos
                let mut q_pos = 0; // 0-based query pos

                for &(op, len) in record.cigar_ops() {
                    match op {
                        CigarKind::Match
                        | CigarKind::SequenceMatch
                        | CigarKind::SequenceMismatch => {
                            for i in 0..len {
                                let curr_t = t_pos + i;
                                let s_usize = *s as usize;
                                let e_usize = *e as usize;
                                if curr_t >= s_usize && curr_t <= e_usize {
                                    let exon_offset = curr_t - s_usize;
                                    let mut final_cds_pos = offset + exon_offset;

                                    if strand == Strand::Reverse {
                                        final_cds_pos = cds_len - final_cds_pos - 1;
                                    }

                                    // Base
                                    let base_byte = seq.get(q_pos + i).copied().unwrap_or(b'N');
                                    let char_base = match base_byte {
                                        b'A' => 'A',
                                        b'C' => 'C',
                                        b'G' => 'G',
                                        b'T' => 'T',
                                        _ => 'N',
                                    };

                                    let final_base = if strand == Strand::Reverse {
                                        match char_base {
                                            'A' => 'T',
                                            'T' => 'A',
                                            'C' => 'G',
                                            'G' => 'C',
                                            x => x,
                                        }
                                    } else {
                                        char_base
                                    };

                                    if final_cds_pos < read_vec.len() {
                                        read_vec[final_cds_pos] = final_base;
                                    }

                                    if final_cds_pos < covg.len() {
                                        covg[final_cds_pos][0] += 1;
                                        let idx = match final_base {
                                            'A' => 1,
                                            'C' => 2,
                                            'G' => 3,
                                            'T' => 4,
                                            'N' => 5,
                                            '-' => 6,
                                            _ => 5,
                                        };
                                        covg[final_cds_pos][idx] += 1;
                                    }
                                }
                            }
                            t_pos += len;
                            q_pos += len;
                        }
                        CigarKind::Deletion => {
                            for i in 0..len {
                                let curr_t = t_pos + i;
                                let s_usize = *s as usize;
                                let e_usize = *e as usize;
                                if curr_t >= s_usize && curr_t <= e_usize {
                                    let exon_offset = curr_t - s_usize;
                                    let mut final_cds_pos = offset + exon_offset;
                                    if strand == Strand::Reverse {
                                        final_cds_pos = cds_len - final_cds_pos - 1;
                                    }

                                    let final_base = '-';

                                    if final_cds_pos < read_vec.len() {
                                        read_vec[final_cds_pos] = final_base;
                                    }

                                    if final_cds_pos < covg.len() {
                                        covg[final_cds_pos][0] += 1;
                                        covg[final_cds_pos][6] += 1;
                                    }
                                }
                            }
                            t_pos += len;
                        }
                        CigarKind::Insertion => {
                            q_pos += len;
                        }
                        CigarKind::SoftClip => {
                            q_pos += len;
                        }
                        CigarKind::Skip => {
                            t_pos += len;
                        }
                        _ => {}
                    }
                }
            }

            offset += (e - s + 1) as usize;
        }

        // Process Ref Seq if Reverse
        let final_ref_seq = if strand == Strand::Reverse {
            rev_comp(&ref_cds_seq)
        } else {
            ref_cds_seq.clone()
        };

        // Analyze Variations
        let mut muts: Vec<(usize, char, char, f32)> = Vec::new(); // pos(1-based), ref, alt, af

        for (i, cov) in covg.iter().enumerate().take(cds_len) {
            let depth = cov[0] as f32;
            if (depth as f64) < min_coverage {
                continue;
            }

            // Find max alt
            let ref_char = final_ref_seq.chars().nth(i).unwrap_or('N');
            let mut best_alt = 'N';
            let mut best_cnt = 0;

            for (j, c) in [(1, 'A'), (2, 'C'), (3, 'G'), (4, 'T'), (6, '-')] {
                if c != ref_char && cov[j] > best_cnt {
                    best_cnt = cov[j];
                    best_alt = c;
                }
            }

            let af = best_cnt as f32 / depth;
            if (af as f64) > min_allele_freq {
                muts.push((i + 1, ref_char, best_alt, af));
            }
        }

        // Filter mutations to only target pharmacogenomic positions (matching reference behavior)
        let muts: Vec<(usize, char, char, f32)> = if let Some(targets) = gene_targets {
            let target_positions: std::collections::HashSet<usize> = targets
                .values()
                .flat_map(|loci| loci.iter().copied())
                .collect();
            muts.into_iter()
                .filter(|(pos, _, _, _)| target_positions.contains(pos))
                .collect()
        } else {
            muts
        };

        // Calculate average depth
        let total_depth: u32 = covg.iter().map(|c| c[0]).sum();
        let avg_depth = total_depth as f64 / cds_len as f64;

        // Perform Phasing
        let mut phasing_map: HashMap<String, Option<bool>> = HashMap::new();

        for i in 0..muts.len() {
            let (p1, r1, a1, _) = muts[i];
            let idx1 = p1 - 1;

            for &(p2, r2, a2, _) in muts.iter().skip(i + 1) {
                let idx2 = p2 - 1;

                let key = format!("{}-{}", p1, p2);

                let mut cis_count = 0;
                let mut trans_count = 0;

                for seq in reads.values() {
                    if idx1 >= seq.len() || idx2 >= seq.len() {
                        continue;
                    }
                    let b1 = seq[idx1];
                    let b2 = seq[idx2];

                    if b1 == ' ' || b2 == ' ' {
                        continue;
                    }

                    if (b1 == r1 && b2 == r2) || (b1 == a1 && b2 == a2) {
                        cis_count += 1;
                    } else if (b1 == r1 && b2 == a2) || (b1 == a1 && b2 == r2) {
                        trans_count += 1;
                    }
                }

                let phase_status = if cis_count > 1 && cis_count > trans_count * phasing_ratio {
                    Some(true) // Cis
                } else if trans_count > 1 && trans_count > cis_count * phasing_ratio {
                    Some(false) // Trans
                } else {
                    None
                };

                phasing_map.insert(key, phase_status);
            }
        }

        // Genotype
        let mut matches = Vec::new();
        let mut used_muts = std::collections::HashSet::new();

        if let Some(target_muts) = gene_targets {
            for (name, loci) in target_muts {
                let matching_indices: Vec<usize> = muts
                    .iter()
                    .enumerate()
                    .filter(|(i, (p, _, a, _))| {
                        loci.contains(p) && !used_muts.contains(i) && *a != '-'
                    })
                    .map(|(i, _)| i)
                    .collect();

                if matching_indices.len() == loci.len() {
                    let mut phased_ok = true;
                    if loci.len() > 1 {
                        for i in 0..matching_indices.len() {
                            for j in (i + 1)..matching_indices.len() {
                                let m1 = muts[matching_indices[i]];
                                let m2 = muts[matching_indices[j]];

                                let key = format!("{}-{}", m1.0, m2.0);
                                if let Some(status) = phasing_map.get(&key) {
                                    if *status != Some(true) {
                                        phased_ok = false;
                                        break;
                                    }
                                } else {
                                    phased_ok = false;
                                    break;
                                }
                            }
                            if !phased_ok {
                                break;
                            }
                        }
                    }

                    if phased_ok {
                        matches.push(name.as_str());

                        if loci.len() == 1 {
                            let vaf = muts[matching_indices[0]].3;
                            if (vaf as f64) > homozygosity_vaf {
                                matches.push(name.as_str());
                            }
                        } else if loci.len() > 1 {
                            let avg_vaf: f32 =
                                matching_indices.iter().map(|&i| muts[i].3).sum::<f32>()
                                    / matching_indices.len() as f32;
                            if (avg_vaf as f64) > homozygosity_vaf {
                                matches.push(name.as_str());
                            }
                        }

                        for idx in matching_indices {
                            used_muts.insert(idx);
                        }
                    }
                }
            }
        }

        let remaining_muts: Vec<String> = muts
            .iter()
            .enumerate()
            .filter(|(i, _)| !used_muts.contains(i))
            .map(|(_, (p, r, a, af))| format!("{}{}>{} ({:.2})", r, p, a, af))
            .collect();

        if matches.is_empty() {
            matches.push("*1");
            matches.push("*1");
        } else if matches.len() == 1 {
            matches.push("*1");
        }

        let genotype = matches.join("/");

        // Translation check
        let mut cons_seq = final_ref_seq.chars().collect::<Vec<char>>();
        for (p, _r, a, _af) in &muts {
            if *p <= cons_seq.len() {
                cons_seq[(*p) - 1] = *a;
            }
        }
        let cons_str: String = cons_seq.into_iter().collect();
        let aa = translate(&cons_str, &codon_table);
        let ref_aa = translate(&final_ref_seq, &codon_table);

        // Collect AA muts
        let mut aa_muts = Vec::new();
        for i in 0..aa.len() {
            let r = ref_aa.chars().nth(i).unwrap_or(' ');
            let a = aa.chars().nth(i).unwrap_or(' ');
            if r != a && r != ' ' && a != ' ' {
                aa_muts.push(format!("{}{}{}", r, i + 1, a));
            }
        }

        // Output
        let mut phase_output = HashMap::new();
        for (k, v) in phasing_map {
            let val = match v {
                Some(true) => "alt/alt",
                Some(false) => "ref/alt",
                None => "?",
            };
            phase_output.insert(k, val);
        }

        let mutations_dict: HashMap<String, String> = muts
            .iter()
            .map(|(p, r, a, af)| (format!("{}{}>{}", p, r, a), format!("{:.2}", af)))
            .collect();

        genes.insert(gene.clone(), SnvGeneResult {
            genotype,
            coverage: Some(avg_depth),
            mutations: mutations_dict,
            aa_changes: aa_muts,
            phase: Some(serde_json::to_value(&phase_output)?),
            unassigned_muts: if remaining_muts.is_empty() { None } else { Some(remaining_muts) },
        });
    }

    Ok(SnvOutput { genes })
}

#[cfg(test)]
mod tests {
    use crate::config::SnvTranscriptConfig;
    use std::collections::HashMap;

    #[test]
    fn test_config_parsing() {
        let data = r#"{
    "TPMT": {
      "transcript": "rna-NM_001346817.1",
      "variants": {
        "*5": [146],
        "*3A": [460, 719]
      }
    },
    "NUDT15": {
        "variants": {
            "*5": [52]
        }
    }
}"#;
        let config: HashMap<String, SnvTranscriptConfig> = serde_json::from_str(data).unwrap();

        assert!(config.contains_key("TPMT"));
        let tpmt = config.get("TPMT").unwrap();
        assert_eq!(tpmt.transcript.as_ref().unwrap(), "rna-NM_001346817.1");

        assert_eq!(tpmt.variants.get("*5").unwrap(), &vec![146]);
        assert_eq!(tpmt.variants.get("*3A").unwrap(), &vec![460, 719]);
    }
}
