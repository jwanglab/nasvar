use crate::input::{AlignmentInput, CigarKind};
use crate::config::ItdRegion;
use crate::output::{ItdOutput, ItdEvent};
use log::{info, error};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

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

    // create a hashmap to store indels
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

pub fn call_itds(
    bam: &mut AlignmentInput,
    config: &ItdConfig,
) -> Result<ItdOutput, Box<dyn std::error::Error>>
{
    let mut genes: HashMap<String, Vec<ItdEvent>> = HashMap::new();

    for (name, region) in config {
        info!("Processing ITDs for {}...", name);
        let (chr, s, e) = (&region.chrom, region.start, region.end);
        // Translate config chromosome name to BAM naming convention
        let bam_chr = bam.contig_mapper.to_bam_name(chr);
        let mut covg = vec![0; e - s + 1];
        let indels = match get_indels(bam, &mut covg, bam_chr, s, e) {
            Ok(x) => x,
            Err(e) => {
                error!("Error processing indels: {}", e);
                return Err(e);
            }
        };

        // collect and merge indel calls
        let mut events: Vec<ItdEvent> = Vec::new();
        let min_length = region.min_length;
        let min_frequency = region.min_frequency;
        let mut sorted_ins_positions: Vec<&usize> = indels.keys().collect();
        sorted_ins_positions.sort_unstable();
        for t in &sorted_ins_positions {
            let cts: Vec<(&isize, &usize)> = indels[t]
                .iter()
                .filter(|(l, count)| **l >= min_length && **count > 1)
                .collect();
            if cts.is_empty() {
                continue;
            }
            let mut sorted_ins: Vec<&isize> = indels[t].keys().collect();
            sorted_ins.sort_unstable();
            for l in sorted_ins {
                if *l >= min_length && indels[t][l] > 1 {
                    events.push(ItdEvent {
                        position: **t as i64,
                        length: *l as i64,
                        merged: indels[t][l] as i64,
                        coverage: covg[**t - s] as i64,
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
        // filter by minimum allelic frequency
        if min_frequency > 0.0 {
            events.retain(|e| {
                let merged = e.merged as f64;
                let coverage = e.coverage as f64;
                let denom = coverage - merged;
                denom > 0.0 && merged / denom >= min_frequency
            });
        }
        genes.insert(name.clone(), events);
    }

    Ok(ItdOutput { genes })
}
