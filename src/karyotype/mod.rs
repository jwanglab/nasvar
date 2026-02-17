mod plot;

use log::{info, warn, debug};
use crate::bam::ContigMapper;
use crate::config::{KaryotypeThresholds, ReferenceConfig};
use crate::output::KaryotypeOutput;
use crate::var::maf::Site;
use crate::utils::bed::BedRegion;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

// GC bias correction method

/// Method used for GC bias correction during karyotype inference.
#[derive(Debug, Clone, Copy, PartialEq, Default, clap::ValueEnum)]
pub enum GcCorrectionMethod {
    /// No GC bias correction
    None,
    /// Simple linear regression (default)
    #[default]
    Linear,
    /// LOESS local regression (handles non-linear bias)
    Loess,
}

// Data structures

#[derive(Debug, Clone)]
pub struct CoverageBin {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub coverage: f64,         // Normalized coverage
    pub maf_data: Option<f64>, // MAF at this bin location (if any)
    pub gc_content: Option<f64>, // GC content fraction for this bin
}

// --------------------------------------------------------------------------------
// Constants & Lookups
// --------------------------------------------------------------------------------

fn get_segment(chrom: &str, start: u32, end: u32, ref_config: &ReferenceConfig) -> Option<String> {
    // Use ContigMapper to normalize accession IDs to chr names
    let mapper = ContigMapper::new();
    let normalized_chrom = mapper.to_chr_name(chrom);

    if !normalized_chrom.starts_with("chr") {
        return None;
    }

    if let Some(&(cen_start, cen_end)) = ref_config.centromeres.get(&normalized_chrom) {
        if end < cen_start {
            return Some(format!("{}p", normalized_chrom.trim_start_matches("chr")));
        } else if start > cen_end {
            return Some(format!("{}q", normalized_chrom.trim_start_matches("chr")));
        } else if start < cen_end && end > cen_start {
            // Overlaps centromere, skip
            return None;
        }
    } else if let Some(&(p_start, p_end)) = ref_config.par1.get(&normalized_chrom) {
        if start < p_end && end > p_start {
            if start < p_start || end > p_end {
                return None; // Partial overlap
            } else {
                return Some("PAR1".to_string());
            }
        }
    } else if let Some(&(p_start, p_end)) = ref_config.par2.get(&normalized_chrom)
        && end > p_start && start < p_end {
            return None; // Overlaps PAR2, skip
        }

    // Check if it is a canonical chromosome that we want to keep even if not split
    let suffix = normalized_chrom.trim_start_matches("chr");
    // Canonical: 1-22, X, Y
    let is_canonical = match suffix {
        "X" | "Y" => true,
        n => n
            .parse::<usize>()
            .map(|x| (1..=22).contains(&x))
            .unwrap_or(false),
    };

    if is_canonical {
        Some(suffix.to_string())
    } else {
        None
    }
}

// --------------------------------------------------------------------------------
// String Generation
// --------------------------------------------------------------------------------

fn generate_karyotype_string(cns: &HashMap<String, usize>) -> String {
    let mut counts: HashMap<usize, Vec<String>> = HashMap::new();
    for (k, v) in cns {
        counts.entry(*v).or_default().push(k.clone());
    }

    let mut chrom_ct = 0;
    let mut mods = Vec::new();

    // Autosomes: if ANY arm is diploid, count as 2; else if ANY is haploid, count as 1; else min(p,q)
    for i in 1..=22 {
        let s = i.to_string();
        let sp = format!("{}p", i);
        let sq = format!("{}q", i);

        // Check if whole chromosome or any arm is in counts[2]
        let in_2 = counts.get(&2).is_some_and(|v| {
            v.contains(&s) || v.contains(&sp) || v.contains(&sq)
        });
        // Check if whole chromosome or any arm is in counts[1]
        let in_1 = counts.get(&1).is_some_and(|v| {
            v.contains(&s) || v.contains(&sp) || v.contains(&sq)
        });

        let cn = if in_2 {
            2
        } else if in_1 {
            1
        } else if let Some(&c) = cns.get(&s) {
            c
        } else {
            let cp = cns.get(&sp).cloned().unwrap_or(2);
            let cq = cns.get(&sq).cloned().unwrap_or(2);
            std::cmp::min(cp, cq)
        };
        chrom_ct += cn;
    }

    // Baseline ploidy
    let baseline = if chrom_ct <= 30 { 1 } else { 2 };

    // Generate mod strings
    for i in 1..=22 {
        let s = i.to_string();
        let sp = format!("{}p", i);
        let sq = format!("{}q", i);

        if let Some(&c) = cns.get(&s) {
            // Whole chromosome called
            if c < baseline {
                for _ in 0..(baseline - c) {
                    mods.push(format!("-{}", i));
                }
            } else {
                for _ in 0..(c - baseline) {
                    mods.push(format!("+{}", i));
                }
            }
        } else {
            // Arms
            let cp = cns.get(&sp).cloned().unwrap_or(baseline);
            let cq = cns.get(&sq).cloned().unwrap_or(baseline);

            // If both arms have the same CN, treat as whole chromosome
            if cp == cq {
                if cp < baseline {
                    for _ in 0..(baseline - cp) {
                        mods.push(format!("-{}", i));
                    }
                } else {
                    for _ in 0..(cp - baseline) {
                        mods.push(format!("+{}", i));
                    }
                }
            } else {
                // Split
                if cp < baseline {
                    for _ in 0..(baseline - cp) {
                        mods.push(format!("-{}p", i));
                    }
                } else {
                    for _ in 0..(cp - baseline) {
                        mods.push(format!("+{}p", i));
                    }
                }

                if cq < baseline {
                    for _ in 0..(baseline - cq) {
                        mods.push(format!("-{}q", i));
                    }
                } else {
                    for _ in 0..(cq - baseline) {
                        mods.push(format!("+{}q", i));
                    }
                }
            }
        }
    }

    // Sex chromosomes
    let cx = *cns.get("X").unwrap_or(&0); // Default 0 if missing? 
    let cy = *cns.get("Y").unwrap_or(&0);

    chrom_ct += cx + cy;

    let xy_str = if cy >= 1 {
        "XY"
    } else if cx >= 2 {
        "XX"
    } else if cx == 1 {
        "X"
    } else {
        "??"
    };
    let expected_x = xy_str.chars().filter(|&c| c == 'X').count();
    let expected_y = xy_str.chars().filter(|&c| c == 'Y').count();

    if cx < expected_x {
        for _ in 0..(expected_x - cx) {
            mods.push("-X".to_string());
        }
    } else {
        for _ in 0..(cx - expected_x) {
            mods.push("+X".to_string());
        }
    }

    if cy < expected_y {
        for _ in 0..(expected_y - cy) {
            mods.push("-Y".to_string());
        }
    } else {
        for _ in 0..(cy - expected_y) {
            mods.push("+Y".to_string());
        }
    }

    let ploidy_prefix = if baseline == 1 { "(1n) " } else { "" };
    if mods.is_empty() {
        format!("{}{}; {}", chrom_ct, xy_str, ploidy_prefix)
    } else {
        format!(
            "{}{}; {}{}",
            chrom_ct,
            xy_str,
            ploidy_prefix,
            mods.join(", ")
        )
    }
}

fn generate_iscn_string(cns: &HashMap<String, usize>) -> String {
    let mut groups: Vec<String> = Vec::new();
    let mut cts: HashMap<usize, Vec<String>> = HashMap::new();

    // Reverse map
    for (k, v) in cns {
        cts.entry(*v).or_default().push(k.clone());
    }

    // Merge p and q if same
    for (_, list) in cts.iter_mut() {
        // We need to iterate 1..22 and check if both p and q are in list and if so replace with whole
        let metacentrics = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20];

        for i in metacentrics {
            let sp = format!("{}p", i);
            let sq = format!("{}q", i);
            let s = i.to_string();

            if list.contains(&sp) && list.contains(&sq) {
                list.retain(|x| x != &sp && x != &sq);
                list.push(s);
            }
        }
        list.sort_by(|a, b| human_sort(a, b));
    }

    if let Some(list) = cts.get(&1)
        && !list.is_empty() {
            groups.push(format!("({})x1", list.join(", ")));
        }

    // Check if autosomes are ALL diploid
    let dip_list = if let Some(l) = cts.get(&2) {
        l.clone()
    } else {
        Vec::new()
    };

    // Check coverage of 1-22
    let mut aut_diploid = true;
    for i in 1..=22 {
        if !dip_list.contains(&i.to_string()) {
            aut_diploid = false;
            break;
        }
    }

    if aut_diploid || dip_list.contains(&"X".to_string()) || dip_list.contains(&"Y".to_string()) {
        let mut dips = Vec::new();
        if dip_list.contains(&"X".to_string()) {
            dips.push("X");
        }
        if dip_list.contains(&"Y".to_string()) {
            dips.push("Y");
        }
        if aut_diploid {
            dips.push("1-22");
        }

        groups.push(format!("({})x2", dips.join(", ")));
    }

    // > 2
    let max_cn = cts.keys().max().cloned().unwrap_or(0);
    for n in 3..=max_cn {
        if let Some(list) = cts.get(&n)
            && !list.is_empty() {
                groups.push(format!("({})x{}", list.join(", "), n));
            }
    }

    format!("seq {}", groups.join(", "))
}

// sort by numeric part first, then string part
fn human_sort(a: &str, b: &str) -> std::cmp::Ordering {
    let split = |s: &str| -> (u32, String) {
        let num_end = s.find(|c: char| !c.is_numeric()).unwrap_or(s.len());
        let (n, rest) = s.split_at(num_end);
        let num = n.parse::<u32>().unwrap_or(999); // X, Y go to end
        // Map X->100, Y->101, PAR->102 etc if we want specific order
        let num = if s.starts_with('X') {
            100
        } else if s.starts_with('Y') {
            101
        } else {
            num
        };
        (num, rest.to_string())
    };
    split(a).cmp(&split(b))
}

// --------------------------------------------------------------------------------
// Parsing Logic
// --------------------------------------------------------------------------------

pub fn parse_coverage(cov_path: &str) -> Result<Vec<CoverageBin>, Box<dyn std::error::Error>> {
    let file = File::open(cov_path)?;
    let reader = BufReader::new(file);
    let mut bins = Vec::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with("chromosome") {
            continue;
        }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }

        let chrom = parts[0].to_string();
        let start: u32 = parts[1].parse()?;
        let end: u32 = parts[2].parse()?;
        let cov: f64 = parts[3].parse()?;
        let gc = if parts.len() >= 5 {
            parts[4].parse::<f64>().ok().filter(|v| v.is_finite())
        } else {
            None
        };

        bins.push(CoverageBin {
            chrom,
            start,
            end,
            coverage: cov,
            maf_data: None,
            gc_content: gc,
        });
    }
    Ok(bins)
}

// --------------------------------------------------------------------------------
// MAF Helpers
// --------------------------------------------------------------------------------

fn get_segment_from_pos(chrom: &str, pos: u32, ref_config: &ReferenceConfig) -> Option<String> {
    get_segment(chrom, pos, pos, ref_config)
}

pub fn parse_maf(maf_path: &str, ref_config: &ReferenceConfig) -> Result<HashMap<String, Vec<f64>>, Box<dyn std::error::Error>> {
    let file = File::open(maf_path)?;
    let reader = BufReader::new(file);
    let mut mafs: HashMap<String, Vec<f64>> = HashMap::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with("chrom") || l.starts_with("WARNING") {
            continue;
        }
        let parts: Vec<&str> = l.split('\t').collect();
        // chrom pos major minor (count) -> or ref alt major minor ?
        // file format: chrom pos major minor ...
        if parts.len() < 4 {
            continue;
        }

        let chrom = parts[0].to_string();
        let pos: u32 = parts[1].parse().unwrap_or(0);
        let ref_ct: f64 = parts[2].parse().unwrap_or(0.0);
        let alt_ct: f64 = parts[3].parse().unwrap_or(0.0);

        let total = ref_ct + alt_ct;
        if total == 0.0 {
            continue;
        }
        let maf = ref_ct.min(alt_ct) / total;

        if maf > 0.1 {
            // MAF > 0.1 threshold filters out homozygous sites
            if let Some(segment) = get_segment_from_pos(&chrom, pos, ref_config) {
                mafs.entry(segment).or_default().push(maf);
            }
        }
    }
    Ok(mafs)
}

/// Parse MAF file returning BAF (alt / (ref+alt)) values grouped by segment for plotting.
/// No threshold filter — all sites that passed depth in calc_maf are included.
pub fn parse_maf_for_plot(maf_path: &str, ref_config: &ReferenceConfig) -> Result<HashMap<String, Vec<f64>>, Box<dyn std::error::Error>> {
    let file = File::open(maf_path)?;
    let reader = BufReader::new(file);
    let mut bafs: HashMap<String, Vec<f64>> = HashMap::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with("chrom") || l.starts_with("WARNING") {
            continue;
        }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }

        let chrom = parts[0].to_string();
        let pos: u32 = parts[1].parse().unwrap_or(0);
        let ref_ct: f64 = parts[2].parse().unwrap_or(0.0);
        let alt_ct: f64 = parts[3].parse().unwrap_or(0.0);

        let total = ref_ct + alt_ct;
        if total == 0.0 {
            continue;
        }
        let baf = alt_ct / total;

        if let Some(segment) = get_segment_from_pos(&chrom, pos, ref_config) {
            bafs.entry(segment).or_default().push(baf);
        }
    }
    Ok(bafs)
}

fn quantile(data: &[f64], q: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let pos = (data.len() - 1) as f64 * q;
    let base = pos.floor() as usize;
    let rest = pos - base as f64;
    if (base + 1) < data.len() {
        data[base] + rest * (data[base + 1] - data[base])
    } else {
        data[base]
    }
}

pub fn find_bin_width(data: &[f64]) -> f64 {
    // defaults: alpha=0.12, trim=95, deviation=35
    let alpha = 0.12;
    let deviation = 0.35;

    let iqr = quantile(data, 0.5 + deviation) - quantile(data, 0.5 - deviation);
    let n = data.len() as f64;
    let maxi = quantile(data, 0.95);

    let binsize = ((1.0 - alpha) * iqr + alpha * maxi / 2.0) * 2.0 * n.powf(-1.0 / 3.0);
    if binsize <= 0.0 { 0.01 } else { binsize }
}

pub fn find_levels(chrom_bins: &HashMap<String, Vec<f64>>) -> Vec<(f64, usize)> {
    // 1. Concatenate all data
    let mut all_values: Vec<f64> = chrom_bins.values().flatten().cloned().collect();
    if all_values.is_empty() {
        return Vec::new();
    }
    all_values.sort_by(|a, b| a.total_cmp(b));

    // 2. Bin width - round to nearest integer
    let bin_size_raw = find_bin_width(&all_values);
    let bin_size = bin_size_raw.round().max(1.0); // Round to integer, minimum 1
    debug!("Calculated histogram bin size: {:.0}", bin_size);

    // 3. Histogram - range determined by max of 95th percentile per segment
    let trim_percentile = 0.95;
    let max_val = chrom_bins
        .values()
        .filter(|v| !v.is_empty())
        .map(|vals| {
            let mut sorted = vals.clone();
            sorted.sort_by(|a, b| a.total_cmp(b));
            quantile(&sorted, trim_percentile)
        })
        .fold(0.0_f64, |a, b| a.max(b));

    let range_end = max_val + 2.0 * bin_size;
    // Integer bins: range(0, max + 2*binsize, binsize)
    let num_bins = ((range_end / bin_size).floor() as usize) + 1;
    let mut histogram = vec![0; num_bins];

    for v in &all_values {
        let bin_idx = (*v / bin_size).floor() as usize;
        if bin_idx < histogram.len() {
            histogram[bin_idx] += 1;
        }
    }

    // 4. Find peaks: local maxima with plateau handling + distance post-filter
    // For plateaus: find left edge, verify it's a local max, then take midpoint
    let mut local_maxima: Vec<(usize, usize)> = Vec::new(); // (bin_idx, count)
    let mut i = 1;
    while i < histogram.len().saturating_sub(1) {
        let v = histogram[i];
        // Check if this is the start of a plateau or local max
        if v >= histogram[i - 1] {
            // Find the extent of the plateau (consecutive equal values)
            let mut j = i;
            while j + 1 < histogram.len() && histogram[j + 1] == v {
                j += 1;
            }
            // Now i..=j is the plateau. Check if it's a local maximum
            // (left side < plateau AND right side < plateau)
            let is_left_lower = histogram[i - 1] < v;
            let is_right_lower = j + 1 >= histogram.len() || histogram[j + 1] < v;

            if is_left_lower && is_right_lower {
                // Take the midpoint of the plateau as the peak index
                let mid = (i + j) / 2;
                local_maxima.push((mid, v));
            }
            i = j + 1; // Skip past the plateau
        } else {
            i += 1;
        }
    }

    // Sort by height descending
    local_maxima.sort_by(|a, b| b.1.cmp(&a.1));

    // Apply distance=3 post-filter: keep tallest, remove others within distance
    let distance = 3usize;
    let mut kept_indices: Vec<usize> = Vec::new();
    for (idx, _) in &local_maxima {
        let dominated = kept_indices
            .iter()
            .any(|&kept| (*idx as isize - kept as isize).unsigned_abs() < distance);
        if !dominated {
            kept_indices.push(*idx);
        }
    }

    // Build peaks from kept indices
    let mut peaks: Vec<(f64, usize)> = Vec::new();
    for idx in &kept_indices {
        let center = (*idx as f64 * bin_size) + (bin_size / 2.0);
        peaks.push((center, histogram[*idx]));
    }

    // Sort by height descending (already sorted, but kept_indices order may differ)
    peaks.sort_by(|a, b| b.1.cmp(&a.1));

    // 5. Fine tune peaks
    // Median of points under peak +- 1.5 bins
    let mut refined_peaks = Vec::new();
    for (center, _count) in peaks {
        let low = center - bin_size * 1.5;
        let high = center + bin_size * 1.5;
        let mut sub: Vec<f64> = all_values
            .iter()
            .filter(|&v| *v > low && *v < high)
            .cloned()
            .collect();
        if sub.is_empty() {
            continue;
        }
        sub.sort_by(|a, b| a.total_cmp(b));
        let new_peak = sub[sub.len() / 2];

        let new_low = new_peak - bin_size;
        let new_high = new_peak + bin_size;
        let new_count = all_values
            .iter()
            .filter(|&v| *v > new_low && *v < new_high)
            .count();
        refined_peaks.push((new_peak, new_count));
    }

    if refined_peaks.is_empty() {
        return Vec::new();
    }

    // 6. Filter by height (> 5% of max)
    // Safety: refined_peaks is non-empty (checked above)
    let max_height = refined_peaks.iter().map(|x| x.1).max().unwrap();
    let threshold = max_height as f64 * 0.05;
    refined_peaks.retain(|x| x.1 as f64 > threshold);

    // 7. Filter: must have segments assigned
    // Calculate medians for each segment
    let mut medians: HashMap<String, f64> = HashMap::new();
    for (chr, vals) in chrom_bins {
        let mut v = vals.clone();
        v.sort_by(|a, b| a.total_cmp(b));
        if !v.is_empty() {
            medians.insert(chr.clone(), v[v.len() / 2]);
        }
    }

    // Assign each segment to its closest level (argmin of ratio difference)
    // Then keep only levels that have at least one AUTOSOMAL segment assigned
    // (levels with only sex chromosomes shouldn't define tumor ploidy)
    let level_values: Vec<f64> = refined_peaks.iter().map(|(v, _)| *v).collect();
    let mut level_autosomal_counts: HashMap<usize, usize> = HashMap::new();

    for (chr, med) in &medians {
        // Find closest level (argmin of |level/median - 1|)
        let mut best_idx = 0;
        let mut best_ratio_diff = f64::MAX;
        for (i, &lvl) in level_values.iter().enumerate() {
            let ratio = lvl / med;
            let diff = (ratio - 1.0).abs();
            if diff < best_ratio_diff {
                best_ratio_diff = diff;
                best_idx = i;
            }
        }
        // Only assign if within 0.3 tolerance AND is an autosome
        let is_autosome = chr != "X" && chr != "Y";
        if best_ratio_diff < 0.3 && is_autosome {
            *level_autosomal_counts.entry(best_idx).or_insert(0) += 1;
        }
    }

    // Keep only levels with at least one autosomal segment assigned
    let mut active_peaks = Vec::new();
    for (i, &(p_val, p_count)) in refined_peaks.iter().enumerate() {
        if level_autosomal_counts.get(&i).copied().unwrap_or(0) > 0 {
            active_peaks.push((p_val, p_count));
        }
    }

    // 8. Filter duplicates (close peaks)
    // Tolerance 9.5%
    // Peaks are sorted by height; keep the tallest, discard others within tolerance.

    let mut final_peaks = Vec::new();
    let mut discarded_indices = std::collections::HashSet::new();
    let tolerance = 0.095;

    for i in 0..active_peaks.len() {
        if discarded_indices.contains(&i) {
            continue;
        }

        // Keep i
        let (val_i, _) = active_peaks[i];
        final_peaks.push(active_peaks[i]);

        let lower = val_i * (1.0 - tolerance);
        let upper = val_i * (1.0 + tolerance);

        for (j, &(val_j, _)) in active_peaks.iter().enumerate().skip(i + 1) {
            if discarded_indices.contains(&j) {
                continue;
            }
            if val_j >= lower && val_j <= upper {
                discarded_indices.insert(j);
            }
        }
    }

    final_peaks
}

fn find_levels_maf(
    maf_data: &HashMap<String, Vec<f64>>,
    levels: &[(f64, usize)],
    medians: &HashMap<String, f64>,
) -> HashMap<usize, Vec<f64>> {
    // Only look at first 3 levels
    let top_levels: Vec<f64> = levels.iter().take(3).map(|x| x.0).collect();
    let mut levels_maf: HashMap<usize, Vec<f64>> = HashMap::new();

    // Assign segments to levels
    let mut lvls_seg: HashMap<usize, Vec<String>> = HashMap::new(); // level_idx -> vec<segment>

    for (seg, &med) in medians {
        // Find closest level
        let mut best_idx = 0;
        let mut min_diff = f64::MAX;

        for (i, &lvl) in top_levels.iter().enumerate() {
            let ratio = med / lvl;
            let diff = (ratio - 1.0).abs();
            if diff < min_diff {
                min_diff = diff;
                best_idx = i;
            }
        }

        if min_diff < 0.3 {
            lvls_seg
                .entry(best_idx)
                .or_default()
                .push(seg.clone());
        }
    }

    // Combine MAF
    for (idx, _lvl) in top_levels.iter().enumerate() {
        if let Some(segs) = lvls_seg.get(&idx) {
            let mut all_mafs = Vec::new();
            for seg in segs {
                if let Some(v) = maf_data.get(seg) {
                    all_mafs.extend(v);
                }
            }
            if all_mafs.len() >= 50 {
                // filter low count
                levels_maf.insert(idx, all_mafs);
            }
        }
    }

    levels_maf
}

fn find_maf_peak(levels_maf: &HashMap<usize, Vec<f64>>) -> HashMap<usize, f64> {
    let mut levels_maf_peak = HashMap::new();

    for (&idx, mafs) in levels_maf {
        if mafs.len() < 10 {
            continue;
        }

        let binsize = find_bin_width(mafs).min(0.05); // cap at 0.05

        // Histogram 0.0 to 0.6
        let num_bins = (0.6 / binsize).ceil() as usize;
        let mut hist = vec![0; num_bins + 1];

        for &m in mafs {
            if (0.0..=0.6).contains(&m) {
                let b = (m / binsize).floor() as usize;
                if b < hist.len() {
                    hist[b] += 1;
                }
            }
        }

        let mut peaks = Vec::new(); // (center, count)

        for i in 2..num_bins.saturating_sub(2) {
            let v = hist[i];
            if v > hist[i - 1] && v > hist[i - 2] && v > hist[i + 1] && v > hist[i + 2] {
                let center = (i as f64 * binsize) + (binsize / 2.0);
                peaks.push((center, v));
            }
        }

        if peaks.is_empty() {
            continue;
        }

        // Sort by height desc
        peaks.sort_by(|a, b| b.1.cmp(&a.1));
        let max_height = peaks[0].1;

        peaks.sort_by(|a, b| b.0.total_cmp(&a.0)); // desc pos

        for p in peaks {
            if p.1 as f64 > 0.5 * max_height as f64 {
                levels_maf_peak.insert(idx, p.0);
                break;
            }
        }
    }

    levels_maf_peak
}

/// Compute the number of variant sites per chromosome arm segment that fall
/// within enriched BED regions. Used to assess MAF data sufficiency.
pub fn compute_seg_bases(
    sites: &[Vec<Site>],
    enriched: &[BedRegion],
    ref_config: &ReferenceConfig,
) -> HashMap<String, usize> {
    let mapper = ContigMapper::new();

    // Group enriched intervals by chromosome index
    let mut bed_by_idx: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 24];
    for region in enriched {
        if let Some(idx) = mapper.get_chr_index(&region.segment) {
            if idx < 24 {
                bed_by_idx[idx].push((region.start, region.end));
            }
        }
    }

    let mut seg_counts: HashMap<String, usize> = HashMap::new();

    for (idx, chr_sites) in sites.iter().enumerate() {
        if idx >= 24 { break; }
        let chrom = match ContigMapper::chr_name_from_index(idx) {
            Some(name) => name,
            None => continue,
        };
        let intervals = &bed_by_idx[idx];
        if intervals.is_empty() { continue; }

        for site in chr_sites {
            let pos = site.pos as u32;
            let in_enriched = intervals.iter().any(|&(start, end)| pos >= start && pos <= end);
            if in_enriched {
                if let Some(seg) = get_segment_from_pos(chrom, pos, ref_config) {
                    *seg_counts.entry(seg).or_insert(0) += 1;
                }
            }
        }
    }

    seg_counts
}

pub fn call_karyotype(
    cov_path: &str,
    maf_path: Option<&str>,
    reads_aligned: Option<u64>,
    ref_config: &ReferenceConfig,
    thresholds: &KaryotypeThresholds,
    seg_bases: Option<&HashMap<String, usize>>,
) -> Result<KaryotypeOutput, Box<dyn std::error::Error>> {
    let bins = parse_coverage(cov_path)?;
    let mut warnings: Vec<String> = Vec::new();

    // Group by chrom (segment) using get_segment
    let mut chrom_bins: HashMap<String, Vec<f64>> = HashMap::new();
    for b in &bins {
        if let Some(segment) = get_segment(&b.chrom, b.start, b.end, ref_config) {
            chrom_bins
                .entry(segment)
                .or_default()
                .push(b.coverage);
        }
    }

    // Check spread early
    let spread = within_segment_spread(&chrom_bins);
    debug!("Within-segment spread: {:.4}", spread);
    if spread > thresholds.spread_warning {
        let w = format!(
            "High within-segment spread (>{:.2}): {:.4}. This may indicate poor quality data or waviness.",
            thresholds.spread_warning, spread
        );
        warn!("{}", w);
        warnings.push(w);
    }

    let levels = find_levels(&chrom_bins);
    info!("Found {} valid coverage levels.", levels.len());
    for (l, c) in &levels {
        debug!("Level: {:.4} (count {})", l, c);
    }

    // Calculate Medians
    let mut medians: HashMap<String, f64> = HashMap::new();
    for (chr, vals) in &chrom_bins {
        let mut v = vals.clone();
        v.sort_by(|a, b| a.total_cmp(b));
        if !v.is_empty() {
            medians.insert(chr.clone(), v[v.len() / 2]);
        }
    }

    // Process MAF
    let mut levels_maf_peaks: Option<HashMap<usize, f64>> = None;
    if let Some(mp) = maf_path {
        let maf_data = parse_maf(mp, ref_config)?;
        if !maf_data.is_empty() {
            // Check sufficiency: count segments with >1% MAF density
            let maf_sufficient = if let Some(sb) = seg_bases {
                let mut sufficient_segs = 0;
                for (seg, data) in &maf_data {
                    if let Some(&bases) = sb.get(seg) {
                        if bases > 0 {
                            let ratio = data.len() as f64 / bases as f64;
                            if ratio > 0.01 {
                                sufficient_segs += 1;
                            }
                        }
                    }
                }
                if sufficient_segs < 20 {
                    let w = "Not enough MAF data for karyotyping!".to_string();
                    warn!("{}", w);
                    warnings.push(w);
                    false
                } else {
                    true
                }
            } else {
                debug!("seg_bases not provided, skipping MAF sufficiency check");
                true
            };

            if maf_sufficient {
                let lvls_maf = find_levels_maf(&maf_data, &levels, &medians);
                let peaks = find_maf_peak(&lvls_maf);
                if !peaks.is_empty() {
                    levels_maf_peaks = Some(peaks);
                } else {
                    let w = "Not enough minor allele frequency data for confident automated karyotype estimation.".to_string();
                    warn!("{}", w);
                    warnings.push(w);
                }
            }
        }
    }

    // Infer CN states (predict_karyo_v2 + MAF logic)
    let (cn1, cn2, cn3) = resolve_cn_states(&levels, levels_maf_peaks.as_ref());

    debug!(
        "Estimated levels: CN1={:.2}, CN2={:.2}, CN3={:.2}",
        cn1, cn2, cn3
    );

    let delta = cn3 - cn2;

    // Assign CN to chromosomes
    let mut karyotype: HashMap<String, usize> = HashMap::new();

    for (chr, &med) in &medians {
        let cn = if med > cn2 - delta / 2.0 && med < cn2 + delta / 2.0 {
            2
        } else if med > cn1 - delta / 2.0 && med < cn1 + delta / 2.0 {
            1
        } else if med > cn3 - delta / 2.0 && med < cn3 + delta / 2.0 {
            3
        } else if med > cn3 + delta / 2.0 {
            // Est > 3
            ((med - cn2) / delta + 2.0).round() as usize
        } else if med < delta / 2.0 {
            0
        } else {
            // Ambiguous
            debug!(
                "Warning: Ambiguous copy number for {} (cov {:.2})",
                chr, med
            );
            // Closest?
            if (med - cn2).abs() < (med - cn1).abs() {
                2
            } else {
                1
            } // Naive fallback
        };

        karyotype.insert(chr.clone(), cn);
    }

    // ============ POST-ASSIGNMENT ADJUSTMENTS ============

    // Collect cn_depths for adjustment logic
    let mut cn_depths: HashMap<usize, Vec<f64>> = HashMap::new();
    for (chr, &cn) in &karyotype {
        let med = medians.get(chr).unwrap_or(&0.0);
        cn_depths.entry(cn).or_default().push(*med);
    }

    // Adjustment 0: X/Y adjustment for ambiguous sex chromosome calls
    // This handles cases where Y (and possibly X) fall outside clear CN bins,
    // especially in low-blast samples where germline 1n differs from tumor 1n
    let x_cn = karyotype.get("X").cloned();
    let y_cn_initial = karyotype.get("Y").cloned();
    let x_med = medians.get("X").cloned();
    let y_med = medians.get("Y").cloned();

    if let (Some(_y_cn_val), Some(y_med_val)) = (y_cn_initial, y_med) {
        // Get average diploid coverage from autosomes
        let cn2_depths: Vec<f64> = karyotype
            .iter()
            .filter(|(chr, cn)| {
                **cn == 2 && !chr.contains("X") && !chr.contains("Y") && !chr.contains("PAR")
            })
            .filter_map(|(chr, _)| medians.get(chr).cloned())
            .collect();

        if !cn2_depths.is_empty() {
            let cn2_avg = cn2_depths.iter().sum::<f64>() / cn2_depths.len() as f64;
            let germline_cn1_est = cn2_avg / 2.0;

            // Check if Y is around germline 1n (cn2/2)
            if y_med_val > germline_cn1_est * 0.9 && y_med_val < germline_cn1_est * 1.1 {
                debug!(
                    "Y coverage ({:.2}) ~ cn2/2 ({:.2}); setting Y to 1n",
                    y_med_val, germline_cn1_est
                );
                karyotype.insert("Y".to_string(), 1);
                let y_cn1 = y_med_val;

                // Now adjust X if it was ambiguous (assigned 1 or is ambiguous)
                if let (Some(x_cn_val), Some(x_med_val)) = (x_cn, x_med)
                    && (x_cn_val == 1 || x_cn_val == 2) {
                        // Calculate preliminary blast ratio from autosomes for X adjustment
                        let mut cn_depths_auto: HashMap<usize, Vec<f64>> = HashMap::new();
                        for (chr, &cn) in &karyotype {
                            if chr.contains("X") || chr.contains("Y") || chr.contains("PAR") {
                                continue;
                            }
                            if let Some(&med) = medians.get(chr) {
                                cn_depths_auto.entry(cn).or_default().push(med);
                            }
                        }

                        if let Some(blast_ratio) = estimate_blast_ratio(&cn_depths_auto) {
                            // Estimate expected X coverage if it were 2n
                            // For XX: 2 copies in both tumor and germline
                            // For XY with low blast: blend of tumor (2n in blast) and germline (1n)
                            let estimated_x_cn2 =
                                2.0 * y_cn1 * blast_ratio + (1.0 - blast_ratio) * y_cn1;

                            if x_med_val > estimated_x_cn2 * 0.9
                                && x_med_val < estimated_x_cn2 * 1.1
                            {
                                debug!(
                                    "X coverage ({:.2}) ~ expected 2n ({:.2}); setting X to 2n (with blast ratio {:.2})",
                                    x_med_val, estimated_x_cn2, blast_ratio
                                );
                                karyotype.insert("X".to_string(), 2);
                            } else {
                                debug!(
                                    "X coverage ({:.2}) not matching expected 2n ({:.2}); setting X to 1n (with blast ratio {:.2})",
                                    x_med_val, estimated_x_cn2, blast_ratio
                                );
                                karyotype.insert("X".to_string(), 1);
                            }
                        } else {
                            // No blast ratio available, use simple distance check
                            if (cn2_avg - x_med_val).abs() > (x_med_val - germline_cn1_est).abs() {
                                debug!(
                                    "X coverage ({:.2}) closer to 1n; setting X to 1n (no blast ratio)",
                                    x_med_val
                                );
                                karyotype.insert("X".to_string(), 1);
                            } else {
                                debug!(
                                    "X coverage ({:.2}) closer to 2n; setting X to 2n (no blast ratio)",
                                    x_med_val
                                );
                                karyotype.insert("X".to_string(), 2);
                            }
                        }
                    }
            }
        }
    }

    // Rebuild cn_depths after X/Y adjustment
    cn_depths.clear();
    for (chr, &cn) in &karyotype {
        let med = medians.get(chr).unwrap_or(&0.0);
        cn_depths.entry(cn).or_default().push(*med);
    }

    // Adjustment 1: Low-penetrance hypodiploid
    // When Y is the only 1n and "3n" ≈ 2× "1n", we misidentified levels
    // The "diploid" is actually monosomic loss, and "triploid" is actually diploid
    let y_cn = karyotype.get("Y").cloned().unwrap_or(0);
    let non_sex_1n_count = karyotype
        .iter()
        .filter(|(chr, _)| !chr.contains("X") && !chr.contains("Y"))
        .filter(|(_, cn)| **cn == 1)
        .count();

    // Critical: Only trigger if Y genuinely fits in the TUMOR cn1 bin, not just assigned
    // via distance fallback or X/Y adjustment. This prevents false triggering on
    // hyperdiploid samples where Y is at germline 1n (cn2/2) but not tumor 1n.
    let y_med = medians.get("Y").cloned().unwrap_or(0.0);
    let y_in_tumor_cn1_bin = y_med > cn1 - delta / 2.0 && y_med < cn1 + delta / 2.0;

    // Also skip if MAF data confirms the lower level is genuinely diploid (MAF peak > 0.4)
    // MAF-based karyotype prediction doesn't use hypodiploid adjustment
    let maf_confirms_diploid = levels_maf_peaks
        .as_ref()
        .and_then(|peaks| peaks.get(&0)) // Level 0 is the cn2 level (highest count)
        .map(|&peak| peak > 0.4)
        .unwrap_or(false);

    if y_cn == 1 && non_sex_1n_count == 0 && y_in_tumor_cn1_bin && !maf_confirms_diploid {
        // Check if "3n" average ≈ 2× "1n" average
        let cn1_avg = cn_depths
            .get(&1)
            .map(|v| v.iter().sum::<f64>() / v.len() as f64)
            .unwrap_or(0.0);
        let cn3_avg = cn_depths
            .get(&3)
            .map(|v| v.iter().sum::<f64>() / v.len() as f64)
            .unwrap_or(0.0);

        if cn3_avg > 0.0 && cn1_avg > 0.0 {
            let ratio = cn3_avg / 2.0;
            if ratio > cn1_avg * 0.9 && ratio < cn1_avg * 1.1 {
                debug!("---------- Adjusting to low-penetrance hypodiploid! ----------");
                debug!(
                    "CN1 avg: {:.2}, CN3 avg: {:.2}, CN3/2: {:.2}",
                    cn1_avg, cn3_avg, ratio
                );
                // Shift: 2→1 (diploid becomes haploid loss), 3→2 (triploid becomes diploid)
                for (_, cn) in karyotype.iter_mut() {
                    if *cn == 2 {
                        *cn = 1;
                    } else if *cn == 3 {
                        *cn = 2;
                    }
                }
            }
        }
    }

    // Adjustment 2: High 3n count (>25 chromosomes at 3n, none at 1n)
    // Suggests levels should be shifted down
    let cn1_count = karyotype.values().filter(|&&cn| cn == 1).count();
    let cn3_count = karyotype.values().filter(|&&cn| cn == 3).count();

    if cn1_count == 0 && cn3_count > 25 {
        debug!("---------- Adjusting: >25 at 3n with none at 1n, shifting down ----------");
        for (_, cn) in karyotype.iter_mut() {
            if *cn > 0 {
                *cn -= 1;
            }
        }
    }

    // ============ END ADJUSTMENTS ============

    // Estimate blast ratio
    // We need to collect medians by assigned CN, excluding X and Y
    let mut cn_depths_auto: HashMap<usize, Vec<f64>> = HashMap::new();

    for (chr, &cn) in &karyotype {
        // Exclude X, Y, PAR, M...
        if chr.contains("X") || chr.contains("Y") || chr.contains("M") || chr.contains("PAR") {
            continue;
        }

        let med = medians.get(chr).unwrap_or(&0.0);
        cn_depths_auto.entry(cn).or_default().push(*med);
    }

    let blast_ratio = estimate_blast_ratio(&cn_depths_auto);
    if let Some(br) = blast_ratio {
        info!("Estimated blast ratio: {:.4}", br);
    } else {
        warn!("Could not estimate blast ratio (insufficient data).");
    }

    // Generate Karyotype String
    let karyo_str = generate_karyotype_string(&karyotype);
    info!("Karyotype: {}", karyo_str);

    // Generate ISCN String
    let iscn_str = generate_iscn_string(&karyotype);
    info!("ISCN: {}", iscn_str);

    // Build output struct
    Ok(KaryotypeOutput {
        karyotype,
        medians,
        levels_found: Some(levels),
        karyotype_string: karyo_str,
        iscn_string: iscn_str,
        blast_ratio,
        reads_aligned,
        within_segment_spread: Some(spread),
        warnings: if warnings.is_empty() { None } else { Some(warnings.join("\n")) },
        maf_peaks: levels_maf_peaks,
    })
}

fn within_segment_spread(chrom_bins: &HashMap<String, Vec<f64>>) -> f64 {
    // 1. Calculate IQR for each segment
    let mut keys: Vec<&String> = chrom_bins.keys().collect();
    keys.sort_by(|a, b| human_sort(a, b));

    let mut seg_iqrs = Vec::new();
    for k in &keys {
        let vals = &chrom_bins[*k];
        let mut v = vals.clone();
        v.sort_by(|a, b| a.total_cmp(b));
        if !v.is_empty() {
            let iqr = quantile(&v, 0.75) - quantile(&v, 0.25);
            seg_iqrs.push(iqr);
        }
    }

    // 2. Remove last one (Y approx) if we have enough
    if seg_iqrs.len() > 1 {
        seg_iqrs.pop();
    }

    // 3. IQR of IQRs
    seg_iqrs.sort_by(|a, b| a.total_cmp(b));
    let spread_val = quantile(&seg_iqrs, 0.75) - quantile(&seg_iqrs, 0.25);

    // 4. Global median
    let mut all_values: Vec<f64> = chrom_bins.values().flatten().cloned().collect();
    all_values.sort_by(|a, b| a.total_cmp(b));
    let median = if !all_values.is_empty() {
        all_values[all_values.len() / 2]
    } else {
        1.0 // Avoid div by zero
    };

    if median == 0.0 {
        0.0
    } else {
        spread_val / median
    }
}

fn resolve_cn_states(
    levels: &[(f64, usize)],
    maf_peaks: Option<&HashMap<usize, f64>>,
) -> (f64, f64, f64) {
    if levels.is_empty() {
        return (1.0, 2.0, 3.0);
    } // Fallback

    // find_levels returns peaks sorted by height/count desc
    // We need to sort them by value to distinguish low vs high coverage levels
    let mut sorted_levels: Vec<f64> = levels.iter().map(|x| x.0).collect();

    // If we have more than 2 levels, take the top 2 by count (which are the first 2 in `levels`)
    if levels.len() >= 2 {
        sorted_levels = vec![levels[0].0, levels[1].0];
    }
    sorted_levels.sort_by(|a, b| a.total_cmp(b));

    let v1 = sorted_levels[0];
    let v2 = if sorted_levels.len() > 1 {
        sorted_levels[1]
    } else {
        v1
    };

    // If only 1 level (or very close), assume 2n unless MAF indicates otherwise
    if sorted_levels.len() < 2 || (v2 - v1).abs() < 0.05 {
        // Check MAF
        if let Some(peaks) = maf_peaks {
            // We need to map back to the index in `levels`.
            // If len < 2, v1 corresponds to levels[0] (which has index 0 in our map).
            if let Some(&peak) = peaks.get(&0) {
                if peak > 0.4 {
                    debug!("Single level, MAF peak {:.2} confirms 2n", peak);
                } else {
                    debug!(
                        "Single level, MAF peak {:.2} suggests NOT 2n (likely 1n)",
                        peak
                    );
                }
            }
        }
        // Default to 2n behavior for single peak
        return (v1 * 0.5, v1, v1 * 1.5);
    }

    // Two levels: v1 (lower), v2 (higher)
    // Find closest original indices
    let idx_v1 = levels
        .iter()
        .position(|x| (x.0 - v1).abs() < 0.001)
        .unwrap_or(0);
    let idx_v2 = levels
        .iter()
        .position(|x| (x.0 - v2).abs() < 0.001)
        .unwrap_or(0);

    if let Some(peaks) = maf_peaks {
        // We have MAF data
        let p1 = peaks.get(&idx_v1).cloned().unwrap_or(0.0);
        let p2 = peaks.get(&idx_v2).cloned().unwrap_or(0.0);

        debug!("Level 1 ({:.2}, idx {}): MAF peak {:.2}", v1, idx_v1, p1);
        debug!("Level 2 ({:.2}, idx {}): MAF peak {:.2}", v2, idx_v2, p2);

        // Case 0: BOTH levels have high MAF (>0.4) - compare them to decide
        // This is critical for hypodiploid cases where monosomic regions also have MAF ~0.5
        // BUT only use MAF comparison if the difference is meaningful (> 0.03)
        if p1 > 0.4 && p2 > 0.4 {
            let maf_diff = (p2 - p1).abs();
            if maf_diff < 0.03 {
                // MAF difference too small to be meaningful - fall back to ratio logic
                debug!(
                    "Both MAF > 0.4 but diff too small ({:.3}). Falling back to ratio logic.",
                    maf_diff
                );
                // Don't return here - let it fall through to the ratio-based logic below
            } else if p2 > p1 {
                // Higher level has higher MAF → v2 is more likely diploid, v1 is monosomic
                // BUT add ratio sanity check: if ratio < 1.6, levels are likely 2n/3n not 1n/2n
                if v2 / v1 >= 1.6 {
                    debug!(
                        "Both MAF > 0.4, p2 ({:.2}) > p1 ({:.2}), ratio >= 1.6. v1=1n, v2=2n. cn2/cn1 = {:.2}",
                        p2,
                        p1,
                        v2 / v1
                    );
                    return (v1, v2, v2 + (v2 - v1));
                } else {
                    // Ratio too low for 1n/2n - MAF difference is likely noise
                    debug!(
                        "Both MAF > 0.4, p2 ({:.2}) > p1 ({:.2}), but ratio < 1.6. Likely 2n/3n. cn3/cn2 = {:.2}",
                        p2,
                        p1,
                        v2 / v1
                    );
                    return (2.0 * v1 - v2, v1, v2);
                }
            } else {
                // Lower level has higher/equal MAF → v1 is diploid (with ratio sanity check)
                if v2 / v1 >= 1.6 {
                    debug!(
                        "Both MAF > 0.4, p1 >= p2, but ratio >= 1.6. v1=1n, v2=2n. cn2/cn1 = {:.2}",
                        v2 / v1
                    );
                    return (v1, v2, v2 + (v2 - v1));
                } else {
                    debug!(
                        "Both MAF > 0.4, p1 ({:.2}) >= p2 ({:.2}). v1=2n, v2=3n. cn3/cn2 = {:.2}",
                        p1,
                        p2,
                        v2 / v1
                    );
                    return (2.0 * v1 - v2, v1, v2);
                }
            }
        }

        // Case 1: Only lower level has high MAF (>0.4)
        if p1 > 0.4 {
            // Coverage ratio sanity check: if v2/v1 >= 1.6, the levels are actually 1n/2n
            // (the high MAF on lower level is misleading, possibly due to low blast fraction)
            if v2 / v1 >= 1.6 {
                debug!(
                    "MAF suggests v1 is 2n, but ratio >= 1.6. Actually v1=1n, v2=2n. cn2/cn1 = {:.2}",
                    v2 / v1
                );
                return (v1, v2, v2 + (v2 - v1));
            } else {
                debug!(
                    "MAF suggests v1 is 2n (peak > 0.4). v2 likely 3n. cn3/cn2 = {:.2}",
                    v2 / v1
                );
                // cn1 = cn2 - (cn3 - cn2) = 2*v1 - v2
                return (2.0 * v1 - v2, v1, v2);
            }
        }

        // Case 2: Only higher level has high MAF (>0.4)
        if p2 > 0.4 {
            debug!(
                "MAF suggests v2 is 2n (peak > 0.4). v1 likely 1n. cn2/cn1 = {:.2}",
                v2 / v1
            );
            return (v1, v2, v2 + (v2 - v1));
        }

        // Ambiguous MAF? Fallback to ratio logic
        debug!("MAF ambiguous. Falling back to ratio logic.");
    }

    // heuristic fallback logic
    if v2 < v1 * 1.51 {
        (2.0 * v1 - v2, v1, v2)
    } else {
        (v1, v2, v2 + (v2 - v1))
    }
}

fn estimate_blast_ratio(cn_depths: &HashMap<usize, Vec<f64>>) -> Option<f64> {
    let get_mean = |k: usize| -> Option<f64> {
        if let Some(vals) = cn_depths.get(&k)
            && !vals.is_empty() {
                let sum: f64 = vals.iter().sum();
                return Some(sum / vals.len() as f64);
            }
        None
    };

    let get_count = |k: usize| -> usize { cn_depths.get(&k).map(|v| v.len()).unwrap_or(0) };

    let cn1 = get_mean(1);
    let cn2 = get_mean(2);
    let cn3 = get_mean(3);

    let cn1_count = get_count(1);
    let cn2_count = get_count(2);
    let cn3_count = get_count(3);

    debug!("{} segments are 2n", cn2_count);
    debug!("{} segments are 1n", cn1_count);
    debug!("{} segments are 3n", cn3_count);

    // If all three are available, prefer the formula with more supporting segments
    if let (Some(v1), Some(v2), Some(v3)) = (cn1, cn2, cn3) {
        if cn1_count >= cn3_count {
            // Use 1n and 2n: (1 - CN1/CN2) * 2
            let ratio = (1.0 - (v1 / v2)) * 2.0;
            debug!("blast_ratio from 1n 2n: {:.4}", ratio);
            return Some(ratio);
        } else {
            // Use 2n and 3n: (CN3/CN2 * 2) - 2
            let ratio = (v3 / v2 * 2.0) - 2.0;
            debug!("blast_ratio from 2n 3n: {:.4}", ratio);
            return Some(ratio);
        }
    }

    // Fallback: use whichever pair is available
    if let (Some(v1), Some(v2)) = (cn1, cn2) {
        // (1 - (CN1/CN2)) * 2
        let ratio = (1.0 - (v1 / v2)) * 2.0;
        debug!("blast_ratio from 1n 2n: {:.4}", ratio);
        return Some(ratio);
    }

    if let (Some(v2), Some(v3)) = (cn2, cn3) {
        // (CN3/CN2 * 2) - 2
        let ratio = (v3 / v2 * 2.0) - 2.0;
        debug!("blast_ratio from 2n 3n: {:.4}", ratio);
        return Some(ratio);
    }

    debug!("There is only autosomal diploid, we cannot estimate blast ratio");
    None
}

// --------------------------------------------------------------------------------
// GC Bias Correction
// --------------------------------------------------------------------------------

/// Fit a simple linear regression y = m*x + b, returning (m, b).
fn linear_regression(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let sum_x: f64 = xs.iter().sum();
    let sum_y: f64 = ys.iter().sum();
    let sum_xy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();
    let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
    let denom = n * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-12 {
        return (0.0, sum_y / n);
    }
    let m = (n * sum_xy - sum_x * sum_y) / denom;
    let b = (sum_y - m * sum_x) / n;
    (m, b)
}

/// Predict a single value using LOESS (locally weighted linear regression).
///
/// At each query point, selects the nearest `bandwidth` fraction of training data,
/// applies tricube kernel weights, and fits a local weighted linear regression.
fn loess_predict_single(xs: &[f64], ys: &[f64], query: f64, bandwidth: f64) -> f64 {
    let n = xs.len();
    let k = ((n as f64 * bandwidth).ceil() as usize).max(3).min(n);

    // Find k nearest neighbors by distance to query
    let mut dists: Vec<(usize, f64)> = xs.iter().enumerate()
        .map(|(i, &x)| (i, (x - query).abs()))
        .collect();
    dists.sort_by(|a, b| a.1.total_cmp(&b.1));
    let max_dist = dists[k - 1].1;

    if max_dist < 1e-12 {
        // All neighbors at same point, return weighted mean
        let sum: f64 = dists[..k].iter().map(|(i, _)| ys[*i]).sum();
        return sum / k as f64;
    }

    // Tricube kernel weights: w(u) = (1 - |u|^3)^3 for |u| < 1
    let weights: Vec<f64> = dists[..k].iter().map(|(_, d)| {
        let u = d / max_dist;
        let t = 1.0 - u * u * u;
        t * t * t
    }).collect();

    // Weighted linear regression
    let sum_w: f64 = weights.iter().sum();
    let sum_wx: f64 = dists[..k].iter().zip(&weights).map(|((i, _), w)| w * xs[*i]).sum();
    let sum_wy: f64 = dists[..k].iter().zip(&weights).map(|((i, _), w)| w * ys[*i]).sum();
    let sum_wxx: f64 = dists[..k].iter().zip(&weights).map(|((i, _), w)| w * xs[*i] * xs[*i]).sum();
    let sum_wxy: f64 = dists[..k].iter().zip(&weights).map(|((i, _), w)| w * xs[*i] * ys[*i]).sum();

    let denom = sum_w * sum_wxx - sum_wx * sum_wx;
    if denom.abs() < 1e-12 {
        return sum_wy / sum_w;
    }
    let slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denom;
    let intercept = (sum_wy - slope * sum_wx) / sum_w;

    slope * query + intercept
}

/// Generate LOESS fitted curve across the typical GC range for plotting.
fn loess_fit_curve(xs: &[f64], ys: &[f64], bandwidth: f64) -> Vec<(f64, f64)> {
    let n_points = 50;
    let gc_min = 0.3;
    let gc_max = 0.7;
    let step = (gc_max - gc_min) / (n_points - 1) as f64;
    (0..n_points)
        .map(|i| {
            let gc = gc_min + i as f64 * step;
            let predicted = loess_predict_single(xs, ys, gc, bandwidth);
            (gc, predicted)
        })
        .collect()
}

const LOESS_BANDWIDTH: f64 = 0.3;
const GC_REFERENCE: f64 = 0.41;

/// Apply GC bias correction to coverage bins.
///
/// 1. Identify majority ploidy from the initial karyotype
/// 2. Select bins belonging to majority-ploidy autosomal segments
/// 3. Fit regression model (linear or LOESS) on those bins
/// 4. Correct each bin: adjusted = observed * reference_at_0.41 / predicted_at_gc
fn gc_correct_coverage(
    bins: &[CoverageBin],
    karyotype: &HashMap<String, usize>,
    ref_config: &ReferenceConfig,
    method: GcCorrectionMethod,
) -> (Vec<CoverageBin>, Vec<(f64, f64)>) {
    // Find majority CN among autosomal segments
    let mut cn_counts: HashMap<usize, usize> = HashMap::new();
    for (seg, &cn) in karyotype {
        if seg.contains("X") || seg.contains("Y") || seg.contains("M") || seg.contains("PAR") {
            continue;
        }
        *cn_counts.entry(cn).or_default() += 1;
    }
    let majority_cn = cn_counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(cn, _)| cn)
        .unwrap_or(2);
    debug!("GC correction: majority ploidy = {}", majority_cn);

    // Collect majority-ploidy segments
    let majority_segments: std::collections::HashSet<String> = karyotype
        .iter()
        .filter(|(seg, cn)| {
            **cn == majority_cn
                && !seg.contains("X")
                && !seg.contains("Y")
                && !seg.contains("M")
                && !seg.contains("PAR")
        })
        .map(|(seg, _)| seg.clone())
        .collect();

    // Collect (gc, cov) pairs from bins in majority-ploidy segments
    let mut gc_vals = Vec::new();
    let mut cov_vals = Vec::new();
    for b in bins {
        if let (Some(gc), Some(segment)) =
            (b.gc_content, get_segment(&b.chrom, b.start, b.end, ref_config))
            && majority_segments.contains(&segment) {
                gc_vals.push(gc);
                cov_vals.push(b.coverage);
            }
    }

    if gc_vals.len() < 3 {
        warn!("GC correction: insufficient data points ({}), skipping", gc_vals.len());
        return (bins.to_vec(), Vec::new());
    }

    // Fit model and build prediction function + curve for plotting
    let (reference_cov, predict_fn, curve_points): (f64, Box<dyn Fn(f64) -> f64>, Vec<(f64, f64)>) =
        match method {
            GcCorrectionMethod::Linear => {
                let (m, b) = linear_regression(&gc_vals, &cov_vals);
                debug!("GC correction (linear): m={:.4}, b={:.4} from {} bins", m, b, gc_vals.len());
                let reference = m * GC_REFERENCE + b;
                let curve = vec![(0.3, m * 0.3 + b), (0.7, m * 0.7 + b)];
                (reference, Box::new(move |gc| m * gc + b), curve)
            }
            GcCorrectionMethod::Loess => {
                let reference = loess_predict_single(&gc_vals, &cov_vals, GC_REFERENCE, LOESS_BANDWIDTH);
                let curve = loess_fit_curve(&gc_vals, &cov_vals, LOESS_BANDWIDTH);
                debug!("GC correction (LOESS): reference={:.4} at GC={}, bandwidth={}, {} bins",
                    reference, GC_REFERENCE, LOESS_BANDWIDTH, gc_vals.len());
                let xs = gc_vals.clone();
                let ys = cov_vals.clone();
                (reference, Box::new(move |gc| loess_predict_single(&xs, &ys, gc, LOESS_BANDWIDTH)), curve)
            }
            GcCorrectionMethod::None => unreachable!(),
        };

    // Apply correction: adjusted = observed * reference / predicted
    let corrected = bins.iter()
        .map(|bin| {
            let mut corrected = bin.clone();
            if let Some(gc) = bin.gc_content
                && bin.coverage > 0.0 {
                    let predicted = predict_fn(gc);
                    if predicted > 0.0 {
                        corrected.coverage = bin.coverage * reference_cov / predicted;
                    }
                    if corrected.coverage < 0.0 {
                        corrected.coverage = 0.0;
                    }
                }
            corrected
        })
        .collect();
    (corrected, curve_points)
}

/// Write coverage bins to a TSV file.
fn write_coverage_tsv(
    path: &str,
    bins: &[CoverageBin],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    writeln!(file, "chromosome\tbin_start\tbin_end\tn_reads\tgc_content")?;
    for b in bins {
        let gc_str = match b.gc_content {
            Some(gc) => format!("{:.4}", gc),
            None => "NaN".to_string(),
        };
        writeln!(file, "{}\t{}\t{}\t{:.4}\t{}", b.chrom, b.start, b.end, b.coverage, gc_str)?;
    }
    Ok(())
}

/// Two-pass karyotype calling with GC bias correction.
///
/// Pass 1: Initial karyotype estimate from raw coverage
/// Pass 2: GC-correct coverage using majority-ploidy bins, then re-run karyotype
///
/// When `gc_correction` is `None`, skips GC correction entirely (returns Pass 1 result).
pub fn call_karyotype_gc_corrected(
    cov_path: &str,
    maf_path: Option<&str>,
    out_prefix: &str,
    reads_aligned: Option<u64>,
    ref_config: &ReferenceConfig,
    thresholds: &KaryotypeThresholds,
    gc_correction: GcCorrectionMethod,
    seg_bases: Option<&HashMap<String, usize>>,
) -> Result<KaryotypeOutput, Box<dyn std::error::Error>> {
    // Pass 1: initial karyotype
    info!("=== GC Correction Pass 1: Initial karyotype ===");
    let initial_output = call_karyotype(cov_path, maf_path, reads_aligned, ref_config, thresholds, seg_bases)?;

    // If no GC correction requested, return Pass 1 directly
    if gc_correction == GcCorrectionMethod::None {
        info!("GC correction disabled, returning Pass 1 result");
        return Ok(initial_output);
    }

    let karyotype = &initial_output.karyotype;

    if karyotype.is_empty() {
        warn!("GC correction: no karyotype data, skipping correction");
        return Ok(initial_output);
    }

    // Parse bins (with GC content)
    let bins = parse_coverage(cov_path)?;

    // Plot pre-correction karyotype
    let pre_plot_path = format!("{}.karyotype.svg", out_prefix);
    plot::plot_karyotype(&bins, ref_config, &pre_plot_path, "Karyotype (raw)");

    // Apply GC correction
    info!("=== GC Correction ({:?}): Adjusting coverage ===", gc_correction);
    let (corrected_bins, curve_points) = gc_correct_coverage(&bins, karyotype, ref_config, gc_correction);

    // Plot GC vs coverage (raw)
    let gc_plot_path = format!("{}.gc_vs_coverage.svg", out_prefix);
    plot::plot_gc_vs_coverage(&bins, ref_config, karyotype, &curve_points, &gc_plot_path, "Coverage vs GC Content");

    // Plot GC vs coverage (adjusted)
    let gc_adj_plot_path = format!("{}.gc_vs_coverage.gc_corrected.svg", out_prefix);
    plot::plot_gc_vs_coverage(&corrected_bins, ref_config, karyotype, &curve_points, &gc_adj_plot_path, "Coverage vs GC Content (GC Corrected)");

    // Write corrected coverage
    let corrected_cov_path = format!("{}.coverage.gc_adjusted.tsv", out_prefix);
    write_coverage_tsv(&corrected_cov_path, &corrected_bins)?;
    info!("Wrote GC-adjusted coverage to {}", corrected_cov_path);

    // Pass 2: re-run karyotype on corrected coverage
    info!("=== GC Correction Pass 2: Corrected karyotype ===");
    let result = call_karyotype(
        &corrected_cov_path,
        maf_path,
        reads_aligned,
        ref_config,
        thresholds,
        seg_bases,
    );

    // Plot post-correction karyotype
    let post_plot_path = format!("{}.karyotype.gc_corrected.svg", out_prefix);
    plot::plot_karyotype(&corrected_bins, ref_config, &post_plot_path, "Karyotype (GC corrected)");

    // Plot combined karyotype + BAF
    if let Some(maf_p) = maf_path {
        match parse_maf_for_plot(maf_p, ref_config) {
            Ok(baf_data) => {
                let combined_path = format!("{}.karyotype_baf.gc_corrected.svg", out_prefix);
                plot::plot_karyotype_with_baf(
                    &corrected_bins, &baf_data, ref_config,
                    &combined_path, "Karyotype + BAF (GC corrected)",
                );
            }
            Err(e) => warn!("Could not parse MAF for BAF plot: {}", e),
        }
    }

    result
}

