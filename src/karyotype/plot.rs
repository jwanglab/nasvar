use crate::config::ReferenceConfig;
use log::{info, warn, error};
use crate::plotting::prelude::*;
use std::collections::HashMap;

use super::{get_segment, CoverageBin};

/// Canonical segment order for karyotype x-axis layout.
const SEGMENT_ORDER: &[&str] = &[
    "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
    "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q",
    "13", "14", "15", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q",
    "20p", "20q", "21", "22", "X", "Y",
];

const GAP: f64 = 3.0;

/// Downsample scatter data by binning into pixel-width buckets and keeping at most
/// `max_per_bucket` evenly-spaced points per bucket. Small datasets pass through unchanged.
fn downsample_scatter(
    x: &[f64],
    y: &[f64],
    plot_width_px: f64,
    x_min: f64,
    x_max: f64,
    max_per_bucket: usize,
) -> (Vec<f64>, Vec<f64>) {
    let target = plot_width_px as usize * max_per_bucket;
    if x.len() <= target {
        return (x.to_vec(), y.to_vec());
    }

    let x_range = x_max - x_min;
    if x_range <= 0.0 {
        return (x.to_vec(), y.to_vec());
    }

    // Bucket points by pixel column
    let n_buckets = plot_width_px as usize;
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n_buckets + 1];
    for (i, &xi) in x.iter().enumerate() {
        let col = (((xi - x_min) / x_range) * plot_width_px) as usize;
        let col = col.min(n_buckets);
        buckets[col].push(i);
    }

    // Stratified sampling: take evenly spaced indices from each bucket
    let mut out_x = Vec::with_capacity(target);
    let mut out_y = Vec::with_capacity(target);
    for bucket in &buckets {
        let n = bucket.len();
        if n <= max_per_bucket {
            for &idx in bucket {
                out_x.push(x[idx]);
                out_y.push(y[idx]);
            }
        } else {
            let stride = n as f64 / max_per_bucket as f64;
            for j in 0..max_per_bucket {
                let idx = bucket[(j as f64 * stride) as usize];
                out_x.push(x[idx]);
                out_y.push(y[idx]);
            }
        }
    }

    (out_x, out_y)
}

/// Plot a karyotype coverage scatter with segment medians.
pub fn plot_karyotype(
    bins: &[CoverageBin],
    ref_config: &ReferenceConfig,
    out_path: &str,
    title: &str,
) {
    // Group bins by segment
    let mut seg_bins: HashMap<String, Vec<f64>> = HashMap::new();
    for b in bins {
        if let Some(seg) = get_segment(&b.chrom, b.start, b.end, ref_config) {
            seg_bins.entry(seg).or_default().push(b.coverage);
        }
    }

    // Lay out segments on x-axis
    let mut all_x: Vec<f64> = Vec::new();
    let mut all_y: Vec<f64> = Vec::new();
    let mut median_x_start: Vec<f64> = Vec::new();
    let mut median_x_end: Vec<f64> = Vec::new();
    let mut median_y: Vec<f64> = Vec::new();
    let mut separator_x: Vec<f64> = Vec::new();
    let mut tick_positions: Vec<f64> = Vec::new();
    let mut tick_labels: Vec<String> = Vec::new();

    let mut x = 0.0;
    for &seg in SEGMENT_ORDER {
        let vals = match seg_bins.get(seg) {
            Some(v) if !v.is_empty() => v,
            _ => continue,
        };

        let seg_start = x;
        for (i, &cov) in vals.iter().enumerate() {
            all_x.push(x + i as f64);
            all_y.push(cov);
        }
        let seg_end = x + (vals.len() - 1) as f64;

        // Median
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = sorted[sorted.len() / 2];
        median_x_start.push(seg_start);
        median_x_end.push(seg_end);
        median_y.push(med);

        // Tick at center
        tick_positions.push((seg_start + seg_end) / 2.0);
        tick_labels.push(seg.to_string());

        x = seg_end + GAP;
        separator_x.push(x - GAP / 2.0);
    }

    if all_y.is_empty() {
        warn!("plot_karyotype: no data to plot");
        return;
    }

    // Y-limit: 95th percentile
    let mut sorted_y = all_y.clone();
    sorted_y.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let y_max = sorted_y[(sorted_y.len() as f64 * 0.99) as usize].max(1.0);
    let x_max = x;

    let mut fig = Figure::new(1600.0, 500.0);
    let ax = fig.gca();

    // Scatter: coverage bins
    ax.scatter(&all_x, &all_y)
        .color("blue")
        .size(3.0)
        .alpha(0.5)
        .edge_width(0.0)
        .build();

    // Median lines (draw as short line segments)
    for i in 0..median_y.len() {
        ax.plot(
            [median_x_start[i], median_x_end[i]],
            [median_y[i], median_y[i]],
        )
        .color("red")
        .linewidth(2.5)
        .build();
    }

    // Vertical separator lines
    for &sx in &separator_x {
        ax.plot([sx, sx], [0.0, y_max])
            .color("black")
            .linewidth(0.5)
            .build();
    }

    ax.set_xlim(-(GAP / 2.0), x_max - GAP / 2.0);
    ax.set_ylim(0.0, y_max);
    ax.set_title(title);
    ax.set_ylabel("Reads per Mbp");
    ax.x_axis.tick_positions = Some(tick_positions);
    ax.x_axis.tick_labels = Some(tick_labels);

    if let Err(e) = fig.save(out_path) {
        error!("Error saving karyotype plot to {}: {}", out_path, e);
    } else {
        info!("Saved karyotype plot: {}", out_path);
    }
}

/// Plot GC content vs coverage scatter with fitted curve.
///
/// `fit_curve` contains (gc, predicted_coverage) points describing the regression fit.
/// For linear regression this is 2 points (a line); for LOESS it is ~50 points (a smooth curve).
pub fn plot_gc_vs_coverage(
    bins: &[CoverageBin],
    ref_config: &ReferenceConfig,
    karyotype: &HashMap<String, usize>,
    fit_curve: &[(f64, f64)],
    out_path: &str,
    title: &str,
) {
    // Find majority CN
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

    let mut maj_gc = Vec::new();
    let mut maj_cov = Vec::new();
    let mut other_gc = Vec::new();
    let mut other_cov = Vec::new();

    for bin in bins {
        if let Some(gc) = bin.gc_content
            && let Some(seg) = get_segment(&bin.chrom, bin.start, bin.end, ref_config) {
                if majority_segments.contains(&seg) {
                    maj_gc.push(gc);
                    maj_cov.push(bin.coverage);
                } else {
                    other_gc.push(gc);
                    other_cov.push(bin.coverage);
                }
            }
    }

    let mut fig = Figure::new(800.0, 600.0);
    let ax = fig.gca();

    // Other-ploidy bins in gray
    if !other_gc.is_empty() {
        ax.scatter(&other_gc, &other_cov)
            .color("gray")
            .size(3.0)
            .alpha(0.5)
            .edge_width(0.0)
            .label("Other CN")
            .build();
    }

    // Majority-ploidy bins in blue
    if !maj_gc.is_empty() {
        ax.scatter(&maj_gc, &maj_cov)
            .color("blue")
            .size(3.0)
            .alpha(0.5)
            .edge_width(0.0)
            .label(format!("CN{} (fit)", majority_cn))
            .build();
    }

    // Fitted curve (line for linear, smooth curve for LOESS)
    if !fit_curve.is_empty() {
        let line_x: Vec<f64> = fit_curve.iter().map(|(x, _)| *x).collect();
        let line_y: Vec<f64> = fit_curve.iter().map(|(_, y)| *y).collect();
        let label = if fit_curve.len() <= 2 { "Linear fit" } else { "LOESS fit" };
        ax.plot(&line_x, &line_y)
            .color("red")
            .linewidth(2.0)
            .label(label)
            .build();
    }

    ax.set_title(title);
    ax.set_xlabel("GC Content");
    ax.set_ylabel("Coverage (reads/Mbp)");
    ax.legend();

    if let Err(e) = fig.save(out_path) {
        error!("Error saving GC vs coverage plot to {}: {}", out_path, e);
    } else {
        info!("Saved GC vs coverage plot: {}", out_path);
    }
}

/// Shared segment x-axis layout used by coverage and BAF panels.
struct SegmentLayout {
    tick_positions: Vec<f64>,
    tick_labels: Vec<String>,
    separator_x: Vec<f64>,
    x_max: f64,
    /// Per-segment: (start_x, number of bins)
    seg_ranges: Vec<(String, f64, usize)>,
}

fn layout_segments(
    bins: &[CoverageBin],
    ref_config: &ReferenceConfig,
) -> (SegmentLayout, HashMap<String, Vec<f64>>) {
    let mut seg_bins: HashMap<String, Vec<f64>> = HashMap::new();
    for b in bins {
        if let Some(seg) = get_segment(&b.chrom, b.start, b.end, ref_config) {
            seg_bins.entry(seg).or_default().push(b.coverage);
        }
    }

    let mut tick_positions = Vec::new();
    let mut tick_labels = Vec::new();
    let mut separator_x = Vec::new();
    let mut seg_ranges = Vec::new();
    let mut x = 0.0;

    for &seg in SEGMENT_ORDER {
        let vals = match seg_bins.get(seg) {
            Some(v) if !v.is_empty() => v,
            _ => continue,
        };
        let seg_start = x;
        let seg_end = x + (vals.len() - 1) as f64;
        tick_positions.push((seg_start + seg_end) / 2.0);
        tick_labels.push(seg.to_string());
        seg_ranges.push((seg.to_string(), seg_start, vals.len()));
        x = seg_end + GAP;
        separator_x.push(x - GAP / 2.0);
    }

    let layout = SegmentLayout {
        tick_positions,
        tick_labels,
        separator_x,
        x_max: x,
        seg_ranges,
    };
    (layout, seg_bins)
}

/// Two-panel plot: top = GC-corrected coverage, bottom = B-allele frequency.
pub fn plot_karyotype_with_baf(
    bins: &[CoverageBin],
    baf_by_segment: &HashMap<String, Vec<f64>>,
    ref_config: &ReferenceConfig,
    out_path: &str,
    title: &str,
) {
    let (layout, seg_bins) = layout_segments(bins, ref_config);

    if layout.seg_ranges.is_empty() {
        warn!("plot_karyotype_with_baf: no data to plot");
        return;
    }

    // Collect coverage scatter data
    let mut cov_x: Vec<f64> = Vec::new();
    let mut cov_y: Vec<f64> = Vec::new();
    let mut median_x_start: Vec<f64> = Vec::new();
    let mut median_x_end: Vec<f64> = Vec::new();
    let mut median_y: Vec<f64> = Vec::new();

    for (seg, seg_start, n) in &layout.seg_ranges {
        let vals = &seg_bins[seg];
        for (i, &cov) in vals.iter().enumerate() {
            cov_x.push(seg_start + i as f64);
            cov_y.push(cov);
        }
        let seg_end = seg_start + (*n as f64 - 1.0);
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = sorted[sorted.len() / 2];
        median_x_start.push(*seg_start);
        median_x_end.push(seg_end);
        median_y.push(med);
    }

    // Collect BAF scatter data using same segment layout
    let mut baf_x: Vec<f64> = Vec::new();
    let mut baf_y: Vec<f64> = Vec::new();

    for (seg, seg_start, _n) in &layout.seg_ranges {
        if let Some(bafs) = baf_by_segment.get(seg) {
            for (i, &b) in bafs.iter().enumerate() {
                // Spread BAF points evenly across the segment width
                let n_sites = bafs.len();
                let seg_width = seg_bins.get(seg).map(|v| v.len() as f64 - 1.0).unwrap_or(1.0).max(1.0);
                let x_pos = if n_sites > 1 {
                    seg_start + (i as f64 / (n_sites - 1) as f64) * seg_width
                } else {
                    seg_start + seg_width / 2.0
                };
                baf_x.push(x_pos);
                baf_y.push(b);
            }
        }
    }

    // Y-limit for coverage: 95th percentile
    let mut sorted_y = cov_y.clone();
    sorted_y.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let y_max = sorted_y[(sorted_y.len() as f64 * 0.99) as usize].max(1.0);

    let mut fig = Figure::new(1600.0, 700.0);
    fig = fig.suptitle(title);

    // --- Top panel: coverage ---
    {
        let ax = fig.add_subplot(2, 1, 1);
        ax.scatter(&cov_x, &cov_y)
            .color("blue")
            .size(3.0)
            .alpha(0.5)
            .edge_width(0.0)
            .build();

        for i in 0..median_y.len() {
            ax.plot(
                [median_x_start[i], median_x_end[i]],
                [median_y[i], median_y[i]],
            )
            .color("red")
            .linewidth(2.5)
            .build();
        }

        for &sx in &layout.separator_x {
            ax.plot([sx, sx], [0.0, y_max])
                .color("black")
                .linewidth(0.5)
                .build();
        }

        ax.set_xlim(-(GAP / 2.0), layout.x_max - GAP / 2.0);
        ax.set_ylim(0.0, y_max);
        ax.set_ylabel("Reads per Mbp");
        ax.x_axis.tick_positions = Some(layout.tick_positions.clone());
        ax.x_axis.tick_labels = Some(layout.tick_labels.clone());
    }

    // --- Bottom panel: BAF ---
    {
        let ax = fig.add_subplot(2, 1, 2);

        if !baf_x.is_empty() {
            let (ds_x, ds_y) = downsample_scatter(
                &baf_x, &baf_y,
                1600.0,
                -(GAP / 2.0),
                layout.x_max - GAP / 2.0,
                10,
            );
            ax.scatter(&ds_x, &ds_y)
                .color("blue")
                .size(3.0)
                .alpha(0.5)
                .edge_width(0.0)
                .build();
        }

        for &sx in &layout.separator_x {
            ax.plot([sx, sx], [0.0, 1.0])
                .color("black")
                .linewidth(0.5)
                .build();
        }

        ax.set_xlim(-(GAP / 2.0), layout.x_max - GAP / 2.0);
        ax.set_ylim(0.0, 1.0);
        ax.set_ylabel("B-Allele Frequency");
        ax.x_axis.tick_positions = Some(layout.tick_positions.clone());
        ax.x_axis.tick_labels = Some(layout.tick_labels.clone());
    }

    if let Err(e) = fig.save(out_path) {
        error!("Error saving karyotype+BAF plot to {}: {}", out_path, e);
    } else {
        info!("Saved karyotype+BAF plot: {}", out_path);
    }
}
