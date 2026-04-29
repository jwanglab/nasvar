use anyhow::Result;
use log::{info, warn};

use crate::config::FocalDepthPlot;
use crate::input::AlignmentInput;
use crate::plotting::prelude::*;
use crate::plotting::style::Color;

/// Parse "chrom:start-end" into (chrom, start_bp, end_bp) (1-based, inclusive).
fn parse_region(region: &str) -> Result<(String, u64, u64)> {
    let colon = region.find(':')
        .ok_or_else(|| anyhow::anyhow!("Invalid region '{}': missing ':'", region))?;
    let chrom = region[..colon].to_string();
    let rest = &region[colon + 1..];
    let dash = rest.find('-')
        .ok_or_else(|| anyhow::anyhow!("Invalid region '{}': missing '-'", region))?;
    let start: u64 = rest[..dash].parse()?;
    let end: u64 = rest[dash + 1..].parse()?;
    Ok((chrom, start, end))
}

/// Parse "start-end/label" into (start_bp, end_bp, label).
fn parse_annotation(anno: &str) -> Result<(u64, u64, String)> {
    let slash = anno.find('/')
        .ok_or_else(|| anyhow::anyhow!("Invalid annotation '{}': missing '/'", anno))?;
    let label = anno[slash + 1..].to_string();
    let coords = &anno[..slash];
    let dash = coords.find('-')
        .ok_or_else(|| anyhow::anyhow!("Invalid annotation '{}': missing '-'", anno))?;
    let start: u64 = coords[..dash].parse()?;
    let end: u64 = coords[dash + 1..].parse()?;
    Ok((start, end, label))
}

fn fmt_mbp(pos: f64) -> String {
    let mbp = pos / 1_000_000.0;
    if (mbp.fract() * 10.0).round() == 0.0 {
        format!("{:.0} Mbp", mbp)
    } else {
        format!("{:.1} Mbp", mbp)
    }
}

fn format_binsize(binsize: u64) -> String {
    if binsize >= 1_000_000 {
        let v = binsize / 1_000_000;
        format!("{} Mbp", v)
    } else if binsize >= 1_000 {
        let v = binsize / 1_000;
        format!("{} Kbp", v)
    } else {
        format!("{} bp", binsize)
    }
}

/// Map data coordinates to figure pixel coordinates.
/// `xlim` / `ylim` are the explicit limits passed to set_xlim / set_ylim.
/// The axes render_svg pads the data bounds by 5%, so we replicate that here.
fn data_to_pixel(
    data_x: f64,
    data_y: f64,
    fig_w: f64,
    fig_h: f64,
    xlim: (f64, f64),
    ylim: (f64, f64),
    axes_pos: (f64, f64, f64, f64), // (left, right, bottom, top) normalized
) -> (f64, f64) {
    let (ax_left, ax_right, ax_bottom, ax_top) = axes_pos;
    let px_left = ax_left * fig_w;
    let px_right = ax_right * fig_w;
    let py_top = (1.0 - ax_top) * fig_h;
    let py_bottom = (1.0 - ax_bottom) * fig_h;

    // Replicate the 5% padding applied by axes.render_svg
    let x_pad = (xlim.1 - xlim.0) * 0.05;
    let y_pad = (ylim.1 - ylim.0) * 0.05;
    let dx_min = xlim.0 - x_pad;
    let dx_max = xlim.1 + x_pad;
    let dy_min = ylim.0 - y_pad;
    let dy_max = ylim.1 + y_pad;

    let x_norm = (data_x - dx_min) / (dx_max - dx_min);
    let y_norm = (data_y - dy_min) / (dy_max - dy_min);

    let px = px_left + x_norm * (px_right - px_left);
    let py = py_bottom - y_norm * (py_bottom - py_top);

    (px, py)
}

/// Plot a focal coverage depth track for one FocalDepthPlot configuration entry.
/// Outputs an SVG at `{out_prefix}/{plot_cfg.output}` (`.png` extension is replaced with `.svg`).
pub fn plot_focal_depth(
    bam: &mut AlignmentInput,
    plot_cfg: &FocalDepthPlot,
    out_prefix: &str,
) -> Result<()> {
    let (chrom, region_start, region_end) = parse_region(&plot_cfg.region)?;

    let bam_chrom = bam.contig_mapper.to_bam_name(&chrom);
    let bam_region = format!("{}:{}-{}", bam_chrom, region_start, region_end);

    let binsize = if plot_cfg.binsize > 0 { plot_cfg.binsize as u64 } else { 1 };
    let n_bins = ((region_end - region_start) / binsize + 1) as usize;
    let mut counts = vec![0u32; n_bins];

    info!("Counting reads in {} with binsize {} bp ({} bins)", bam_region, binsize, n_bins);

    let query = bam.query(&bam_region)?;
    for result in query {
        let rec = result?;
        if rec.flag & 0x900 != 0 { continue; } // skip secondary/supplementary
        if rec.pos < 0 { continue; }
        let pos = rec.pos as u64;
        if pos < region_start || pos >= region_end { continue; }
        let bin = ((pos - region_start) / binsize) as usize;
        if bin < counts.len() {
            counts[bin] += 1;
        }
    }

    // Bin center positions in bp
    let x_vals: Vec<f64> = (0..n_bins)
        .map(|i| (region_start + i as u64 * binsize + binsize / 2) as f64)
        .collect();
    let y_vals: Vec<f64> = counts.iter().map(|&c| c as f64).collect();

    let chrom_display = chrom.strip_prefix("chr").unwrap_or(&chrom).to_string();
    let y_data_max = y_vals.iter().cloned().fold(0.0_f64, f64::max);
    let y_max = y_data_max.max(1.0);
    let ylim = (0.0_f64, y_max * 1.05);
    let xlim = (region_start as f64, region_end as f64);

    // Parse annotations
    let annotations: Vec<(u64, u64, String)> = plot_cfg.annotations.iter().filter_map(|a| {
        match parse_annotation(a) {
            Ok(t) => Some(t),
            Err(e) => { warn!("Skipping annotation '{}': {}", a, e); None }
        }
    }).collect();

    let fig_w = 1000.0_f64;
    let fig_h = 600.0_f64;

    let mut fig = Figure::new(fig_w, fig_h);
    let ax = fig.gca();

    // Line (semi-transparent)
    ax.plot(&x_vals, &y_vals)
        .color(Color::rgba(70, 130, 180, 0.4)) // steelblue 40%
        .linewidth(1.5)
        .build();

    // Scatter dots
    ax.scatter(&x_vals, &y_vals)
        .color(Color::rgb(70, 130, 180))
        .size(4.0)
        .alpha(1.0)
        .edge_width(0.0)
        .build();

    // Annotation lines in data space
    let anno_y_line = y_max * 0.93;
    let anno_y_tick = y_max * 0.96;
    for (as_, ae, _) in &annotations {
        ax.plot([*as_ as f64, *ae as f64], [anno_y_line, anno_y_line])
            .color("red")
            .linewidth(2.0)
            .build();
        let center = (*as_ + *ae) as f64 / 2.0;
        ax.plot([center, center], [anno_y_line, anno_y_tick])
            .color("red")
            .linewidth(2.0)
            .build();
    }

    ax.set_xlim(xlim.0, xlim.1);
    ax.set_ylim(ylim.0, ylim.1);

    // X-axis ticks at round Mbp intervals
    let range_bp = (region_end - region_start) as f64;
    let approx_step_mbp = range_bp / 1_000_000.0 / 6.0;
    let tick_step_mbp = if approx_step_mbp <= 1.0 { 1.0 }
        else if approx_step_mbp <= 2.0 { 2.0 }
        else if approx_step_mbp <= 5.0 { 5.0 }
        else { (approx_step_mbp / 5.0).ceil() * 5.0 };
    let tick_step_bp = (tick_step_mbp * 1_000_000.0) as u64;
    let first_tick = (region_start / tick_step_bp + 1) * tick_step_bp;
    let mut tick_positions: Vec<f64> = Vec::new();
    let mut t = first_tick;
    while t < region_end {
        tick_positions.push(t as f64);
        t += tick_step_bp;
    }
    let tick_labels: Vec<String> = tick_positions.iter().map(|&p| fmt_mbp(p)).collect();
    ax.x_axis.tick_positions = Some(tick_positions);
    ax.x_axis.tick_labels = Some(tick_labels);
    ax.x_axis.num_ticks = 6;

    ax.set_xlabel(format!("Chromosome {}", chrom_display));
    ax.set_ylabel(if plot_cfg.binsize > 0 {
        format!("Reads / {}", format_binsize(binsize))
    } else {
        "Sequencing depth (bp)".to_string()
    });
    ax.grid(true);

    let mut svg = fig.render();

    // Inject annotation text labels by computing pixel coordinates
    if !annotations.is_empty() {
        // axes position from add_subplot(1,1,1) with margin=0.06
        let axes_pos = (0.10_f64, 0.93, 0.11, 0.92);
        let anno_text_y = y_max * 0.98;

        let mut text_svg = String::new();
        for (as_, ae, label) in &annotations {
            let center = (*as_ + *ae) as f64 / 2.0;
            let (px, py) = data_to_pixel(
                center, anno_text_y,
                fig_w, fig_h,
                xlim, ylim,
                axes_pos,
            );
            text_svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" \
                 font-size=\"13\" fill=\"red\" font-style=\"italic\">{}</text>\n",
                px, py,
                crate::plotting::element::escape_xml(label)
            ));
        }

        if let Some(pos) = svg.rfind("</svg>") {
            svg.insert_str(pos, &text_svg);
        }
    }

    // Output: {out_prefix}.{filename}.svg  (flat prefix convention, .png -> .svg)
    let stem = if plot_cfg.output.ends_with(".png") {
        &plot_cfg.output[..plot_cfg.output.len() - 4]
    } else if plot_cfg.output.ends_with(".svg") {
        &plot_cfg.output[..plot_cfg.output.len() - 4]
    } else {
        &plot_cfg.output
    };
    let out_path = format!("{}.{}.svg", out_prefix, stem);

    std::fs::write(&out_path, svg.as_bytes())?;
    info!("Saved focal depth plot: {}", out_path);

    Ok(())
}
