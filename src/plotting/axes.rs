//! Axes (subplot) implementation.

use crate::plotting::data::IntoPlotData;
use crate::plotting::element::{
    AxisConfig, Bounds, GridConfig, Legend, LegendPosition, Text,
};
use crate::plotting::plot::{LinePlot, Plot, ScatterPlot};
use crate::plotting::scale::{LinearScale, Scale};
use crate::plotting::style::{cycle_color, Color, DashPattern, Marker, ThemeConfig};

/// An axes object representing a single plot area.
pub struct Axes {
    /// Position within figure (normalized coordinates)
    pub position: Bounds,
    /// X-axis scale
    pub x_scale: Box<dyn Scale>,
    /// Y-axis scale
    pub y_scale: Box<dyn Scale>,
    /// Plots contained in this axes
    plots: Vec<Box<dyn Plot>>,
    /// Title
    pub title: Option<Text>,
    /// X-axis label
    pub x_label: Option<Text>,
    /// Y-axis label
    pub y_label: Option<Text>,
    /// Legend configuration
    pub legend: Option<Legend>,
    /// Grid configuration
    pub grid: GridConfig,
    /// X-axis configuration
    pub x_axis: AxisConfig,
    /// Y-axis configuration
    pub y_axis: AxisConfig,
    /// Data bounds (computed from plots)
    data_bounds: Option<Bounds>,
    /// Manual x-axis limits
    x_lim: Option<(f64, f64)>,
    /// Manual y-axis limits
    y_lim: Option<(f64, f64)>,
    /// Current color cycle index
    color_index: usize,
    /// Theme configuration
    pub theme: ThemeConfig,
}

impl Axes {
    /// Create a new axes with default settings.
    pub fn new() -> Self {
        Axes {
            position: Bounds::new(0.1, 0.9, 0.1, 0.9),
            x_scale: Box::new(LinearScale::auto()),
            y_scale: Box::new(LinearScale::auto()),
            plots: Vec::new(),
            title: None,
            x_label: None,
            y_label: None,
            legend: None,
            grid: GridConfig::default(),
            x_axis: AxisConfig::default(),
            y_axis: AxisConfig::default(),
            data_bounds: None,
            x_lim: None,
            y_lim: None,
            color_index: 0,
            theme: ThemeConfig::default(),
        }
    }

    /// Set the position within the figure.
    pub fn position(mut self, left: f64, right: f64, bottom: f64, top: f64) -> Self {
        self.position = Bounds::new(left, right, bottom, top);
        self
    }

    /// Add a line plot.
    pub fn plot(
        &mut self,
        x: impl IntoPlotData,
        y: impl IntoPlotData,
    ) -> LinePlotBuilder<'_> {
        LinePlotBuilder {
            axes: self,
            x: x.into_plot_data(),
            y: y.into_plot_data(),
            color: None,
            linewidth: None,
            linestyle: None,
            marker: None,
            markersize: None,
            label: None,
            format: None,
        }
    }

    /// Add a scatter plot.
    pub fn scatter(
        &mut self,
        x: impl IntoPlotData,
        y: impl IntoPlotData,
    ) -> ScatterPlotBuilder<'_> {
        ScatterPlotBuilder {
            axes: self,
            x: x.into_plot_data(),
            y: y.into_plot_data(),
            color: None,
            size: None,
            marker: None,
            alpha: None,
            label: None,
            edge_width: None,
        }
    }

    /// Set the title.
    pub fn set_title(&mut self, title: impl Into<String>) -> &mut Self {
        self.title = Some(Text::new(title, 0.0, 0.0).style(self.theme.title_style.clone()));
        self
    }

    /// Set the x-axis label.
    pub fn set_xlabel(&mut self, label: impl Into<String>) -> &mut Self {
        self.x_label = Some(Text::new(label, 0.0, 0.0).style(self.theme.label_style.clone()));
        self
    }

    /// Set the y-axis label.
    pub fn set_ylabel(&mut self, label: impl Into<String>) -> &mut Self {
        self.y_label = Some(Text::new(label, 0.0, 0.0).style(self.theme.label_style.clone()));
        self
    }

    /// Set the x-axis limits.
    pub fn set_xlim(&mut self, min: f64, max: f64) -> &mut Self {
        self.x_lim = Some((min, max));
        self
    }

    /// Set the y-axis limits.
    pub fn set_ylim(&mut self, min: f64, max: f64) -> &mut Self {
        self.y_lim = Some((min, max));
        self
    }

    /// Enable or disable the grid.
    pub fn grid(&mut self, visible: bool) -> &mut Self {
        self.grid.visible = visible;
        self
    }

    /// Show the legend.
    pub fn legend(&mut self) -> &mut Self {
        if let Some(ref mut legend) = self.legend {
            legend.visible = true;
        } else {
            self.legend = Some(Legend::new());
        }
        self
    }

    /// Show the legend at a specific position.
    pub fn legend_at(&mut self, position: LegendPosition) -> &mut Self {
        if let Some(ref mut legend) = self.legend {
            legend.visible = true;
            legend.position = position;
        } else {
            self.legend = Some(Legend::new().position(position));
        }
        self
    }

    /// Add a plot to this axes.
    pub(crate) fn add_plot(&mut self, plot: Box<dyn Plot>) {
        // Update data bounds
        if let Some(plot_bounds) = plot.bounds() {
            if let Some(ref mut bounds) = self.data_bounds {
                bounds.include_bounds(&plot_bounds);
            } else {
                self.data_bounds = Some(plot_bounds);
            }
        }

        // Add legend entry if plot has a label
        if let Some(entry) = plot.legend_entry() {
            if self.legend.is_none() {
                self.legend = Some(Legend::new().visible(false));
            }
            if let Some(ref mut legend) = self.legend {
                legend.add_entry(entry);
            }
        }

        self.plots.push(plot);
    }

    /// Get the next color from the cycle.
    pub(crate) fn next_color(&mut self) -> Color {
        let color = cycle_color(self.color_index);
        self.color_index += 1;
        color
    }

    /// Get the effective x range.
    fn get_x_range(&self) -> (f64, f64) {
        if let Some((min, max)) = self.x_lim {
            (min, max)
        } else if let Some(ref bounds) = self.data_bounds {
            (bounds.x_min, bounds.x_max)
        } else {
            (0.0, 1.0)
        }
    }

    /// Get the effective y range.
    fn get_y_range(&self) -> (f64, f64) {
        if let Some((min, max)) = self.y_lim {
            (min, max)
        } else if let Some(ref bounds) = self.data_bounds {
            (bounds.y_min, bounds.y_max)
        } else {
            (0.0, 1.0)
        }
    }

    /// Render the axes to SVG.
    pub fn render_svg(&mut self, figure_width: f64, figure_height: f64) -> String {
        let mut svg = String::new();

        // Calculate pixel bounds for the plot area
        let pixel_bounds = Bounds::new(
            self.position.x_min * figure_width,
            self.position.x_max * figure_width,
            (1.0 - self.position.y_max) * figure_height,
            (1.0 - self.position.y_min) * figure_height,
        );

        // Get data bounds with padding
        let (x_min, x_max) = self.get_x_range();
        let (y_min, y_max) = self.get_y_range();

        let mut data_bounds = Bounds::new(x_min, x_max, y_min, y_max);

        // Add padding if range is zero
        if data_bounds.width() == 0.0 {
            data_bounds.x_min -= 0.5;
            data_bounds.x_max += 0.5;
        }
        if data_bounds.height() == 0.0 {
            data_bounds.y_min -= 0.5;
            data_bounds.y_max += 0.5;
        }

        // Add small padding
        let data_bounds = data_bounds.pad(0.05);

        // Update scales
        let _ = self.x_scale.set_range(data_bounds.x_min, data_bounds.x_max);
        let _ = self.y_scale.set_range(data_bounds.y_min, data_bounds.y_max);

        // Render background
        svg.push_str(&format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>\n",
            pixel_bounds.x_min,
            pixel_bounds.y_min,
            pixel_bounds.width(),
            pixel_bounds.height(),
            self.theme.plot_background.to_svg_string()
        ));

        // Render grid
        if self.grid.visible {
            svg.push_str(&self.render_grid(&data_bounds, &pixel_bounds));
        }

        // Render plots (clipped to plot area)
        let clip_id = format!(
            "plot-clip-{:.0}-{:.0}",
            pixel_bounds.x_min, pixel_bounds.y_min
        );
        svg.push_str(&format!(
            "<defs><clipPath id=\"{}\"><rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\"/></clipPath></defs>\n",
            clip_id, pixel_bounds.x_min, pixel_bounds.y_min, pixel_bounds.width(), pixel_bounds.height()
        ));
        svg.push_str(&format!("<g clip-path=\"url(#{})\">", clip_id));
        for plot in &self.plots {
            svg.push_str(&plot.render_svg(&data_bounds, &pixel_bounds));
        }
        svg.push_str("</g>\n");

        // Render axes
        svg.push_str(&self.render_axes(&data_bounds, &pixel_bounds));

        // Render title
        if let Some(ref title) = self.title {
            let x = (pixel_bounds.x_min + pixel_bounds.x_max) / 2.0;
            let y = pixel_bounds.y_min - 10.0;
            let mut t = title.clone();
            t.x = x;
            t.y = y;
            t.style.anchor = crate::plotting::style::TextAnchor::Middle;
            svg.push_str(&t.to_svg());
            svg.push('\n');
        }

        // Render x label
        if let Some(ref label) = self.x_label {
            let x = (pixel_bounds.x_min + pixel_bounds.x_max) / 2.0;
            let y = pixel_bounds.y_max + 40.0;
            let mut t = label.clone();
            t.x = x;
            t.y = y;
            t.style.anchor = crate::plotting::style::TextAnchor::Middle;
            svg.push_str(&t.to_svg());
            svg.push('\n');
        }

        // Render y label
        if let Some(ref label) = self.y_label {
            let x = pixel_bounds.x_min - 45.0;
            let y = (pixel_bounds.y_min + pixel_bounds.y_max) / 2.0;
            let mut t = label.clone();
            t.x = x;
            t.y = y;
            t.style.anchor = crate::plotting::style::TextAnchor::Middle;
            t.style.rotation = -90.0;
            svg.push_str(&t.to_svg());
            svg.push('\n');
        }

        // Render legend
        if let Some(ref legend) = self.legend {
            if legend.visible && !legend.entries.is_empty() {
                svg.push_str(&self.render_legend(legend, &pixel_bounds));
            }
        }

        svg
    }

    fn render_grid(&self, _data_bounds: &Bounds, pixel_bounds: &Bounds) -> String {
        let mut svg = String::new();
        let style = &self.grid.major_style;

        // X grid lines
        if self.grid.show_x {
            let ticks = self.x_scale.ticks(self.x_axis.num_ticks);
            for tick in ticks {
                let x_norm = self.x_scale.transform(tick);
                let px = pixel_bounds.x_min + x_norm * pixel_bounds.width();
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{}\" stroke-opacity=\"{}\"/>\n",
                    px, pixel_bounds.y_min, px, pixel_bounds.y_max,
                    style.color.to_svg_string(), style.width, self.grid.grid_opacity
                ));
            }
        }

        // Y grid lines
        if self.grid.show_y {
            let ticks = self.y_scale.ticks(self.y_axis.num_ticks);
            for tick in ticks {
                let y_norm = self.y_scale.transform(tick);
                let py = pixel_bounds.y_max - y_norm * pixel_bounds.height();
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{}\" stroke-opacity=\"{}\"/>\n",
                    pixel_bounds.x_min, py, pixel_bounds.x_max, py,
                    style.color.to_svg_string(), style.width, self.grid.grid_opacity
                ));
            }
        }

        svg
    }

    fn render_axes(&self, _data_bounds: &Bounds, pixel_bounds: &Bounds) -> String {
        let mut svg = String::new();

        // X axis
        if self.x_axis.visible {
            // Axis line
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{}\"/>\n",
                pixel_bounds.x_min, pixel_bounds.y_max,
                pixel_bounds.x_max, pixel_bounds.y_max,
                self.x_axis.line_style.color.to_svg_string(),
                self.x_axis.line_style.width
            ));

            // Ticks and labels
            if self.x_axis.show_ticks {
                // Use custom tick positions/labels if provided, otherwise auto-generate
                let (tick_values, tick_labels): (Vec<f64>, Vec<String>) =
                    if let Some(ref positions) = self.x_axis.tick_positions {
                        let labels = if let Some(ref custom_labels) = self.x_axis.tick_labels {
                            custom_labels.clone()
                        } else {
                            positions.iter().map(|t| self.x_axis.format_tick(*t)).collect()
                        };
                        (positions.clone(), labels)
                    } else {
                        let ticks = self.x_scale.ticks(self.x_axis.num_ticks);
                        let labels = ticks.iter().map(|t| self.x_axis.format_tick(*t)).collect();
                        (ticks, labels)
                    };

                for (tick, label) in tick_values.iter().zip(tick_labels.iter()) {
                    let x_norm = self.x_scale.transform(*tick);
                    let px = pixel_bounds.x_min + x_norm * pixel_bounds.width();

                    // Tick mark
                    svg.push_str(&format!(
                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>\n",
                        px, pixel_bounds.y_max,
                        px, pixel_bounds.y_max + self.x_axis.tick_length,
                        self.x_axis.line_style.color.to_svg_string()
                    ));

                    // Label
                    svg.push_str(&format!(
                        "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\" {} >{}</text>\n",
                        px,
                        pixel_bounds.y_max + self.x_axis.tick_length + self.x_axis.tick_padding + 10.0,
                        self.x_axis.tick_label_style.to_svg_attrs(),
                        label
                    ));
                }
            }
        }

        // Y axis
        if self.y_axis.visible {
            // Axis line
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{}\"/>\n",
                pixel_bounds.x_min, pixel_bounds.y_min,
                pixel_bounds.x_min, pixel_bounds.y_max,
                self.y_axis.line_style.color.to_svg_string(),
                self.y_axis.line_style.width
            ));

            // Ticks and labels
            if self.y_axis.show_ticks {
                let (tick_values, tick_labels): (Vec<f64>, Vec<String>) =
                    if let Some(ref positions) = self.y_axis.tick_positions {
                        let labels = if let Some(ref custom_labels) = self.y_axis.tick_labels {
                            custom_labels.clone()
                        } else {
                            positions.iter().map(|t| self.y_axis.format_tick(*t)).collect()
                        };
                        (positions.clone(), labels)
                    } else {
                        let ticks = self.y_scale.ticks(self.y_axis.num_ticks);
                        let labels = ticks.iter().map(|t| self.y_axis.format_tick(*t)).collect();
                        (ticks, labels)
                    };

                for (tick, label) in tick_values.iter().zip(tick_labels.iter()) {
                    let y_norm = self.y_scale.transform(*tick);
                    let py = pixel_bounds.y_max - y_norm * pixel_bounds.height();

                    // Tick mark
                    svg.push_str(&format!(
                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>\n",
                        pixel_bounds.x_min - self.y_axis.tick_length, py,
                        pixel_bounds.x_min, py,
                        self.y_axis.line_style.color.to_svg_string()
                    ));

                    // Label
                    svg.push_str(&format!(
                        "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\" dominant-baseline=\"middle\" {}>{}</text>\n",
                        pixel_bounds.x_min - self.y_axis.tick_length - self.y_axis.tick_padding,
                        py,
                        self.y_axis.tick_label_style.to_svg_attrs(),
                        label
                    ));
                }
            }
        }

        svg
    }

    fn render_legend(&self, legend: &Legend, pixel_bounds: &Bounds) -> String {
        let mut svg = String::new();

        let (anchor_x, anchor_y) = legend.position.anchor();

        // Calculate legend dimensions
        let line_height = legend.text_style.font_size * 1.5;
        let legend_height = legend.entries.len() as f64 * line_height + legend.padding * 2.0;
        let legend_width = 100.0 + legend.padding * 2.0; // Approximate width

        // Calculate legend position
        let (lx, ly) = match legend.position {
            LegendPosition::TopRight | LegendPosition::Right | LegendPosition::BottomRight => {
                (
                    pixel_bounds.x_max - legend_width - 5.0,
                    pixel_bounds.y_min + (pixel_bounds.height() * (1.0 - anchor_y)),
                )
            }
            LegendPosition::TopLeft | LegendPosition::Left | LegendPosition::BottomLeft => {
                (
                    pixel_bounds.x_min + 5.0,
                    pixel_bounds.y_min + (pixel_bounds.height() * (1.0 - anchor_y)),
                )
            }
            _ => {
                (
                    pixel_bounds.x_min + anchor_x * pixel_bounds.width() - legend_width / 2.0,
                    pixel_bounds.y_min + (1.0 - anchor_y) * pixel_bounds.height(),
                )
            }
        };

        // Background
        svg.push_str(&format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" {}/>\n",
            lx, ly, legend_width, legend_height,
            legend.background.to_svg_style()
        ));

        // Entries
        for (i, entry) in legend.entries.iter().enumerate() {
            let ey = ly + legend.padding + (i as f64 + 0.5) * line_height;
            let line_x1 = lx + legend.padding;
            let line_x2 = lx + legend.padding + legend.line_length;
            let line_mid = (line_x1 + line_x2) / 2.0;

            // Determine what to draw based on available styles
            let has_line = entry.line_style.is_some();
            let has_marker = entry.marker_style.is_some();
            let has_fill = entry.fill_style.is_some();

            // Draw line if available (for line plots)
            if let Some(ref line_style) = entry.line_style {
                let dash_attr = line_style.dash.to_svg_dasharray()
                    .map(|d| format!(" stroke-dasharray=\"{}\"", d))
                    .unwrap_or_default();
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{}\"{}/>\n",
                    line_x1, ey, line_x2, ey,
                    line_style.color.to_svg_string(),
                    line_style.width,
                    dash_attr
                ));
            }

            // Draw fill rectangle if available (for bar/histogram) and no line
            if let Some(ref fill_style) = entry.fill_style {
                if !has_line && !has_marker {
                    let rect_size = 10.0;
                    svg.push_str(&format!(
                        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" {}/>\n",
                        line_mid - rect_size / 2.0,
                        ey - rect_size / 2.0,
                        rect_size,
                        rect_size,
                        fill_style.to_svg_style()
                    ));
                }
            }

            // Draw marker if available (for scatter plots or line plots with markers)
            if let Some(ref marker_style) = entry.marker_style {
                let marker_svg = marker_style.render_at(line_mid, ey);
                svg.push_str(&marker_svg);
            }

            // Fallback: if nothing was drawn, draw a simple colored line
            if !has_line && !has_marker && !has_fill {
                let fallback_color = cycle_color(i);
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"2\"/>\n",
                    line_x1, ey, line_x2, ey,
                    fallback_color.to_svg_string()
                ));
            }

            // Label
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" dominant-baseline=\"middle\" {}>{}</text>\n",
                lx + legend.padding + legend.line_length + legend.label_gap,
                ey,
                legend.text_style.to_svg_attrs(),
                crate::plotting::element::escape_xml(&entry.label)
            ));
        }

        svg
    }
}

impl Default for Axes {
    fn default() -> Self {
        Self::new()
    }
}

// Builder types for fluent API

/// Builder for line plots.
pub struct LinePlotBuilder<'a> {
    axes: &'a mut Axes,
    x: Vec<f64>,
    y: Vec<f64>,
    color: Option<Color>,
    linewidth: Option<f64>,
    linestyle: Option<DashPattern>,
    marker: Option<Marker>,
    markersize: Option<f64>,
    label: Option<String>,
    format: Option<String>,
}

impl<'a> LinePlotBuilder<'a> {
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.color = Some(color.into());
        self
    }

    pub fn linewidth(mut self, width: f64) -> Self {
        self.linewidth = Some(width);
        self
    }

    pub fn linestyle(mut self, style: DashPattern) -> Self {
        self.linestyle = Some(style);
        self
    }

    pub fn marker(mut self, marker: Marker) -> Self {
        self.marker = Some(marker);
        self
    }

    pub fn markersize(mut self, size: f64) -> Self {
        self.markersize = Some(size);
        self
    }

    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    pub fn format(mut self, fmt: impl Into<String>) -> Self {
        self.format = Some(fmt.into());
        self
    }

    pub fn build(self) -> &'a mut Axes {
        let color = self.color.unwrap_or_else(|| self.axes.next_color());

        let mut plot = LinePlot::new(self.x, self.y).color(color);

        if let Some(fmt) = self.format {
            plot = plot.format(&fmt);
        }
        if let Some(width) = self.linewidth {
            plot = plot.linewidth(width);
        }
        if let Some(style) = self.linestyle {
            plot = plot.linestyle(style);
        }
        if let Some(marker) = self.marker {
            plot = plot.marker(marker);
        }
        if let Some(size) = self.markersize {
            plot = plot.markersize(size);
        }
        if let Some(label) = self.label {
            plot = plot.label(label);
        }

        self.axes.add_plot(Box::new(plot));
        self.axes
    }
}

/// Builder for scatter plots.
pub struct ScatterPlotBuilder<'a> {
    axes: &'a mut Axes,
    x: Vec<f64>,
    y: Vec<f64>,
    color: Option<Color>,
    size: Option<f64>,
    marker: Option<Marker>,
    alpha: Option<f64>,
    label: Option<String>,
    edge_width: Option<f64>,
}

impl<'a> ScatterPlotBuilder<'a> {
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.color = Some(color.into());
        self
    }

    pub fn size(mut self, size: f64) -> Self {
        self.size = Some(size);
        self
    }

    pub fn marker(mut self, marker: Marker) -> Self {
        self.marker = Some(marker);
        self
    }

    pub fn alpha(mut self, alpha: f64) -> Self {
        self.alpha = Some(alpha);
        self
    }

    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    pub fn edge_width(mut self, width: f64) -> Self {
        self.edge_width = Some(width);
        self
    }

    pub fn build(self) -> &'a mut Axes {
        let color = self.color.unwrap_or_else(|| self.axes.next_color());

        let mut plot = ScatterPlot::new(self.x, self.y).color(color);

        if let Some(size) = self.size {
            plot = plot.size(size);
        }
        if let Some(marker) = self.marker {
            plot = plot.marker(marker);
        }
        if let Some(alpha) = self.alpha {
            plot = plot.alpha(alpha);
        }
        if let Some(label) = self.label {
            plot = plot.label(label);
        }
        if let Some(ew) = self.edge_width {
            plot = plot.edge_width(ew);
        }

        self.axes.add_plot(Box::new(plot));
        self.axes
    }
}


