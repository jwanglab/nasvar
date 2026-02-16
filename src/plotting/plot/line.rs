//! Line plot implementation.

use crate::plotting::element::Bounds;
use crate::plotting::plot::Plot;
use crate::plotting::style::{Color, DashPattern, LineStyle, Marker, MarkerStyle};

/// A line plot connecting data points.
#[derive(Debug, Clone)]
pub struct LinePlot {
    /// X coordinates
    pub x: Vec<f64>,
    /// Y coordinates
    pub y: Vec<f64>,
    /// Line style
    pub line_style: LineStyle,
    /// Marker style (optional)
    pub marker_style: Option<MarkerStyle>,
    /// Label for legend
    pub label: Option<String>,
}

impl LinePlot {
    /// Create a new line plot from x and y data.
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        LinePlot {
            x,
            y,
            line_style: LineStyle::default(),
            marker_style: None,
            label: None,
        }
    }

    /// Set the line color.
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.line_style.color = color.into();
        if let Some(ref mut marker) = self.marker_style {
            marker.fill = self.line_style.color.clone();
            marker.edge_color = self.line_style.color.clone();
        }
        self
    }

    /// Set the line width.
    pub fn linewidth(mut self, width: f64) -> Self {
        self.line_style.width = width;
        self
    }

    /// Set the line style (dash pattern).
    pub fn linestyle(mut self, dash: DashPattern) -> Self {
        self.line_style.dash = dash;
        self
    }

    /// Add markers to the line plot.
    pub fn marker(mut self, marker: Marker) -> Self {
        let color = self.line_style.color.clone();
        self.marker_style = Some(MarkerStyle {
            marker,
            fill: color.clone(),
            edge_color: color,
            ..Default::default()
        });
        self
    }

    /// Set marker size.
    pub fn markersize(mut self, size: f64) -> Self {
        if let Some(ref mut marker) = self.marker_style {
            marker.size = size;
        }
        self
    }

    /// Set the label for the legend.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Set the complete line style.
    pub fn line_style(mut self, style: LineStyle) -> Self {
        self.line_style = style;
        self
    }

    /// Set the complete marker style.
    pub fn marker_style(mut self, style: MarkerStyle) -> Self {
        self.marker_style = Some(style);
        self
    }

    /// Parse matplotlib-style format string (e.g., "r-", "b--o", "g:s")
    pub fn format(mut self, fmt: &str) -> Self {
        let mut chars = fmt.chars().peekable();

        // Parse color (single letter)
        if let Some(&c) = chars.peek() {
            if let Some(color) = parse_color_char(c) {
                self.line_style.color = color.clone();
                if let Some(ref mut marker) = self.marker_style {
                    marker.fill = color.clone();
                    marker.edge_color = color;
                }
                chars.next();
            }
        }

        // Parse line style
        let remaining: String = chars.clone().collect();
        if let Some(dash) = DashPattern::from_format_str(&remaining) {
            self.line_style.dash = dash;
            // Skip the line style characters
            if remaining.starts_with("--") || remaining.starts_with("-.") {
                chars.next();
                chars.next();
            } else if remaining.starts_with('-') || remaining.starts_with(':') {
                chars.next();
            }
        }

        // Parse marker
        for c in chars {
            if let Some(marker) = Marker::from_format_char(c) {
                let color = self.line_style.color.clone();
                self.marker_style = Some(MarkerStyle {
                    marker,
                    fill: color.clone(),
                    edge_color: color,
                    ..Default::default()
                });
                break;
            }
        }

        self
    }
}

fn parse_color_char(c: char) -> Option<Color> {
    match c {
        'b' => Some(Color::BLUE),
        'g' => Some(Color::GREEN),
        'r' => Some(Color::RED),
        'c' => Some(Color::CYAN),
        'm' => Some(Color::MAGENTA),
        'y' => Some(Color::YELLOW),
        'k' => Some(Color::BLACK),
        'w' => Some(Color::WHITE),
        _ => None,
    }
}

impl Plot for LinePlot {
    fn bounds(&self) -> Option<Bounds> {
        if self.x.is_empty() || self.y.is_empty() {
            return None;
        }

        let mut bounds = Bounds::new(
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );

        for (&x, &y) in self.x.iter().zip(self.y.iter()) {
            if x.is_finite() && y.is_finite() {
                bounds.include_point(x, y);
            }
        }

        if bounds.x_min.is_finite() && bounds.x_max.is_finite() {
            Some(bounds)
        } else {
            None
        }
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn line_style(&self) -> Option<LineStyle> {
        Some(self.line_style.clone())
    }

    fn marker_style(&self) -> Option<MarkerStyle> {
        self.marker_style.clone()
    }

    fn render_svg(&self, data_bounds: &Bounds, pixel_bounds: &Bounds) -> String {
        let mut svg = String::new();

        if self.x.len() < 2 {
            return svg;
        }

        // Transform points
        let points: Vec<(f64, f64)> = self
            .x
            .iter()
            .zip(self.y.iter())
            .filter(|(x, y)| x.is_finite() && y.is_finite())
            .map(|(&x, &y)| transform_point(x, y, data_bounds, pixel_bounds))
            .collect();

        if points.is_empty() {
            return svg;
        }

        // Render line
        let path_data: String = points
            .iter()
            .enumerate()
            .map(|(i, (x, y))| {
                if i == 0 {
                    format!("M{:.2},{:.2}", x, y)
                } else {
                    format!(" L{:.2},{:.2}", x, y)
                }
            })
            .collect();

        svg.push_str(&format!(
            "<path d=\"{}\" {}/>\n",
            path_data,
            self.line_style.to_svg_style()
        ));

        // Render markers
        if let Some(ref marker_style) = self.marker_style {
            let marker = &marker_style.marker;
            let size = marker_style.size / 2.0; // radius

            for (x, y) in &points {
                if marker.is_circle() {
                    svg.push_str(&format!(
                        "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" {}/>\n",
                        x, y, size, marker_style.to_svg_style()
                    ));
                } else if let Some(path) = marker.to_svg_path(size) {
                    svg.push_str(&format!(
                        "<path d=\"{}\" transform=\"translate({:.2},{:.2})\" {}/>\n",
                        path, x, y, marker_style.to_svg_style()
                    ));
                }
            }
        }

        svg
    }
}

fn transform_point(x: f64, y: f64, data: &Bounds, pixel: &Bounds) -> (f64, f64) {
    let x_norm = (x - data.x_min) / data.width();
    let y_norm = (y - data.y_min) / data.height();

    let px = pixel.x_min + x_norm * pixel.width();
    let py = pixel.y_max - y_norm * pixel.height(); // Flip Y

    (px, py)
}

