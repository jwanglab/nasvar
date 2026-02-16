//! Scatter plot implementation.

use std::fmt::Write as _;

use crate::plotting::element::Bounds;
use crate::plotting::plot::Plot;
use crate::plotting::style::{Color, Marker, MarkerStyle};

/// A scatter plot showing individual data points.
#[derive(Debug, Clone)]
pub struct ScatterPlot {
    /// X coordinates
    pub x: Vec<f64>,
    /// Y coordinates
    pub y: Vec<f64>,
    /// Marker style
    pub marker_style: MarkerStyle,
    /// Label for legend
    pub label: Option<String>,
    /// Optional sizes for each point (for bubble charts)
    pub sizes: Option<Vec<f64>>,
    /// Optional colors for each point
    pub colors: Option<Vec<Color>>,
    /// Alpha/opacity value
    pub alpha: f64,
}

impl ScatterPlot {
    /// Create a new scatter plot from x and y data.
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        ScatterPlot {
            x,
            y,
            marker_style: MarkerStyle::default(),
            label: None,
            sizes: None,
            colors: None,
            alpha: 1.0,
        }
    }

    /// Set the marker color.
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        let c = color.into();
        self.marker_style.fill = c.clone();
        self.marker_style.edge_color = c;
        self
    }

    /// Set the marker type.
    pub fn marker(mut self, marker: Marker) -> Self {
        self.marker_style.marker = marker;
        self
    }

    /// Set the marker size.
    pub fn size(mut self, size: f64) -> Self {
        self.marker_style.size = size;
        self
    }

    /// Set individual sizes for each point (bubble chart).
    pub fn sizes(mut self, sizes: Vec<f64>) -> Self {
        self.sizes = Some(sizes);
        self
    }

    /// Set individual colors for each point.
    pub fn colors(mut self, colors: Vec<Color>) -> Self {
        self.colors = Some(colors);
        self
    }

    /// Set the alpha/opacity value.
    pub fn alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha.clamp(0.0, 1.0);
        self.marker_style.fill_opacity = self.alpha;
        self
    }

    /// Set the edge color.
    pub fn edge_color(mut self, color: impl Into<Color>) -> Self {
        self.marker_style.edge_color = color.into();
        self
    }

    /// Set the edge width.
    pub fn edge_width(mut self, width: f64) -> Self {
        self.marker_style.edge_width = width;
        self
    }

    /// Set the label for the legend.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Set the complete marker style.
    pub fn marker_style(mut self, style: MarkerStyle) -> Self {
        self.marker_style = style;
        self
    }
}

impl Plot for ScatterPlot {
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

    fn marker_style(&self) -> Option<MarkerStyle> {
        Some(self.marker_style.clone())
    }

    fn render_svg(&self, data_bounds: &Bounds, pixel_bounds: &Bounds) -> String {
        let marker = &self.marker_style.marker;
        let default_size = self.marker_style.size;
        let uniform_style = self.colors.is_none() && self.sizes.is_none();

        // Tier 1: Compound path for large uniform-style circle datasets (1 DOM element)
        if uniform_style && marker.is_circle() && self.x.len() > 500 {
            let radius = default_size / 2.0;
            let style = self.marker_style.to_svg_style();
            let mut path_d = String::with_capacity(self.x.len() * 25);
            for (&x, &y) in self.x.iter().zip(self.y.iter()) {
                if !x.is_finite() || !y.is_finite() {
                    continue;
                }
                let (px, py) = transform_point(x, y, data_bounds, pixel_bounds);
                // Small rect subpath (visually identical to circle at r≤2px)
                let d = radius * 2.0;
                let _ = write!(path_d, "M{:.0},{:.0}h{:.0}v{:.0}h{:.0}z", px - radius, py - radius, d, d, -d);
            }
            return format!("<path d=\"{}\" {}/>\n", path_d, style);
        }

        // Tier 2: Group wrapper for uniform-style circles (shared attributes)
        if uniform_style && marker.is_circle() {
            let style = self.marker_style.to_svg_style();
            let radius = default_size / 2.0;
            let mut svg = format!("<g {}>\n", style);
            for (&x, &y) in self.x.iter().zip(self.y.iter()) {
                if !x.is_finite() || !y.is_finite() {
                    continue;
                }
                let (px, py) = transform_point(x, y, data_bounds, pixel_bounds);
                svg.push_str(&format!(
                    "<circle cx=\"{:.0}\" cy=\"{:.0}\" r=\"{:.0}\"/>\n",
                    px, py, radius
                ));
            }
            svg.push_str("</g>\n");
            return svg;
        }

        // Tier 3: Per-element rendering (variable colors/sizes or non-circle markers)
        let mut svg = String::new();
        for (i, (&x, &y)) in self.x.iter().zip(self.y.iter()).enumerate() {
            if !x.is_finite() || !y.is_finite() {
                continue;
            }

            let (px, py) = transform_point(x, y, data_bounds, pixel_bounds);

            let size = if let Some(ref sizes) = self.sizes {
                sizes.get(i).copied().unwrap_or(default_size)
            } else {
                default_size
            } / 2.0;

            let style = if let Some(ref colors) = self.colors {
                if let Some(color) = colors.get(i) {
                    let mut s = self.marker_style.clone();
                    s.fill = color.clone();
                    s.edge_color = color.clone();
                    s.to_svg_style()
                } else {
                    self.marker_style.to_svg_style()
                }
            } else {
                self.marker_style.to_svg_style()
            };

            if marker.is_circle() {
                svg.push_str(&format!(
                    "<circle cx=\"{:.0}\" cy=\"{:.0}\" r=\"{:.0}\" {}/>\n",
                    px, py, size, style
                ));
            } else if let Some(path) = marker.to_svg_path(size) {
                svg.push_str(&format!(
                    "<path d=\"{}\" transform=\"translate({:.0},{:.0})\" {}/>\n",
                    path, px, py, style
                ));
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

