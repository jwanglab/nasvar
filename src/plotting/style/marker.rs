//! Marker styles for scatter plots and line plot points.

use super::color::Color;

/// Marker shapes for data points.
#[derive(Debug, Clone, PartialEq)]
pub enum Marker {
    /// No marker
    None,
    /// Circle marker
    Circle,
    /// Square marker
    Square,
    /// Upward-pointing triangle
    Triangle,
    /// Downward-pointing triangle
    TriangleDown,
    /// Diamond marker
    Diamond,
    /// Plus sign
    Plus,
    /// X/Cross marker
    Cross,
    /// Star marker
    Star,
    /// Pentagon marker
    Pentagon,
    /// Hexagon marker
    Hexagon,
    /// Custom SVG path
    Custom(String),
}

impl Marker {
    /// Parse from matplotlib-style format character.
    pub fn from_format_char(c: char) -> Option<Self> {
        match c {
            'o' => Some(Marker::Circle),
            's' => Some(Marker::Square),
            '^' => Some(Marker::Triangle),
            'v' => Some(Marker::TriangleDown),
            'D' | 'd' => Some(Marker::Diamond),
            '+' => Some(Marker::Plus),
            'x' | 'X' => Some(Marker::Cross),
            '*' => Some(Marker::Star),
            'p' => Some(Marker::Pentagon),
            'h' | 'H' => Some(Marker::Hexagon),
            _ => None,
        }
    }

    /// Generate SVG path data for the marker centered at origin.
    /// The size parameter is the radius (half the marker size).
    pub fn to_svg_path(&self, size: f64) -> Option<String> {
        match self {
            Marker::None => None,
            Marker::Circle => None, // Use <circle> element instead
            Marker::Square => {
                let s = size;
                Some(format!("M{},{} L{},{} L{},{} L{},{} Z",
                    -s, -s, s, -s, s, s, -s, s))
            }
            Marker::Triangle => {
                let h = size * 1.1547; // Equilateral triangle height factor
                Some(format!("M0,{} L{},{} L{},{} Z",
                    -h, -size, h * 0.5, size, h * 0.5))
            }
            Marker::TriangleDown => {
                let h = size * 1.1547;
                Some(format!("M0,{} L{},{} L{},{} Z",
                    h, -size, -h * 0.5, size, -h * 0.5))
            }
            Marker::Diamond => {
                let s = size * 1.2;
                Some(format!("M0,{} L{},0 L0,{} L{},0 Z",
                    -s, s, s, -s))
            }
            Marker::Plus => {
                let s = size;
                let w = size * 0.3;
                Some(format!(
                    "M{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} L{},{} Z",
                    -w, -s, w, -s, w, -w, s, -w, s, w, w, w, w, s, -w, s, -w, w, -s, w, -s, -w, -w, -w
                ))
            }
            Marker::Cross => {
                let s = size * 0.707; // 1/sqrt(2)
                // Simple X shape using two crossed lines
                Some(format!(
                    "M{:.2},{:.2} L{:.2},{:.2} M{:.2},{:.2} L{:.2},{:.2}",
                    -s, -s, s, s, -s, s, s, -s
                ))
            }
            Marker::Star => {
                // 5-pointed star
                let outer = size;
                let inner = size * 0.4;
                let mut path = String::new();
                for i in 0..10 {
                    let r = if i % 2 == 0 { outer } else { inner };
                    let angle = std::f64::consts::PI * (i as f64) / 5.0 - std::f64::consts::FRAC_PI_2;
                    let x = r * angle.cos();
                    let y = r * angle.sin();
                    if i == 0 {
                        path.push_str(&format!("M{:.2},{:.2}", x, y));
                    } else {
                        path.push_str(&format!(" L{:.2},{:.2}", x, y));
                    }
                }
                path.push_str(" Z");
                Some(path)
            }
            Marker::Pentagon => {
                let r = size;
                let mut path = String::new();
                for i in 0..5 {
                    let angle = std::f64::consts::PI * 2.0 * (i as f64) / 5.0 - std::f64::consts::FRAC_PI_2;
                    let x = r * angle.cos();
                    let y = r * angle.sin();
                    if i == 0 {
                        path.push_str(&format!("M{:.2},{:.2}", x, y));
                    } else {
                        path.push_str(&format!(" L{:.2},{:.2}", x, y));
                    }
                }
                path.push_str(" Z");
                Some(path)
            }
            Marker::Hexagon => {
                let r = size;
                let mut path = String::new();
                for i in 0..6 {
                    let angle = std::f64::consts::PI * (i as f64) / 3.0;
                    let x = r * angle.cos();
                    let y = r * angle.sin();
                    if i == 0 {
                        path.push_str(&format!("M{:.2},{:.2}", x, y));
                    } else {
                        path.push_str(&format!(" L{:.2},{:.2}", x, y));
                    }
                }
                path.push_str(" Z");
                Some(path)
            }
            Marker::Custom(path) => Some(path.clone()),
        }
    }

    /// Check if this marker should be rendered as a circle element.
    pub fn is_circle(&self) -> bool {
        matches!(self, Marker::Circle)
    }
}

impl Default for Marker {
    fn default() -> Self {
        Marker::None
    }
}

/// Style configuration for markers.
#[derive(Debug, Clone)]
pub struct MarkerStyle {
    /// The marker shape
    pub marker: Marker,
    /// Marker size (diameter in pixels)
    pub size: f64,
    /// Fill color
    pub fill: Color,
    /// Edge/stroke color
    pub edge_color: Color,
    /// Edge/stroke width
    pub edge_width: f64,
    /// Fill opacity
    pub fill_opacity: f64,
}

impl MarkerStyle {
    /// Create a new marker style.
    pub fn new(marker: Marker) -> Self {
        MarkerStyle {
            marker,
            ..Default::default()
        }
    }

    /// Set the marker shape.
    pub fn marker(mut self, marker: Marker) -> Self {
        self.marker = marker;
        self
    }

    /// Set the marker size.
    pub fn size(mut self, size: f64) -> Self {
        self.size = size;
        self
    }

    /// Set the fill color.
    pub fn fill(mut self, color: impl Into<Color>) -> Self {
        self.fill = color.into();
        self
    }

    /// Set the edge color.
    pub fn edge_color(mut self, color: impl Into<Color>) -> Self {
        self.edge_color = color.into();
        self
    }

    /// Set the edge width.
    pub fn edge_width(mut self, width: f64) -> Self {
        self.edge_width = width;
        self
    }

    /// Set the fill opacity.
    pub fn fill_opacity(mut self, opacity: f64) -> Self {
        self.fill_opacity = opacity.clamp(0.0, 1.0);
        self
    }

    /// Generate SVG style attributes.
    pub fn to_svg_style(&self) -> String {
        let mut attrs = vec![
            format!("fill=\"{}\"", self.fill.to_svg_string()),
            format!("stroke=\"{}\"", self.edge_color.to_svg_string()),
            format!("stroke-width=\"{}\"", self.edge_width),
        ];

        if self.fill_opacity < 1.0 {
            attrs.push(format!("fill-opacity=\"{}\"", self.fill_opacity));
        }

        attrs.join(" ")
    }

    /// Render the marker at a specific position, returning SVG elements.
    pub fn render_at(&self, x: f64, y: f64) -> String {
        let style = self.to_svg_style();
        let radius = self.size / 2.0;

        if self.marker.is_circle() {
            format!(
                "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" {}/>\n",
                x, y, radius, style
            )
        } else if let Some(path) = self.marker.to_svg_path(radius) {
            format!(
                "<path d=\"{}\" transform=\"translate({:.2},{:.2})\" {}/>\n",
                path, x, y, style
            )
        } else {
            String::new()
        }
    }
}

impl Default for MarkerStyle {
    fn default() -> Self {
        MarkerStyle {
            marker: Marker::Circle,
            size: 6.0,
            fill: Color::default(),
            edge_color: Color::default(),
            edge_width: 1.0,
            fill_opacity: 1.0,
        }
    }
}

