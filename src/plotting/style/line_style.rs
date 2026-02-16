//! Line styling options.

use super::color::Color;

/// Dash pattern for lines.
#[derive(Debug, Clone, PartialEq)]
pub enum DashPattern {
    /// Solid line
    Solid,
    /// Dashed line (default dash length)
    Dashed,
    /// Dotted line
    Dotted,
    /// Alternating dash-dot pattern
    DashDot,
    /// Alternating dash-dot-dot pattern
    DashDotDot,
    /// Custom dash array [dash_length, gap_length, ...]
    Custom(Vec<f64>),
}

impl DashPattern {
    /// Convert to SVG stroke-dasharray value.
    pub fn to_svg_dasharray(&self) -> Option<String> {
        match self {
            DashPattern::Solid => None,
            DashPattern::Dashed => Some("8,4".to_string()),
            DashPattern::Dotted => Some("2,2".to_string()),
            DashPattern::DashDot => Some("8,4,2,4".to_string()),
            DashPattern::DashDotDot => Some("8,4,2,4,2,4".to_string()),
            DashPattern::Custom(arr) => {
                if arr.is_empty() {
                    None
                } else {
                    Some(
                        arr.iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    )
                }
            }
        }
    }

    /// Parse from matplotlib-style format string.
    pub fn from_format_char(c: char) -> Option<Self> {
        match c {
            '-' => Some(DashPattern::Solid),
            ':' => Some(DashPattern::Dotted),
            _ => None,
        }
    }

    /// Parse from matplotlib-style format string (prefix match).
    pub fn from_format_str(s: &str) -> Option<Self> {
        // Check longest patterns first
        if s.starts_with("--") {
            Some(DashPattern::Dashed)
        } else if s.starts_with("-.") {
            Some(DashPattern::DashDot)
        } else if s.starts_with(':') {
            Some(DashPattern::Dotted)
        } else if s.starts_with('-') {
            Some(DashPattern::Solid)
        } else {
            None
        }
    }
}

impl Default for DashPattern {
    fn default() -> Self {
        DashPattern::Solid
    }
}

/// Style configuration for lines.
#[derive(Debug, Clone)]
pub struct LineStyle {
    /// Line color
    pub color: Color,
    /// Line width in pixels
    pub width: f64,
    /// Dash pattern
    pub dash: DashPattern,
    /// Line cap style
    pub cap: LineCap,
    /// Line join style
    pub join: LineJoin,
    /// Opacity (0.0 - 1.0)
    pub opacity: f64,
}

impl LineStyle {
    /// Create a new line style with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the line color.
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.color = color.into();
        self
    }

    /// Set the line width.
    pub fn width(mut self, width: f64) -> Self {
        self.width = width;
        self
    }

    /// Set the dash pattern.
    pub fn dash(mut self, dash: DashPattern) -> Self {
        self.dash = dash;
        self
    }

    /// Set the line cap style.
    pub fn cap(mut self, cap: LineCap) -> Self {
        self.cap = cap;
        self
    }

    /// Set the line join style.
    pub fn join(mut self, join: LineJoin) -> Self {
        self.join = join;
        self
    }

    /// Set the opacity.
    pub fn opacity(mut self, opacity: f64) -> Self {
        self.opacity = opacity.clamp(0.0, 1.0);
        self
    }

    /// Generate SVG style attributes.
    pub fn to_svg_style(&self) -> String {
        let mut attrs = vec![
            format!("stroke=\"{}\"", self.color.to_svg_string()),
            format!("stroke-width=\"{}\"", self.width),
            format!("stroke-linecap=\"{}\"", self.cap.to_svg_string()),
            format!("stroke-linejoin=\"{}\"", self.join.to_svg_string()),
            "fill=\"none\"".to_string(),
        ];

        if self.opacity < 1.0 {
            attrs.push(format!("stroke-opacity=\"{}\"", self.opacity));
        }

        if let Some(dasharray) = self.dash.to_svg_dasharray() {
            attrs.push(format!("stroke-dasharray=\"{}\"", dasharray));
        }

        attrs.join(" ")
    }
}

impl Default for LineStyle {
    fn default() -> Self {
        LineStyle {
            color: Color::default(),
            width: 1.5,
            dash: DashPattern::Solid,
            cap: LineCap::Round,
            join: LineJoin::Round,
            opacity: 1.0,
        }
    }
}

/// Line cap styles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LineCap {
    /// Flat end at the exact endpoint
    Butt,
    /// Rounded end
    #[default]
    Round,
    /// Square end extending past the endpoint
    Square,
}

impl LineCap {
    pub fn to_svg_string(&self) -> &'static str {
        match self {
            LineCap::Butt => "butt",
            LineCap::Round => "round",
            LineCap::Square => "square",
        }
    }
}

/// Line join styles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LineJoin {
    /// Sharp corner
    Miter,
    /// Rounded corner
    #[default]
    Round,
    /// Beveled corner
    Bevel,
}

impl LineJoin {
    pub fn to_svg_string(&self) -> &'static str {
        match self {
            LineJoin::Miter => "miter",
            LineJoin::Round => "round",
            LineJoin::Bevel => "bevel",
        }
    }
}

