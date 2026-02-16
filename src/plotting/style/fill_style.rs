//! Fill styling options for shapes.

use super::color::Color;

/// Style configuration for filled shapes.
#[derive(Debug, Clone)]
pub struct FillStyle {
    /// Fill color
    pub color: Color,
    /// Fill opacity (0.0 - 1.0)
    pub opacity: f64,
    /// Stroke/border color (None for no stroke)
    pub stroke: Option<Color>,
    /// Stroke width
    pub stroke_width: f64,
}

impl FillStyle {
    /// Create a new fill style with the given color.
    pub fn new(color: impl Into<Color>) -> Self {
        FillStyle {
            color: color.into(),
            ..Default::default()
        }
    }

    /// Create a fill style with no fill (transparent).
    pub fn none() -> Self {
        FillStyle {
            color: Color::TRANSPARENT,
            opacity: 0.0,
            stroke: None,
            stroke_width: 0.0,
        }
    }

    /// Set the fill color.
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.color = color.into();
        self
    }

    /// Set the fill opacity.
    pub fn opacity(mut self, opacity: f64) -> Self {
        self.opacity = opacity.clamp(0.0, 1.0);
        self
    }

    /// Set the stroke color.
    pub fn stroke(mut self, color: impl Into<Color>) -> Self {
        self.stroke = Some(color.into());
        self
    }

    /// Set the stroke width.
    pub fn stroke_width(mut self, width: f64) -> Self {
        self.stroke_width = width;
        self
    }

    /// Generate SVG style attributes.
    pub fn to_svg_style(&self) -> String {
        let mut attrs = vec![format!("fill=\"{}\"", self.color.to_svg_string())];

        if self.opacity < 1.0 {
            attrs.push(format!("fill-opacity=\"{}\"", self.opacity));
        }

        if let Some(ref stroke) = self.stroke {
            attrs.push(format!("stroke=\"{}\"", stroke.to_svg_string()));
            attrs.push(format!("stroke-width=\"{}\"", self.stroke_width));
        } else {
            attrs.push("stroke=\"none\"".to_string());
        }

        attrs.join(" ")
    }
}

impl Default for FillStyle {
    fn default() -> Self {
        FillStyle {
            color: Color::default(),
            opacity: 1.0,
            stroke: None,
            stroke_width: 1.0,
        }
    }
}

impl From<Color> for FillStyle {
    fn from(color: Color) -> Self {
        FillStyle::new(color)
    }
}

impl From<&str> for FillStyle {
    fn from(s: &str) -> Self {
        FillStyle::new(Color::from(s))
    }
}

