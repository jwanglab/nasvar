//! Text styling options.

use super::color::Color;

/// Font weight options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FontWeight {
    /// Normal weight
    #[default]
    Normal,
    /// Bold weight
    Bold,
    /// Light weight
    Light,
    /// Numeric weight (100-900)
    Numeric(u16),
}

impl FontWeight {
    pub fn to_svg_string(&self) -> String {
        match self {
            FontWeight::Normal => "normal".to_string(),
            FontWeight::Bold => "bold".to_string(),
            FontWeight::Light => "300".to_string(),
            FontWeight::Numeric(n) => n.to_string(),
        }
    }
}

/// Font style options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FontStyle {
    /// Normal style
    #[default]
    Normal,
    /// Italic style
    Italic,
    /// Oblique style
    Oblique,
}

impl FontStyle {
    pub fn to_svg_string(&self) -> &'static str {
        match self {
            FontStyle::Normal => "normal",
            FontStyle::Italic => "italic",
            FontStyle::Oblique => "oblique",
        }
    }
}

/// Text anchor position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TextAnchor {
    /// Anchor at the start (left for LTR text)
    #[default]
    Start,
    /// Anchor at the middle
    Middle,
    /// Anchor at the end (right for LTR text)
    End,
}

impl TextAnchor {
    pub fn to_svg_string(&self) -> &'static str {
        match self {
            TextAnchor::Start => "start",
            TextAnchor::Middle => "middle",
            TextAnchor::End => "end",
        }
    }
}

/// Vertical alignment for text.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DominantBaseline {
    /// Align to the baseline
    #[default]
    Auto,
    /// Align to the middle
    Middle,
    /// Align to the top (hanging)
    Hanging,
    /// Align to the alphabetic baseline
    Alphabetic,
    /// Align to the ideographic baseline
    Ideographic,
}

impl DominantBaseline {
    pub fn to_svg_string(&self) -> &'static str {
        match self {
            DominantBaseline::Auto => "auto",
            DominantBaseline::Middle => "middle",
            DominantBaseline::Hanging => "hanging",
            DominantBaseline::Alphabetic => "alphabetic",
            DominantBaseline::Ideographic => "ideographic",
        }
    }
}

/// Style configuration for text elements.
#[derive(Debug, Clone)]
pub struct TextStyle {
    /// Font family (e.g., "Arial", "sans-serif")
    pub font_family: String,
    /// Font size in pixels
    pub font_size: f64,
    /// Font weight
    pub weight: FontWeight,
    /// Font style
    pub style: FontStyle,
    /// Text color
    pub color: Color,
    /// Horizontal anchor
    pub anchor: TextAnchor,
    /// Vertical alignment
    pub baseline: DominantBaseline,
    /// Rotation angle in degrees
    pub rotation: f64,
    /// Opacity
    pub opacity: f64,
}

impl TextStyle {
    /// Create a new text style with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the font family.
    pub fn font_family(mut self, family: impl Into<String>) -> Self {
        self.font_family = family.into();
        self
    }

    /// Set the font size.
    pub fn font_size(mut self, size: f64) -> Self {
        self.font_size = size;
        self
    }

    /// Set the font weight.
    pub fn weight(mut self, weight: FontWeight) -> Self {
        self.weight = weight;
        self
    }

    /// Set bold weight.
    pub fn bold(mut self) -> Self {
        self.weight = FontWeight::Bold;
        self
    }

    /// Set the font style.
    pub fn style(mut self, style: FontStyle) -> Self {
        self.style = style;
        self
    }

    /// Set italic style.
    pub fn italic(mut self) -> Self {
        self.style = FontStyle::Italic;
        self
    }

    /// Set the text color.
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        self.color = color.into();
        self
    }

    /// Set the text anchor.
    pub fn anchor(mut self, anchor: TextAnchor) -> Self {
        self.anchor = anchor;
        self
    }

    /// Set the dominant baseline.
    pub fn baseline(mut self, baseline: DominantBaseline) -> Self {
        self.baseline = baseline;
        self
    }

    /// Set the rotation angle in degrees.
    pub fn rotation(mut self, degrees: f64) -> Self {
        self.rotation = degrees;
        self
    }

    /// Set the opacity.
    pub fn opacity(mut self, opacity: f64) -> Self {
        self.opacity = opacity.clamp(0.0, 1.0);
        self
    }

    /// Generate SVG style attributes (excluding positioning attributes).
    pub fn to_svg_attrs(&self) -> String {
        let mut attrs = vec![
            format!("font-family=\"{}\"", self.font_family),
            format!("font-size=\"{}\"", self.font_size),
            format!("font-weight=\"{}\"", self.weight.to_svg_string()),
            format!("font-style=\"{}\"", self.style.to_svg_string()),
            format!("fill=\"{}\"", self.color.to_svg_string()),
        ];

        if self.opacity < 1.0 {
            attrs.push(format!("fill-opacity=\"{}\"", self.opacity));
        }

        attrs.join(" ")
    }

    /// Generate full SVG style attributes including positioning.
    pub fn to_svg_attrs_full(&self) -> String {
        let mut attrs = vec![
            format!("font-family=\"{}\"", self.font_family),
            format!("font-size=\"{}\"", self.font_size),
            format!("font-weight=\"{}\"", self.weight.to_svg_string()),
            format!("font-style=\"{}\"", self.style.to_svg_string()),
            format!("fill=\"{}\"", self.color.to_svg_string()),
            format!("text-anchor=\"{}\"", self.anchor.to_svg_string()),
            format!("dominant-baseline=\"{}\"", self.baseline.to_svg_string()),
        ];

        if self.opacity < 1.0 {
            attrs.push(format!("fill-opacity=\"{}\"", self.opacity));
        }

        attrs.join(" ")
    }
}

impl Default for TextStyle {
    fn default() -> Self {
        TextStyle {
            font_family: "sans-serif".to_string(),
            font_size: 12.0,
            weight: FontWeight::Normal,
            style: FontStyle::Normal,
            color: Color::BLACK,
            anchor: TextAnchor::Start,
            baseline: DominantBaseline::Auto,
            rotation: 0.0,
            opacity: 1.0,
        }
    }
}

