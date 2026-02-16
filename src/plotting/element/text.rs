//! Text element for labels and annotations.

use crate::plotting::style::TextStyle;

/// A text element that can be rendered on a plot.
#[derive(Debug, Clone)]
pub struct Text {
    /// The text content
    pub content: String,
    /// X position
    pub x: f64,
    /// Y position
    pub y: f64,
    /// Style configuration
    pub style: TextStyle,
}

impl Text {
    /// Create a new text element.
    pub fn new(content: impl Into<String>, x: f64, y: f64) -> Self {
        Text {
            content: content.into(),
            x,
            y,
            style: TextStyle::default(),
        }
    }

    /// Set the style for this text.
    pub fn style(mut self, style: TextStyle) -> Self {
        self.style = style;
        self
    }

    /// Set the font size.
    pub fn font_size(mut self, size: f64) -> Self {
        self.style.font_size = size;
        self
    }

    /// Make the text bold.
    pub fn bold(mut self) -> Self {
        self.style = self.style.bold();
        self
    }

    /// Set the text color.
    pub fn color(mut self, color: impl Into<crate::plotting::style::Color>) -> Self {
        self.style.color = color.into();
        self
    }

    /// Set the rotation angle in degrees.
    pub fn rotation(mut self, degrees: f64) -> Self {
        self.style.rotation = degrees;
        self
    }

    /// Generate SVG for this text element.
    pub fn to_svg(&self) -> String {
        let transform = if self.style.rotation != 0.0 {
            format!(" transform=\"rotate({},{},{})\"", self.style.rotation, self.x, self.y)
        } else {
            String::new()
        };

        format!(
            "<text x=\"{}\" y=\"{}\" {}{}>{}</text>",
            self.x,
            self.y,
            self.style.to_svg_attrs(),
            transform,
            escape_xml(&self.content)
        )
    }
}

/// Escape special XML characters.
pub fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

