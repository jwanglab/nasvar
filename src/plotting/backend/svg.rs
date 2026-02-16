//! SVG rendering backend.

use crate::plotting::element::text::escape_xml;
use crate::plotting::style::{Color, FillStyle, LineStyle, TextStyle};

/// SVG rendering backend.
#[derive(Debug)]
pub struct SvgBackend {
    /// Image width in pixels
    pub width: f64,
    /// Image height in pixels
    pub height: f64,
    /// SVG content accumulated during rendering
    content: Vec<String>,
    /// SVG defs section (for gradients, patterns, etc.)
    defs: Vec<String>,
    /// Whether to include XML declaration
    include_declaration: bool,
}

impl SvgBackend {
    /// Create a new SVG backend with the given dimensions.
    pub fn new(width: f64, height: f64) -> Self {
        SvgBackend {
            width,
            height,
            content: Vec::new(),
            defs: Vec::new(),
            include_declaration: true,
        }
    }

    /// Set whether to include XML declaration.
    pub fn include_declaration(mut self, include: bool) -> Self {
        self.include_declaration = include;
        self
    }

    /// Add raw SVG content.
    pub fn add_content(&mut self, content: String) {
        self.content.push(content);
    }

    /// Add content to the defs section.
    pub fn add_def(&mut self, def: String) {
        self.defs.push(def);
    }

    /// Draw a line between two points.
    pub fn draw_line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, style: &LineStyle) {
        self.content.push(format!(
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" {}/>",
            x1, y1, x2, y2,
            style.to_svg_style()
        ));
    }

    /// Draw a polyline.
    pub fn draw_polyline(&mut self, points: &[(f64, f64)], style: &LineStyle) {
        if points.is_empty() {
            return;
        }

        let points_str: String = points
            .iter()
            .map(|(x, y)| format!("{:.2},{:.2}", x, y))
            .collect::<Vec<_>>()
            .join(" ");

        self.content.push(format!(
            "<polyline points=\"{}\" {}/>",
            points_str,
            style.to_svg_style()
        ));
    }

    /// Draw a path.
    pub fn draw_path(&mut self, path_data: &str, style: &LineStyle) {
        self.content.push(format!(
            "<path d=\"{}\" {}/>",
            path_data,
            style.to_svg_style()
        ));
    }

    /// Draw a rectangle.
    pub fn draw_rect(&mut self, x: f64, y: f64, width: f64, height: f64, style: &FillStyle) {
        self.content.push(format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" {}/>",
            x, y, width, height,
            style.to_svg_style()
        ));
    }

    /// Draw a circle.
    pub fn draw_circle(&mut self, cx: f64, cy: f64, r: f64, style: &FillStyle) {
        self.content.push(format!(
            "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" {}/>",
            cx, cy, r,
            style.to_svg_style()
        ));
    }

    /// Draw text.
    pub fn draw_text(&mut self, x: f64, y: f64, text: &str, style: &TextStyle) {
        let transform = if style.rotation != 0.0 {
            format!(" transform=\"rotate({},{:.2},{:.2})\"", style.rotation, x, y)
        } else {
            String::new()
        };

        self.content.push(format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" {}{}>{}</text>",
            x, y,
            style.to_svg_attrs(),
            transform,
            escape_xml(text)
        ));
    }

    /// Start a group with optional attributes.
    pub fn start_group(&mut self, attrs: &str) {
        if attrs.is_empty() {
            self.content.push("<g>".to_string());
        } else {
            self.content.push(format!("<g {}>", attrs));
        }
    }

    /// End the current group.
    pub fn end_group(&mut self) {
        self.content.push("</g>".to_string());
    }

    /// Start a clip path.
    pub fn start_clip(&mut self, id: &str, x: f64, y: f64, width: f64, height: f64) {
        self.defs.push(format!(
            "<clipPath id=\"{}\"><rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\"/></clipPath>",
            id, x, y, width, height
        ));
        self.content
            .push(format!("<g clip-path=\"url(#{})\">", id));
    }

    /// End the current clip.
    pub fn end_clip(&mut self) {
        self.content.push("</g>".to_string());
    }

    /// Render the final SVG string.
    pub fn render(self) -> String {
        let declaration = if self.include_declaration {
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        } else {
            ""
        };

        let defs_section = if self.defs.is_empty() {
            String::new()
        } else {
            format!("  <defs>\n    {}\n  </defs>\n", self.defs.join("\n    "))
        };

        format!(
            r#"{}<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">
{}{}
</svg>"#,
            declaration,
            self.width,
            self.height,
            self.width,
            self.height,
            defs_section,
            self.content.join("\n  ")
        )
    }

    /// Get the current content as a string (for debugging).
    pub fn content_preview(&self) -> String {
        self.content.join("\n")
    }
}

impl Default for SvgBackend {
    fn default() -> Self {
        Self::new(800.0, 600.0)
    }
}

/// Helper to create a linear gradient definition.
#[allow(dead_code)]
pub fn linear_gradient(
    id: &str,
    x1: &str,
    y1: &str,
    x2: &str,
    y2: &str,
    stops: &[(f64, Color)],
) -> String {
    let stops_str: String = stops
        .iter()
        .map(|(offset, color)| {
            format!(
                "<stop offset=\"{}%\" stop-color=\"{}\"/>",
                offset * 100.0,
                color.to_svg_string()
            )
        })
        .collect::<Vec<_>>()
        .join("");

    format!(
        "<linearGradient id=\"{}\" x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\">{}</linearGradient>",
        id, x1, y1, x2, y2, stops_str
    )
}

