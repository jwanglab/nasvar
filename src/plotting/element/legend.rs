//! Legend configuration and rendering.

use crate::plotting::style::{Color, FillStyle, LineStyle, MarkerStyle, TextStyle};

/// Position of the legend.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum LegendPosition {
    /// Top-left corner
    TopLeft,
    /// Top-right corner
    #[default]
    TopRight,
    /// Bottom-left corner
    BottomLeft,
    /// Bottom-right corner
    BottomRight,
    /// Center top
    Top,
    /// Center bottom
    Bottom,
    /// Center left
    Left,
    /// Center right
    Right,
    /// Center
    Center,
    /// Custom position (x, y in normalized axes coordinates)
    Custom(f64, f64),
}

impl LegendPosition {
    /// Get the anchor point for this position (in normalized coordinates).
    pub fn anchor(&self) -> (f64, f64) {
        match self {
            LegendPosition::TopLeft => (0.02, 0.98),
            LegendPosition::TopRight => (0.98, 0.98),
            LegendPosition::BottomLeft => (0.02, 0.02),
            LegendPosition::BottomRight => (0.98, 0.02),
            LegendPosition::Top => (0.5, 0.98),
            LegendPosition::Bottom => (0.5, 0.02),
            LegendPosition::Left => (0.02, 0.5),
            LegendPosition::Right => (0.98, 0.5),
            LegendPosition::Center => (0.5, 0.5),
            LegendPosition::Custom(x, y) => (*x, *y),
        }
    }

    /// Get the text anchor for this position.
    pub fn text_anchor(&self) -> &'static str {
        match self {
            LegendPosition::TopLeft
            | LegendPosition::BottomLeft
            | LegendPosition::Left => "start",
            LegendPosition::TopRight
            | LegendPosition::BottomRight
            | LegendPosition::Right => "end",
            _ => "middle",
        }
    }
}

/// A single entry in the legend.
#[derive(Debug, Clone)]
pub struct LegendEntry {
    /// Label text
    pub label: String,
    /// Line style (if applicable)
    pub line_style: Option<LineStyle>,
    /// Marker style (if applicable)
    pub marker_style: Option<MarkerStyle>,
    /// Fill style (for bar charts, etc.)
    pub fill_style: Option<FillStyle>,
}

impl LegendEntry {
    /// Create a new legend entry with just a label.
    pub fn new(label: impl Into<String>) -> Self {
        LegendEntry {
            label: label.into(),
            line_style: None,
            marker_style: None,
            fill_style: None,
        }
    }

    /// Set the line style.
    pub fn line_style(mut self, style: LineStyle) -> Self {
        self.line_style = Some(style);
        self
    }

    /// Set the marker style.
    pub fn marker_style(mut self, style: MarkerStyle) -> Self {
        self.marker_style = Some(style);
        self
    }

    /// Set the fill style.
    pub fn fill_style(mut self, style: FillStyle) -> Self {
        self.fill_style = Some(style);
        self
    }
}

/// Legend configuration.
#[derive(Debug, Clone)]
pub struct Legend {
    /// Legend entries
    pub entries: Vec<LegendEntry>,
    /// Position of the legend
    pub position: LegendPosition,
    /// Whether the legend is visible
    pub visible: bool,
    /// Background fill style
    pub background: FillStyle,
    /// Text style for labels
    pub text_style: TextStyle,
    /// Padding inside the legend box
    pub padding: f64,
    /// Spacing between entries
    pub entry_spacing: f64,
    /// Length of the line sample in the legend
    pub line_length: f64,
    /// Gap between line/marker and label
    pub label_gap: f64,
}

impl Legend {
    /// Create a new legend.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an entry to the legend.
    pub fn add_entry(&mut self, entry: LegendEntry) {
        self.entries.push(entry);
    }

    /// Set the position.
    pub fn position(mut self, position: LegendPosition) -> Self {
        self.position = position;
        self
    }

    /// Set visibility.
    pub fn visible(mut self, visible: bool) -> Self {
        self.visible = visible;
        self
    }

    /// Set the background style.
    pub fn background(mut self, style: FillStyle) -> Self {
        self.background = style;
        self
    }

    /// Set the text style.
    pub fn text_style(mut self, style: TextStyle) -> Self {
        self.text_style = style;
        self
    }
}

impl Default for Legend {
    fn default() -> Self {
        Legend {
            entries: Vec::new(),
            position: LegendPosition::TopRight,
            visible: true,
            background: FillStyle::new(Color::WHITE)
                .opacity(0.9)
                .stroke(Color::GRAY)
                .stroke_width(0.5),
            text_style: TextStyle::new().font_size(10.0),
            padding: 8.0,
            entry_spacing: 4.0,
            line_length: 20.0,
            label_gap: 8.0,
        }
    }
}

