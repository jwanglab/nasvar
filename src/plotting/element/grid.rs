//! Grid configuration and rendering.

use crate::plotting::style::{Color, DashPattern, LineStyle};

/// Configuration for grid lines.
#[derive(Debug, Clone)]
pub struct GridConfig {
    /// Whether the grid is visible
    pub visible: bool,
    /// Style for major grid lines
    pub major_style: LineStyle,
    /// Whether to show minor grid lines
    pub show_minor: bool,
    /// Style for minor grid lines
    pub minor_style: LineStyle,
    /// Number of minor grid lines between major lines
    pub minor_divisions: usize,
    /// Whether to show X grid lines
    pub show_x: bool,
    /// Whether to show Y grid lines
    pub show_y: bool,
    /// Grid line opacity
    pub grid_opacity: f64,
}

impl GridConfig {
    /// Create a new grid configuration.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set grid visibility.
    pub fn visible(mut self, visible: bool) -> Self {
        self.visible = visible;
        self
    }

    /// Set the major grid line style.
    pub fn major_style(mut self, style: LineStyle) -> Self {
        self.major_style = style;
        self
    }

    /// Set the minor grid line style.
    pub fn minor_style(mut self, style: LineStyle) -> Self {
        self.minor_style = style;
        self
    }

    /// Set whether to show minor grid lines.
    pub fn show_minor(mut self, show: bool) -> Self {
        self.show_minor = show;
        self
    }

    /// Set the number of minor divisions.
    pub fn minor_divisions(mut self, divisions: usize) -> Self {
        self.minor_divisions = divisions;
        self
    }

    /// Set the grid color (affects both major and minor).
    pub fn color(mut self, color: impl Into<Color>) -> Self {
        let c = color.into();
        self.major_style.color = c.clone();
        self.minor_style.color = c;
        self
    }

    /// Set which axes to show grid for.
    pub fn axes(mut self, show_x: bool, show_y: bool) -> Self {
        self.show_x = show_x;
        self.show_y = show_y;
        self
    }
}

impl Default for GridConfig {
    fn default() -> Self {
        GridConfig {
            visible: true,
            major_style: LineStyle::new()
                .color(Color::LIGHT_GRAY)
                .width(0.5)
                .dash(DashPattern::Solid),
            show_minor: false,
            minor_style: LineStyle::new()
                .color(Color::LIGHT_GRAY)
                .width(0.25)
                .dash(DashPattern::Dotted),
            minor_divisions: 2,
            show_x: true,
            show_y: true,
            grid_opacity: 0.8,
        }
    }
}

