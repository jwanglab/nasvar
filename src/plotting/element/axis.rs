//! Axis configuration and rendering.

use crate::plotting::style::{Color, LineStyle, TextStyle};

/// Position of an axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisPosition {
    Bottom,
    Top,
    Left,
    Right,
}

/// Configuration for an axis.
#[derive(Debug, Clone)]
pub struct AxisConfig {
    /// Whether the axis is visible
    pub visible: bool,
    /// Axis line style
    pub line_style: LineStyle,
    /// Whether to show tick marks
    pub show_ticks: bool,
    /// Length of tick marks in pixels
    pub tick_length: f64,
    /// Style for tick labels
    pub tick_label_style: TextStyle,
    /// Number of ticks to generate
    pub num_ticks: usize,
    /// Padding between tick marks and labels
    pub tick_padding: f64,
    /// Custom tick positions (overrides automatic generation)
    pub tick_positions: Option<Vec<f64>>,
    /// Custom tick labels (must match tick_positions length)
    pub tick_labels: Option<Vec<String>>,
    /// Format string for tick labels (printf-style)
    pub tick_format: Option<String>,
}

impl AxisConfig {
    /// Create a new axis configuration with defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set axis visibility.
    pub fn visible(mut self, visible: bool) -> Self {
        self.visible = visible;
        self
    }

    /// Set the line style.
    pub fn line_style(mut self, style: LineStyle) -> Self {
        self.line_style = style;
        self
    }

    /// Set whether to show tick marks.
    pub fn show_ticks(mut self, show: bool) -> Self {
        self.show_ticks = show;
        self
    }

    /// Set the tick length.
    pub fn tick_length(mut self, length: f64) -> Self {
        self.tick_length = length;
        self
    }

    /// Set the number of ticks.
    pub fn num_ticks(mut self, num: usize) -> Self {
        self.num_ticks = num;
        self
    }

    /// Set custom tick positions.
    pub fn tick_positions(mut self, positions: Vec<f64>) -> Self {
        self.tick_positions = Some(positions);
        self
    }

    /// Set custom tick labels.
    pub fn tick_labels(mut self, labels: Vec<String>) -> Self {
        self.tick_labels = Some(labels);
        self
    }

    /// Set the tick label format.
    pub fn tick_format(mut self, format: impl Into<String>) -> Self {
        self.tick_format = Some(format.into());
        self
    }

    /// Format a tick value as a label.
    pub fn format_tick(&self, value: f64) -> String {
        if let Some(ref format) = self.tick_format {
            // Simple format parsing (supports basic %f, %e, %g patterns)
            if format.contains("%e") || format.contains("%E") {
                format!("{:e}", value)
            } else if format.contains("%g") || format.contains("%G") {
                // General format - use exponential for very large/small numbers
                if value.abs() >= 1e6 || (value != 0.0 && value.abs() < 1e-4) {
                    format!("{:e}", value)
                } else {
                    format!("{}", value)
                }
            } else {
                // Default precision
                let precision = format
                    .find('.')
                    .and_then(|i| format[i + 1..].chars().next())
                    .and_then(|c| c.to_digit(10))
                    .unwrap_or(2) as usize;
                format!("{:.prec$}", value, prec = precision)
            }
        } else {
            // Auto-format based on value magnitude
            if value == 0.0 {
                "0".to_string()
            } else if value.abs() >= 1000.0 || value.abs() < 0.01 {
                format!("{:.2e}", value)
            } else if value.fract() == 0.0 {
                format!("{:.0}", value)
            } else {
                format!("{:.2}", value)
            }
        }
    }
}

impl Default for AxisConfig {
    fn default() -> Self {
        AxisConfig {
            visible: true,
            line_style: LineStyle::new().color(Color::BLACK).width(1.0),
            show_ticks: true,
            tick_length: 5.0,
            tick_label_style: TextStyle::new().font_size(10.0),
            num_ticks: 5,
            tick_padding: 3.0,
            tick_positions: None,
            tick_labels: None,
            tick_format: None,
        }
    }
}

