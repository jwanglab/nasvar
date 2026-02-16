//! Color definitions and utilities.

use std::fmt;

/// Represents a color for plotting elements.
#[derive(Debug, Clone, PartialEq)]
pub enum Color {
    /// RGB color with values 0-255
    Rgb(u8, u8, u8),
    /// RGBA color with alpha 0.0-1.0
    Rgba(u8, u8, u8, f64),
    /// Named color (e.g., "red", "blue", "C0")
    Named(String),
}

impl Color {
    /// Create a new RGB color.
    pub fn rgb(r: u8, g: u8, b: u8) -> Self {
        Color::Rgb(r, g, b)
    }

    /// Create a new RGBA color.
    pub fn rgba(r: u8, g: u8, b: u8, a: f64) -> Self {
        Color::Rgba(r, g, b, a.clamp(0.0, 1.0))
    }

    /// Create a color from a hex string (e.g., "#FF0000" or "FF0000").
    pub fn from_hex(hex: &str) -> Option<Self> {
        let hex = hex.trim_start_matches('#');
        if hex.len() == 6 {
            let r = u8::from_str_radix(&hex[0..2], 16).ok()?;
            let g = u8::from_str_radix(&hex[2..4], 16).ok()?;
            let b = u8::from_str_radix(&hex[4..6], 16).ok()?;
            Some(Color::Rgb(r, g, b))
        } else if hex.len() == 8 {
            let r = u8::from_str_radix(&hex[0..2], 16).ok()?;
            let g = u8::from_str_radix(&hex[2..4], 16).ok()?;
            let b = u8::from_str_radix(&hex[4..6], 16).ok()?;
            let a = u8::from_str_radix(&hex[6..8], 16).ok()?;
            Some(Color::Rgba(r, g, b, a as f64 / 255.0))
        } else {
            None
        }
    }

    /// Convert the color to an SVG-compatible string.
    pub fn to_svg_string(&self) -> String {
        match self {
            Color::Rgb(r, g, b) => format!("rgb({},{},{})", r, g, b),
            Color::Rgba(r, g, b, a) => format!("rgba({},{},{},{})", r, g, b, a),
            Color::Named(name) => resolve_named_color(name),
        }
    }

    /// Get the alpha value (opacity) of the color.
    pub fn alpha(&self) -> f64 {
        match self {
            Color::Rgb(_, _, _) => 1.0,
            Color::Rgba(_, _, _, a) => *a,
            Color::Named(_) => 1.0,
        }
    }

    /// Convert to RGB tuple, resolving named colors.
    pub fn to_rgb(&self) -> (u8, u8, u8) {
        match self {
            Color::Rgb(r, g, b) => (*r, *g, *b),
            Color::Rgba(r, g, b, _) => (*r, *g, *b),
            Color::Named(name) => {
                let resolved = resolve_named_color(name);
                if let Some(color) = Color::from_hex(&resolved) {
                    color.to_rgb()
                } else {
                    (0, 0, 0) // fallback to black
                }
            }
        }
    }

    // Predefined colors
    pub const BLACK: Color = Color::Rgb(0, 0, 0);
    pub const WHITE: Color = Color::Rgb(255, 255, 255);
    pub const RED: Color = Color::Rgb(255, 0, 0);
    pub const GREEN: Color = Color::Rgb(0, 128, 0);
    pub const BLUE: Color = Color::Rgb(0, 0, 255);
    pub const YELLOW: Color = Color::Rgb(255, 255, 0);
    pub const CYAN: Color = Color::Rgb(0, 255, 255);
    pub const MAGENTA: Color = Color::Rgb(255, 0, 255);
    pub const ORANGE: Color = Color::Rgb(255, 165, 0);
    pub const PURPLE: Color = Color::Rgb(128, 0, 128);
    pub const GRAY: Color = Color::Rgb(128, 128, 128);
    pub const LIGHT_GRAY: Color = Color::Rgb(211, 211, 211);
    pub const DARK_GRAY: Color = Color::Rgb(64, 64, 64);
    pub const TRANSPARENT: Color = Color::Rgba(0, 0, 0, 0.0);
}

impl Default for Color {
    fn default() -> Self {
        Color::BLUE
    }
}

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_svg_string())
    }
}

impl From<&str> for Color {
    fn from(s: &str) -> Self {
        if s.starts_with('#') || s.chars().all(|c| c.is_ascii_hexdigit()) && s.len() >= 6 {
            Color::from_hex(s).unwrap_or_else(|| Color::Named(s.to_string()))
        } else {
            Color::Named(s.to_string())
        }
    }
}

impl From<String> for Color {
    fn from(s: String) -> Self {
        Color::from(s.as_str())
    }
}

impl From<(u8, u8, u8)> for Color {
    fn from((r, g, b): (u8, u8, u8)) -> Self {
        Color::Rgb(r, g, b)
    }
}

impl From<(u8, u8, u8, f64)> for Color {
    fn from((r, g, b, a): (u8, u8, u8, f64)) -> Self {
        Color::Rgba(r, g, b, a)
    }
}

/// Resolve a named color to its hex value.
fn resolve_named_color(name: &str) -> String {
    match name.to_lowercase().as_str() {
        // Basic colors
        "black" => "#000000".to_string(),
        "white" => "#FFFFFF".to_string(),
        "red" => "#FF0000".to_string(),
        "green" => "#008000".to_string(),
        "blue" => "#0000FF".to_string(),
        "yellow" => "#FFFF00".to_string(),
        "cyan" | "aqua" => "#00FFFF".to_string(),
        "magenta" | "fuchsia" => "#FF00FF".to_string(),
        "orange" => "#FFA500".to_string(),
        "purple" => "#800080".to_string(),
        "gray" | "grey" => "#808080".to_string(),
        "lightgray" | "lightgrey" => "#D3D3D3".to_string(),
        "darkgray" | "darkgrey" => "#404040".to_string(),
        "pink" => "#FFC0CB".to_string(),
        "brown" => "#A52A2A".to_string(),
        "navy" => "#000080".to_string(),
        "teal" => "#008080".to_string(),
        "olive" => "#808000".to_string(),
        "maroon" => "#800000".to_string(),
        "lime" => "#00FF00".to_string(),
        "silver" => "#C0C0C0".to_string(),

        // Matplotlib-style cycle colors (C0-C9)
        "c0" => "#1f77b4".to_string(),
        "c1" => "#ff7f0e".to_string(),
        "c2" => "#2ca02c".to_string(),
        "c3" => "#d62728".to_string(),
        "c4" => "#9467bd".to_string(),
        "c5" => "#8c564b".to_string(),
        "c6" => "#e377c2".to_string(),
        "c7" => "#7f7f7f".to_string(),
        "c8" => "#bcbd22".to_string(),
        "c9" => "#17becf".to_string(),

        // Single-letter shortcuts (matplotlib style)
        "b" => "#0000FF".to_string(),
        "g" => "#008000".to_string(),
        "r" => "#FF0000".to_string(),
        "c" => "#00FFFF".to_string(),
        "m" => "#FF00FF".to_string(),
        "y" => "#FFFF00".to_string(),
        "k" => "#000000".to_string(),
        "w" => "#FFFFFF".to_string(),

        // Default: return as-is (might be a valid CSS color)
        _ => name.to_string(),
    }
}

/// Get a color from the default color cycle by index.
pub fn cycle_color(index: usize) -> Color {
    let colors = [
        Color::from_hex("#1f77b4").unwrap(),
        Color::from_hex("#ff7f0e").unwrap(),
        Color::from_hex("#2ca02c").unwrap(),
        Color::from_hex("#d62728").unwrap(),
        Color::from_hex("#9467bd").unwrap(),
        Color::from_hex("#8c564b").unwrap(),
        Color::from_hex("#e377c2").unwrap(),
        Color::from_hex("#7f7f7f").unwrap(),
        Color::from_hex("#bcbd22").unwrap(),
        Color::from_hex("#17becf").unwrap(),
    ];
    colors[index % colors.len()].clone()
}

