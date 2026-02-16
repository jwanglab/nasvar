//! Predefined themes for plots.

use super::color::Color;
use super::text_style::TextStyle;

/// Predefined plot themes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Theme {
    /// Default theme with white background
    #[default]
    Default,
    /// Dark theme with dark background
    Dark,
    /// Minimal theme with reduced visual elements
    Minimal,
    /// Seaborn-inspired theme
    Seaborn,
    /// High contrast theme for accessibility
    HighContrast,
}

/// Theme configuration containing all style settings.
#[derive(Debug, Clone)]
pub struct ThemeConfig {
    /// Background color for the figure
    pub background: Color,
    /// Background color for the plot area
    pub plot_background: Color,
    /// Color for axis lines
    pub axis_color: Color,
    /// Color for grid lines
    pub grid_color: Color,
    /// Grid line opacity
    pub grid_opacity: f64,
    /// Whether to show grid by default
    pub show_grid: bool,
    /// Color for text
    pub text_color: Color,
    /// Default title style
    pub title_style: TextStyle,
    /// Default label style
    pub label_style: TextStyle,
    /// Default tick label style
    pub tick_style: TextStyle,
    /// Color cycle for plot series
    pub color_cycle: Vec<Color>,
    /// Default line width
    pub line_width: f64,
    /// Default marker size
    pub marker_size: f64,
    /// Axis line width
    pub axis_width: f64,
    /// Grid line width
    pub grid_width: f64,
}

impl Theme {
    /// Get the configuration for this theme.
    pub fn config(&self) -> ThemeConfig {
        match self {
            Theme::Default => ThemeConfig::default_theme(),
            Theme::Dark => ThemeConfig::dark_theme(),
            Theme::Minimal => ThemeConfig::minimal_theme(),
            Theme::Seaborn => ThemeConfig::seaborn_theme(),
            Theme::HighContrast => ThemeConfig::high_contrast_theme(),
        }
    }
}

impl ThemeConfig {
    fn default_theme() -> Self {
        ThemeConfig {
            background: Color::WHITE,
            plot_background: Color::WHITE,
            axis_color: Color::BLACK,
            grid_color: Color::LIGHT_GRAY,
            grid_opacity: 0.8,
            show_grid: true,
            text_color: Color::BLACK,
            title_style: TextStyle::new()
                .font_size(14.0)
                .bold()
                .color(Color::BLACK),
            label_style: TextStyle::new()
                .font_size(12.0)
                .color(Color::BLACK),
            tick_style: TextStyle::new()
                .font_size(10.0)
                .color(Color::BLACK),
            color_cycle: default_color_cycle(),
            line_width: 1.5,
            marker_size: 6.0,
            axis_width: 1.0,
            grid_width: 0.5,
        }
    }

    fn dark_theme() -> Self {
        let text_color = Color::rgb(220, 220, 220);
        ThemeConfig {
            background: Color::rgb(30, 30, 30),
            plot_background: Color::rgb(40, 40, 40),
            axis_color: Color::rgb(180, 180, 180),
            grid_color: Color::rgb(80, 80, 80),
            grid_opacity: 0.5,
            show_grid: true,
            text_color: text_color.clone(),
            title_style: TextStyle::new()
                .font_size(14.0)
                .bold()
                .color(text_color.clone()),
            label_style: TextStyle::new()
                .font_size(12.0)
                .color(text_color.clone()),
            tick_style: TextStyle::new()
                .font_size(10.0)
                .color(text_color),
            color_cycle: dark_color_cycle(),
            line_width: 1.5,
            marker_size: 6.0,
            axis_width: 1.0,
            grid_width: 0.5,
        }
    }

    fn minimal_theme() -> Self {
        ThemeConfig {
            background: Color::WHITE,
            plot_background: Color::WHITE,
            axis_color: Color::GRAY,
            grid_color: Color::LIGHT_GRAY,
            grid_opacity: 0.3,
            show_grid: false,
            text_color: Color::DARK_GRAY,
            title_style: TextStyle::new()
                .font_size(13.0)
                .color(Color::DARK_GRAY),
            label_style: TextStyle::new()
                .font_size(11.0)
                .color(Color::GRAY),
            tick_style: TextStyle::new()
                .font_size(9.0)
                .color(Color::GRAY),
            color_cycle: default_color_cycle(),
            line_width: 1.0,
            marker_size: 5.0,
            axis_width: 0.5,
            grid_width: 0.25,
        }
    }

    fn seaborn_theme() -> Self {
        ThemeConfig {
            background: Color::WHITE,
            plot_background: Color::rgb(234, 234, 242),
            axis_color: Color::rgb(100, 100, 100),
            grid_color: Color::WHITE,
            grid_opacity: 1.0,
            show_grid: true,
            text_color: Color::rgb(50, 50, 50),
            title_style: TextStyle::new()
                .font_size(14.0)
                .bold()
                .color(Color::rgb(50, 50, 50)),
            label_style: TextStyle::new()
                .font_size(12.0)
                .color(Color::rgb(50, 50, 50)),
            tick_style: TextStyle::new()
                .font_size(10.0)
                .color(Color::rgb(100, 100, 100)),
            color_cycle: seaborn_color_cycle(),
            line_width: 1.75,
            marker_size: 6.0,
            axis_width: 1.0,
            grid_width: 1.0,
        }
    }

    fn high_contrast_theme() -> Self {
        ThemeConfig {
            background: Color::WHITE,
            plot_background: Color::WHITE,
            axis_color: Color::BLACK,
            grid_color: Color::BLACK,
            grid_opacity: 0.2,
            show_grid: true,
            text_color: Color::BLACK,
            title_style: TextStyle::new()
                .font_size(16.0)
                .bold()
                .color(Color::BLACK),
            label_style: TextStyle::new()
                .font_size(14.0)
                .bold()
                .color(Color::BLACK),
            tick_style: TextStyle::new()
                .font_size(12.0)
                .color(Color::BLACK),
            color_cycle: high_contrast_color_cycle(),
            line_width: 2.5,
            marker_size: 8.0,
            axis_width: 2.0,
            grid_width: 0.5,
        }
    }
}

impl Default for ThemeConfig {
    fn default() -> Self {
        Self::default_theme()
    }
}

fn default_color_cycle() -> Vec<Color> {
    vec![
        Color::from_hex("#1f77b4").unwrap(), // blue
        Color::from_hex("#ff7f0e").unwrap(), // orange
        Color::from_hex("#2ca02c").unwrap(), // green
        Color::from_hex("#d62728").unwrap(), // red
        Color::from_hex("#9467bd").unwrap(), // purple
        Color::from_hex("#8c564b").unwrap(), // brown
        Color::from_hex("#e377c2").unwrap(), // pink
        Color::from_hex("#7f7f7f").unwrap(), // gray
        Color::from_hex("#bcbd22").unwrap(), // olive
        Color::from_hex("#17becf").unwrap(), // cyan
    ]
}

fn dark_color_cycle() -> Vec<Color> {
    vec![
        Color::from_hex("#58a6ff").unwrap(), // bright blue
        Color::from_hex("#f0883e").unwrap(), // orange
        Color::from_hex("#3fb950").unwrap(), // green
        Color::from_hex("#f85149").unwrap(), // red
        Color::from_hex("#a371f7").unwrap(), // purple
        Color::from_hex("#db6d28").unwrap(), // brown
        Color::from_hex("#ff7b72").unwrap(), // pink
        Color::from_hex("#8b949e").unwrap(), // gray
        Color::from_hex("#d29922").unwrap(), // yellow
        Color::from_hex("#56d4dd").unwrap(), // cyan
    ]
}

fn seaborn_color_cycle() -> Vec<Color> {
    vec![
        Color::from_hex("#4c72b0").unwrap(),
        Color::from_hex("#dd8452").unwrap(),
        Color::from_hex("#55a868").unwrap(),
        Color::from_hex("#c44e52").unwrap(),
        Color::from_hex("#8172b3").unwrap(),
        Color::from_hex("#937860").unwrap(),
        Color::from_hex("#da8bc3").unwrap(),
        Color::from_hex("#8c8c8c").unwrap(),
        Color::from_hex("#ccb974").unwrap(),
        Color::from_hex("#64b5cd").unwrap(),
    ]
}

fn high_contrast_color_cycle() -> Vec<Color> {
    vec![
        Color::from_hex("#0000FF").unwrap(), // blue
        Color::from_hex("#FF0000").unwrap(), // red
        Color::from_hex("#008000").unwrap(), // green
        Color::from_hex("#FF8C00").unwrap(), // dark orange
        Color::from_hex("#800080").unwrap(), // purple
        Color::from_hex("#000000").unwrap(), // black
    ]
}

