//! Figure (canvas) implementation.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::plotting::axes::Axes;
use crate::plotting::backend::SvgBackend;
use crate::plotting::error::PlotResult;
use crate::plotting::style::{Color, Theme, ThemeConfig};

/// A figure containing one or more axes (subplots).
pub struct Figure {
    /// Figure width in pixels
    pub width: f64,
    /// Figure height in pixels
    pub height: f64,
    /// Background color
    pub background: Color,
    /// Axes (subplots) in this figure
    axes: Vec<Axes>,
    /// Theme configuration
    pub theme: ThemeConfig,
    /// Figure title
    pub title: Option<String>,
}

impl Figure {
    /// Create a new figure with the given dimensions.
    pub fn new(width: f64, height: f64) -> Self {
        let theme = Theme::Default.config();
        Figure {
            width,
            height,
            background: theme.background.clone(),
            axes: Vec::new(),
            theme,
            title: None,
        }
    }

    /// Create a figure with default dimensions (800x600).
    pub fn default_size() -> Self {
        Self::new(800.0, 600.0)
    }

    /// Set the figure size.
    pub fn size(mut self, width: f64, height: f64) -> Self {
        self.width = width;
        self.height = height;
        self
    }

    /// Set the background color.
    pub fn background(mut self, color: impl Into<Color>) -> Self {
        self.background = color.into();
        self
    }

    /// Set the theme.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme.config();
        self.background = self.theme.background.clone();
        // Update existing axes themes
        for ax in &mut self.axes {
            ax.theme = self.theme.clone();
        }
        self
    }

    /// Set the figure title.
    pub fn suptitle(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Add a subplot at the given position.
    /// Uses matplotlib-style indexing: (rows, cols, index) where index is 1-based.
    pub fn add_subplot(&mut self, rows: usize, cols: usize, index: usize) -> &mut Axes {
        let index = index.saturating_sub(1); // Convert to 0-based
        let row = index / cols;
        let col = index % cols;

        // Calculate position with margins
        let margin = 0.06;
        let subplot_width = (1.0 - 2.0 * margin) / cols as f64;
        let subplot_height = (1.0 - 2.0 * margin) / rows as f64;

        let left = margin + col as f64 * subplot_width + 0.04;
        let right = margin + (col + 1) as f64 * subplot_width - 0.01;
        let bottom = margin + (rows - 1 - row) as f64 * subplot_height + 0.05;
        let top = margin + (rows - row) as f64 * subplot_height - 0.02;

        let mut axes = Axes::new();
        axes.position = crate::plotting::element::Bounds::new(left, right, bottom, top);
        axes.theme = self.theme.clone();
        axes.grid.visible = self.theme.show_grid;

        self.axes.push(axes);
        self.axes.last_mut().unwrap()
    }

    /// Create a grid of subplots and return references to all of them.
    pub fn subplots(&mut self, rows: usize, cols: usize) -> Vec<&mut Axes> {
        for i in 1..=(rows * cols) {
            self.add_subplot(rows, cols, i);
        }

        // Return mutable references to all axes
        self.axes.iter_mut().collect()
    }

    /// Get the first (or only) axes, creating it if necessary.
    pub fn gca(&mut self) -> &mut Axes {
        if self.axes.is_empty() {
            self.add_subplot(1, 1, 1);
        }
        self.axes.last_mut().unwrap()
    }

    /// Get all axes.
    pub fn get_axes(&mut self) -> &mut [Axes] {
        &mut self.axes
    }

    /// Adjust subplot parameters to give specified padding.
    pub fn tight_layout(&mut self) {
        // Simple implementation - just ensure reasonable margins
        let n = self.axes.len();
        if n == 0 {
            return;
        }

        // For now, just leave the positions as set by add_subplot
        // A more sophisticated implementation would measure text extents
    }

    /// Render the figure to an SVG string.
    pub fn render(&mut self) -> String {
        let mut backend = SvgBackend::new(self.width, self.height);

        // Background
        backend.add_content(format!(
            "<rect width=\"{}\" height=\"{}\" fill=\"{}\"/>",
            self.width,
            self.height,
            self.background.to_svg_string()
        ));

        // Figure title
        if let Some(ref title) = self.title {
            backend.add_content(format!(
                "<text x=\"{}\" y=\"30\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">{}</text>",
                self.width / 2.0,
                crate::plotting::element::escape_xml(title)
            ));
        }

        // Render each axes
        for axes in &mut self.axes {
            let svg = axes.render_svg(self.width, self.height);
            backend.add_content(svg);
        }

        backend.render()
    }

    /// Save the figure to a file.
    pub fn save(&mut self, path: impl AsRef<Path>) -> PlotResult<()> {
        let svg = self.render();
        let mut file = File::create(path)?;
        file.write_all(svg.as_bytes())?;
        Ok(())
    }

    /// Save the figure to a file and return self for chaining.
    pub fn savefig(mut self, path: impl AsRef<Path>) -> PlotResult<Self> {
        self.save(path)?;
        Ok(self)
    }
}

impl Default for Figure {
    fn default() -> Self {
        Self::default_size()
    }
}

