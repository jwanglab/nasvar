//! Plot types for visualizing data.

mod line;
mod scatter;

pub use line::LinePlot;
pub use scatter::ScatterPlot;

use crate::plotting::element::{Bounds, LegendEntry};
use crate::plotting::style::{FillStyle, LineStyle, MarkerStyle};

/// Trait for plot types that can be rendered.
pub trait Plot {
    /// Get the data bounds for this plot.
    fn bounds(&self) -> Option<Bounds>;

    /// Get the label for this plot (for legend).
    fn label(&self) -> Option<&str>;

    /// Get the line style for legend display.
    fn line_style(&self) -> Option<LineStyle> {
        None
    }

    /// Get the marker style for legend display.
    fn marker_style(&self) -> Option<MarkerStyle> {
        None
    }

    /// Get the fill style for legend display.
    fn fill_style(&self) -> Option<FillStyle> {
        None
    }

    /// Create a legend entry for this plot.
    fn legend_entry(&self) -> Option<LegendEntry> {
        self.label().map(|label| {
            let mut entry = LegendEntry::new(label);
            if let Some(style) = self.line_style() {
                entry = entry.line_style(style);
            }
            if let Some(style) = self.marker_style() {
                entry = entry.marker_style(style);
            }
            if let Some(style) = self.fill_style() {
                entry = entry.fill_style(style);
            }
            entry
        })
    }

    /// Render this plot to SVG, returning the SVG elements as a string.
    fn render_svg(
        &self,
        data_bounds: &Bounds,
        pixel_bounds: &Bounds,
    ) -> String;
}
