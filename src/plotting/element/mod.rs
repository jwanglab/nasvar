//! Visual elements for plots.

mod axis;
mod grid;
mod legend;
pub mod text;

pub use axis::{AxisConfig, AxisPosition};
pub use grid::GridConfig;
pub use legend::{Legend, LegendEntry, LegendPosition};
pub use text::{escape_xml, Text};

use crate::plotting::error::PlotResult;

/// Bounding box for elements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

impl Bounds {
    /// Create a new bounds with the given values.
    pub fn new(x_min: f64, x_max: f64, y_min: f64, y_max: f64) -> Self {
        Bounds {
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    /// Create bounds from corner points.
    pub fn from_points(x1: f64, y1: f64, x2: f64, y2: f64) -> Self {
        Bounds {
            x_min: x1.min(x2),
            x_max: x1.max(x2),
            y_min: y1.min(y2),
            y_max: y1.max(y2),
        }
    }

    /// Create a unit bounds (0 to 1).
    pub fn unit() -> Self {
        Bounds::new(0.0, 1.0, 0.0, 1.0)
    }

    /// Get the width of the bounds.
    pub fn width(&self) -> f64 {
        self.x_max - self.x_min
    }

    /// Get the height of the bounds.
    pub fn height(&self) -> f64 {
        self.y_max - self.y_min
    }

    /// Get the center point.
    pub fn center(&self) -> (f64, f64) {
        (
            (self.x_min + self.x_max) / 2.0,
            (self.y_min + self.y_max) / 2.0,
        )
    }

    /// Expand bounds to include a point.
    pub fn include_point(&mut self, x: f64, y: f64) {
        self.x_min = self.x_min.min(x);
        self.x_max = self.x_max.max(x);
        self.y_min = self.y_min.min(y);
        self.y_max = self.y_max.max(y);
    }

    /// Expand bounds to include another bounds.
    pub fn include_bounds(&mut self, other: &Bounds) {
        self.x_min = self.x_min.min(other.x_min);
        self.x_max = self.x_max.max(other.x_max);
        self.y_min = self.y_min.min(other.y_min);
        self.y_max = self.y_max.max(other.y_max);
    }

    /// Add padding as a fraction of the range.
    pub fn pad(&self, fraction: f64) -> Bounds {
        let x_pad = self.width() * fraction;
        let y_pad = self.height() * fraction;
        Bounds {
            x_min: self.x_min - x_pad,
            x_max: self.x_max + x_pad,
            y_min: self.y_min - y_pad,
            y_max: self.y_max + y_pad,
        }
    }

    /// Check if a point is inside the bounds.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.x_min && x <= self.x_max && y >= self.y_min && y <= self.y_max
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Bounds::unit()
    }
}

/// Trait for renderable elements.
pub trait Renderable {
    /// Render this element to the given context.
    fn render(&self, ctx: &mut RenderContext) -> PlotResult<()>;

    /// Get the data bounds for this element.
    fn bounds(&self) -> Option<Bounds>;
}

/// Context for rendering operations.
pub struct RenderContext {
    /// The data bounds being rendered.
    pub data_bounds: Bounds,
    /// The pixel bounds to render into.
    pub pixel_bounds: Bounds,
    /// Current color cycle index.
    pub color_index: usize,
}

impl RenderContext {
    /// Transform a data point to pixel coordinates.
    pub fn transform(&self, x: f64, y: f64) -> (f64, f64) {
        let x_norm = (x - self.data_bounds.x_min) / self.data_bounds.width();
        let y_norm = (y - self.data_bounds.y_min) / self.data_bounds.height();

        let px = self.pixel_bounds.x_min + x_norm * self.pixel_bounds.width();
        // Flip Y axis since SVG has Y increasing downward
        let py = self.pixel_bounds.y_max - y_norm * self.pixel_bounds.height();

        (px, py)
    }
}

