//! Linear scale transformation.

use super::{nice_ticks, Scale};
use crate::plotting::error::{PlotError, PlotResult};

/// A linear scale for axis transformation.
#[derive(Debug, Clone)]
pub struct LinearScale {
    min: f64,
    max: f64,
}

impl LinearScale {
    /// Create a new linear scale with the given range.
    pub fn new(min: f64, max: f64) -> PlotResult<Self> {
        if min >= max {
            return Err(PlotError::InvalidConfig(
                "Scale minimum must be less than maximum".to_string(),
            ));
        }
        Ok(LinearScale { min, max })
    }

    /// Create a linear scale with automatic range (to be set later).
    pub fn auto() -> Self {
        LinearScale { min: 0.0, max: 1.0 }
    }
}

impl Default for LinearScale {
    fn default() -> Self {
        LinearScale::auto()
    }
}

impl Scale for LinearScale {
    fn transform(&self, value: f64) -> f64 {
        let range = self.max - self.min;
        if range == 0.0 {
            return 0.5;
        }
        (value - self.min) / range
    }

    fn inverse(&self, normalized: f64) -> f64 {
        let range = self.max - self.min;
        self.min + normalized * range
    }

    fn set_range(&mut self, min: f64, max: f64) -> PlotResult<()> {
        if min >= max {
            // Add a small padding if min equals max
            let padding = if min == 0.0 { 1.0 } else { min.abs() * 0.1 };
            self.min = min - padding;
            self.max = max + padding;
        } else {
            self.min = min;
            self.max = max;
        }
        Ok(())
    }

    fn range(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    fn ticks(&self, num_ticks: usize) -> Vec<f64> {
        nice_ticks(self.min, self.max, num_ticks)
    }

    fn clone_box(&self) -> Box<dyn Scale> {
        Box::new(self.clone())
    }
}

