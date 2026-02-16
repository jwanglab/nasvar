//! Axis scaling and transformation.

mod linear;

pub use linear::LinearScale;

use crate::plotting::error::PlotResult;

/// Trait for axis scale transformations.
pub trait Scale: std::fmt::Debug {
    /// Transform a data value to normalized coordinates (0.0 - 1.0).
    fn transform(&self, value: f64) -> f64;

    /// Inverse transform from normalized coordinates to data value.
    fn inverse(&self, normalized: f64) -> f64;

    /// Set the data range for this scale.
    fn set_range(&mut self, min: f64, max: f64) -> PlotResult<()>;

    /// Get the current data range.
    fn range(&self) -> (f64, f64);

    /// Generate nice tick values for this scale.
    fn ticks(&self, num_ticks: usize) -> Vec<f64>;

    /// Clone the scale into a boxed trait object.
    fn clone_box(&self) -> Box<dyn Scale>;
}

impl Clone for Box<dyn Scale> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Compute "nice" numbers for axis ticks.
pub fn nice_number(range: f64, round: bool) -> f64 {
    let exponent = range.log10().floor();
    let fraction = range / 10_f64.powf(exponent);

    let nice_fraction = if round {
        if fraction < 1.5 {
            1.0
        } else if fraction < 3.0 {
            2.0
        } else if fraction < 7.0 {
            5.0
        } else {
            10.0
        }
    } else if fraction <= 1.0 {
        1.0
    } else if fraction <= 2.0 {
        2.0
    } else if fraction <= 5.0 {
        5.0
    } else {
        10.0
    };

    nice_fraction * 10_f64.powf(exponent)
}

/// Generate nice tick positions for a given range.
pub fn nice_ticks(min: f64, max: f64, num_ticks: usize) -> Vec<f64> {
    if num_ticks < 2 {
        return vec![(min + max) / 2.0];
    }

    let range = nice_number(max - min, false);
    let tick_spacing = nice_number(range / (num_ticks - 1) as f64, true);
    let nice_min = (min / tick_spacing).floor() * tick_spacing;
    let nice_max = (max / tick_spacing).ceil() * tick_spacing;

    let mut ticks = Vec::new();
    let mut tick = nice_min;
    while tick <= nice_max + tick_spacing * 0.5 {
        if tick >= min - tick_spacing * 0.001 && tick <= max + tick_spacing * 0.001 {
            ticks.push(tick);
        }
        tick += tick_spacing;
    }

    ticks
}

