#![allow(dead_code)]
//! Internal minimal plotting module
//! Provides matplotlib-like scatter/line plot API with SVG output.

pub mod axes;
pub mod backend;
pub mod data;
pub mod element;
pub mod error;
pub mod figure;
pub mod plot;
pub mod scale;
pub mod style;

pub use data::IntoPlotData;
pub use element::Bounds;
pub use error::{PlotError, PlotResult};
pub use figure::Figure;
pub use plot::{LinePlot, Plot, ScatterPlot};
pub use scale::{LinearScale, Scale};
pub use style::{
    Color, DashPattern, FillStyle, LineStyle, Marker, MarkerStyle, TextStyle, Theme, ThemeConfig,
};

/// Prelude module for convenient imports.
pub mod prelude {
    pub use crate::plotting::axes::Axes;
    pub use crate::plotting::data::IntoPlotData;
    pub use crate::plotting::element::Bounds;
    pub use crate::plotting::error::{PlotError, PlotResult};
    pub use crate::plotting::figure::Figure;
    pub use crate::plotting::plot::{LinePlot, ScatterPlot};
    pub use crate::plotting::style::{
        Color, DashPattern, FillStyle, LineStyle, Marker, MarkerStyle, TextStyle, Theme,
    };
}

