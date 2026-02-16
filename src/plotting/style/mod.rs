//! Styling module for rustplot.
//!
//! This module contains all style-related types including colors,
//! line styles, markers, text styles, and themes.

pub mod color;
pub mod fill_style;
pub mod line_style;
pub mod marker;
pub mod text_style;
pub mod theme;

pub use color::{cycle_color, Color};
pub use fill_style::FillStyle;
pub use line_style::{DashPattern, LineCap, LineJoin, LineStyle};
pub use marker::{Marker, MarkerStyle};
pub use text_style::{DominantBaseline, FontStyle, FontWeight, TextAnchor, TextStyle};
pub use theme::{Theme, ThemeConfig};
