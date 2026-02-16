//! Error types for the rustplot library.

use std::fmt;
use std::io;

/// The main error type for rustplot operations.
#[derive(Debug)]
pub enum PlotError {
    /// Error during IO operations (file writing, etc.)
    Io(io::Error),
    /// Invalid data provided for plotting
    InvalidData(String),
    /// Invalid configuration or parameters
    InvalidConfig(String),
    /// Rendering error
    RenderError(String),
    /// Empty data provided where non-empty data is required
    EmptyData,
}

impl fmt::Display for PlotError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PlotError::Io(err) => write!(f, "IO error: {}", err),
            PlotError::InvalidData(msg) => write!(f, "Invalid data: {}", msg),
            PlotError::InvalidConfig(msg) => write!(f, "Invalid configuration: {}", msg),
            PlotError::RenderError(msg) => write!(f, "Render error: {}", msg),
            PlotError::EmptyData => write!(f, "Empty data provided"),
        }
    }
}

impl std::error::Error for PlotError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            PlotError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for PlotError {
    fn from(err: io::Error) -> Self {
        PlotError::Io(err)
    }
}

/// Result type alias for rustplot operations.
pub type PlotResult<T> = Result<T, PlotError>;
