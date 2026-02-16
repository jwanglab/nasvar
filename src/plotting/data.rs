//! Data handling traits and utilities.

/// Trait for types that can be converted into plot data.
pub trait IntoPlotData {
    /// Convert into a vector of f64 values.
    fn into_plot_data(self) -> Vec<f64>;
}

impl IntoPlotData for Vec<f64> {
    fn into_plot_data(self) -> Vec<f64> {
        self
    }
}

impl IntoPlotData for &Vec<f64> {
    fn into_plot_data(self) -> Vec<f64> {
        self.clone()
    }
}

impl IntoPlotData for &[f64] {
    fn into_plot_data(self) -> Vec<f64> {
        self.to_vec()
    }
}

impl IntoPlotData for Vec<f32> {
    fn into_plot_data(self) -> Vec<f64> {
        self.into_iter().map(|x| x as f64).collect()
    }
}

impl IntoPlotData for &[f32] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

impl IntoPlotData for Vec<i32> {
    fn into_plot_data(self) -> Vec<f64> {
        self.into_iter().map(|x| x as f64).collect()
    }
}

impl IntoPlotData for &[i32] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

impl IntoPlotData for Vec<i64> {
    fn into_plot_data(self) -> Vec<f64> {
        self.into_iter().map(|x| x as f64).collect()
    }
}

impl IntoPlotData for &[i64] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

impl IntoPlotData for Vec<usize> {
    fn into_plot_data(self) -> Vec<f64> {
        self.into_iter().map(|x| x as f64).collect()
    }
}

impl IntoPlotData for &[usize] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

impl<const N: usize> IntoPlotData for [f64; N] {
    fn into_plot_data(self) -> Vec<f64> {
        self.to_vec()
    }
}

impl<const N: usize> IntoPlotData for &[f64; N] {
    fn into_plot_data(self) -> Vec<f64> {
        self.to_vec()
    }
}

impl<const N: usize> IntoPlotData for [f32; N] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

impl<const N: usize> IntoPlotData for [i32; N] {
    fn into_plot_data(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

