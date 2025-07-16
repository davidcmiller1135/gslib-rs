//! Empirical variogram calculation and model fitting

pub fn compute_empirical(
    coords: &[(f64, f64, f64)],
    values: &[f64],
    max_dist: f64,
    n_lags: usize,
) -> (Vec<f64>, Vec<f64>) {
    // returns (lag distances, semivariance)
    unimplemented!()
}

pub fn fit_model(lags: &[f64], semivars: &[f64], model: &str) -> Vec<f64> {
    unimplemented!()
}
