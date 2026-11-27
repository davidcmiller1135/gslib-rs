// gslib_rs/src/variogram.rs

//! Empirical variogram calculation and simple model evaluation.

/// Compute the empirical semivariogram:
/// - coords: list of (x,y,z) locations
/// - values: matching variable values
/// - max_dist: maximum separation distance to consider
/// - n_lags: number of lag classes
///
/// Returns (lag_centers, semivariances).
pub fn compute_empirical(
    coords: &[(f64, f64, f64)],
    values: &[f64],
    max_dist: f64,
    n_lags: usize,
) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(
        coords.len(),
        values.len(),
        "coords and values length mismatch"
    );
    let n = coords.len();
    let hmax = max_dist;
    let bin_width = hmax / (n_lags as f64);

    // lag‐class centers
    let lag_centers: Vec<f64> = (0..n_lags).map(|i| bin_width * (i as f64 + 0.5)).collect();

    let mut sum = vec![0.0; n_lags];
    let mut count = vec![0_usize; n_lags];

    // accumulate semivariances
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i].0 - coords[j].0;
            let dy = coords[i].1 - coords[j].1;
            let dz = coords[i].2 - coords[j].2;
            let h = (dx * dx + dy * dy + dz * dz).sqrt();
            if h <= hmax {
                let mut bin = (h / bin_width).floor() as usize;
                if bin >= n_lags {
                    bin = n_lags - 1;
                }
                let semi_var = (values[i] - values[j]).powi(2) * 0.5;
                sum[bin] += semi_var;
                count[bin] += 1;
            }
        }
    }

    // compute averages, zero if no pairs
    let semivars: Vec<f64> = sum
        .into_iter()
        .zip(count)
        .map(|(s, c)| if c > 0 { s / (c as f64) } else { 0.0 })
        .collect();

    (lag_centers, semivars)
}

/// Spherical variogram model:
/// γ(h) = sill * [1.5*(h/a) - 0.5*(h/a)^3] for h < a, else sill.
pub fn spherical_model(lags: &[f64], range: f64, sill: f64) -> Vec<f64> {
    lags.iter()
        .map(|&h| {
            if h >= range {
                sill
            } else {
                let hr = h / range;
                sill * (1.5 * hr - 0.5 * hr.powi(3))
            }
        })
        .collect()
}

/// Exponential variogram model:
/// γ(h) = sill * [1 - exp(-3h/a)].
pub fn exponential_model(lags: &[f64], range: f64, sill: f64) -> Vec<f64> {
    lags.iter()
        .map(|&h| sill * (1.0 - (-3.0 * h / range).exp()))
        .collect()
}

/// Simple “fit” wrapper: takes sill = last semivar, range = last lag
/// and evaluates the chosen model (“spherical” or “exponential”).
pub fn fit_model(lags: &[f64], semivars: &[f64], model: &str) -> Vec<f64> {
    let sill = *semivars.last().unwrap_or(&0.0);
    let range = *lags.last().unwrap_or(&0.0);
    match model {
        "spherical" => spherical_model(lags, range, sill),
        "exponential" => exponential_model(lags, range, sill),
        _ => panic!("Unknown model type: {}", model),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    #[test]
    fn test_compute_empirical_two_points() {
        let pts = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let vals = vec![0.0, 1.0];
        let (lags, semis) = compute_empirical(&pts, &vals, 1.0, 1);
        assert_eq!(lags.len(), 1);
        // center = 0.5
        assert!((lags[0] - 0.5).abs() < EPS);
        // semivar = (1^2)/2 = 0.5
        assert!((semis[0] - 0.5).abs() < EPS);
    }

    #[test]
    fn test_spherical_model_behavior() {
        let lags = vec![0.0, 1.0, 2.0];
        let out = spherical_model(&lags, 2.0, 4.0);
        // h=0
        assert!((out[0] - 0.0).abs() < EPS);
        // h=1: 4*[1.5*(1/2) - 0.5*(1/2)^3] = 4*(0.75 - 0.0625) = 2.75
        assert!((out[1] - 2.75).abs() < EPS);
        // h>=range
        assert!((out[2] - 4.0).abs() < EPS);
    }

    #[test]
    fn test_exponential_model_behavior() {
        let lags = vec![0.0, 1.0];
        let out = exponential_model(&lags, 1.0, 2.0);
        // h=0 => 2*(1-exp(0)) = 0
        assert!((out[0] - 0.0).abs() < EPS);
        // h=1 => 2*(1-exp(-3))
        let expected = 2.0f64 * (1.0f64 - (-3.0f64).exp());
        assert!((out[1] - expected).abs() < EPS);
    }

    #[test]
    fn test_fit_model_dispatch() {
        let lags = vec![0.0, 2.0];
        let semis = vec![0.0, 5.0];
        let sph = fit_model(&lags, &semis, "spherical");
        assert_eq!(sph[1], 5.0);
        let exp = fit_model(&lags, &semis, "exponential");
        let exp_expected = 5.0f64 * (1.0f64 - (-3.0f64 * 2.0f64 / 2.0f64).exp());
        assert!((exp[1] - exp_expected).abs() < EPS);
    }
}
