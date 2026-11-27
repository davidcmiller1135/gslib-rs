// gslib_rs/src/kriging.rs

//! Kriging estimators: ordinary and simple kriging (cokriging stub).

use nalgebra::{DMatrix, DVector};
use crate::utils::{solve_kriging_system, solve_kriging_system_pivot};

/// Euclidean distance in 3D.
fn distance(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    let dz = a.2 - b.2;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Ordinary Kriging:
/// Builds the (n+1)x(n+1) semivariance matrix with Lagrange multiplier,
/// solves for weights, and returns the estimated value.
pub fn ordinary_kriging<F>(
    coords: &[(f64, f64, f64)],
    values: &[f64],
    target: (f64, f64, f64),
    variogram: F,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let n = coords.len();
    assert_eq!(values.len(), n, "coords/values length mismatch");

    // Build the (n+1)x(n+1) kriging matrix
    let mut mat = DMatrix::zeros(n + 1, n + 1);
    for i in 0..n {
        for j in 0..n {
            mat[(i, j)] = variogram(distance(coords[i], coords[j]));
        }
        mat[(i, n)] = 1.0;
        mat[(n, i)] = 1.0;
    }
    mat[(n, n)] = 0.0;

    // Build right‐hand side: variogram between each sample and the target
    let mut rhs = DVector::zeros(n + 1);
    for i in 0..n {
        rhs[i] = variogram(distance(coords[i], target));
    }
    rhs[n] = 1.0;

    // Solve with pivoting for stability
    let sol = solve_kriging_system_pivot(&mat, &rhs)
        .expect("Ordinary kriging system could not be solved");

    // Estimate = Σ w_i * value_i
    sol.iter()
        .take(n)
        .zip(values.iter())
        .map(|(w, &v)| w * v)
        .sum()
}

/// Simple Kriging with a known stationary `mean`.
///
/// `covariance` is the covariance function C(h), **not** a variogram: simple
/// kriging has no unbiasedness constraint, so its system must be written in
/// covariances. For a variogram model with sill `s`, pass `|h| s - gamma(h)`.
///
/// Solves C·w = c and returns mean + Σ w_i·(v_i − mean).
pub fn simple_kriging<F>(
    coords: &[(f64, f64, f64)],
    values: &[f64],
    target: (f64, f64, f64),
    covariance: F,
    mean: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let n = coords.len();
    assert_eq!(values.len(), n, "coords/values length mismatch");

    // Build the n×n data-to-data covariance matrix
    let mut mat = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            mat[(i, j)] = covariance(distance(coords[i], coords[j]));
        }
    }

    // Right-hand side: covariance between each sample and the target
    let mut rhs = DVector::zeros(n);
    for i in 0..n {
        rhs[i] = covariance(distance(coords[i], target));
    }

    // Solve (no pivoting)
    let weights = solve_kriging_system(&mat, &rhs)
        .expect("Simple kriging system could not be solved");

    // Estimate = mean + Σ w_i * (v_i - mean)
    weights
        .iter()
        .zip(values.iter())
        .fold(mean, |acc, (&w, &v)| acc + w * (v - mean))
}



/// Simple Cokriging (no sum‐to‐one constraints):
/// - `z_coords`, `z_vals`: primary-data locations & values
/// - `y_coords`, `y_vals`: secondary-data locations & values
/// - `cov_zz`: Z–Z covariance function CZZ(h)
/// - `cov_zy`: Z–Y cross‐covariance function CZY(h)
/// - `cov_yy`: Y–Y covariance function CYY(h)
/// - `m_z`, `m_y`: stationary means of Z and Y
/// 
/// Returns the cokriged estimate at `target`.
#[allow(clippy::too_many_arguments)]
pub fn simple_cokriging<Z, Y, H>(
    z_coords: &[(f64, f64, f64)],
    z_vals: &[f64],
    y_coords: &[(f64, f64, f64)],
    y_vals: &[f64],
    target: (f64, f64, f64),
    cov_zz: Z,
    cov_zy: Y,
    cov_yy: H,
    m_z: f64,
    m_y: f64,
) -> f64
where
    Z: Fn(f64) -> f64,
    Y: Fn(f64) -> f64,
    H: Fn(f64) -> f64,
{
    let n = z_coords.len();
    let m = y_coords.len();
    assert_eq!(z_vals.len(), n, "primary lengths mismatch");
    assert_eq!(y_vals.len(), m, "secondary lengths mismatch");

    let tot = n + m;
    // Build the (n+m)x(n+m) covariance matrix
    let mut mat = DMatrix::zeros(tot, tot);
    // Z–Z block
    for i in 0..n {
        for j in 0..n {
            let dx = z_coords[i].0 - z_coords[j].0;
            let dy = z_coords[i].1 - z_coords[j].1;
            let dz = z_coords[i].2 - z_coords[j].2;
            let h = (dx*dx + dy*dy + dz*dz).sqrt();
            mat[(i, j)] = cov_zz(h);
        }
    }
    // Z–Y and Y–Z blocks
    for i in 0..n {
        for j in 0..m {
            let dx = z_coords[i].0 - y_coords[j].0;
            let dy = z_coords[i].1 - y_coords[j].1;
            let dz = z_coords[i].2 - y_coords[j].2;
            let h = (dx*dx + dy*dy + dz*dz).sqrt();
            let czy = cov_zy(h);
            mat[(i, n + j)] = czy;
            mat[(n + j, i)] = czy;
        }
    }
    // Y–Y block
    for i in 0..m {
        for j in 0..m {
            let dx = y_coords[i].0 - y_coords[j].0;
            let dy = y_coords[i].1 - y_coords[j].1;
            let dz = y_coords[i].2 - y_coords[j].2;
            let h = (dx*dx + dy*dy + dz*dz).sqrt();
            mat[(n + i, n + j)] = cov_yy(h);
        }
    }

    // Build RHS vector [CZZ(h_i,u); CZY(h_j,u)]
    let mut rhs = DVector::zeros(tot);
    for i in 0..n {
        let dx = z_coords[i].0 - target.0;
        let dy = z_coords[i].1 - target.1;
        let dz = z_coords[i].2 - target.2;
        rhs[i] = cov_zz((dx*dx + dy*dy + dz*dz).sqrt());
    }
    for j in 0..m {
        let dx = y_coords[j].0 - target.0;
        let dy = y_coords[j].1 - target.1;
        let dz = y_coords[j].2 - target.2;
        rhs[n + j] = cov_zy((dx*dx + dy*dy + dz*dz).sqrt());
    }


    // Solve for weights (unconstrained), but if singular, return the primary mean:
    let w = match solve_kriging_system(&mat, &rhs) {
        Some(weights) => weights,
        None => return m_z,
    };



    // Build estimate: m_z + Σ w_i (Z_i - m_z) + Σ w_{n+j} (Y_j - m_y)
    let mut est = m_z;
    for i in 0..n {
        est += w[i] * (z_vals[i] - m_z);
    }
    for j in 0..m {
        est += w[n + j] * (y_vals[j] - m_y);
    }
    est
}

#[cfg(test)]
mod cokrig_tests {
    use super::*;

    /// Trivial check: if all data equal their means,
    /// estimate should equal m_z.
    #[test]
    fn test_simple_cokriging_constant() {
        let zc = vec![(0.0,0.0,0.0)];
        let yc = vec![(1.0,1.0,1.0)];
        let zv = vec![10.0];
        let yv = vec![20.0];
        let mz = 10.0;
        let my = 20.0;
        // constant covariances for any distance
        let cz = |_| 1.0;
        let cy = |_| 1.0;
        let yy = |_| 1.0;

        let est = simple_cokriging(
            &zc, &zv,
            &yc, &yv,
            (2.0,2.0,2.0),
            cz, cy, yy,
            mz, my,
        );
        assert!((est - mz).abs() < 1e-12);
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    /// If all known values are equal, OK returns that constant.
    #[test]
    fn test_ok_constant_values() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)];
        let values = vec![5.0, 5.0, 5.0];
        let est = ordinary_kriging(&coords, &values, (0.2, 0.3, 0.0), |h| h);
        assert!((est - 5.0).abs() < EPS);
    }

    /// Two‐point OK with variogram = distance should give the midpoint average.
    #[test]
    fn test_ok_two_points_identity_variogram() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let values = vec![1.0, 2.0];
        // target at 0.5 => equal weights
        let est = ordinary_kriging(&coords, &values, (0.5, 0.0, 0.0), |h| h);
        assert!((est - 1.5).abs() < EPS);
    }

    /// Exponential covariance with unit sill and practical range 3.
    fn exp_cov(h: f64) -> f64 {
        (-h).exp()
    }

    /// Simple kriging with both values equal to the mean returns the mean.
    #[test]
    fn test_simple_kriging_equal_to_mean() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let values = vec![7.0, 7.0];
        let mean = 7.0;
        let est = simple_kriging(&coords, &values, (0.5, 0.0, 0.0), exp_cov, mean);
        assert!((est - mean).abs() < EPS);
    }

    /// SK is an exact interpolator: at a data location it returns the datum.
    #[test]
    fn test_simple_kriging_exact_at_data() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let values = vec![12.0, 14.0];
        let est = simple_kriging(&coords, &values, (1.0, 0.0, 0.0), exp_cov, 10.0);
        assert!((est - 14.0).abs() < 1e-10, "est = {est}");
    }

    /// Far beyond the correlation range, SK reverts to the known mean.
    #[test]
    fn test_simple_kriging_reverts_to_mean_far_away() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let values = vec![12.0, 14.0];
        let est = simple_kriging(&coords, &values, (100.0, 0.0, 0.0), exp_cov, 10.0);
        assert!((est - 10.0).abs() < 1e-10, "est = {est}");
    }

    /// One datum: weight is C(h)/C(0), so est = m + C(h)·(v − m) for unit sill.
    #[test]
    fn test_simple_kriging_single_datum_weight() {
        let est = simple_kriging(&[(0.0, 0.0, 0.0)], &[12.0], (1.0, 0.0, 0.0), exp_cov, 10.0);
        let expected = 10.0 + (-1.0f64).exp() * 2.0;
        assert!((est - expected).abs() < EPS, "est = {est}, expected = {expected}");
    }
}
