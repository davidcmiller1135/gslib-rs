// gslib_rs/src/simulation.rs

//! Geostatistical simulation: Sequential Gaussian Simulation (SGS).

use nalgebra::{DMatrix, DVector};
use crate::utils::{solve_kriging_system, rand_f, randn};
use std::cmp::Ordering;

use crate::variogram::compute_empirical;
use std::collections::HashSet;
use std::f64;



/// Euclidean distance in 3D.
fn distance(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    let dz = a.2 - b.2;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Perform **Sequential Gaussian Simulation** on `sim_coords`, conditioned on
/// `data_coords`/`data_vals` (which must already be normal‐score values),
/// using covariance model `cov_func` and known `mean`.  
///
/// Returns a Vec of length `n_real`, each entry a Vec<f64> of simulated values
/// (normal scores) at each `sim_coords` node.
///
/// Assumes SK means/variances at each step (no neighborhood truncation).
pub fn sgsim<F>(
    data_coords: &[(f64, f64, f64)],
    data_vals: &[f64],
    sim_coords: &[(f64, f64, f64)],
    cov_func: F,
    mean: f64,
    n_real: usize,
) -> Vec<Vec<f64>>
where
    F: Fn(f64) -> f64 + Copy,
{
    let n_data = data_coords.len();
    let n_sim = sim_coords.len();
    assert_eq!(
        data_vals.len(),
        n_data,
        "data_coords/data_vals length mismatch"
    );

    let mut realizations = Vec::with_capacity(n_real);

    for _ in 0..n_real {
        // 1) Initialize conditioning set
        let mut cond_coords = data_coords.to_vec();
        let mut cond_vals = data_vals.to_vec();

        // 2) Generate a random path over the sim nodes
        let mut order: Vec<(f64, usize)> = (0..n_sim)
            .map(|i| (rand_f(), i))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
        let path: Vec<usize> = order.into_iter().map(|(_, idx)| idx).collect();

        // 3) Simulate along the path
        let mut sim_vals = vec![0.0; n_sim];
        for &i in &path {
            let target = sim_coords[i];
            let m = cond_coords.len();

            // Build SK system A·w = c
            let mut a = DMatrix::zeros(m, m);
            for ii in 0..m {
                for jj in 0..m {
                    a[(ii, jj)] = cov_func(distance(cond_coords[ii], cond_coords[jj]));
                }
            }
            let mut c = DVector::zeros(m);
            for ii in 0..m {
                c[ii] = cov_func(distance(cond_coords[ii], target));
            }

            // Solve for weights; if singular, fall back to mean (zero variance)
            let (sk_mean, sk_var) = match solve_kriging_system(&a, &c) {
                Some(w) => {
                    // SK estimate = mean + Σ w_i (v_i - mean)
                    let est = w
                        .iter()
                        .zip(cond_vals.iter())
                        .fold(mean, |acc, (&wi, &vi)| acc + wi * (vi - mean));
                    // SK variance = C(0) - Σ w_i * c_i
                    let mut var = cov_func(0.0);
                    for (wi, ci) in w.iter().zip(c.iter()) {
                        var -= wi * ci;
                    }
                    (est, var.max(0.0))
                }
                None => (mean, 0.0),
            };

            // Draw a Gaussian deviate and form the simulated normal score
            let sim = if sk_var > 0.0 {
                sk_mean + randn() * sk_var.sqrt()
            } else {
                sk_mean
            };
            sim_vals[i] = sim;

            // Append simulated node to conditioning set
            cond_coords.push(target);
            cond_vals.push(sim);
        }

        realizations.push(sim_vals);
    }

    realizations
}

#[cfg(test)]
mod tests {
    use super::*;

    const MEAN: f64 = 5.0;

    /// If covariance is constant and data equal the mean, each sim node
    /// should also equal the mean (zero variance).
    #[test]
    fn test_sgsim_constant_covariance() {
        let data_coords = vec![(0.0, 0.0, 0.0)];
        let data_vals = vec![MEAN];
        let sim_coords = vec![(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)];
        let cov_const = |_| 1.0;
        let realizations = sgsim(&data_coords, &data_vals, &sim_coords, cov_const, MEAN, 3);
        assert_eq!(realizations.len(), 3);
        for sim in realizations {
            assert_eq!(sim, vec![MEAN, MEAN]);
        }
    }
}




/// Sequential Indicator Simulation (single threshold).
///
/// - `data_coords` / `data_vals`: conditioning (real) data.
/// - `sim_coords`             : locations to simulate.
/// - `threshold`              : the cutoff (values < threshold ⇒ indicator=1).
/// - `cov_func`               : covariance model C(h).
/// - `n_real`                 : number of realizations.
///
/// Returns a Vec of length `n_real`, each being a Vec<u8> of 0/1 indicators.
pub fn sisim<F>(
    data_coords: &[(f64, f64, f64)],
    data_vals: &[f64],
    sim_coords: &[(f64, f64, f64)],
    threshold: f64,
    cov_func: F,
    n_real: usize,
) -> Vec<Vec<u8>>
where
    F: Fn(f64) -> f64 + Copy,
{
    let n_data = data_coords.len();
    let n_sim = sim_coords.len();
    assert_eq!(data_vals.len(), n_data);

    // Precompute the indicator coding of the conditioning data
    let init_inds: Vec<f64> = data_vals
        .iter()
        .map(|&v| if v < threshold { 1.0 } else { 0.0 })
        .collect();

    let mut all_reals = Vec::with_capacity(n_real);

    for _ in 0..n_real {
        // 1) start with the original data
        let mut cond_coords = data_coords.to_vec();
        let mut cond_inds = init_inds.clone();

        // 2) random path
        let mut order: Vec<(f64, usize)> = (0..n_sim)
            .map(|i| (rand_f(), i))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let path: Vec<usize> = order.into_iter().map(|(_, i)| i).collect();

        // 3) simulate
        let mut sim_inds = vec![0u8; n_sim];
        for &i in &path {
            let target = sim_coords[i];
            let m = cond_coords.len();

            // build the SK system A·w = c (no mean here, indicators are zero‐mean residuals)
            let mut a = nalgebra::DMatrix::zeros(m, m);
            for ii in 0..m {
                for jj in 0..m {
                    a[(ii, jj)] = cov_func(distance(cond_coords[ii], cond_coords[jj]));
                }
            }
            let mut c = nalgebra::DVector::zeros(m);
            for ii in 0..m {
                c[ii] = cov_func(distance(cond_coords[ii], target));
            }

            // solve for weights; if singular, use marginal probability = mean of cond_inds
            let p = match solve_kriging_system(&a, &c) {
                Some(w) => {
                    let mut prob = 0.0;
                    for (wi, &ind) in w.iter().zip(cond_inds.iter()) {
                        prob += wi * ind;
                    }
                    // clamp into [0,1]
                    prob.max(0.0).min(1.0)
                }
                None => {
                    let sum: f64 = cond_inds.iter().sum();
                    sum / (m as f64)
                }
            };

            // draw the indicator: 1 with probability p
            let draw = if rand_f() < p { 1u8 } else { 0u8 };
            sim_inds[i] = draw;

            // append to conditioning set
            cond_coords.push(target);
            cond_inds.push(draw as f64);
        }

        all_reals.push(sim_inds);
    }

    all_reals
}

#[cfg(test)]
mod sisim_tests {
    use super::*;

    /// If all data < threshold and covariance is constant,
    /// indicators should all be 1.
    #[test]
    fn test_sisim_all_ones() {
        let data_coords = vec![(0.0,0.0,0.0)];
        let data_vals   = vec![0.5];
        let sim_coords  = vec![(1.0,0.0,0.0),(2.0,0.0,0.0)];
        let cov_const = |_| 1.0;
        let reps = sisim(&data_coords, &data_vals, &sim_coords, 1.0, cov_const, 5);
        assert_eq!(reps.len(), 5);
        for sim in reps {
            assert_eq!(sim, vec![1,1]);
        }
    }

    /// If all data ≥ threshold, indicators should all be 0.
    #[test]
    fn test_sisim_all_zeros() {
        let data_coords = vec![(0.0,0.0,0.0)];
        let data_vals   = vec![2.0];
        let sim_coords  = vec![(1.0,0.0,0.0),(2.0,0.0,0.0)];
        let cov_const = |_| 1.0;
        let reps = sisim(&data_coords, &data_vals, &sim_coords, 1.0, cov_const, 3);
        assert_eq!(reps.len(), 3);
        for sim in reps {
            assert_eq!(sim, vec![0,0]);
        }
    }
}

// gslib_rs/src/simulation.rs

// Geostatistical simulation: SASIM (simulated-annealing)


// ----------------------------------
// SASIM PARAMS
// ----------------------------------

/// Parameters controlling a SASIM run.
#[derive(Debug, Clone)]
pub struct SasimParams {
    /// Which objective components to include: [hist, vario, ivario, corr, cpdf]
    pub obj_components: [bool; 5],
    /// Weights for each component (same ordering).
    pub weights: [f64; 5],
    /// Histogram bins, and optional explicit (min,max).  If `hist_range` is None,
    /// range is taken from the primary data.
    pub hist_bins: usize,
    pub hist_range: Option<(f64, f64)>,

    /// Variogram settings: max distance and number of lags.
    pub vario_max_dist: f64,
    pub vario_n_lags: usize,

    /// Threshold for indicator semivariogram.
    pub threshold: f64,

    /// Secondary variables: each entry is a vector of simulated‐grid values,
    /// aligned one‐to‐one with the simulation grid.
    pub secondary_vals: Vec<Vec<f64>>,

    /// If `true`, do not perturb conditioning data nodes.
    pub conditioning: bool,

    /// Sim grid definition: (nx, x0, dx, ny, y0, dy, nz, z0, dz)
    pub grid: (usize, f64, f64, usize, f64, f64, usize, f64, f64),

    /// Number of realizations to generate.
    pub nsim: usize,

    /// Random seed.
    pub seed: u64,

    /// Annealing schedule: (T0, reduction_factor, max_accepts_per_level,
    /// max_trials_per_level, n_levels, T_min)
    pub schedule: (f64, f64, usize, usize, usize, f64),
}

// ----------------------------------
// SASIM IMPLEMENTATION
// ----------------------------------

/// Main entrypoint for SASIM.
pub fn sasim(
    data_coords: &[(f64, f64, f64)],
    data_vals: &[f64],
    params: &SasimParams,
) -> Vec<Vec<f64>> {
    // 0) Seed the RNG
    crate::utils::set_seed(params.seed);

    // 1) Build the simulation grid coords
    let sim_coords = generate_sim_coords(params.grid);
    let n_sim = sim_coords.len();

    // 2) Identify which grid‐indices are conditioning data
    let mut cond_idx = HashSet::new();
    for &dc in data_coords.iter() {
        if let Some(i) = sim_coords.iter().position(|&c| c == dc) {
            cond_idx.insert(i);
        }
    }

    // 3) Precompute targets for each active component

    // 3a) Histogram target (as proportions)
    let mut target_hist = Vec::new();
    if params.obj_components[0] {
        let (hmin, hmax) = params
            .hist_range
            .unwrap_or_else(|| {
                let mn = data_vals
                    .iter()
                    .cloned()
                    .fold(f64::INFINITY, f64::min);
                let mx = data_vals
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max);
                (mn, mx)
            });
        let hist = compute_histogram(data_vals, params.hist_bins, hmin, hmax);
        target_hist = hist
            .into_iter()
            .map(|c| c as f64 / data_vals.len() as f64)
            .collect();
    }

    // 3b) Variogram target
    let mut target_vario = Vec::new();
    if params.obj_components[1] {
        let (_lags, semis) = compute_empirical(
            data_coords,
            data_vals,
            params.vario_max_dist,
            params.vario_n_lags,
        );
        target_vario = semis;
    }

    // 3c) Indicator‐variogram target
    let mut target_ivario = Vec::new();
    if params.obj_components[2] {
        let inds: Vec<f64> = data_vals
            .iter()
            .map(|&v| if v < params.threshold { 1.0 } else { 0.0 })
            .collect();
        let (_lags_i, semis_i) = compute_empirical(
            data_coords,
            &inds,
            params.vario_max_dist,
            params.vario_n_lags,
        );
        target_ivario = semis_i;
    }

    // 3d) Correlation targets, one per secondary variable
    let mut target_corr = Vec::new();
    if params.obj_components[3] {
        // build a sub‐vector of sec‐values at the conditioning indices
        for sec in &params.secondary_vals {
            let sec_sub: Vec<f64> = cond_idx.iter().map(|&i| sec[i]).collect();
            let corr = compute_correlation(data_vals, &sec_sub);
            target_corr.push(corr);
        }
    }

    // 3e) Conditional‐PDF target would go here…

    // 4) Now run `nsim` realizations
    let mut reals = Vec::with_capacity(params.nsim);
    for _ in 0..params.nsim {
        // 4a) Initialize the field by sampling data_vals at random
        let mut field: Vec<f64> = (0..n_sim)
            .map(|_| {
                let idx = ((rand_f() * data_vals.len() as f64).floor() as usize)
                    .min(data_vals.len() - 1);
                data_vals[idx]
            })
            .collect();

        // 4b) Compute its initial “energy”
        let mut current_e = compute_total_energy(
            &field,
            &sim_coords,
            params,
            &target_hist,
            &target_vario,
            &target_ivario,
            &target_corr,
        );

        // 4c) Annealing
        let (mut T, redfac, ka, kmax, n_levels, Tmin) = params.schedule;
        'levels: for _ in 0..n_levels {
            let mut accepted = 0;
            let mut trials = 0;
            while accepted < ka && trials < kmax {
                trials += 1;
                // pick a random node
                let i = ((rand_f() * n_sim as f64).floor() as usize).min(n_sim - 1);
                if params.conditioning && cond_idx.contains(&i) {
                    continue;
                }
                // propose a new value
                let old = field[i];
                let j = ((rand_f() * data_vals.len() as f64).floor() as usize)
                    .min(data_vals.len() - 1);
                let new = data_vals[j];
                field[i] = new;

                // compute new energy
                let new_e = compute_total_energy(
                    &field,
                    &sim_coords,
                    params,
                    &target_hist,
                    &target_vario,
                    &target_ivario,
                    &target_corr,
                );
                let dE = new_e - current_e;
                if dE <= 0.0 || rand_f() < (-dE / T).exp() {
                    current_e = new_e;
                    accepted += 1;
                } else {
                    field[i] = old;
                }
            }
            T *= redfac;
            if T < Tmin {
                break 'levels;
            }
        }

        reals.push(field);
    }

    reals
}

// ----------------------------------
// HELPER FUNCTIONS
// ----------------------------------

/// Generate the full list of grid coordinates from params.grid.
fn generate_sim_coords(
    grid: (usize, f64, f64, usize, f64, f64, usize, f64, f64),
) -> Vec<(f64, f64, f64)> {
    let (nx, x0, dx, ny, y0, dy, nz, z0, dz) = grid;
    let mut coords = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                coords.push((x0 + ix as f64 * dx, y0 + iy as f64 * dy, z0 + iz as f64 * dz));
            }
        }
    }
    coords
}

/// Build a simple histogram (counts per bin) over `[min,max]`, clamping edges.
fn compute_histogram(data: &[f64], bins: usize, min: f64, max: f64) -> Vec<usize> {
    let mut hist = vec![0; bins];
    let range = max - min;
    for &v in data {
        let idx = if v <= min {
            0
        } else if v >= max {
            bins - 1
        } else {
            (((v - min) / range) * bins as f64).floor() as usize
        };
        hist[idx] += 1;
    }
    hist
}

/// Compute Pearson correlation between two slices (must have same len).
fn compute_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    num / (var_x.sqrt() * var_y.sqrt())
}

/// Compute the total weighted energy (objective) given the active targets.
fn compute_total_energy(
    field: &[f64],
    sim_coords: &[(f64, f64, f64)],
    params: &SasimParams,
    target_hist: &[f64],
    target_vario: &[f64],
    target_ivario: &[f64],
    target_corr: &[f64],
) -> f64 {
    let mut E = 0.0;

    // 1) Histogram term
    if params.obj_components[0] {
        let (hmin, hmax) = params.hist_range.unwrap();
        let hist = compute_histogram(field, params.hist_bins, hmin, hmax);
        let prop: Vec<f64> = hist
            .into_iter()
            .map(|c| c as f64 / field.len() as f64)
            .collect();
        let mut sum = 0.0;
        for (i, &t) in target_hist.iter().enumerate() {
            let d = prop[i] - t;
            sum += d * d;
        }
        E += params.weights[0] * sum;
    }

    // 2) Semivariogram term
    if params.obj_components[1] {
        let (_, sim_semis) = compute_empirical(
            sim_coords,
            field,
            params.vario_max_dist,
            params.vario_n_lags,
        );
        let mut sum = 0.0;
        for (i, &t) in target_vario.iter().enumerate() {
            let d = sim_semis[i] - t;
            sum += d * d;
        }
        E += params.weights[1] * sum;
    }

    // 3) Indicator‐semivariogram term
    if params.obj_components[2] {
        let inds: Vec<f64> =
            field.iter().map(|&v| if v < params.threshold { 1.0 } else { 0.0 }).collect();
        let (_, sim_ii) = compute_empirical(
            sim_coords,
            &inds,
            params.vario_max_dist,
            params.vario_n_lags,
        );
        let mut sum = 0.0;
        for (i, &t) in target_ivario.iter().enumerate() {
            let d = sim_ii[i] - t;
            sum += d * d;
        }
        E += params.weights[2] * sum;
    }

    // 4) Correlation term (one or more secondaries)
    if params.obj_components[3] {
        for (j, sec) in params.secondary_vals.iter().enumerate() {
            let r_sim = compute_correlation(field, sec);
            let d = r_sim - target_corr[j];
            E += params.weights[3] * d * d;
        }
    }

    // 5) Conditional‐PDF term
    if params.obj_components[4] {
        // You can drop in your own conditional‐PDF distance here.
        // For now we leave it at zero.
    }

    E
}

// ----------------------------------
// SASIM UNIT TEST (stub)
// ----------------------------------

#[cfg(test)]
mod sasim_tests {
    use super::*;

    #[test]
    fn test_sasim_constant_histogram() {
        // If data_vals are constant, histogram objective alone
        // should produce exactly that constant at every node.
        let data_coords = vec![(0.0, 0.0, 0.0)];
        let data_vals = vec![3.1415];
        let params = SasimParams {
            obj_components: [true, false, false, false, false],
            weights: [1.0, 0.0, 0.0, 0.0, 0.0],
            hist_bins: 1,
            hist_range: Some((3.1415, 3.1415)),
            vario_max_dist: 0.0,
            vario_n_lags: 1,
            threshold: 0.0,
            secondary_vals: vec![],
            conditioning: false,
            grid: (2, 0.0, 1.0, 1, 0.0, 1.0, 1, 0.0, 1.0),
            nsim: 2,
            seed: 42,
            schedule: (1.0, 0.5, 10, 20, 5, 1e-3),
        };

        let reals = sasim(&data_coords, &data_vals, &params);
        assert_eq!(reals.len(), 2);
        for field in reals {
            assert!(field.iter().all(|&v| (v - 3.1415).abs() < 1e-12));
        }
    }
}
