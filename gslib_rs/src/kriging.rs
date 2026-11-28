//! Kriging estimators: simple, ordinary, and simple cokriging.
//!
//! Each function estimates one target location from the samples it is given;
//! choosing a search neighbourhood is the caller's job.

use nalgebra::{DMatrix, DVector};

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::linalg::{solve_general, solve_spd, symmetric_matrix};
use crate::model::{CovarianceModel, VariogramModel};
use crate::sample::Sample;

/// A kriging estimate and its kriging (error) variance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Estimate {
    /// The estimated value.
    pub value: f64,
    /// The kriging variance, clamped to be non-negative.
    pub variance: f64,
}

/// A kriging solution: the estimate and the weights that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// The estimate and its variance.
    pub estimate: Estimate,
    /// One weight per sample, in the order given.
    pub weights: Vec<f64>,
    /// Ordinary kriging's Lagrange multiplier μ (in variogram form; GSLIB's
    /// covariance form reports −μ scaled by its `unbias` constant), or
    /// `None` for simple kriging.
    pub lagrange: Option<f64>,
}

fn dot(a: &DVector<f64>, b: impl IntoIterator<Item = f64>) -> f64 {
    a.iter().zip(b).map(|(w, v)| w * v).sum()
}

/// Simple kriging around a known stationary `mean`.
///
/// Solves C·w = c and returns `mean + Σ wᵢ(zᵢ − mean)` with variance
/// `C(0) − wᵀc`. With no samples, returns the mean and the sill.
///
/// # Errors
/// [`Error::SingularSystem`] if the covariance matrix is not positive definite
/// (typically two samples at the same location).
pub fn simple_kriging(
    samples: &[Sample],
    target: Point3,
    model: &CovarianceModel,
    mean: f64,
) -> Result<Estimate> {
    simple_kriging_solution(samples, target, model, mean).map(|s| s.estimate)
}

/// [`simple_kriging`], also returning the weights.
///
/// # Errors
/// As [`simple_kriging`].
pub fn simple_kriging_solution(
    samples: &[Sample],
    target: Point3,
    model: &CovarianceModel,
    mean: f64,
) -> Result<Solution> {
    let a = symmetric_matrix(samples, |p, q| model.covariance(p.location(), q.location()));
    let c = DVector::from_iterator(
        samples.len(),
        samples
            .iter()
            .map(|s| model.covariance(s.location(), target)),
    );
    let w = solve_spd(a, &c, "simple")?;
    Ok(Solution {
        estimate: Estimate {
            value: mean + dot(&w, samples.iter().map(|s| s.value() - mean)),
            variance: (model.sill() - w.dot(&c)).max(0.0),
        },
        weights: w.iter().copied().collect(),
        lagrange: None,
    })
}

/// Ordinary kriging: unknown constant mean, weights constrained to sum to one.
///
/// Written in variograms, so it also works for models without a sill (power
/// model). Solves `Σⱼ wⱼ γᵢⱼ + μ = γᵢ₀`, `Σ wᵢ = 1`; variance is
/// `Σ wᵢ γᵢ₀ + μ`.
///
/// # Errors
/// [`Error::InvalidInput`] if `samples` is empty (the mean is unknown), or
/// [`Error::SingularSystem`] if the system cannot be solved.
pub fn ordinary_kriging(
    samples: &[Sample],
    target: Point3,
    model: &VariogramModel,
) -> Result<Estimate> {
    ordinary_kriging_solution(samples, target, model).map(|s| s.estimate)
}

/// [`ordinary_kriging`], also returning the weights and μ.
///
/// # Errors
/// As [`ordinary_kriging`].
pub fn ordinary_kriging_solution(
    samples: &[Sample],
    target: Point3,
    model: &VariogramModel,
) -> Result<Solution> {
    if samples.is_empty() {
        return Err(Error::invalid(
            "samples",
            "ordinary kriging needs at least one sample",
        ));
    }
    let n = samples.len();
    // [Γ 1; 1ᵀ 0]: the variogram matrix bordered by the unbiasedness constraint.
    let a = DMatrix::from_fn(n + 1, n + 1, |i, j| {
        match (samples.get(i), samples.get(j)) {
            (Some(p), Some(q)) => model.gamma(p.location(), q.location()),
            (None, None) => 0.0,
            (Some(_), None) | (None, Some(_)) => 1.0,
        }
    });
    let rhs_gamma: Vec<f64> = samples
        .iter()
        .map(|s| model.gamma(s.location(), target))
        .collect();
    let b = DVector::from_iterator(n + 1, rhs_gamma.iter().copied().chain(std::iter::once(1.0)));
    let x = solve_general(a, &b, "ordinary")?;
    let weights = x.rows(0, n).into_owned();
    let mu: f64 = x.iter().skip(n).sum();
    Ok(Solution {
        estimate: Estimate {
            value: dot(&weights, samples.iter().map(Sample::value)),
            variance: (dot(&weights, rhs_gamma) + mu).max(0.0),
        },
        weights: weights.iter().copied().collect(),
        lagrange: Some(mu),
    })
}

/// The three models simple cokriging needs.
#[derive(Debug, Clone, PartialEq)]
pub struct CokrigingModel {
    /// Primary–primary covariance.
    pub primary: CovarianceModel,
    /// Primary–secondary cross-covariance.
    pub cross: CovarianceModel,
    /// Secondary–secondary covariance.
    pub secondary: CovarianceModel,
}

/// Stationary means of the primary and secondary variables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Means {
    /// Mean of the primary variable.
    pub primary: f64,
    /// Mean of the secondary variable.
    pub secondary: f64,
}

/// Simple cokriging of the primary variable from primary and secondary
/// samples (no unbiasedness constraints).
///
/// # Errors
/// [`Error::SingularSystem`] if the joint covariance matrix is not positive
/// definite, e.g. the cross-covariance is too strong for the direct models.
pub fn simple_cokriging(
    primary: &[Sample],
    secondary: &[Sample],
    target: Point3,
    model: &CokrigingModel,
    means: Means,
) -> Result<Estimate> {
    #[derive(Clone, Copy)]
    enum Var {
        Primary,
        Secondary,
    }
    let all: Vec<(Var, Sample)> = primary
        .iter()
        .map(|&s| (Var::Primary, s))
        .chain(secondary.iter().map(|&s| (Var::Secondary, s)))
        .collect();
    let cov = |a: Var, b: Var| match (a, b) {
        (Var::Primary, Var::Primary) => &model.primary,
        (Var::Secondary, Var::Secondary) => &model.secondary,
        (Var::Primary, Var::Secondary) | (Var::Secondary, Var::Primary) => &model.cross,
    };
    let a = symmetric_matrix(&all, |&(va, sa), &(vb, sb)| {
        cov(va, vb).covariance(sa.location(), sb.location())
    });
    let c = DVector::from_iterator(
        all.len(),
        all.iter()
            .map(|&(v, s)| cov(Var::Primary, v).covariance(s.location(), target)),
    );
    let w = solve_spd(a, &c, "cokriging")?;
    let residuals = all.iter().map(|&(v, s)| match v {
        Var::Primary => s.value() - means.primary,
        Var::Secondary => s.value() - means.secondary,
    });
    Ok(Estimate {
        value: means.primary + dot(&w, residuals),
        variance: (model.primary.sill() - w.dot(&c)).max(0.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Kind, Structure};
    use approx::assert_abs_diff_eq;

    const EPS: f64 = 1e-10;

    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }
    fn at(x: f64) -> Point3 {
        Point3::new(x, 0.0, 0.0)
    }
    fn exp_model(sill: f64, range: f64) -> VariogramModel {
        VariogramModel::new(
            0.0,
            vec![Structure::isotropic(Kind::Exponential, sill, range).unwrap()],
        )
        .unwrap()
    }
    fn exp_cov(sill: f64, range: f64) -> CovarianceModel {
        CovarianceModel::try_from(exp_model(sill, range)).unwrap()
    }
    fn zero_cov() -> CovarianceModel {
        CovarianceModel::try_from(VariogramModel::cross(0.0, vec![]).unwrap()).unwrap()
    }

    #[test]
    fn sk_is_exact_at_data() {
        let e = simple_kriging(
            &[s(0.0, 12.0), s(1.0, 14.0)],
            at(1.0),
            &exp_cov(1.0, 3.0),
            10.0,
        )
        .unwrap();
        assert_abs_diff_eq!(e.value, 14.0, epsilon = EPS);
        assert_abs_diff_eq!(e.variance, 0.0, epsilon = EPS);
    }

    #[test]
    fn sk_reverts_to_mean_and_sill_far_away() {
        let e = simple_kriging(
            &[s(0.0, 12.0), s(1.0, 14.0)],
            at(1000.0),
            &exp_cov(2.0, 3.0),
            10.0,
        )
        .unwrap();
        assert_abs_diff_eq!(e.value, 10.0, epsilon = EPS);
        assert_abs_diff_eq!(e.variance, 2.0, epsilon = EPS);
    }

    #[test]
    fn sk_single_datum_closed_form() {
        // w = C(h)/C(0); variance = C(0) − C(h)²/C(0).
        let model = exp_cov(1.0, 3.0);
        let ch = model.covariance(at(0.0), at(1.0));
        let e = simple_kriging(&[s(0.0, 12.0)], at(1.0), &model, 10.0).unwrap();
        assert_abs_diff_eq!(e.value, 10.0 + ch * 2.0, epsilon = EPS);
        assert_abs_diff_eq!(e.variance, 1.0 - ch * ch, epsilon = EPS);
    }

    #[test]
    fn sk_with_no_data_is_the_prior() {
        let e = simple_kriging(&[], at(0.0), &exp_cov(3.0, 1.0), 5.0).unwrap();
        assert_abs_diff_eq!(e.value, 5.0);
        assert_abs_diff_eq!(e.variance, 3.0);
    }

    #[test]
    fn duplicate_locations_are_an_error_not_a_fallback() {
        let r = simple_kriging(
            &[s(0.0, 1.0), s(0.0, 2.0)],
            at(0.5),
            &exp_cov(1.0, 3.0),
            0.0,
        );
        assert!(matches!(r, Err(Error::SingularSystem { .. })), "{r:?}");
    }

    #[test]
    fn ok_constant_data() {
        let samples = [s(0.0, 5.0), s(1.0, 5.0), s(3.0, 5.0)];
        let e = ordinary_kriging(&samples, at(0.4), &exp_model(1.0, 3.0)).unwrap();
        assert_abs_diff_eq!(e.value, 5.0, epsilon = EPS);
    }

    #[test]
    fn ok_midpoint_gets_equal_weights() {
        let e =
            ordinary_kriging(&[s(0.0, 1.0), s(1.0, 2.0)], at(0.5), &exp_model(1.0, 3.0)).unwrap();
        assert_abs_diff_eq!(e.value, 1.5, epsilon = EPS);
    }

    #[test]
    fn ok_single_datum_variance_is_twice_the_sill_far_away() {
        // Var(Z0 − Z1) = 2σ² when Z0 and Z1 are uncorrelated.
        let e = ordinary_kriging(&[s(0.0, 7.0)], at(1e6), &exp_model(1.5, 3.0)).unwrap();
        assert_abs_diff_eq!(e.value, 7.0, epsilon = EPS);
        assert_abs_diff_eq!(e.variance, 3.0, epsilon = 1e-9);
    }

    #[test]
    fn ok_works_with_power_model() {
        let m = VariogramModel::new(
            0.0,
            vec![Structure::isotropic(Kind::Power { exponent: 1.0 }, 1.0, 1.0).unwrap()],
        )
        .unwrap();
        let e = ordinary_kriging(&[s(0.0, 1.0), s(2.0, 3.0)], at(1.0), &m).unwrap();
        assert_abs_diff_eq!(e.value, 2.0, epsilon = EPS);
    }

    #[test]
    fn ok_needs_data() {
        assert!(
            ordinary_kriging(&[], at(0.0), &exp_model(1.0, 1.0)).is_err(),
            "no data"
        );
    }

    #[test]
    fn cokriging_without_secondary_equals_sk() {
        let primary = [s(0.0, 12.0), s(2.0, 9.0)];
        let model = CokrigingModel {
            primary: exp_cov(1.0, 3.0),
            cross: zero_cov(),
            secondary: exp_cov(1.0, 3.0),
        };
        let means = Means {
            primary: 10.0,
            secondary: 0.0,
        };
        let ck = simple_cokriging(&primary, &[], at(1.0), &model, means).unwrap();
        let sk = simple_kriging(&primary, at(1.0), &model.primary, 10.0).unwrap();
        assert_abs_diff_eq!(ck.value, sk.value, epsilon = EPS);
        assert_abs_diff_eq!(ck.variance, sk.variance, epsilon = EPS);
    }

    #[test]
    fn uncorrelated_secondary_gets_zero_weight() {
        let primary = [s(0.0, 12.0)];
        let secondary = [s(1.0, 100.0)];
        let model = CokrigingModel {
            primary: exp_cov(1.0, 3.0),
            cross: zero_cov(),
            secondary: exp_cov(1.0, 3.0),
        };
        let means = Means {
            primary: 10.0,
            secondary: 50.0,
        };
        let ck = simple_cokriging(&primary, &secondary, at(1.0), &model, means).unwrap();
        let sk = simple_kriging(&primary, at(1.0), &model.primary, 10.0).unwrap();
        assert_abs_diff_eq!(ck.value, sk.value, epsilon = EPS);
    }

    #[test]
    fn correlated_secondary_pulls_estimate() {
        let model = CokrigingModel {
            primary: exp_cov(1.0, 3.0),
            cross: CovarianceModel::try_from(
                VariogramModel::cross(
                    0.0,
                    vec![Structure::isotropic(Kind::Exponential, 0.8, 3.0).unwrap()],
                )
                .unwrap(),
            )
            .unwrap(),
            secondary: exp_cov(1.0, 3.0),
        };
        let means = Means {
            primary: 0.0,
            secondary: 0.0,
        };
        let e = simple_cokriging(&[], &[s(0.0, 1.0)], at(0.0), &model, means).unwrap();
        // Secondary at the target with correlation 0.8 ⇒ estimate 0.8, variance 1 − 0.64.
        assert_abs_diff_eq!(e.value, 0.8, epsilon = EPS);
        assert_abs_diff_eq!(e.variance, 0.36, epsilon = EPS);
    }

    #[test]
    fn solutions_report_their_weights() {
        let samples = [s(0.0, 1.0), s(10.0, 3.0)];
        let ok = ordinary_kriging_solution(&samples, at(5.0), &exp_model(1.0, 30.0)).unwrap();
        assert_abs_diff_eq!(ok.weights.iter().sum::<f64>(), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(ok.weights[0], 0.5, epsilon = 1e-12);
        assert!(ok.lagrange.is_some());
        assert_eq!(
            ok.estimate,
            ordinary_kriging(&samples, at(5.0), &exp_model(1.0, 30.0)).unwrap()
        );
        let sk = simple_kriging_solution(&samples, at(5.0), &exp_cov(1.0, 30.0), 2.0).unwrap();
        assert_eq!(sk.weights.len(), 2);
        assert_eq!(sk.lagrange, None);
        let value: f64 = sk.weights[1].mul_add(1.0, sk.weights[0].mul_add(-1.0, 2.0));
        assert_abs_diff_eq!(sk.estimate.value, value, epsilon = 1e-12);
    }
}
