//! Indicator variables for indicator variograms (GSLIB `gamv` types 9 and
//! 10; `docs/spec/gamv.md`, D-050).
//!
//! GSLIB builds an indicator column per cut-off and computes its
//! semivariogram. Here the same is done explicitly: derive the indicator with
//! [`indicator_continuous`] or [`indicator_categorical`], add it to
//! [`Data`](super::Data) as a column, and request a
//! [`VariogramKind::Semivariogram`](super::VariogramKind::Semivariogram) of
//! it, standardised by [`Indicator::sill`] if wanted.
//!
//! # Conventions, and how they differ from `gamv`
//!
//! - **Continuous:** `I = 1` if `v ≤ c`, the convention of IK3D and of
//!   indicator kriging and simulation. `gamv` uses `I = 1` if `v ≥ c`
//!   (`GAMV.FOR` l. 318). γ is the same unless data equal the cut-off, which
//!   the two put on opposite sides.
//! - **Categorical:** `I = 1` if `v = k` exactly, and `v` and `k` must be
//!   integers. `gamv` compares `int(v + 0.5)` with `int(k + 0.5)`, which
//!   truncates toward zero, so it misclassifies negative values (l. 327).
//! - **Sill:** `p(1 − p)` for both, with p the proportion of ones. `gamv`
//!   leaves the categorical sill at 0 (ll. 317–333).
//!
//! GSLIB's conventions are a legacy policy of the parameter-file adapter
//! (`Indicators::GslibGamv`), not of these functions (D-050).

use crate::error::{Error, Result};
use crate::num::to_f64;

/// An indicator variable and its sill.
#[derive(Debug, Clone, PartialEq)]
pub struct Indicator {
    /// 0 or 1 per datum; `None` where the source value is missing.
    pub values: Vec<Option<f64>>,
    /// `p(1 − p)`, with p the proportion of ones among the present values;
    /// `None` if no value is present.
    pub sill: Option<f64>,
}

fn indicator(values: &[Option<f64>], is_one: impl Fn(f64) -> bool) -> Indicator {
    let values: Vec<Option<f64>> = values
        .iter()
        .map(|v| v.map(|v| if is_one(v) { 1.0 } else { 0.0 }))
        .collect();
    let (ones, present) = values
        .iter()
        .flatten()
        .fold((0.0, 0_usize), |(ones, n), &i| (ones + i, n + 1));
    let sill = (present > 0).then(|| {
        let p = ones / to_f64(present);
        p * (1.0 - p)
    });
    Indicator { values, sill }
}

/// The continuous indicator `I = 1` if `v ≤ cutoff`, else 0.
///
/// # Errors
/// [`Error::InvalidInput`] if `cutoff` or a present value is not finite.
pub fn indicator_continuous(values: &[Option<f64>], cutoff: f64) -> Result<Indicator> {
    if !cutoff.is_finite() {
        return Err(Error::invalid(
            "indicator cut-off",
            format!("must be finite, got {cutoff}"),
        ));
    }
    if values.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::invalid("indicator data", "values must be finite"));
    }
    Ok(indicator(values, |v| v <= cutoff))
}

/// The categorical indicator `I = 1` if `v = code`, else 0.
///
/// # Errors
/// [`Error::InvalidInput`] if `code` or a present value is not a finite
/// integer.
pub fn indicator_categorical(values: &[Option<f64>], code: f64) -> Result<Indicator> {
    let integral = |v: f64| v.is_finite() && v.fract() == 0.0;
    if !integral(code) {
        return Err(Error::invalid(
            "category code",
            format!("must be an integer, got {code}"),
        ));
    }
    if let Some(bad) = values.iter().flatten().find(|v| !integral(**v)) {
        return Err(Error::invalid(
            "categorical data",
            format!("values must be integer codes, got {bad}"),
        ));
    }
    #[expect(
        clippy::float_cmp,
        reason = "both sides are checked to be integers, so equality is exact"
    )]
    let is_one = move |v: f64| v == code;
    Ok(indicator(values, is_one))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::variogram::{
        CentredLags, Data, Direction, GamvPolicies, GamvRequest, VariogramKind, VariogramSpec, gamv,
    };
    use approx::assert_abs_diff_eq;

    fn column(values: &[f64]) -> Vec<Option<f64>> {
        values.iter().copied().map(Some).collect()
    }

    #[test]
    fn continuous_is_one_at_and_below_the_cutoff() {
        let ind = indicator_continuous(&column(&[1.0, 2.0, 2.5, 3.0]), 2.0).unwrap();
        assert_eq!(ind.values, column(&[1.0, 1.0, 0.0, 0.0]), "v = c gives 1");
        assert_abs_diff_eq!(ind.sill.unwrap(), 0.25, epsilon = 1e-15);
    }

    #[test]
    fn categorical_matches_the_code_exactly() {
        let ind = indicator_categorical(&column(&[0.0, 1.0, 1.0, 2.0, -1.0]), 1.0).unwrap();
        assert_eq!(ind.values, column(&[0.0, 1.0, 1.0, 0.0, 0.0]));
        // p = 2/5: the sill is p(1 − p), where gamv leaves it at 0.
        assert_abs_diff_eq!(ind.sill.unwrap(), 0.4 * 0.6, epsilon = 1e-15);
        // Negative codes are categories like any other (gamv's int(v + 0.5)
        // puts −1 into category 0).
        let neg = indicator_categorical(&column(&[-1.0, 0.0]), -1.0).unwrap();
        assert_eq!(neg.values, column(&[1.0, 0.0]));
    }

    #[test]
    fn non_integral_categories_are_errors() {
        assert!(indicator_categorical(&column(&[1.0]), 1.5).is_err(), "code");
        assert!(
            indicator_categorical(&column(&[1.0, 0.7]), 1.0).is_err(),
            "value"
        );
        assert!(
            indicator_categorical(&column(&[1.0]), f64::NAN).is_err(),
            "NaN code"
        );
        assert!(
            indicator_continuous(&column(&[1.0]), f64::INFINITY).is_err(),
            "cut-off"
        );
        assert!(
            indicator_continuous(&[Some(f64::NAN)], 1.0).is_err(),
            "NaN value"
        );
    }

    #[test]
    fn missing_stays_missing() {
        let values = vec![Some(1.0), None, Some(3.0)];
        let ind = indicator_continuous(&values, 2.0).unwrap();
        assert_eq!(ind.values, vec![Some(1.0), None, Some(0.0)]);
        assert_abs_diff_eq!(ind.sill.unwrap(), 0.25, epsilon = 1e-15);
        let none = indicator_categorical(&[None, None], 1.0).unwrap();
        assert_eq!(none.sill, None, "no present values: no sill");
    }

    #[test]
    fn constant_indicator_has_zero_sill() {
        let ind = indicator_continuous(&column(&[1.0, 2.0]), 5.0).unwrap();
        assert_eq!(ind.sill, Some(0.0));
    }

    /// End to end: the standardised indicator semivariogram, as GSLIB's type 9
    /// with `isill = 1`. Four data along a line, indicator [1, 1, 0, 0]: at
    /// distance 1 the pairs differ once in three, so γ = (1/3)/2 = 1/6, and
    /// standardised by p(1 − p) = 1/4 it is 2/3.
    #[test]
    fn standardised_indicator_semivariogram() {
        let ind = indicator_continuous(&column(&[0.5, 1.0, 3.0, 4.0]), 1.0).unwrap();
        let points = (0..4)
            .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
            .collect();
        let data = Data::new(points, vec![ind.values]).unwrap();
        let request = GamvRequest {
            lags: CentredLags::new(2, 1.0, 0.25).unwrap(),
            directions: vec![Direction::new(90.0, 10.0, 1e9, 0.0, 90.0, 1e9).unwrap()],
            variograms: vec![VariogramSpec {
                tail: 0,
                head: 0,
                kind: VariogramKind::Semivariogram,
                standardize: true,
            }],
            sills: vec![ind.sill],
            policies: GamvPolicies::default(),
        };
        let lag = gamv(&data, &request).unwrap()[0].lags[2];
        assert_eq!(lag.pairs, 3);
        assert_abs_diff_eq!(lag.gamma.unwrap(), (1.0 / 6.0) / 0.25, epsilon = 1e-12);
    }
}
