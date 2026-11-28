//! Normal-score transform and back-transform (GSLIB `nscore`, `backtr`).

use std::f64::consts::SQRT_2;

use statrs::function::erf::erfc_inv;

use crate::error::{Error, Result, check_lengths};

/// Standard normal quantile Φ⁻¹(p) for 0 < p < 1 (statrs `erfc_inv`,
/// accurate to ~1e-16 relative).
fn normal_quantile(p: f64) -> f64 {
    -SQRT_2 * erfc_inv(2.0 * p)
}

/// Standard normal CDF Φ(y).
///
/// Uses `libm::erfc` (FreeBSD msun, < 1 ulp): statrs's `erfc` is only
/// accurate to ~5e-11 relative, which would leak into every back-transform.
fn normal_cdf(y: f64) -> f64 {
    0.5 * libm::erfc(-y / SQRT_2)
}

/// The standard normal quantile Φ⁻¹(p), to about 1e-16.
///
/// GSLIB's `gauinv` is a rational approximation accurate to about 1.5e-8
/// (F-007 in `docs/spec/decisions.md`); this is the exact function it
/// approximates.
///
/// # Errors
/// [`Error::InvalidInput`] unless `0 < p < 1`.
pub fn standard_normal_quantile(p: f64) -> Result<f64> {
    if p > 0.0 && p < 1.0 {
        Ok(normal_quantile(p))
    } else {
        Err(Error::invalid(
            "probability",
            format!("must be strictly between 0 and 1, got {p}"),
        ))
    }
}

/// The standard normal CDF Φ(y), to within an ulp or so.
///
/// GSLIB's `gcum` is a polynomial approximation accurate to about 7.5e-8
/// absolute (F-007); this is the exact function it approximates. NaN maps to
/// NaN.
#[must_use]
pub fn standard_normal_cdf(y: f64) -> f64 {
    normal_cdf(y)
}

/// One row of a transformation table.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Entry {
    value: f64,
    score: f64,
    cdf: f64,
}

/// The mapping between original values and normal scores, sorted ascending
/// with strictly increasing values and scores.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalScoreTable {
    entries: Vec<Entry>,
}

impl NormalScoreTable {
    /// A table from `(value, score)` pairs, e.g. read back from a file.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the table is empty, has non-finite entries,
    /// or values or scores are not strictly increasing.
    pub fn from_pairs(pairs: &[(f64, f64)]) -> Result<Self> {
        if pairs.is_empty() {
            return Err(Error::invalid("transformation table", "is empty"));
        }
        if pairs
            .iter()
            .any(|&(v, s)| !(v.is_finite() && s.is_finite()))
        {
            return Err(Error::invalid(
                "transformation table",
                "has non-finite entries",
            ));
        }
        if pairs
            .windows(2)
            .any(|w| matches!(w, [(v0, s0), (v1, s1)] if !(v0 < v1 && s0 < s1)))
        {
            return Err(Error::invalid(
                "transformation table",
                "values and scores must be strictly increasing",
            ));
        }
        Ok(Self {
            entries: pairs
                .iter()
                .map(|&(value, score)| Entry {
                    value,
                    score,
                    cdf: normal_cdf(score),
                })
                .collect(),
        })
    }

    /// The `(value, score)` rows, ascending.
    pub fn pairs(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.entries.iter().map(|e| (e.value, e.score))
    }
}

/// Result of [`normal_scores`].
#[derive(Debug, Clone, PartialEq)]
pub struct NormalScores {
    /// The score of each input value, in input order.
    pub scores: Vec<f64>,
    /// The transformation table, for back-transforming.
    pub table: NormalScoreTable,
}

/// Transform values to standard normal scores by rank (GSLIB `nscore`).
///
/// Each distinct value gets the normal quantile of the midpoint of its
/// cumulative weight interval. **Ties share one score** (the midpoint of their
/// combined interval) instead of being ranked arbitrarily. With `None`
/// weights, every value has weight 1.
///
/// # Errors
/// [`Error::InvalidInput`] if `values` is empty or contains non-finite values,
/// or a weight is not finite and positive; [`Error::LengthMismatch`] if
/// `weights` has the wrong length.
pub fn normal_scores(values: &[f64], weights: Option<&[f64]>) -> Result<NormalScores> {
    if values.is_empty() {
        return Err(Error::invalid(
            "values",
            "cannot normal-score an empty data set",
        ));
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(Error::invalid("values", "contain NaN or infinite entries"));
    }
    let weights: Vec<f64> = match weights {
        Some(w) => {
            check_lengths("values/weights", values.len(), w.len())?;
            if w.iter().any(|&x| !(x.is_finite() && x > 0.0)) {
                return Err(Error::invalid("weights", "must all be finite and positive"));
            }
            w.to_vec()
        }
        None => vec![1.0; values.len()],
    };
    let total: f64 = weights.iter().sum();

    let mut sorted: Vec<(usize, f64, f64)> = values
        .iter()
        .zip(&weights)
        .enumerate()
        .map(|(i, (&v, &w))| (i, v, w))
        .collect();
    sorted.sort_by(|a, b| a.1.total_cmp(&b.1));

    let mut cumulative = 0.0;
    let mut entries = Vec::new();
    let mut scored: Vec<(usize, f64)> = Vec::with_capacity(values.len());
    for group in sorted.chunk_by(|a, b| a.1.total_cmp(&b.1).is_eq()) {
        let group_weight: f64 = group.iter().map(|g| g.2).sum();
        let cdf = group_weight.mul_add(0.5, cumulative) / total;
        let score = normal_quantile(cdf);
        cumulative += group_weight;
        scored.extend(group.iter().map(|g| (g.0, score)));
        if let Some(&(_, value, _)) = group.first() {
            entries.push(Entry { value, score, cdf });
        }
    }
    scored.sort_by_key(|&(i, _)| i);
    Ok(NormalScores {
        scores: scored.into_iter().map(|(_, s)| s).collect(),
        table: NormalScoreTable { entries },
    })
}

/// How to back-transform scores beyond the ends of the table.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Tail {
    /// Return the nearest table value (no extrapolation).
    Clamp,
    /// Extrapolate linearly in probability to this bound, which is reached at
    /// probability 0 (lower tail) or 1 (upper tail). GSLIB tail option 1
    /// (verified against `BACKTR.FOR`; GSLIB computes the probabilities with
    /// its `gcum` approximation, F-007).
    Linear(f64),
}

/// A validated back-transform: a table plus tail models.
#[derive(Debug, Clone, PartialEq)]
pub struct BackTransform {
    table: NormalScoreTable,
    first: Entry,
    last: Entry,
    lower: Tail,
    upper: Tail,
}

impl BackTransform {
    /// # Errors
    /// [`Error::InvalidInput`] if a linear tail bound is not finite or lies
    /// inside the table's value range.
    pub fn new(table: NormalScoreTable, lower: Tail, upper: Tail) -> Result<Self> {
        let (Some(&first), Some(&last)) = (table.entries.first(), table.entries.last()) else {
            return Err(Error::invalid("transformation table", "is empty"));
        };
        if let Tail::Linear(zmin) = lower
            && !(zmin.is_finite() && zmin <= first.value)
        {
            return Err(Error::invalid(
                "lower tail bound",
                format!("must be finite and <= the minimum value {}", first.value),
            ));
        }
        if let Tail::Linear(zmax) = upper
            && !(zmax.is_finite() && zmax >= last.value)
        {
            return Err(Error::invalid(
                "upper tail bound",
                format!("must be finite and >= the maximum value {}", last.value),
            ));
        }
        Ok(Self {
            table,
            first,
            last,
            lower,
            upper,
        })
    }

    /// Back-transform one normal score to the original scale.
    #[must_use]
    pub fn apply(&self, score: f64) -> f64 {
        let (first, last) = (&self.first, &self.last);
        if score <= first.score {
            return match self.lower {
                Tail::Clamp => first.value,
                Tail::Linear(zmin) => {
                    (first.value - zmin).mul_add(normal_cdf(score) / first.cdf, zmin)
                }
            };
        }
        if score >= last.score {
            return match self.upper {
                Tail::Clamp => last.value,
                Tail::Linear(zmax) => {
                    let frac = (normal_cdf(score) - last.cdf) / (1.0 - last.cdf);
                    (zmax - last.value).mul_add(frac, last.value)
                }
            };
        }
        // Interior: linear interpolation in score space between the
        // bracketing table rows, as GSLIB's `backtr` does (F-007).
        let entries = &self.table.entries;
        let above = entries.partition_point(|e| e.score <= score);
        match (
            above.checked_sub(1).and_then(|i| entries.get(i)),
            entries.get(above),
        ) {
            (Some(lo), Some(hi)) => {
                let frac = (score - lo.score) / (hi.score - lo.score);
                (hi.value - lo.value).mul_add(frac, lo.value)
            }
            _ => last.value,
        }
    }
}

/// The measured bound on `|gslib_gauinv(p) − Φ⁻¹(p)|` for
/// `1e-10 ≤ p ≤ 1 − 1e-10`: the worst error is 1.501e-8, near p = 0.033
/// (F-007, D-016; pinned by a test).
pub const GSLIB_GAUINV_MAX_ERROR: f64 = 2e-8;

/// GSLIB's `gauinv`, reproduced exactly (`vendor/gslib-2.0/src/GSLIB/GAUINV.FOR`;
/// `docs/spec/rng.md` RN6).
///
/// A rational approximation to Φ⁻¹ (Kennedy & Gentle 1980, as its comment
/// cites), within [`GSLIB_GAUINV_MAX_ERROR`] of
/// [`standard_normal_quantile`]. It exists for exact-stream parity, where
/// GSLIB's deviates must be reproduced bit for bit; use
/// [`standard_normal_quantile`] otherwise.
///
/// **GSLIB's contract, kept as is:** below `1e-10` it returns `−1e10`, above
/// `1 − 1e-10` it returns `1e10`, and exactly at 0.5 it returns 0. The
/// arithmetic follows the Fortran's evaluation order without fused
/// multiply-adds, as the reference build is compiled
/// (`-ffp-contract=off`).
#[must_use]
#[expect(
    clippy::suboptimal_flops,
    reason = "reproduce GSLIB's unfused arithmetic bit for bit"
)]
pub fn gslib_gauinv(p: f64) -> f64 {
    const LIM: f64 = 1.0e-10;
    let (p0, p1, p2, p3, p4) = (
        -0.322_232_431_088,
        -1.0,
        -0.342_242_088_547,
        -0.020_423_121_024_5,
        -0.000_045_364_221_014_8,
    );
    let (q0, q1, q2, q3, q4) = (
        0.099_348_462_606_0,
        0.588_581_570_495,
        0.531_103_462_366,
        0.103_537_752_850,
        0.003_856_070_063_4,
    );
    if p < LIM {
        return -1.0e10;
    }
    if p > 1.0 - LIM {
        return 1.0e10;
    }
    let pp = if p > 0.5 { 1.0 - p } else { p };
    if p.to_bits() == 0.5_f64.to_bits() {
        return 0.0;
    }
    let y = (1.0 / (pp * pp)).ln().sqrt();
    let xp = y
        + ((((y * p4 + p3) * y + p2) * y + p1) * y + p0)
            / ((((y * q4 + q3) * y + q2) * y + q1) * y + q0);
    // GSLIB: `if(real(p).eq.real(pp)) xp = -xp`.
    if p.to_bits() == pp.to_bits() { -xp } else { xp }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn quantile_and_cdf_are_inverse() {
        assert_abs_diff_eq!(normal_quantile(0.5), 0.0, epsilon = EPS);
        // Standard values to 16 digits; they agree within 4.4e-16 with
        // Python's statistics.NormalDist (CPython 3.13), an independent
        // implementation (F-007).
        assert_abs_diff_eq!(
            normal_quantile(0.975),
            1.959_963_984_540_054,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            normal_quantile(0.1),
            -1.281_551_565_544_600_4,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            normal_quantile(1e-10),
            -6.361_340_902_404_056,
            epsilon = 1e-14
        );
        assert_abs_diff_eq!(normal_cdf(-1.281_551_565_544_600_4), 0.1, epsilon = 1e-16);
        assert_abs_diff_eq!(normal_cdf(normal_quantile(0.1)), 0.1, epsilon = 1e-16);
    }

    #[test]
    fn public_quantile_rejects_out_of_range() {
        assert!(standard_normal_quantile(0.0).is_err(), "p = 0");
        assert!(standard_normal_quantile(1.0).is_err(), "p = 1");
        assert!(standard_normal_quantile(f64::NAN).is_err(), "NaN");
        assert_abs_diff_eq!(standard_normal_quantile(0.5).unwrap(), 0.0, epsilon = EPS);
        assert_abs_diff_eq!(standard_normal_cdf(0.0), 0.5, epsilon = EPS);
    }

    #[test]
    fn scores_are_symmetric_and_centered() {
        let out = normal_scores(&[3.0, 1.0, 2.0], None).unwrap();
        assert_abs_diff_eq!(out.scores[2], 0.0, epsilon = EPS);
        assert_abs_diff_eq!(out.scores[0], -out.scores[1], epsilon = EPS);
        assert_abs_diff_eq!(out.scores[0], normal_quantile(5.0 / 6.0), epsilon = EPS);
    }

    #[test]
    fn ties_share_a_score() {
        let out = normal_scores(&[1.0, 2.0, 2.0, 3.0], None).unwrap();
        assert_abs_diff_eq!(out.scores[1], out.scores[2]);
        assert_abs_diff_eq!(out.scores[1], 0.0, epsilon = EPS);
        assert_eq!(out.table.pairs().count(), 3);
    }

    #[test]
    fn weights_shift_quantiles() {
        // Weight 3 on the low value: its interval is [0, 0.75], midpoint 0.375.
        let out = normal_scores(&[1.0, 2.0], Some(&[3.0, 1.0])).unwrap();
        assert_abs_diff_eq!(out.scores[0], normal_quantile(0.375), epsilon = EPS);
        assert_abs_diff_eq!(out.scores[1], normal_quantile(0.875), epsilon = EPS);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(normal_scores(&[], None).is_err(), "empty");
        assert!(normal_scores(&[1.0, f64::NAN], None).is_err(), "NaN");
        assert!(
            normal_scores(&[1.0, 2.0], Some(&[1.0, 0.0])).is_err(),
            "zero weight"
        );
        assert!(normal_scores(&[1.0, 2.0], Some(&[1.0])).is_err(), "length");
    }

    #[test]
    fn back_transform_interior_and_clamped_tails() {
        let table =
            NormalScoreTable::from_pairs(&[(10.0, -1.0), (20.0, 0.0), (30.0, 1.0)]).unwrap();
        let bt = BackTransform::new(table, Tail::Clamp, Tail::Clamp).unwrap();
        assert_abs_diff_eq!(bt.apply(0.0), 20.0);
        assert_abs_diff_eq!(bt.apply(0.5), 25.0, epsilon = EPS);
        assert_abs_diff_eq!(bt.apply(-5.0), 10.0);
        assert_abs_diff_eq!(bt.apply(5.0), 30.0);
    }

    #[test]
    fn linear_tails_reach_bounds_at_extreme_probabilities() {
        let table = NormalScoreTable::from_pairs(&[(10.0, -1.0), (30.0, 1.0)]).unwrap();
        let bt = BackTransform::new(table, Tail::Linear(0.0), Tail::Linear(50.0)).unwrap();
        assert_abs_diff_eq!(bt.apply(-40.0), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(bt.apply(40.0), 50.0, epsilon = 1e-9);
        assert_abs_diff_eq!(bt.apply(-1.0), 10.0, epsilon = EPS);
        let lower_mid = bt.apply(-1.5);
        assert!(lower_mid > 0.0 && lower_mid < 10.0, "{lower_mid}");
    }

    #[test]
    fn tail_bounds_must_lie_outside_data() {
        let table = NormalScoreTable::from_pairs(&[(10.0, -1.0), (30.0, 1.0)]).unwrap();
        assert!(
            BackTransform::new(table.clone(), Tail::Linear(15.0), Tail::Clamp).is_err(),
            "zmin inside"
        );
        assert!(
            BackTransform::new(table, Tail::Clamp, Tail::Linear(20.0)).is_err(),
            "zmax inside"
        );
    }

    #[test]
    fn table_must_be_increasing() {
        assert!(
            NormalScoreTable::from_pairs(&[(1.0, 0.0), (1.0, 1.0)]).is_err(),
            "repeated value"
        );
        assert!(
            NormalScoreTable::from_pairs(&[(1.0, 1.0), (2.0, 0.0)]).is_err(),
            "decreasing score"
        );
    }

    proptest! {
        #[test]
        fn back_transform_inverts_scores(values in prop::collection::vec(-1e3..1e3_f64, 1..60)) {
            let out = normal_scores(&values, None).unwrap();
            let bt = BackTransform::new(out.table, Tail::Clamp, Tail::Clamp).unwrap();
            for (&v, &y) in values.iter().zip(&out.scores) {
                prop_assert!((bt.apply(y) - v).abs() <= 1e-9 * v.abs().max(1.0));
            }
        }

        #[test]
        fn scores_preserve_order(values in prop::collection::vec(-1e3..1e3_f64, 2..60)) {
            let out = normal_scores(&values, None).unwrap();
            for (i, (&a, &ya)) in values.iter().zip(&out.scores).enumerate() {
                for (&b, &yb) in values.iter().zip(&out.scores).skip(i + 1) {
                    prop_assert_eq!(a.total_cmp(&b), ya.total_cmp(&yb));
                }
            }
        }
    }

    /// Probabilities on a log grid from 1e-10 to 0.5 plus a linear grid,
    /// mirrored above 0.5.
    fn probabilities() -> Vec<f64> {
        let log = (0..=20_000).map(|k| 10_f64.powf(-10.0 + 10.0 * f64::from(k) / 20_000.0));
        let lin = (1..100_000).map(|k| f64::from(k) / 200_000.0);
        let lower: Vec<f64> = log
            .chain(lin)
            .filter(|&p| (1e-10..=0.5).contains(&p))
            .collect();
        let upper: Vec<f64> = lower.iter().map(|p| 1.0 - p).collect();
        lower.into_iter().chain(upper).collect()
    }

    #[test]
    fn gslib_gauinv_error_bound() {
        // RN6: the bound D-016 measured, pinned in both directions.
        let worst = probabilities()
            .into_iter()
            .map(|p| (gslib_gauinv(p) - standard_normal_quantile(p).unwrap()).abs())
            .fold(0.0, f64::max);
        assert!(
            worst <= GSLIB_GAUINV_MAX_ERROR,
            "error {worst:e} exceeds the bound"
        );
        assert!(
            worst > GSLIB_GAUINV_MAX_ERROR / 2.0,
            "bound is loose: {worst:e}"
        );
    }

    #[test]
    fn gslib_gauinv_keeps_gslibs_contract() {
        assert_eq!(gslib_gauinv(0.5).to_bits(), 0.0_f64.to_bits());
        assert_eq!(gslib_gauinv(0.0).to_bits(), (-1.0e10_f64).to_bits());
        assert_eq!(gslib_gauinv(1.0 - 1e-11).to_bits(), 1.0e10_f64.to_bits());
        assert!((gslib_gauinv(0.3) + gslib_gauinv(0.7)).abs() < 1e-15, "odd");
        assert!(gslib_gauinv(0.25) < 0.0 && gslib_gauinv(0.75) > 0.0);
    }
}
