//! Experimental (empirical) semivariograms.
//!
//! Two ways to group pairs by separation, both behind [`LagScheme`]:
//! - [`LagClasses`]: equal-width, non-overlapping classes over `[0, max]`;
//! - [`CentredLags`]: GSLIB `gamv`'s classes, centred on multiples of a lag
//!   spacing with a tolerance (they may overlap), plus a zero-distance class.

use std::ops::Range;

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::num::{floor_index, to_f64};
use crate::sample::Sample;

/// A way of assigning pair separations to lag classes.
pub trait LagScheme {
    /// Number of classes.
    fn len(&self) -> usize;

    /// Whether there are no classes (never true for a validated scheme).
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The nominal distance of a class.
    fn center(&self, class: usize) -> f64;

    /// The classes a separation `h` belongs to: a possibly empty range,
    /// since classes may overlap.
    fn classes(&self, h: f64) -> Range<usize>;
}

/// Equal-width distance classes covering `[0, max_distance]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LagClasses {
    width: f64,
    count: usize,
    max_distance: f64,
}

impl LagClasses {
    /// `count` classes of width `max_distance / count`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `count` is zero or `max_distance` is not
    /// finite and positive.
    pub fn new(max_distance: f64, count: usize) -> Result<Self> {
        if count == 0 {
            return Err(Error::invalid("lag count", "must be at least 1"));
        }
        if !(max_distance.is_finite() && max_distance > 0.0) {
            return Err(Error::invalid(
                "maximum lag distance",
                format!("must be finite and positive, got {max_distance}"),
            ));
        }
        Ok(Self {
            width: max_distance / to_f64(count),
            count,
            max_distance,
        })
    }

    /// Number of classes.
    #[must_use]
    pub const fn count(&self) -> usize {
        self.count
    }
}

impl LagScheme for LagClasses {
    fn len(&self) -> usize {
        self.count
    }

    fn center(&self, class: usize) -> f64 {
        self.width * (to_f64(class) + 0.5)
    }

    fn classes(&self, h: f64) -> Range<usize> {
        if h > self.max_distance {
            return 0..0;
        }
        // h == max_distance belongs to the last class.
        floor_index(h / self.width).map_or(0..0, |i| {
            let c = i.min(self.count - 1);
            c..c + 1
        })
    }
}

/// GSLIB `gamv`'s lag classes (`vendor/gslib-2.0/src/GAMV/GAMV.FOR`, F-008).
///
/// Class 0 holds coincident pairs (separation ≤ 1e-20, GSLIB's `EPSLON`).
/// Class `k + 1`, for `k = 0..=count`, holds separations within `tolerance`
/// of `k · spacing`, inclusive at both ends, so classes overlap when
/// `tolerance > spacing / 2` and a pair then counts in several. Pairs farther
/// apart than `(count + 0.5) · spacing` are ignored.
///
/// Unlike GSLIB, a datum is never paired with itself: GSLIB's class 0 also
/// counts every datum against itself, diluting γ there towards zero
/// (ROADMAP §2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CentredLags {
    count: usize,
    spacing: f64,
    tolerance: f64,
}

/// GSLIB `gamv`'s `EPSLON`: separations up to this are "zero".
const ZERO_SEPARATION: f64 = 1e-20;

impl CentredLags {
    /// `count` lags of `spacing` (GSLIB `nlag`, `xlag`) with half-width
    /// `tolerance` (`xltol`).
    ///
    /// GSLIB replaces a non-positive tolerance with `spacing / 2`; here it is
    /// an error (explicit failure, `AGENTS.md`).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `count` is zero, or `spacing` or
    /// `tolerance` is not finite and positive.
    pub fn new(count: usize, spacing: f64, tolerance: f64) -> Result<Self> {
        if count == 0 {
            return Err(Error::invalid("lag count", "must be at least 1"));
        }
        for (what, v) in [("lag spacing", spacing), ("lag tolerance", tolerance)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(Error::invalid(
                    what,
                    format!("must be finite and positive, got {v}"),
                ));
            }
        }
        Ok(Self {
            count,
            spacing,
            tolerance,
        })
    }
}

impl LagScheme for CentredLags {
    fn len(&self) -> usize {
        self.count + 2
    }

    fn center(&self, class: usize) -> f64 {
        class
            .checked_sub(1)
            .map_or(0.0, |k| self.spacing * to_f64(k))
    }

    fn classes(&self, h: f64) -> Range<usize> {
        let limit = (to_f64(self.count) + 0.5) * self.spacing;
        if h.is_nan() || h > limit {
            return 0..0;
        }
        if h <= ZERO_SEPARATION {
            return 0..1;
        }
        let within = |k: usize| {
            let c = self.spacing * to_f64(k);
            h >= c - self.tolerance && h <= c + self.tolerance
        };
        let first = (0..=self.count).find(|&k| within(k));
        let last = (0..=self.count).rev().find(|&k| within(k));
        match (first, last) {
            (Some(a), Some(b)) => a + 1..b + 2,
            _ => 0..0,
        }
    }
}

/// One lag of an experimental variogram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lag {
    /// Nominal distance of the class.
    pub center: f64,
    /// Mean separation of the pairs in the class (`None` if there are none).
    pub mean_distance: Option<f64>,
    /// Number of pairs in the class.
    pub pairs: usize,
    /// Semivariance γ, or `None` if the class has no pairs.
    pub gamma: Option<f64>,
    /// Mean of the values in the class's pairs (both ends of every pair);
    /// GSLIB's head and tail means for an omnidirectional variogram.
    pub pair_mean: Option<f64>,
}

/// Omnidirectional experimental semivariogram of `samples`.
#[must_use]
pub fn experimental_variogram<L: LagScheme + ?Sized>(samples: &[Sample], lags: &L) -> Vec<Lag> {
    let points: Vec<Point3> = samples.iter().map(Sample::location).collect();
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    semivariogram(&points, &values, lags)
}

/// Experimental semivariogram of `values` at `points` (zipped; any excess in
/// the longer slice is ignored). Each unordered pair of distinct data counts
/// once in every class its separation belongs to.
pub(crate) fn semivariogram<L: LagScheme + ?Sized>(
    points: &[Point3],
    values: &[f64],
    lags: &L,
) -> Vec<Lag> {
    #[derive(Clone, Copy, Default)]
    struct Acc {
        gamma: f64,
        distance: f64,
        sum: f64,
        pairs: usize,
    }
    let data: Vec<(Point3, f64)> = points.iter().copied().zip(values.iter().copied()).collect();
    let mut acc = vec![Acc::default(); lags.len()];
    for (i, &(p, v)) in data.iter().enumerate() {
        for &(q, w) in data.iter().skip(i + 1) {
            let h = p.distance(q);
            let d = v - w;
            for class in lags.classes(h) {
                if let Some(a) = acc.get_mut(class) {
                    a.gamma = (0.5 * d).mul_add(d, a.gamma);
                    a.distance += h;
                    a.sum += v + w;
                    a.pairs += 1;
                }
            }
        }
    }
    acc.into_iter()
        .enumerate()
        .map(|(class, a)| {
            let n = (a.pairs > 0).then(|| to_f64(a.pairs));
            Lag {
                center: lags.center(class),
                mean_distance: n.map(|n| a.distance / n),
                pairs: a.pairs,
                gamma: n.map(|n| a.gamma / n),
                pair_mean: n.map(|n| a.sum / (2.0 * n)),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }

    #[test]
    fn two_points_one_lag() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(1.0, 1.0)],
            &LagClasses::new(1.0, 1).unwrap(),
        );
        assert_eq!(lags.len(), 1);
        assert_abs_diff_eq!(lags[0].center, 0.5);
        assert_eq!(lags[0].pairs, 1);
        assert_abs_diff_eq!(lags[0].gamma.unwrap(), 0.5);
        assert_abs_diff_eq!(lags[0].mean_distance.unwrap(), 1.0);
    }

    #[test]
    fn empty_lags_are_none_not_zero() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(3.5, 1.0)],
            &LagClasses::new(4.0, 4).unwrap(),
        );
        assert_eq!(
            lags.iter().map(|l| l.pairs).collect::<Vec<_>>(),
            vec![0, 0, 0, 1]
        );
        assert_eq!(lags[0].gamma, None);
        assert!(lags[3].gamma.is_some(), "the populated lag has a value");
    }

    #[test]
    fn pairs_beyond_max_distance_are_excluded() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(5.0, 1.0)],
            &LagClasses::new(4.0, 2).unwrap(),
        );
        assert!(lags.iter().all(|l| l.pairs == 0), "{lags:?}");
    }

    #[test]
    fn linear_trend_gives_quadratic_gamma() {
        // z = x on a unit grid ⇒ γ(h) = h²/2 exactly.
        let samples: Vec<Sample> = (0..10).map(|i| s(f64::from(i), f64::from(i))).collect();
        let lags = experimental_variogram(&samples, &LagClasses::new(3.5, 7).unwrap());
        for lag in lags.iter().filter(|l| l.pairs > 0) {
            let h = lag.mean_distance.unwrap();
            assert_abs_diff_eq!(lag.gamma.unwrap(), h * h / 2.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn invalid_classes_rejected() {
        assert!(LagClasses::new(1.0, 0).is_err(), "zero lags");
        assert!(LagClasses::new(0.0, 3).is_err(), "zero distance");
        assert!(LagClasses::new(f64::NAN, 3).is_err(), "NaN distance");
    }

    fn centred(count: usize, spacing: f64, tolerance: f64) -> CentredLags {
        CentredLags::new(count, spacing, tolerance).unwrap()
    }

    #[test]
    fn centred_lags_follow_gamv() {
        // nlag = 3, xlag = 10, xltol = 5: classes 0 (zero), then centres 0,
        // 10, 20, 30; GAMV.FOR's bounds are inclusive at both ends.
        let lags = centred(3, 10.0, 5.0);
        assert_eq!(lags.len(), 5);
        assert_eq!(lags.classes(0.0), 0..1, "coincident pairs");
        assert_eq!(lags.classes(3.0), 1..2, "centre 0");
        assert_eq!(lags.classes(5.0), 1..3, "on a shared bound: both classes");
        assert_eq!(lags.classes(12.0), 2..3);
        assert_eq!(lags.classes(35.0), 4..5, "up to (nlag + 0.5) * xlag");
        assert_eq!(lags.classes(35.5), 0..0, "beyond the limit");
        assert_eq!(lags.classes(f64::NAN), 0..0, "NaN is never classified");
        assert_eq!(
            (0..5).map(|c| lags.center(c)).collect::<Vec<_>>(),
            vec![0.0, 0.0, 10.0, 20.0, 30.0]
        );
    }

    #[test]
    fn centred_lags_overlap_with_wide_tolerance() {
        // xltol > xlag / 2: 10 lies within 8 of centre 10 only; 14 lies
        // within 8 of both 10 and 20.
        let lags = centred(4, 10.0, 8.0);
        assert_eq!(lags.classes(10.0), 2..3);
        assert_eq!(lags.classes(14.0), 2..4, "counted in two classes");
    }

    #[test]
    fn centred_variogram_counts_each_distinct_pair_once_per_class() {
        let data = [s(0.0, 1.0), s(0.0, 3.0), s(10.0, 2.0), s(14.0, 6.0)];
        let lags = experimental_variogram(&data, &centred(2, 10.0, 8.0));
        // Zero class: one coincident pair (1, 3), not the four self-pairs.
        assert_eq!(lags[0].pairs, 1);
        assert_abs_diff_eq!(lags[0].gamma.unwrap(), 2.0);
        assert_abs_diff_eq!(lags[0].pair_mean.unwrap(), 2.0);
        // Separation 4 (10-14) is within 8 of centre 0 and of 10.
        assert_eq!(lags[1].pairs, 1);
        // Class 2 is centred on 10 with tolerance 8, i.e. [2, 18].
        assert_eq!(lags[2].pairs, 5, "pairs at 10, 10, 14, 14 and 4");
        // Class 3 is centred on 20, i.e. [12, 28].
        assert_eq!(lags[3].pairs, 2, "pairs at 14 and 14");
    }

    #[test]
    fn centred_lags_reject_invalid() {
        assert!(CentredLags::new(0, 1.0, 0.5).is_err(), "zero lags");
        assert!(CentredLags::new(3, 0.0, 0.5).is_err(), "zero spacing");
        assert!(
            CentredLags::new(3, 1.0, 0.0).is_err(),
            "GSLIB's default tolerance must be explicit"
        );
    }
}
