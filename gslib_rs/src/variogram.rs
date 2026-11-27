//! Experimental (empirical) semivariograms.

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::num::{floor_index, to_f64};
use crate::sample::Sample;

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

    /// Which class a separation distance falls in; `None` beyond the maximum.
    fn classify(&self, h: f64) -> Option<usize> {
        if h > self.max_distance {
            return None;
        }
        // h == max_distance belongs to the last class.
        floor_index(h / self.width).map(|i| i.min(self.count - 1))
    }

    fn center(&self, class: usize) -> f64 {
        self.width * (to_f64(class) + 0.5)
    }
}

/// One lag of an experimental variogram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lag {
    /// Midpoint of the distance class.
    pub center: f64,
    /// Mean separation of the pairs in the class (`None` if there are none).
    pub mean_distance: Option<f64>,
    /// Number of pairs in the class.
    pub pairs: usize,
    /// Semivariance γ, or `None` if the class has no pairs.
    pub gamma: Option<f64>,
}

/// Omnidirectional experimental semivariogram of `samples`.
#[must_use]
pub fn experimental_variogram(samples: &[Sample], lags: &LagClasses) -> Vec<Lag> {
    let points: Vec<Point3> = samples.iter().map(Sample::location).collect();
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    semivariogram(&points, &values, lags)
}

/// Experimental semivariogram of `values` at `points` (zipped; any excess in
/// the longer slice is ignored).
pub(crate) fn semivariogram(points: &[Point3], values: &[f64], lags: &LagClasses) -> Vec<Lag> {
    #[derive(Clone, Copy, Default)]
    struct Acc {
        gamma: f64,
        distance: f64,
        pairs: usize,
    }
    let data: Vec<(Point3, f64)> = points.iter().copied().zip(values.iter().copied()).collect();
    let mut acc = vec![Acc::default(); lags.count];
    for (i, &(p, v)) in data.iter().enumerate() {
        for &(q, w) in data.iter().skip(i + 1) {
            let h = p.distance(q);
            if let Some(slot) = lags.classify(h).and_then(|c| acc.get_mut(c)) {
                let d = v - w;
                slot.gamma = (0.5 * d).mul_add(d, slot.gamma);
                slot.distance += h;
                slot.pairs += 1;
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
}
