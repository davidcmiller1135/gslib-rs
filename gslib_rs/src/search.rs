//! Search neighbourhoods: which conditioning data to use at each location.

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::sample::Sample;

/// Use at most `max_samples` of the nearest samples within `radius`
/// (GSLIB `ndmax` and an isotropic `radius`, as in `KT3D.FOR`; F-012).
///
/// Distance is Euclidean. GSLIB's anisotropic search ellipsoid, `ndmin`,
/// octants and super-blocks are specified in `docs/spec/search.md` and not
/// implemented yet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Search {
    max_samples: usize,
    radius: f64,
}

impl Search {
    /// # Errors
    /// [`Error::InvalidInput`] if `max_samples` is zero or `radius` is not
    /// positive (`f64::INFINITY` means unlimited).
    pub fn new(max_samples: usize, radius: f64) -> Result<Self> {
        if max_samples == 0 {
            return Err(Error::invalid("search", "max_samples must be at least 1"));
        }
        if radius.is_nan() || radius <= 0.0 {
            return Err(Error::invalid(
                "search",
                format!("radius must be positive, got {radius}"),
            ));
        }
        Ok(Self {
            max_samples,
            radius,
        })
    }

    /// Every sample, however many and however far (cost grows quickly).
    #[must_use]
    pub const fn unlimited() -> Self {
        Self {
            max_samples: usize::MAX,
            radius: f64::INFINITY,
        }
    }

    /// The nearest samples to `target` allowed by this neighbourhood, in no
    /// particular order.
    #[must_use]
    pub fn select(&self, samples: &[Sample], target: Point3) -> Vec<Sample> {
        let mut within: Vec<(f64, Sample)> = samples
            .iter()
            .map(|&s| (s.location().distance(target), s))
            .filter(|&(d, _)| d <= self.radius)
            .collect();
        if within.len() > self.max_samples {
            within.select_nth_unstable_by(self.max_samples, |a, b| a.0.total_cmp(&b.0));
            within.truncate(self.max_samples);
        }
        within.into_iter().map(|(_, s)| s).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(x: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), x).unwrap()
    }

    fn selected(search: Search, xs: &[f64], at: f64) -> Vec<f64> {
        let samples: Vec<Sample> = xs.iter().map(|&x| s(x)).collect();
        let mut v: Vec<f64> = search
            .select(&samples, Point3::new(at, 0.0, 0.0))
            .iter()
            .map(Sample::value)
            .collect();
        v.sort_by(f64::total_cmp);
        v
    }

    #[test]
    fn keeps_the_nearest() {
        let search = Search::new(2, f64::INFINITY).unwrap();
        assert_eq!(
            selected(search, &[0.0, 10.0, 4.0, 6.0, -3.0], 5.0),
            vec![4.0, 6.0]
        );
    }

    #[test]
    fn respects_radius() {
        let search = Search::new(10, 2.0).unwrap();
        assert_eq!(selected(search, &[0.0, 4.0, 6.0, 7.5], 5.0), vec![4.0, 6.0]);
    }

    #[test]
    fn unlimited_keeps_everything() {
        assert_eq!(
            selected(Search::unlimited(), &[3.0, 1.0, 2.0], 0.0),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn rejects_invalid() {
        assert!(Search::new(0, 1.0).is_err(), "zero samples");
        assert!(Search::new(1, 0.0).is_err(), "zero radius");
        assert!(Search::new(1, f64::NAN).is_err(), "NaN radius");
    }
}
