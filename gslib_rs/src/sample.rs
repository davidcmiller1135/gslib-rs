//! A measured value at a location.

use crate::error::{Error, Result};
use crate::geometry::Point3;

/// One datum: a value measured at a location. Both are guaranteed finite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample {
    location: Point3,
    value: f64,
}

impl Sample {
    /// A sample, rejecting non-finite coordinates or values.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if any coordinate or the value is NaN or infinite.
    pub fn new(location: Point3, value: f64) -> Result<Self> {
        if !location.is_finite() {
            return Err(Error::invalid(
                "sample location",
                format!("{location:?} is not finite"),
            ));
        }
        if !value.is_finite() {
            return Err(Error::invalid(
                "sample value",
                format!("{value} is not finite"),
            ));
        }
        Ok(Self { location, value })
    }

    /// Where the sample was taken.
    #[must_use]
    pub const fn location(&self) -> Point3 {
        self.location
    }

    /// The measured value.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }
}

/// Arithmetic mean of the sample values, or `None` if there are none.
#[must_use]
pub fn mean_value(samples: &[Sample]) -> Option<f64> {
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    crate::num::mean(&values)
}

/// Pair up locations and values into samples.
///
/// # Errors
/// [`Error::LengthMismatch`] if the lengths differ, or [`Error::InvalidInput`]
/// if any entry is not finite.
pub fn zip_samples(locations: &[Point3], values: &[f64]) -> Result<Vec<Sample>> {
    crate::error::check_lengths("locations/values", locations.len(), values.len())?;
    locations
        .iter()
        .zip(values)
        .map(|(&p, &v)| Sample::new(p, v))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_finite() {
        assert!(
            Sample::new(Point3::new(f64::NAN, 0.0, 0.0), 1.0).is_err(),
            "NaN x"
        );
        assert!(
            Sample::new(Point3::default(), f64::INFINITY).is_err(),
            "infinite value"
        );
    }

    #[test]
    fn zip_checks_lengths() {
        assert!(
            matches!(
                zip_samples(&[Point3::default()], &[]),
                Err(Error::LengthMismatch { .. })
            ),
            "length mismatch must be reported"
        );
    }
}
