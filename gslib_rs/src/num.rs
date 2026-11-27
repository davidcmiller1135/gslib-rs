//! Numeric conversions.
//!
//! `as` casts are denied crate-wide (`clippy::as_conversions`). The few
//! conversions numerical code genuinely needs live here, each with its
//! precondition stated, so lossy behaviour is never accidental.

/// Convert a count or index to `f64`.
///
/// Exact for values up to 2^53, far beyond any array this crate can hold.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "counts are < 2^53, so the conversion is exact"
)]
pub(crate) const fn to_f64(n: usize) -> f64 {
    n as f64
}

/// `floor(x)` as an index, or `None` if `x` is negative, NaN, infinite, or too
/// large for `usize`.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "range is checked first, so the cast is exact"
)]
pub(crate) fn floor_index(x: f64) -> Option<usize> {
    let f = x.floor();
    // 2^53: every integer below this is exactly representable and fits usize
    // on 64-bit targets. Bins/indices are far smaller in practice.
    let limit = 9_007_199_254_740_992.0_f64;
    // Float-to-int `as` saturates (never UB), so evaluating it eagerly is safe.
    (f.is_finite() && (0.0..limit).contains(&f)).then_some(f as usize)
}

/// Mean of a slice, or `None` if it is empty.
#[must_use]
pub(crate) fn mean(values: &[f64]) -> Option<f64> {
    (!values.is_empty()).then_some(values.iter().sum::<f64>() / to_f64(values.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_index_accepts_non_negative_finite() {
        assert_eq!(floor_index(0.0), Some(0));
        assert_eq!(floor_index(2.999), Some(2));
        assert_eq!(floor_index(-0.0), Some(0));
    }

    #[test]
    fn floor_index_rejects_negative_nan_and_infinite() {
        assert_eq!(floor_index(-0.5), None);
        assert_eq!(floor_index(f64::NAN), None);
        assert_eq!(floor_index(f64::INFINITY), None);
        assert_eq!(floor_index(1e300), None);
    }

    #[test]
    fn mean_of_empty_is_none() {
        assert_eq!(mean(&[]), None);
        assert_eq!(mean(&[1.0, 2.0, 6.0]), Some(3.0));
    }
}
