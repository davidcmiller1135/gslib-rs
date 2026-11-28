//! GSLIB's normal-distribution approximations: `gcum` ported here, `gauinv`
//! re-exported from the library.
//!
//! They exist only to *justify* tolerances: the tests below measure how far
//! each strays from the exact function in gslib-rs, and pin the bounds the
//! parity tests rely on (F-007, D-016 in `docs/spec/decisions.md`).

/// GSLIB's `gauinv` and its error bound now live in the library
/// (`docs/spec/rng.md` RN6), where the exact-stream mode uses them.
pub use gslib_rs::transform::{GSLIB_GAUINV_MAX_ERROR as GAUINV_MAX_ERROR, gslib_gauinv as gauinv};

/// Bound on |`gcum`(x) − Φ(x)| for all x. Measured maximum: 7.458e-8 near
/// x = 0.717 (see `gcum_error_bound`).
pub const GCUM_MAX_ERROR: f64 = 7.5e-8;

/// GSLIB's `gcum` (`vendor/gslib-2.0/src/GSLIB/GCUM.FOR`): a polynomial
/// approximation to Φ. For |x| > 6 it returns exactly 0 or 1.
#[must_use]
pub fn gcum(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / 0.231_641_9_f64.mul_add(z, 1.0);
    let poly = t * t.mul_add(
        t.mul_add(
            t.mul_add(t.mul_add(1.330_274_429, -1.821_255_978), 1.781_477_937),
            -0.356_563_782,
        ),
        0.319_381_53,
    );
    let e2 = if z <= 6.0 {
        (-z * z / 2.0).exp() * 0.398_942_280_3
    } else {
        0.0
    };
    let g = e2.mul_add(-poly, 1.0);
    if x >= 0.0 { g } else { 1.0 - g }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gslib_rs::transform::standard_normal_cdf;

    #[test]
    fn gcum_error_bound() {
        let worst = (0..=400_000)
            .map(|k| -8.0 + 16.0 * f64::from(k) / 400_000.0)
            .map(|x| (gcum(x) - standard_normal_cdf(x)).abs())
            .fold(0.0, f64::max);
        assert!(
            worst <= GCUM_MAX_ERROR,
            "gcum error {worst:e} exceeds the bound"
        );
        assert!(
            worst > GCUM_MAX_ERROR / 2.0,
            "bound {GCUM_MAX_ERROR:e} is loose: {worst:e}"
        );
    }

    #[test]
    fn ports_are_symmetric() {
        assert!(
            (gcum(1.3) + gcum(-1.3) - 1.0).abs() < 1e-15,
            "gcum symmetric"
        );
        assert!(
            (gcum(7.0) - 1.0).abs() < f64::EPSILON,
            "gcum saturates past 6"
        );
    }
}
