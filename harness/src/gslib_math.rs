//! Ports of GSLIB's normal-distribution approximations.
//!
//! They exist only to *justify* tolerances: the tests below measure how far
//! each strays from the exact function in gslib-rs, and pin the bounds the
//! parity tests rely on (F-007, D-016 in `docs/spec/decisions.md`).

/// Bound on |`gauinv`(p) − Φ⁻¹(p)| for 1e-10 ≤ p ≤ 1 − 1e-10. Measured
/// maximum: 1.501e-8 near p = 0.033 (see `gauinv_error_bound`).
pub const GAUINV_MAX_ERROR: f64 = 2e-8;

/// Bound on |`gcum`(x) − Φ(x)| for all x. Measured maximum: 7.458e-8 near
/// x = 0.717 (see `gcum_error_bound`).
pub const GCUM_MAX_ERROR: f64 = 7.5e-8;

/// GSLIB's `gauinv` (`vendor/gslib-2.0/src/GSLIB/GAUINV.FOR`).
///
/// A rational approximation to Φ⁻¹ (Kennedy & Gentle 1980, as its comment
/// cites), valid for 1e-10 ≤ p ≤ 1 − 1e-10; GSLIB returns ±1e10 outside.
#[must_use]
pub fn gauinv(p: f64) -> f64 {
    const P: [f64; 5] = [
        -0.322_232_431_088,
        -1.0,
        -0.342_242_088_547,
        -0.020_423_121_024_5,
        -0.000_045_364_221_014_8,
    ];
    const Q: [f64; 5] = [
        0.099_348_462_606_0,
        0.588_581_570_495,
        0.531_103_462_366,
        0.103_537_752_850,
        0.003_856_070_063_4,
    ];
    if p < 1e-10 {
        return -1e10;
    }
    if p > 1.0 - 1e-10 {
        return 1e10;
    }
    if (p - 0.5).abs() < f64::EPSILON * 0.5 {
        return 0.0;
    }
    let pp = if p > 0.5 { 1.0 - p } else { p };
    let y = (1.0 / (pp * pp)).ln().sqrt();
    let [p0, p1, p2, p3, p4] = P;
    let [q0, q1, q2, q3, q4] = Q;
    let num = y.mul_add(y.mul_add(y.mul_add(y.mul_add(p4, p3), p2), p1), p0);
    let den = y.mul_add(y.mul_add(y.mul_add(y.mul_add(q4, q3), q2), q1), q0);
    let xp = y + num / den;
    if p > 0.5 { xp } else { -xp }
}

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
    use gslib_rs::transform::{standard_normal_cdf, standard_normal_quantile};

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
    fn gauinv_error_bound() {
        let worst = probabilities()
            .into_iter()
            .map(|p| (gauinv(p) - standard_normal_quantile(p).unwrap()).abs())
            .fold(0.0, f64::max);
        assert!(
            worst <= GAUINV_MAX_ERROR,
            "gauinv error {worst:e} exceeds the bound"
        );
        assert!(
            worst > GAUINV_MAX_ERROR / 2.0,
            "bound {GAUINV_MAX_ERROR:e} is loose: {worst:e}"
        );
    }

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
        assert!((gauinv(0.3) + gauinv(0.7)).abs() < 1e-15, "gauinv odd");
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
