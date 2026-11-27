//! Variogram models: nugget plus nested structures (GSLIB `cova3`).

use std::f64::consts::PI;

use crate::error::{Error, Result};
use crate::geometry::{Angles, Anisotropy, COINCIDENT, Point3, Ranges};

/// The shape of one nested structure, with GSLIB's practical-range scaling.
/// At the range `a`, the exponential reaches 1 − e⁻³ ≈ 95% of its sill and
/// the Gaussian 1 − e⁻⁹ ≈ 99.99%.
///
/// Verified against `vendor/gslib-2.0/src/GSLIB/COVA3.FOR` (F-004 in
/// `docs/spec/decisions.md`).
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Kind {
    /// γ = c·(1.5·h/a − 0.5·(h/a)³) for h < a, else c. GSLIB type 1.
    Spherical,
    /// γ = c·(1 − exp(−3h/a)). GSLIB type 2.
    Exponential,
    /// γ = c·(1 − exp(−(3h/a)²)). GSLIB type 3.
    Gaussian,
    /// γ = c·(h/a)^ω with 0 < ω < 2; no sill. GSLIB type 4.
    ///
    /// GSLIB stores ω in the range slot; here the range `a` is a unit
    /// distance (use 1 to match GSLIB) and its ratios still give anisotropy.
    Power {
        /// The exponent ω, strictly between 0 and 2.
        exponent: f64,
    },
    /// γ = c·(1 − cos(π·h/a)). GSLIB type 5. Only valid in 1-D; GSLIB uses
    /// it along one direction.
    HoleEffect,
}

/// One nested structure: a shape, a sill contribution, and an anisotropy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Structure {
    kind: Kind,
    contribution: f64,
    major: f64,
    anisotropy: Anisotropy,
}

impl Structure {
    /// A structure with the given orientation and ranges.
    ///
    /// The contribution may be negative (needed for cross-variograms);
    /// [`VariogramModel::new`] rejects negative contributions for direct
    /// variograms.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the contribution is not finite, the ranges
    /// or angles are invalid, or a power exponent is outside (0, 2).
    pub fn new(kind: Kind, contribution: f64, angles: Angles, ranges: Ranges) -> Result<Self> {
        if !contribution.is_finite() {
            return Err(Error::invalid("structure contribution", "must be finite"));
        }
        if let Kind::Power { exponent } = kind
            && !(exponent > 0.0 && exponent < 2.0)
        {
            return Err(Error::invalid(
                "power exponent",
                format!("must be strictly between 0 and 2, got {exponent}"),
            ));
        }
        Ok(Self {
            kind,
            contribution,
            major: ranges.major,
            anisotropy: Anisotropy::new(angles, ranges)?,
        })
    }

    /// An isotropic structure.
    ///
    /// # Errors
    /// As for [`Structure::new`].
    pub fn isotropic(kind: Kind, contribution: f64, range: f64) -> Result<Self> {
        Self::new(
            kind,
            contribution,
            Angles::default(),
            Ranges::isotropic(range),
        )
    }

    /// The structure's shape.
    #[must_use]
    pub const fn kind(&self) -> Kind {
        self.kind
    }

    /// The sill contribution `c`.
    #[must_use]
    pub const fn contribution(&self) -> f64 {
        self.contribution
    }

    /// γ for this structure alone at separation `a`–`b`.
    fn gamma(&self, a: Point3, b: Point3) -> f64 {
        let hr = self.anisotropy.scaled_distance(a, b) / self.major;
        let c = self.contribution;
        match self.kind {
            Kind::Spherical => {
                if hr < 1.0 {
                    c * hr * (-0.5 * hr).mul_add(hr, 1.5)
                } else {
                    c
                }
            }
            Kind::Exponential => c * -(-3.0 * hr).exp_m1(),
            Kind::Gaussian => {
                let s = 3.0 * hr;
                c * -(-s * s).exp_m1()
            }
            Kind::Power { exponent } => c * hr.powf(exponent),
            Kind::HoleEffect => c * (1.0 - (PI * hr).cos()),
        }
    }
}

/// A variogram model: a nugget effect plus zero or more nested structures.
#[derive(Debug, Clone, PartialEq)]
pub struct VariogramModel {
    nugget: f64,
    structures: Vec<Structure>,
}

impl VariogramModel {
    /// A model for a direct variogram (one variable with itself).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the nugget or any contribution is negative
    /// or not finite, or the model is identically zero.
    pub fn new(nugget: f64, structures: Vec<Structure>) -> Result<Self> {
        if !(nugget.is_finite() && nugget >= 0.0) {
            return Err(Error::invalid(
                "nugget",
                format!("must be finite and >= 0, got {nugget}"),
            ));
        }
        if let Some(s) = structures.iter().find(|s| s.contribution < 0.0) {
            return Err(Error::invalid(
                "structure contribution",
                format!(
                    "must be >= 0 for a direct variogram, got {}",
                    s.contribution
                ),
            ));
        }
        let total = nugget + structures.iter().map(|s| s.contribution).sum::<f64>();
        if total <= 0.0 {
            return Err(Error::invalid(
                "variogram model",
                "total sill must be positive",
            ));
        }
        Ok(Self { nugget, structures })
    }

    /// A model for a cross-variogram between two variables. Contributions,
    /// including the nugget, may be negative or zero.
    ///
    /// Validity of the joint model (the linear model of coregionalization) is
    /// not checked yet.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the nugget is not finite.
    pub fn cross(nugget: f64, structures: Vec<Structure>) -> Result<Self> {
        if !nugget.is_finite() {
            return Err(Error::invalid("nugget", "must be finite"));
        }
        Ok(Self { nugget, structures })
    }

    /// The nugget effect `c0`.
    #[must_use]
    pub const fn nugget(&self) -> f64 {
        self.nugget
    }

    /// The nested structures.
    #[must_use]
    pub fn structures(&self) -> &[Structure] {
        &self.structures
    }

    /// The semivariogram γ between two locations. γ is exactly 0 for
    /// coincident points; the nugget applies at any non-zero separation.
    #[must_use]
    pub fn gamma(&self, a: Point3, b: Point3) -> f64 {
        if a.distance(b) < COINCIDENT {
            return 0.0;
        }
        self.nugget + self.structures.iter().map(|s| s.gamma(a, b)).sum::<f64>()
    }

    /// The sill, `c0 + Σc`, or `None` if a power structure makes it infinite.
    #[must_use]
    pub fn sill(&self) -> Option<f64> {
        let has_power = self
            .structures
            .iter()
            .any(|s| matches!(s.kind, Kind::Power { .. }));
        (!has_power)
            .then(|| self.nugget + self.structures.iter().map(|s| s.contribution).sum::<f64>())
    }
}

/// A variogram model with a finite sill, and hence a covariance
/// C(h) = sill − γ(h).
///
/// Simple kriging, cokriging, and the sequential simulations need a
/// covariance. GSLIB fakes one for the power model by subtracting γ from a
/// large constant (`PMX = 999.` in `COVA3.FOR`); here that is a type error
/// instead.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceModel {
    model: VariogramModel,
    sill: f64,
}

impl TryFrom<VariogramModel> for CovarianceModel {
    type Error = Error;

    fn try_from(model: VariogramModel) -> Result<Self> {
        let sill = model.sill().ok_or(Error::NoCovariance)?;
        Ok(Self { model, sill })
    }
}

impl CovarianceModel {
    /// The covariance between two locations.
    #[must_use]
    pub fn covariance(&self, a: Point3, b: Point3) -> f64 {
        self.sill - self.model.gamma(a, b)
    }

    /// C(0), the sill.
    #[must_use]
    pub const fn sill(&self) -> f64 {
        self.sill
    }

    /// The underlying variogram model.
    #[must_use]
    pub const fn variogram(&self) -> &VariogramModel {
        &self.model
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const ORIGIN: Point3 = Point3::new(0.0, 0.0, 0.0);
    const EPS: f64 = 1e-12;

    fn at(h: f64) -> Point3 {
        Point3::new(h, 0.0, 0.0)
    }

    fn single(kind: Kind, c: f64, a: f64) -> VariogramModel {
        VariogramModel::new(0.0, vec![Structure::isotropic(kind, c, a).unwrap()]).unwrap()
    }

    #[test]
    fn spherical_known_values() {
        let m = single(Kind::Spherical, 4.0, 2.0);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(1.0)), 2.75, epsilon = EPS);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(2.0)), 4.0, epsilon = EPS);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(10.0)), 4.0, epsilon = EPS);
    }

    #[test]
    fn exponential_reaches_95_percent_at_practical_range() {
        let m = single(Kind::Exponential, 1.0, 10.0);
        assert_abs_diff_eq!(
            m.gamma(ORIGIN, at(10.0)),
            1.0 - (-3.0_f64).exp(),
            epsilon = EPS
        );
    }

    /// Regression: GSLIB's `COVA3.FOR` type 3 is `exp(-(3h/a)*(3h/a))`;
    /// gslib-rs previously used `exp(-3(h/a)²)` (F-004).
    #[test]
    fn gaussian_matches_gslib_cova3() {
        let m = single(Kind::Gaussian, 1.0, 10.0);
        assert_abs_diff_eq!(
            m.gamma(ORIGIN, at(10.0)),
            1.0 - (-9.0_f64).exp(),
            epsilon = EPS
        );
        let h = 4.0_f64;
        let expected = 1.0 - (-(3.0 * h / 10.0) * (3.0 * h / 10.0)).exp();
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(h)), expected, epsilon = EPS);
    }

    #[test]
    fn power_model_has_no_sill_or_covariance() {
        let m = single(Kind::Power { exponent: 1.5 }, 2.0, 1.0);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(4.0)), 2.0 * 8.0, epsilon = EPS);
        assert_eq!(m.sill(), None);
        assert!(
            matches!(CovarianceModel::try_from(m), Err(Error::NoCovariance)),
            "power has no covariance"
        );
    }

    #[test]
    fn hole_effect_peaks_at_range() {
        let m = single(Kind::HoleEffect, 1.0, 5.0);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(5.0)), 2.0, epsilon = EPS);
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(10.0)), 0.0, epsilon = EPS);
    }

    #[test]
    fn nugget_is_discontinuous_at_origin() {
        let m = VariogramModel::new(
            0.3,
            vec![Structure::isotropic(Kind::Spherical, 0.7, 10.0).unwrap()],
        )
        .unwrap();
        assert_abs_diff_eq!(m.gamma(ORIGIN, ORIGIN), 0.0);
        assert!(
            m.gamma(ORIGIN, at(1e-6)) > 0.3,
            "nugget applies at any separation"
        );
        let cov = CovarianceModel::try_from(m).unwrap();
        assert_abs_diff_eq!(cov.sill(), 1.0, epsilon = EPS);
        assert_abs_diff_eq!(cov.covariance(ORIGIN, ORIGIN), 1.0, epsilon = EPS);
        assert_abs_diff_eq!(cov.covariance(ORIGIN, at(20.0)), 0.0, epsilon = EPS);
    }

    #[test]
    fn nested_structures_add() {
        let m = VariogramModel::new(
            0.1,
            vec![
                Structure::isotropic(Kind::Spherical, 0.4, 5.0).unwrap(),
                Structure::isotropic(Kind::Exponential, 0.5, 30.0).unwrap(),
            ],
        )
        .unwrap();
        let h = 3.0;
        let expected = 0.1
            + single(Kind::Spherical, 0.4, 5.0).gamma(ORIGIN, at(h))
            + single(Kind::Exponential, 0.5, 30.0).gamma(ORIGIN, at(h));
        assert_abs_diff_eq!(m.gamma(ORIGIN, at(h)), expected, epsilon = EPS);
    }

    #[test]
    fn anisotropic_structure_uses_scaled_distance() {
        let s = Structure::new(
            Kind::Spherical,
            1.0,
            Angles::default(),
            Ranges {
                major: 10.0,
                minor: 5.0,
                vertical: 1.0,
            },
        )
        .unwrap();
        let m = VariogramModel::new(0.0, vec![s]).unwrap();
        // 5 east (minor axis) is as far as 10 north (major axis): both at range.
        assert_abs_diff_eq!(
            m.gamma(ORIGIN, Point3::new(5.0, 0.0, 0.0)),
            1.0,
            epsilon = EPS
        );
        assert_abs_diff_eq!(
            m.gamma(ORIGIN, Point3::new(0.0, 10.0, 0.0)),
            1.0,
            epsilon = EPS
        );
        assert!(
            m.gamma(ORIGIN, Point3::new(0.0, 5.0, 0.0)) < 1.0,
            "5 north is within range"
        );
    }

    #[test]
    fn rejects_invalid_direct_models() {
        assert!(
            VariogramModel::new(-0.1, vec![]).is_err(),
            "negative nugget"
        );
        assert!(VariogramModel::new(0.0, vec![]).is_err(), "zero model");
        let neg = Structure::isotropic(Kind::Spherical, -1.0, 1.0).unwrap();
        assert!(
            VariogramModel::new(0.5, vec![neg]).is_err(),
            "negative contribution"
        );
        assert!(
            VariogramModel::cross(0.0, vec![neg]).is_ok(),
            "allowed for cross-variograms"
        );
        assert!(
            Structure::isotropic(Kind::Power { exponent: 2.0 }, 1.0, 1.0).is_err(),
            "ω = 2"
        );
    }
}
