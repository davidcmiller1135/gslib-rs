//! Comparing gslib-rs results with GSLIB's (R7, D-010).
//!
//! A value passes when `|ours − gslib| ≤ relative · scale + absolute`, where
//! `scale` is chosen per quantity by the caller (for example, the sill for
//! variogram values, since GSLIB computes γ as `sill − C(h)` and so carries
//! absolute error relative to the sill). Cases may override the defaults per
//! quantity, with a justification, or declare a quantity an expected
//! deviation, which is then not compared.

use oracle_case::Case;

/// A tolerance for one quantity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    /// Maximum error relative to the comparison scale.
    pub relative: f64,
    /// Maximum absolute error, added to the relative bound.
    pub absolute: f64,
}

/// The default for values computed by direct formulas (D-010).
pub const DIRECT: Tolerance = Tolerance {
    relative: 1e-12,
    absolute: 0.0,
};

/// The largest error seen for one quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct Worst {
    /// Largest `|Δ| / scale`.
    pub scaled_error: f64,
    /// Where it occurred.
    pub location: String,
}

/// Accumulates comparisons for one fixture and reports failures.
#[derive(Debug)]
pub struct Comparison<'a> {
    case: &'a Case,
    failures: Vec<String>,
    worst: Vec<(String, Worst)>,
}

impl<'a> Comparison<'a> {
    /// Start comparing values for `case`.
    #[must_use]
    pub const fn new(case: &'a Case) -> Self {
        Self {
            case,
            failures: Vec::new(),
            worst: Vec::new(),
        }
    }

    /// Whether the case declares `quantity` an expected deviation.
    #[must_use]
    pub fn is_expected_deviation(&self, quantity: &str) -> bool {
        self.case
            .expected_deviations
            .iter()
            .any(|d| d.quantity == quantity)
    }

    /// The tolerance for `quantity`: the case's override, or `default`.
    #[must_use]
    pub fn tolerance(&self, quantity: &str, default: Tolerance) -> Tolerance {
        self.case
            .tolerances
            .iter()
            .find(|t| t.quantity == quantity)
            .map_or(default, |t| Tolerance {
                relative: t.relative.unwrap_or(0.0),
                absolute: t.absolute.unwrap_or(0.0),
            })
    }

    /// Compare one value. `scale` must be positive.
    pub fn value(
        &mut self,
        quantity: &str,
        ours: f64,
        gslib: f64,
        scale: f64,
        default: Tolerance,
        location: impl Fn() -> String,
    ) {
        if self.is_expected_deviation(quantity) {
            return;
        }
        let tol = self.tolerance(quantity, default);
        let error = (ours - gslib).abs();
        let scaled = error / scale;
        let ok = error <= tol.relative.mul_add(scale, tol.absolute)
            && ours.is_finite() == gslib.is_finite();
        let index = self
            .worst
            .iter()
            .position(|(q, _)| q == quantity)
            .unwrap_or_else(|| {
                self.worst.push((
                    quantity.to_owned(),
                    Worst {
                        scaled_error: -1.0,
                        location: String::new(),
                    },
                ));
                self.worst.len() - 1
            });
        if let Some((_, slot)) = self.worst.get_mut(index)
            && (scaled > slot.scaled_error || scaled.is_nan())
        {
            slot.scaled_error = scaled;
            slot.location = location();
        }
        if !ok {
            self.failures.push(format!(
                "{quantity} at {}: ours {ours:e}, GSLIB {gslib:e}, |Δ|/scale {scaled:e} > {:e} (+{:e})",
                location(),
                tol.relative,
                tol.absolute
            ));
        }
    }

    /// Record a structural mismatch (shape, counts, labels), which must
    /// match exactly.
    pub fn structure(&mut self, ok: bool, message: impl Fn() -> String) {
        if !ok {
            self.failures.push(format!("structure: {}", message()));
        }
    }

    /// The largest scaled error per quantity, for calibration.
    #[must_use]
    pub fn worst(&self) -> &[(String, Worst)] {
        &self.worst
    }

    /// Every failure, or nothing if the fixture matches.
    #[must_use]
    pub fn failures(&self) -> &[String] {
        &self.failures
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(extra: &str) -> Case {
        let text = format!("schema = 1\nprogram = \"vmodel\"\ndescription = \"test\"\n{extra}");
        oracle_case::parse_case(&text, std::path::Path::new("case.toml")).unwrap()
    }

    #[test]
    fn default_tolerance_is_relative_to_scale() {
        let c = case("");
        let mut cmp = Comparison::new(&c);
        cmp.value("gamma", 1.0 + 5e-13, 1.0, 1.0, DIRECT, String::new);
        assert!(cmp.failures().is_empty(), "within 1e-12");
        cmp.value("gamma", 1.0 + 5e-12, 1.0, 1.0, DIRECT, || {
            "lag 3".to_owned()
        });
        assert_eq!(cmp.failures().len(), 1);
        assert!(cmp.failures()[0].contains("lag 3"), "{:?}", cmp.failures());
        assert_eq!(cmp.worst()[0].1.location, "lag 3");
    }

    #[test]
    fn overrides_and_expected_deviations_apply() {
        let c = case(
            "[[tolerances]]\nquantity = \"gamma\"\nrelative = 1e-9\njustification = \"j\"\n\
             [[expected_deviations]]\nquantity = \"cov\"\ndivergence = \"d\"\nreason = \"r\"\n",
        );
        let mut cmp = Comparison::new(&c);
        cmp.value("gamma", 1.0 + 5e-10, 1.0, 1.0, DIRECT, String::new);
        cmp.value("cov", 0.0, 999.0, 1.0, DIRECT, String::new);
        assert!(cmp.failures().is_empty(), "{:?}", cmp.failures());
    }

    #[test]
    fn nan_never_passes() {
        let c = case("");
        let mut cmp = Comparison::new(&c);
        cmp.value("gamma", f64::NAN, 1.0, 1.0, DIRECT, String::new);
        assert_eq!(cmp.failures().len(), 1);
    }
}
