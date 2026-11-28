//! Trimming limits, non-finite policies and the sentinel guard on the
//! command line (`docs/spec/missing-values.md`, D-020, D-021).

use clap::{Args, ValueEnum};
use gslib_rs::io::{Extraction, GslibData, OnInfinite, OnNan, Trimming};

/// GSLIB's usual missing-value codes, which the guard refuses to treat as
/// data without being told (D-021). `lusim`'s −1 is deliberately absent.
const SENTINELS: [f64; 2] = [-999.0, -99.0];

/// Values this large are GSLIB's "no limit" or missing codes (D-021).
const HUGE: f64 = 1.0e21;

/// How the data file's value column is trimmed.
#[derive(Debug, Clone, Copy, Args)]
pub(crate) struct MissingArgs {
    /// Values below this are missing (GSLIB's `tmin`; e.g. -998 for -999 codes).
    #[arg(long, allow_negative_numbers = true)]
    tmin: Option<f64>,
    /// Values at or above this are missing (GSLIB's `tmax`).
    #[arg(long, allow_negative_numbers = true)]
    tmax: Option<f64>,
    /// No trimming: values such as -999 are real data.
    #[arg(long, conflicts_with_all = ["tmin", "tmax"])]
    no_trimming: bool,
    /// What a NaN value means.
    #[arg(long, value_enum, default_value_t = NanPolicy::Error)]
    nan: NanPolicy,
    /// What an infinite value means.
    #[arg(long, value_enum, default_value_t = InfPolicy::Error)]
    inf: InfPolicy,
}

/// `--nan`.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub(crate) enum NanPolicy {
    /// Fail, naming the row.
    Error,
    /// Treat it as missing.
    Missing,
}

/// `--inf`.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub(crate) enum InfPolicy {
    /// Fail, naming the row.
    Error,
    /// Treat it as missing.
    Missing,
    /// Replace it with the largest (or smallest) finite value kept.
    Clamp,
}

impl MissingArgs {
    /// Whether the user said anything about trimming.
    const fn limits_given(&self) -> bool {
        self.no_trimming || self.tmin.is_some() || self.tmax.is_some()
    }

    /// The library policy. A missing limit is unbounded on that side.
    pub(crate) fn extraction(&self) -> gslib_rs::Result<Extraction> {
        let trimming = if self.tmin.is_none() && self.tmax.is_none() {
            Trimming::NONE
        } else {
            Trimming::limits(
                self.tmin.unwrap_or(f64::NEG_INFINITY),
                self.tmax.unwrap_or(f64::INFINITY),
            )?
        };
        Ok(Extraction {
            trimming,
            nan: match self.nan {
                NanPolicy::Error => OnNan::Error,
                NanPolicy::Missing => OnNan::Missing,
            },
            infinite: match self.inf {
                InfPolicy::Error => OnInfinite::Error,
                InfPolicy::Missing => OnInfinite::Missing,
                InfPolicy::Clamp => OnInfinite::ClampToRange,
            },
        })
    }

    /// The sentinel guard (M6): unless trimming was specified, fail on the
    /// first value of `column` (0-based) that looks like a missing-value
    /// code. Non-finite values are left to the `--nan`/`--inf` policies.
    pub(crate) fn guard(&self, data: &GslibData, column: usize) -> Result<(), String> {
        if self.limits_given() {
            return Ok(());
        }
        let values = data.column(column).map_err(|e| e.to_string())?;
        let suspicious = |v: f64| {
            v.is_finite()
                && (SENTINELS.iter().any(|s| s.to_bits() == v.to_bits()) || v.abs() >= HUGE)
        };
        match values.iter().enumerate().find(|(_, v)| suspicious(**v)) {
            None => Ok(()),
            Some((row, v)) => Err(format!(
                "column {} ({}), row {}: {v} looks like a missing-value code; \
                 give --tmin/--tmax to trim such values (e.g. --tmin -998), \
                 or --no-trimming if they are real data",
                column + 1,
                data.names().get(column).map_or("?", String::as_str),
                row + 1,
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gslib_rs::io::parse_gslib;

    fn args(tmin: Option<f64>, tmax: Option<f64>, no_trimming: bool) -> MissingArgs {
        MissingArgs {
            tmin,
            tmax,
            no_trimming,
            nan: NanPolicy::Error,
            inf: InfPolicy::Error,
        }
    }

    fn data(values: &str) -> GslibData {
        let rows: Vec<&str> = values.split_whitespace().collect();
        parse_gslib(&format!("t\n1\nv\n{}\n", rows.join("\n"))).unwrap()
    }

    #[test]
    fn guard_catches_each_sentinel() {
        let none = args(None, None, false);
        for (text, row) in [("1 -999", 2), ("-99 1", 1), ("1 2 1e21", 3), ("-2e21", 1)] {
            let err = none.guard(&data(text), 0).unwrap_err();
            assert!(err.contains(&format!("row {row}:")), "{text}: {err}");
            assert!(err.contains("column 1 (v)"), "{err}");
        }
    }

    #[test]
    fn guard_passes_clean_data_and_leaves_non_finite_to_the_policies() {
        let none = args(None, None, false);
        assert!(
            none.guard(&data("1 -1 -998.5 9.9e20"), 0).is_ok(),
            "no sentinels"
        );
        assert!(none.guard(&data("1 inf NaN"), 0).is_ok(), "non-finite");
    }

    #[test]
    fn guard_is_off_once_trimming_is_specified() {
        let d = data("1 -999");
        assert!(
            args(Some(-998.0), None, false).guard(&d, 0).is_ok(),
            "--tmin"
        );
        assert!(
            args(None, Some(1.0e21), false).guard(&d, 0).is_ok(),
            "--tmax"
        );
        assert!(args(None, None, true).guard(&d, 0).is_ok(), "--no-trimming");
    }

    #[test]
    fn one_limit_leaves_the_other_side_open() {
        let e = args(Some(-998.0), None, false).extraction().unwrap();
        assert_eq!(
            e.trimming.bounds(),
            Some((-998.0, f64::INFINITY)),
            "tmax unbounded"
        );
        let untrimmed = args(None, None, true).extraction().unwrap();
        assert_eq!(untrimmed.trimming, Trimming::NONE);
        assert!(
            args(Some(5.0), Some(1.0), false).extraction().is_err(),
            "reversed limits"
        );
    }
}
