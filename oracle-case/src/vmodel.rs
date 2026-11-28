//! `vmodel` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout and limits from `vendor/gslib-2.0/src/VMODEL/VMODEL.FOR`: an output
//! file name; `ndir nlag`; one `azm dip lag` line per direction (at most
//! `MDIR = 10`); `nst c0`; then per structure (at most `MAXNST = 4`) an
//! `it cc ang1 ang2 ang3` line and an `a_hmax a_hmin a_vert` line.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// Most directions `vmodel` accepts (`MDIR`).
pub const MAX_DIRECTIONS: usize = 10;
/// Most nested structures `vmodel` accepts (`MAXNST`).
pub const MAX_STRUCTURES: usize = 4;
/// The output file the recipe asks `vmodel` to write.
pub const OUTPUT_FILE: &str = "vmodel.var";

/// One direction along which the model is tabulated.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Direction {
    /// Azimuth, degrees clockwise from north.
    pub azimuth: f64,
    /// Dip, degrees (negative is downward).
    pub dip: f64,
    /// Lag spacing.
    pub lag: f64,
}

/// One nested structure, exactly as GSLIB's parameter file states it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Structure {
    /// GSLIB type: 1 spherical, 2 exponential, 3 Gaussian, 4 power,
    /// 5 hole effect.
    pub kind: u8,
    /// Sill contribution `cc`.
    pub cc: f64,
    /// `ang1` (azimuth), `ang2` (dip), `ang3` (plunge), degrees.
    pub angles: [f64; 3],
    /// `a_hmax a_hmin a_vert`. For the power model, `a_hmax` is the exponent
    /// ω and the others set anisotropy relative to it.
    pub ranges: [f64; 3],
}

/// Parameters of a `vmodel` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Number of lags per direction (`nlag`; `vmodel` writes `nlag + 1`
    /// points, plus a repeat of the first).
    pub nlag: u32,
    /// Directions to tabulate.
    pub directions: Vec<Direction>,
    /// Nugget effect `c0`.
    pub nugget: f64,
    /// Nested structures.
    pub structures: Vec<Structure>,
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema, or a value is outside
    /// what both `vmodel` and gslib-rs accept: counts above GSLIB's limits,
    /// unknown types, non-finite numbers, non-positive lags or ranges, a
    /// negative nugget or contribution, or a power exponent outside (0, 2).
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let params: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        params.validate()?;
        Ok(params)
    }

    fn validate(&self) -> Result<(), String> {
        if self.nlag == 0 {
            return Err("nlag must be at least 1".to_owned());
        }
        if self.directions.is_empty() || self.directions.len() > MAX_DIRECTIONS {
            return Err(format!("need 1 to {MAX_DIRECTIONS} directions"));
        }
        if self.structures.is_empty() || self.structures.len() > MAX_STRUCTURES {
            return Err(format!(
                "need 1 to {MAX_STRUCTURES} structures (vmodel requires nst >= 1)"
            ));
        }
        if !(self.nugget.is_finite() && self.nugget >= 0.0) {
            return Err("nugget must be finite and >= 0".to_owned());
        }
        for d in &self.directions {
            if ![d.azimuth, d.dip, d.lag].iter().all(|v| v.is_finite()) || d.lag <= 0.0 {
                return Err("directions need finite angles and a positive lag".to_owned());
            }
        }
        for s in &self.structures {
            if !(1..=5).contains(&s.kind) {
                return Err(format!("structure type {} is not 1-5", s.kind));
            }
            let numbers = std::iter::once(s.cc).chain(s.angles).chain(s.ranges);
            if numbers.into_iter().any(|v: f64| !v.is_finite()) || s.cc < 0.0 {
                return Err("structures need finite values and cc >= 0".to_owned());
            }
            if s.ranges.iter().any(|&a| a <= 0.0) {
                return Err("ranges must be positive".to_owned());
            }
            if s.kind == 4 && s.ranges[0] >= 2.0 {
                return Err("power exponent (a_hmax) must be below 2".to_owned());
            }
        }
        Ok(())
    }

    /// The GSLIB parameter file for this case. Numbers use Rust's shortest
    /// round-trip form, which gfortran's list-directed READ parses to the
    /// same `f64`.
    #[must_use]
    pub fn to_par(&self) -> String {
        let mut par = String::new();
        // Writing to a String cannot fail.
        let mut line = |text: String| {
            let _ = writeln!(par, "{text}");
        };
        line("                  Parameters for VMODEL".to_owned());
        line("                  *********************".to_owned());
        line(String::new());
        line("START OF PARAMETERS:".to_owned());
        line(format!("{OUTPUT_FILE}   \\file for variogram output"));
        line(format!(
            "{} {}   \\number of directions and lags",
            self.directions.len(),
            self.nlag
        ));
        for d in &self.directions {
            line(format!(
                "{} {} {}   \\azm, dip, lag distance",
                d.azimuth, d.dip, d.lag
            ));
        }
        line(format!(
            "{} {}   \\nst, nugget effect",
            self.structures.len(),
            self.nugget
        ));
        for s in &self.structures {
            let [a1, a2, a3] = s.angles;
            let [r1, r2, r3] = s.ranges;
            line(format!(
                "{} {} {a1} {a2} {a3}   \\it,cc,ang1,ang2,ang3",
                s.kind, s.cc
            ));
            line(format!("{r1} {r2} {r3}   \\a_hmax, a_hmin, a_vert"));
        }
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(text: &str) -> toml::Table {
        text.parse().unwrap()
    }

    const PARAMS: &str = r"
nlag = 5
nugget = 0.1
directions = [{ azimuth = 0.0, dip = 0.0, lag = 1.5 }]
structures = [{ kind = 1, cc = 0.9, angles = [30.0, 0.0, 0.0], ranges = [10.0, 5.0, 2.0] }]
";

    #[test]
    fn par_file_layout() {
        let par = Params::from_table(&table(PARAMS)).unwrap().to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert_eq!(lines[3], "START OF PARAMETERS:");
        assert!(lines[4].starts_with("vmodel.var "), "{par}");
        assert!(lines[5].starts_with("1 5 "), "{par}");
        assert!(lines[6].starts_with("0 0 1.5 "), "{par}");
        assert!(lines[7].starts_with("1 0.1 "), "{par}");
        assert!(lines[8].starts_with("1 0.9 30 0 0 "), "{par}");
        assert!(lines[9].starts_with("10 5 2 "), "{par}");
        assert_eq!(lines.len(), 10);
    }

    #[test]
    fn rejects_values_gslib_or_gslib_rs_cannot_take() {
        let bad =
            |from: &str, to: &str| Params::from_table(&table(&PARAMS.replace(from, to))).is_err();
        assert!(bad("nlag = 5", "nlag = 0"), "nlag 0");
        assert!(bad("kind = 1", "kind = 6"), "type 6");
        assert!(bad("ranges = [10.0", "ranges = [-10.0"), "negative range");
        assert!(bad("lag = 1.5", "lag = 0.0"), "zero lag");
        assert!(
            bad(
                "kind = 1, cc = 0.9, angles = [30.0, 0.0, 0.0], ranges = [10.0",
                "kind = 4, cc = 0.9, angles = [30.0, 0.0, 0.0], ranges = [2.0"
            ),
            "power 2"
        );
        assert!(
            bad("nugget = 0.1", "nugget = 0.1\nextra = 1"),
            "unknown field"
        );
    }
}
