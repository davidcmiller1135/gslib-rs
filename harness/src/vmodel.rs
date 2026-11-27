//! `vmodel` output and the mapping from GSLIB's model to gslib-rs's.
//!
//! Facts from `vendor/gslib-2.0/src/VMODEL/VMODEL.FOR` (F-006):
//! - each direction starts with `Model Variogram for Direction: <id>`;
//! - then rows `il h gam ndir cov ro` for `il = 1..=nlag+1`, at offset
//!   `(il − 1) · (xoff, yoff, zoff)`, with
//!   `xoff = sin(azm)·cos(dip)·lag`, `yoff = cos(azm)·cos(dip)·lag`,
//!   `zoff = sin(dip)·lag` and `DEG2RAD = 3.14159265/180.0`;
//! - the row for `il = 1` is written twice: the second evaluates the model
//!   at `0.0001 · offset` but prints the first row's `h` (it is not
//!   recomputed);
//! - the fourth column is the number of directions, not the direction index;
//! - `gam = maxcov − cov`, with `maxcov` the covariance at zero distance.

use gslib_rs::geometry::{Angles, Point3, Ranges};
use gslib_rs::model::{Kind, Structure, VariogramModel};
use oracle_case::vmodel::Params;

/// `DEG2RAD` as declared in `VMODEL.FOR` line 48. The harness evaluates our
/// model at the points `vmodel` actually used, so it needs GSLIB's value.
#[expect(
    clippy::approx_constant,
    reason = "deliberately GSLIB's truncated pi from VMODEL.FOR, not std's PI"
)]
pub const VMODEL_DEG2RAD: f64 = 3.141_592_65 / 180.0;

/// One row of `vmodel` output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    /// Lag index (1-based).
    pub il: u32,
    /// Separation distance as printed.
    pub h: f64,
    /// Semivariogram value.
    pub gamma: f64,
    /// The number of directions (as `vmodel` prints it).
    pub ndir: u32,
    /// Covariance.
    pub cov: f64,
    /// Correlation `cov / maxcov`.
    pub ro: f64,
}

/// One direction's table.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// The id in the header line.
    pub id: u32,
    /// Its rows, in file order.
    pub rows: Vec<Row>,
}

/// Parse a `vmodel` output file.
///
/// # Errors
/// A message naming the first line that does not fit the format.
pub fn parse_output(text: &str) -> Result<Vec<Table>, String> {
    let mut tables: Vec<Table> = Vec::new();
    for (n, line) in text.lines().enumerate().map(|(i, l)| (i + 1, l.trim())) {
        if line.is_empty() {
            continue;
        }
        if let Some(id) = line.strip_prefix("Model Variogram for Direction:") {
            let id = id
                .trim()
                .parse()
                .map_err(|_| format!("line {n}: bad direction id"))?;
            tables.push(Table {
                id,
                rows: Vec::new(),
            });
            continue;
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let [il, h, gamma, ndir, cov, ro] = tokens.as_slice() else {
            return Err(format!(
                "line {n}: expected 6 columns, got {}",
                tokens.len()
            ));
        };
        let real = |t: &str| {
            t.parse::<f64>()
                .map_err(|_| format!("line {n}: {t:?} is not a number"))
        };
        let int = |t: &str| {
            t.parse::<u32>()
                .map_err(|_| format!("line {n}: {t:?} is not an integer"))
        };
        let row = Row {
            il: int(il)?,
            h: real(h)?,
            gamma: real(gamma)?,
            ndir: int(ndir)?,
            cov: real(cov)?,
            ro: real(ro)?,
        };
        tables
            .last_mut()
            .ok_or_else(|| format!("line {n}: data before the first direction header"))?
            .rows
            .push(row);
    }
    Ok(tables)
}

/// The point `vmodel` evaluates for row `index` (0-based, file order) of a
/// direction: the repeated first row is at `0.0001 · offset`.
#[must_use]
pub fn evaluation_point(
    direction: &oracle_case::vmodel::Direction,
    row_index: usize,
    il: u32,
) -> Point3 {
    let (sa, ca) = (direction.azimuth * VMODEL_DEG2RAD).sin_cos();
    let (sd, cd) = (direction.dip * VMODEL_DEG2RAD).sin_cos();
    let step = if row_index == 1 {
        0.0001
    } else {
        f64::from(il.saturating_sub(1))
    };
    let lag = direction.lag * step;
    Point3::new(sa * cd * lag, ca * cd * lag, sd * lag)
}

/// The gslib-rs model equivalent to a `vmodel` case's model.
///
/// GSLIB types 1–5 map to spherical, exponential, Gaussian, power and hole
/// effect. For the power model GSLIB stores ω in `a_hmax` and does not divide
/// distance by it; gslib-rs takes the exponent separately and a unit range, so
/// the ranges become `(1, a_hmin/ω, a_vert/ω)`.
///
/// # Errors
/// A message if gslib-rs rejects the model.
pub fn model(params: &Params) -> Result<VariogramModel, String> {
    let structures = params
        .structures
        .iter()
        .map(|s| {
            let [a1, a2, a3] = s.angles;
            let [r1, r2, r3] = s.ranges;
            let angles = Angles {
                azimuth: a1,
                dip: a2,
                plunge: a3,
            };
            let (kind, ranges) = match s.kind {
                1 => (
                    Kind::Spherical,
                    Ranges {
                        major: r1,
                        minor: r2,
                        vertical: r3,
                    },
                ),
                2 => (
                    Kind::Exponential,
                    Ranges {
                        major: r1,
                        minor: r2,
                        vertical: r3,
                    },
                ),
                3 => (
                    Kind::Gaussian,
                    Ranges {
                        major: r1,
                        minor: r2,
                        vertical: r3,
                    },
                ),
                4 => (
                    Kind::Power { exponent: r1 },
                    Ranges {
                        major: 1.0,
                        minor: r2 / r1,
                        vertical: r3 / r1,
                    },
                ),
                5 => (
                    Kind::HoleEffect,
                    Ranges {
                        major: r1,
                        minor: r2,
                        vertical: r3,
                    },
                ),
                other => return Err(format!("unknown GSLIB structure type {other}")),
            };
            Structure::new(kind, s.cc, angles, ranges).map_err(|e| e.to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    VariogramModel::new(params.nugget, structures).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_tables() {
        let text = "Model Variogram for Direction:  1\n   1   0.0E+000   0.0E+000      2   1.0E+000   1.0E+000\n   1   0.0E+000   2.5E-001      2   7.5E-001   7.5E-001\nModel Variogram for Direction:  2\n   1   0.0E+000   0.0E+000      2   1.0E+000   1.0E+000\n";
        let t = parse_output(text).unwrap();
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].rows.len(), 2);
        assert!((t[0].rows[1].gamma - 0.25).abs() < 1e-15, "{t:?}");
        assert_eq!(t[1].id, 2);
    }

    #[test]
    fn rejects_malformed_rows() {
        assert!(
            parse_output("   1 0.0 0.0 2 1.0\n").is_err(),
            "data before header"
        );
        assert!(
            parse_output("Model Variogram for Direction:  1\n 1 0.0 x 2 1.0 1.0\n").is_err(),
            "not a number"
        );
    }

    #[test]
    fn evaluation_points_follow_vmodel() {
        let d = oracle_case::vmodel::Direction {
            azimuth: 0.0,
            dip: 0.0,
            lag: 2.0,
        };
        assert_eq!(evaluation_point(&d, 0, 1), Point3::new(0.0, 0.0, 0.0));
        assert_eq!(evaluation_point(&d, 1, 1), Point3::new(0.0, 0.0002, 0.0));
        assert_eq!(evaluation_point(&d, 4, 4), Point3::new(0.0, 6.0, 0.0));
    }
}
