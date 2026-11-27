//! GSLIB (simplified Geo-EAS) ASCII files.
//!
//! [`parse_gslib`] and [`format_gslib`] are pure; [`read_gslib`] and
//! [`write_gslib`] are thin wrappers that touch the filesystem.
//!
//! # Missing values (`docs/spec/missing-values.md`)
//!
//! [`GslibData`] holds a file's values exactly as written, sentinel codes
//! such as −999 included (M1). Missingness is decided only when values are
//! taken out, by an [`Extraction`] policy (D-018):
//!
//! 1. **Non-finite values first** (M9): NaN and ±Inf are errors by default,
//!    or become missing ([`OnNan`], [`OnInfinite`]). Infinities can instead
//!    be clamped to the largest or smallest finite value kept by the
//!    trimming limits (D-020).
//! 2. **Then trimming** (M4): a value is present iff `tmin ≤ v < tmax`, the
//!    half-open test most GSLIB programs use (F-009, D-019).
//!
//! Each value is decided on its own (M7); coordinates are never trimmed
//! (D-022).
//!
//! **Divergence from GSLIB** (ROADMAP §2): GSLIB passes NaN through its
//! trimming tests as data and trims ±Inf by its limits (F-010). Here both are
//! errors unless a policy says otherwise, and the non-finite policy applies
//! even where finite limits would have trimmed the value.

use std::fmt::Write as _;
use std::path::Path;

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::sample::Sample;

/// Trimming limits: a value `v` is present iff `tmin ≤ v < tmax` (D-019).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trimming {
    limits: Option<(f64, f64)>,
}

impl Trimming {
    /// No limits: every finite value is present.
    pub const NONE: Self = Self { limits: None };

    /// Limits `tmin ≤ v < tmax`. Either limit may be infinite.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if a limit is NaN or `tmin ≥ tmax` (an empty
    /// interval would silently discard every value).
    pub fn limits(tmin: f64, tmax: f64) -> Result<Self> {
        if tmin.is_nan() || tmax.is_nan() {
            return Err(Error::invalid(
                "trimming limits",
                format!("tmin {tmin} and tmax {tmax} must not be NaN"),
            ));
        }
        if tmin >= tmax {
            return Err(Error::invalid(
                "trimming limits",
                format!("tmin {tmin} must be less than tmax {tmax}"),
            ));
        }
        Ok(Self {
            limits: Some((tmin, tmax)),
        })
    }

    /// The limits `(tmin, tmax)`, or `None` for no trimming.
    #[must_use]
    pub const fn bounds(&self) -> Option<(f64, f64)> {
        self.limits
    }

    /// Whether `v` lies within the limits. NaN never does, except with
    /// [`Trimming::NONE`]; [`GslibData::extract`] deals with NaN before
    /// trimming.
    #[must_use]
    pub fn contains(&self, v: f64) -> bool {
        self.limits.is_none_or(|(tmin, tmax)| tmin <= v && v < tmax)
    }
}

/// What extraction does with a NaN value (D-020).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OnNan {
    /// Fail with [`Error::NonFinite`].
    #[default]
    Error,
    /// Treat it as missing.
    Missing,
}

/// What extraction does with an infinite value (D-020).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OnInfinite {
    /// Fail with [`Error::NonFinite`].
    #[default]
    Error,
    /// Treat it as missing.
    Missing,
    /// Replace +Inf with the column's largest finite value within the
    /// trimming limits, and −Inf with its smallest. Fails if there is none.
    ClampToRange,
}

/// How values become present or missing when taken out of a [`GslibData`]:
/// the non-finite policies first, then trimming (M9).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extraction {
    /// Trimming limits.
    pub trimming: Trimming,
    /// What to do with NaN.
    pub nan: OnNan,
    /// What to do with ±Inf.
    pub infinite: OnInfinite,
}

impl Extraction {
    /// The given trimming, with non-finite values as errors (the defaults).
    #[must_use]
    pub const fn new(trimming: Trimming) -> Self {
        Self {
            trimming,
            nan: OnNan::Error,
            infinite: OnInfinite::Error,
        }
    }
}

/// The samples of a data set whose values are present, with the rows they
/// came from.
#[derive(Debug, Clone, PartialEq)]
pub struct PresentSamples {
    samples: Vec<Sample>,
    rows: Vec<usize>,
    total_rows: usize,
}

impl PresentSamples {
    /// The present samples, in row order.
    #[must_use]
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    /// The 0-based row of each sample.
    #[must_use]
    pub fn rows(&self) -> &[usize] {
        &self.rows
    }

    /// The number of rows in the data set, present or not.
    #[must_use]
    pub const fn total_rows(&self) -> usize {
        self.total_rows
    }

    /// The samples, without their rows.
    #[must_use]
    pub fn into_samples(self) -> Vec<Sample> {
        self.samples
    }

    /// Spread one result per sample back over the data set's rows: one
    /// entry per row, `None` for rows whose value was missing (M3).
    ///
    /// # Errors
    /// [`Error::LengthMismatch`] unless there is one result per sample.
    pub fn align<T: Copy>(&self, results: &[T]) -> Result<Vec<Option<T>>> {
        crate::error::check_lengths("results/present samples", results.len(), self.samples.len())?;
        let mut aligned = vec![None; self.total_rows];
        for (&row, &r) in self.rows.iter().zip(results) {
            let slot = aligned.get_mut(row).ok_or_else(|| {
                Error::invalid(
                    "present samples",
                    format!("row {row} is outside the data set"),
                )
            })?;
            *slot = Some(r);
        }
        Ok(aligned)
    }
}

/// The value written in place of a missing result (D-023). Always finite,
/// so it reads back as a number that trimming limits can exclude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MissingCode(f64);

impl MissingCode {
    /// −999: GSLIB's `UNEST` in `kt3d`, `kb2d`, `cokb3d` and `postsim`, and
    /// what `nscore` writes for trimmed rows (F-009). gslib-rs also uses it
    /// for `lusim`, where GSLIB writes −1 (D-023, ROADMAP §2).
    pub const MINUS_999: Self = Self(-999.0);
    /// −99: GSLIB's `UNEST` in `sgsim`, the `sisim` family and `postik`
    /// (F-009).
    pub const MINUS_99: Self = Self(-99.0);

    /// A code of the caller's choosing.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `code` is NaN or infinite.
    pub fn new(code: f64) -> Result<Self> {
        if code.is_finite() {
            Ok(Self(code))
        } else {
            Err(Error::invalid(
                "missing-value code",
                format!("{code} is not finite"),
            ))
        }
    }

    /// The code's value.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0
    }
}

/// The contents of a GSLIB file: a title, variable names, and a table of
/// values with one column per variable.
#[derive(Debug, Clone, PartialEq)]
pub struct GslibData {
    title: String,
    names: Vec<String>,
    /// Row-major, `names.len()` values per row.
    values: Vec<f64>,
}

impl GslibData {
    /// A data set from row-major `values`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if there are no variables, or `values.len()`
    /// is not a multiple of the number of variables.
    pub fn new(title: impl Into<String>, names: Vec<String>, values: Vec<f64>) -> Result<Self> {
        if names.is_empty() {
            return Err(Error::invalid("GSLIB data", "needs at least one variable"));
        }
        if !values.len().is_multiple_of(names.len()) {
            return Err(Error::invalid(
                "GSLIB data",
                format!(
                    "{} values do not fill rows of {} variables",
                    values.len(),
                    names.len()
                ),
            ));
        }
        Ok(Self {
            title: title.into(),
            names,
            values,
        })
    }

    /// The title line.
    #[must_use]
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Variable names, one per column.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Number of variables (columns).
    #[must_use]
    pub const fn ncols(&self) -> usize {
        self.names.len()
    }

    /// Number of records (rows).
    #[must_use]
    pub fn nrows(&self) -> usize {
        self.values.len().checked_div(self.ncols()).unwrap_or(0)
    }

    /// Iterate over rows.
    pub fn rows(&self) -> impl Iterator<Item = &[f64]> {
        self.values.chunks_exact(self.ncols())
    }

    /// One column by 0-based index.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the column does not exist.
    pub fn column(&self, index: usize) -> Result<Vec<f64>> {
        self.check_column(index)?;
        Ok(self.rows().filter_map(|r| r.get(index).copied()).collect())
    }

    /// Points from three coordinate columns (0-based).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if a column does not exist.
    pub fn points(&self, x_col: usize, y_col: usize, z_col: usize) -> Result<Vec<Point3>> {
        let (xs, ys, zs) = (
            self.column(x_col)?,
            self.column(y_col)?,
            self.column(z_col)?,
        );
        Ok(xs
            .into_iter()
            .zip(ys)
            .zip(zs)
            .map(|((x, y), z)| Point3::new(x, y, z))
            .collect())
    }

    /// One column by 0-based index, with missingness decided by `policy`:
    /// one entry per row, `None` where the value is missing (M3, M5, M9).
    ///
    /// # Errors
    /// - [`Error::InvalidInput`] if the column does not exist, or if
    ///   [`OnInfinite::ClampToRange`] has no finite value to clamp to.
    /// - [`Error::NonFinite`] for the first NaN or ±Inf that the policy
    ///   rejects, naming its 1-based row and column.
    pub fn extract(&self, index: usize, policy: &Extraction) -> Result<Vec<Option<f64>>> {
        let raw = self.column(index)?;
        let range =
            if policy.infinite == OnInfinite::ClampToRange && raw.iter().any(|v| v.is_infinite()) {
                Some(clamp_range(&raw, policy.trimming).ok_or_else(|| {
                    Error::invalid(
                        "column",
                        format!(
                            "column {}: no finite value within the trimming limits \
                         to clamp infinite values to",
                            index + 1
                        ),
                    )
                })?)
            } else {
                None
            };
        raw.iter()
            .enumerate()
            .map(|(row, &v)| {
                let non_finite = || Error::NonFinite {
                    row: row + 1,
                    column: index + 1,
                    value: v,
                };
                let v = if v.is_nan() {
                    match policy.nan {
                        OnNan::Error => return Err(non_finite()),
                        OnNan::Missing => return Ok(None),
                    }
                } else if v.is_infinite() {
                    match (policy.infinite, range) {
                        (OnInfinite::Missing, _) => return Ok(None),
                        (OnInfinite::ClampToRange, Some((lo, hi))) => {
                            if v > 0.0 {
                                hi
                            } else {
                                lo
                            }
                        }
                        _ => return Err(non_finite()),
                    }
                } else {
                    v
                };
                Ok(policy.trimming.contains(v).then_some(v))
            })
            .collect()
    }

    /// The samples whose values are present under `policy`, from three
    /// coordinate columns and a value column (0-based), with their rows.
    ///
    /// Coordinates are never trimmed (D-022): every coordinate must be
    /// finite, whether or not its row's value is present.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if a column does not exist;
    /// [`Error::NonFinite`] for a non-finite coordinate, or a non-finite value
    /// the policy rejects (see [`GslibData::extract`]).
    pub fn samples(
        &self,
        x_col: usize,
        y_col: usize,
        z_col: usize,
        value_col: usize,
        policy: &Extraction,
    ) -> Result<PresentSamples> {
        let values = self.extract(value_col, policy)?;
        for col in [x_col, y_col, z_col] {
            if let Some((row, &value)) = self
                .column(col)?
                .iter()
                .enumerate()
                .find(|(_, v)| !v.is_finite())
            {
                return Err(Error::NonFinite {
                    row: row + 1,
                    column: col + 1,
                    value,
                });
            }
        }
        let mut samples = Vec::new();
        let mut rows = Vec::new();
        for (row, (point, value)) in self
            .points(x_col, y_col, z_col)?
            .into_iter()
            .zip(values)
            .enumerate()
        {
            if let Some(v) = value {
                samples.push(Sample::new(point, v)?);
                rows.push(row);
            }
        }
        Ok(PresentSamples {
            samples,
            rows,
            total_rows: self.nrows(),
        })
    }

    /// A copy with one more column of results, one per row, `None` written
    /// as `code` (M3, M8). Present values are written unchanged.
    ///
    /// # Errors
    /// - [`Error::LengthMismatch`] unless there is one entry per row.
    /// - [`Error::InvalidInput`] if a present value is not finite, or equals
    ///   `code`: it would read back as missing.
    pub fn with_optional_column(
        &self,
        name: impl Into<String>,
        values: &[Option<f64>],
        code: MissingCode,
    ) -> Result<Self> {
        let name = name.into();
        crate::error::check_lengths("new column/rows", values.len(), self.nrows())?;
        let written = values
            .iter()
            .enumerate()
            .map(|(row, v)| match *v {
                None => Ok(code.value()),
                Some(x) if !x.is_finite() => Err(Error::NonFinite {
                    row: row + 1,
                    column: self.ncols() + 1,
                    value: x,
                }),
                #[expect(
                    clippy::float_cmp,
                    reason = "an exact match is what makes a value indistinguishable from the code"
                )]
                Some(x) if x == code.value() => Err(Error::invalid(
                    "missing-value code",
                    format!(
                        "row {}: the result {x} of {name:?} equals the code; choose another code",
                        row + 1
                    ),
                )),
                Some(x) => Ok(x),
            })
            .collect::<Result<Vec<f64>>>()?;
        self.with_column(name, &written)
    }

    /// A copy with one more column appended.
    ///
    /// # Errors
    /// [`Error::LengthMismatch`] if `values` does not have one entry per row.
    pub fn with_column(&self, name: impl Into<String>, values: &[f64]) -> Result<Self> {
        crate::error::check_lengths("new column/rows", values.len(), self.nrows())?;
        let mut names = self.names.clone();
        names.push(name.into());
        let data = self
            .rows()
            .zip(values)
            .flat_map(|(row, &v)| row.iter().copied().chain(std::iter::once(v)))
            .collect();
        Self::new(self.title.clone(), names, data)
    }

    /// A copy with a different title.
    #[must_use]
    pub fn with_title(mut self, title: impl Into<String>) -> Self {
        self.title = title.into();
        self
    }

    fn check_column(&self, index: usize) -> Result<()> {
        if index < self.ncols() {
            Ok(())
        } else {
            Err(Error::invalid(
                "column",
                format!(
                    "column {} requested but the file has {}",
                    index + 1,
                    self.ncols()
                ),
            ))
        }
    }
}

/// The smallest and largest finite values of `raw` within `trimming`.
fn clamp_range(raw: &[f64], trimming: Trimming) -> Option<(f64, f64)> {
    raw.iter()
        .copied()
        .filter(|&v| v.is_finite() && trimming.contains(v))
        .fold(None, |acc, v| {
            Some(acc.map_or((v, v), |(lo, hi): (f64, f64)| (lo.min(v), hi.max(v))))
        })
}

/// Parse a value, accepting Fortran `D`/`d` exponents (`1.5D+02`).
fn parse_number(token: &str) -> Option<f64> {
    token
        .parse()
        .ok()
        .or_else(|| token.replace(['D', 'd'], "E").parse().ok())
}

/// Parse GSLIB file contents.
///
/// - Line 1: title. Line 2: the number of variables `nvar` as the first token;
///   anything after it (e.g. grid dimensions) is ignored. Next `nvar` lines:
///   variable names.
/// - Remaining non-blank lines: `nvar` whitespace-separated numbers each.
///
/// # Errors
/// [`Error::Parse`] with the 1-based line number of the first problem.
pub fn parse_gslib(text: &str) -> Result<GslibData> {
    let mut lines = text.lines().enumerate().map(|(i, l)| (i + 1, l));
    let parse_err = |line: usize, message: String| Error::Parse { line, message };

    let (_, title) = lines
        .next()
        .ok_or_else(|| parse_err(1, "missing title line".into()))?;
    let (nvar_line, header) = lines
        .next()
        .ok_or_else(|| parse_err(2, "missing number-of-variables line".into()))?;
    let nvar: usize = header
        .split_whitespace()
        .next()
        .and_then(|t| t.parse().ok())
        .filter(|&n| n > 0)
        .ok_or_else(|| {
            parse_err(
                nvar_line,
                format!("expected a positive variable count, got {header:?}"),
            )
        })?;

    let names = (0..nvar)
        .map(|k| {
            lines
                .next()
                .map(|(_, l)| l.trim().to_owned())
                .ok_or_else(|| {
                    parse_err(
                        nvar_line + k + 1,
                        format!("missing name of variable {}", k + 1),
                    )
                })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut values = Vec::new();
    for (line, content) in lines.filter(|(_, l)| !l.trim().is_empty()) {
        let row: Vec<f64> = content
            .split_whitespace()
            .map(|t| {
                parse_number(t).ok_or_else(|| parse_err(line, format!("{t:?} is not a number")))
            })
            .collect::<Result<_>>()?;
        if row.len() != nvar {
            return Err(parse_err(
                line,
                format!("expected {nvar} values, got {}", row.len()),
            ));
        }
        values.extend(row);
    }
    GslibData::new(title.trim_end(), names, values)
}

/// Format data as GSLIB file contents. Values use Rust's shortest
/// round-trip representation, so `parse_gslib(&format_gslib(d)) == d`.
#[must_use]
pub fn format_gslib(data: &GslibData) -> String {
    let mut out = String::new();
    // Writing to a String cannot fail.
    let _ = writeln!(out, "{}", data.title);
    let _ = writeln!(out, "{}", data.ncols());
    for name in &data.names {
        let _ = writeln!(out, "{name}");
    }
    for row in data.rows() {
        let line: Vec<String> = row.iter().map(f64::to_string).collect();
        let _ = writeln!(out, "{}", line.join(" "));
    }
    out
}

/// Read and parse a GSLIB file.
///
/// # Errors
/// [`Error::Io`] if the file cannot be read, or [`Error::Parse`].
pub fn read_gslib(path: impl AsRef<Path>) -> Result<GslibData> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|source| Error::Io {
        path: path.to_owned(),
        source,
    })?;
    parse_gslib(&text)
}

/// Format and write a GSLIB file.
///
/// # Errors
/// [`Error::Io`] if the file cannot be written.
pub fn write_gslib(path: impl AsRef<Path>, data: &GslibData) -> Result<()> {
    let path = path.as_ref();
    std::fs::write(path, format_gslib(data)).map_err(|source| Error::Io {
        path: path.to_owned(),
        source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const SAMPLE: &str = "My Data Title\n3\nX Coordinate\nY Coordinate\nValue\n0.0 1.0 10.0\n2.0 3.0 20.0\n4.0 5.0 30.0\n";

    #[test]
    fn parses_basic_file() {
        let d = parse_gslib(SAMPLE).unwrap();
        assert_eq!(d.title(), "My Data Title");
        assert_eq!(d.names(), ["X Coordinate", "Y Coordinate", "Value"]);
        assert_eq!(d.nrows(), 3);
        assert_eq!(d.column(2).unwrap(), vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn tolerates_blank_lines_grid_headers_and_fortran_exponents() {
        let text = "grid\n1 3 1 1\nv\n1.5D+02\n\n  2.0d0 \n3e0\n\n";
        let d = parse_gslib(text).unwrap();
        assert_eq!(d.column(0).unwrap(), vec![150.0, 2.0, 3.0]);
    }

    fn error_line(text: &str) -> usize {
        match parse_gslib(text) {
            Err(Error::Parse { line, .. }) => line,
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    #[test]
    fn errors_report_line_numbers() {
        assert_eq!(error_line("t\n2\na\nb\n1 2\n3\n"), 6, "short row");
        assert_eq!(error_line("t\n1\na\n1\nabc\n"), 5, "not a number");
        assert_eq!(error_line("t\nx\n"), 2, "bad variable count");
        assert_eq!(error_line("t\n2\na\n"), 4, "missing name");
    }

    #[test]
    fn missing_columns_are_errors() {
        let d = parse_gslib(SAMPLE).unwrap();
        assert!(d.column(3).is_err(), "no fourth column");
        let all = Extraction::new(Trimming::NONE);
        assert!(d.samples(0, 1, 1, 3, &all).is_err(), "no value column");
        assert!(d.extract(3, &all).is_err(), "no fourth column to extract");
    }

    #[test]
    fn with_column_appends() {
        let d = parse_gslib(SAMPLE)
            .unwrap()
            .with_column("New", &[7.0, 8.0, 9.0])
            .unwrap();
        assert_eq!(d.ncols(), 4);
        assert_eq!(d.column(0).unwrap(), vec![0.0, 2.0, 4.0]);
        assert_eq!(d.column(3).unwrap(), vec![7.0, 8.0, 9.0]);
        assert!(
            parse_gslib(SAMPLE)
                .unwrap()
                .with_column("Bad", &[1.0])
                .is_err(),
            "wrong length"
        );
    }

    #[test]
    fn file_roundtrip() {
        let d = parse_gslib(SAMPLE).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.dat");
        write_gslib(&path, &d).unwrap();
        assert_eq!(read_gslib(&path).unwrap(), d);
        assert!(
            matches!(
                read_gslib(dir.path().join("missing")),
                Err(Error::Io { .. })
            ),
            "missing file"
        );
    }

    /// Two variables and coordinates, with a sentinel, a NaN and infinities.
    const MIXED: &str = "mixed\n4\nx\ny\na\nb\n\
        0 0 1.0 5.0\n\
        1 0 -999 6.0\n\
        2 0 3.0 -999\n\
        3 0 NaN 7.0\n\
        4 0 inf 8.0\n\
        5 0 -inf 9.0\n";

    fn values(text: &str, column: usize, policy: &Extraction) -> Result<Vec<Option<f64>>> {
        parse_gslib(text)?.extract(column, policy)
    }

    fn one_column(values: &str) -> String {
        let rows: Vec<&str> = values.split_whitespace().collect();
        format!("t\n1\nv\n{}\n", rows.join("\n"))
    }

    #[test]
    fn trimming_limits_are_validated() {
        assert!(Trimming::limits(f64::NAN, 1.0).is_err(), "NaN tmin");
        assert!(Trimming::limits(0.0, f64::NAN).is_err(), "NaN tmax");
        assert!(Trimming::limits(1.0, 1.0).is_err(), "empty interval");
        assert!(Trimming::limits(2.0, 1.0).is_err(), "reversed");
        let open = Trimming::limits(f64::NEG_INFINITY, f64::INFINITY).unwrap();
        assert!(
            open.contains(-1.0e300) && open.contains(1.0e300),
            "infinite limits"
        );
        assert_eq!(Trimming::NONE.bounds(), None);
    }

    #[test]
    fn trimming_is_half_open() {
        // M4: v = tmin is present, v = tmax is not.
        let t = Trimming::limits(0.0, 10.0).unwrap();
        let got = values(&one_column("-0.5 0 5 10 10.5"), 0, &Extraction::new(t)).unwrap();
        assert_eq!(got, vec![None, Some(0.0), Some(5.0), None, None]);
    }

    #[test]
    fn sentinels_are_data_until_trimmed() {
        // M1/M2: with no limits a sentinel is an ordinary value; the CLI's
        // sentinel guard (M6) is what stops it reaching an algorithm.
        let text = one_column("1 -999 3");
        let none = Extraction::new(Trimming::NONE);
        assert_eq!(
            values(&text, 0, &none).unwrap(),
            vec![Some(1.0), Some(-999.0), Some(3.0)]
        );
        let gslib = Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap());
        assert_eq!(
            values(&text, 0, &gslib).unwrap(),
            vec![Some(1.0), None, Some(3.0)]
        );
    }

    #[test]
    fn non_finite_values_are_errors_by_default() {
        // M5: the error names the 1-based row and column.
        let policy = Extraction::new(Trimming::NONE);
        match values(&one_column("1 2 NaN"), 0, &policy) {
            Err(Error::NonFinite {
                row: 3,
                column: 1,
                value,
            }) => assert!(value.is_nan()),
            other => panic!("expected a NaN error at row 3, got {other:?}"),
        }
        match values(&one_column("1 -inf"), 0, &policy) {
            Err(Error::NonFinite {
                row: 2,
                column: 1,
                value,
            }) => {
                assert!(value.is_infinite() && value.is_sign_negative());
            }
            other => panic!("expected an infinity error at row 2, got {other:?}"),
        }
    }

    #[test]
    fn non_finite_values_can_be_missing() {
        let policy = Extraction {
            nan: OnNan::Missing,
            infinite: OnInfinite::Missing,
            ..Extraction::new(Trimming::NONE)
        };
        assert_eq!(
            values(MIXED, 2, &policy).unwrap(),
            vec![Some(1.0), Some(-999.0), Some(3.0), None, None, None]
        );
        // Each policy covers only its own kind of value.
        let nan_only = Extraction {
            nan: OnNan::Missing,
            ..Extraction::new(Trimming::NONE)
        };
        assert!(
            matches!(
                values(MIXED, 2, &nan_only),
                Err(Error::NonFinite { row: 5, .. })
            ),
            "infinity still an error"
        );
        let inf_only = Extraction {
            infinite: OnInfinite::Missing,
            ..Extraction::new(Trimming::NONE)
        };
        assert!(
            matches!(
                values(MIXED, 2, &inf_only),
                Err(Error::NonFinite { row: 4, .. })
            ),
            "NaN still an error"
        );
    }

    #[test]
    fn infinities_clamp_to_the_trimmed_finite_range() {
        // M5: ±Inf take the largest and smallest finite values that the
        // limits keep, so −999 is not the minimum here.
        let policy = Extraction {
            nan: OnNan::Missing,
            infinite: OnInfinite::ClampToRange,
            ..Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap())
        };
        assert_eq!(
            values(MIXED, 2, &policy).unwrap(),
            vec![Some(1.0), None, Some(3.0), None, Some(3.0), Some(1.0)]
        );
        // Without limits the sentinel is data, and so the minimum.
        let untrimmed = Extraction {
            trimming: Trimming::NONE,
            ..policy
        };
        assert_eq!(values(MIXED, 2, &untrimmed).unwrap()[5], Some(-999.0));
        // Nothing finite to clamp to.
        let clamp = Extraction {
            infinite: OnInfinite::ClampToRange,
            ..Extraction::new(Trimming::NONE)
        };
        assert!(
            matches!(
                values(&one_column("inf -inf"), 0, &clamp),
                Err(Error::InvalidInput { .. })
            ),
            "no finite values"
        );
        let trimmed_away = Extraction {
            trimming: Trimming::limits(10.0, 20.0).unwrap(),
            ..clamp
        };
        assert!(
            values(&one_column("1 inf"), 0, &trimmed_away).is_err(),
            "no finite value within the limits"
        );
        // A column without infinities needs no finite value.
        assert_eq!(
            values(&one_column("1"), 0, &trimmed_away).unwrap(),
            vec![None]
        );
    }

    #[test]
    fn non_finite_policy_applies_before_trimming() {
        // M9: GSLIB trims an infinity by its limits (F-010); here it is an
        // error under the default policy even though the limits exclude it.
        let policy = Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap());
        assert!(
            matches!(
                values(&one_column("1 inf"), 0, &policy),
                Err(Error::NonFinite { row: 2, .. })
            ),
            "+Inf beyond tmax"
        );
        assert!(
            matches!(
                values(&one_column("NaN 1"), 0, &policy),
                Err(Error::NonFinite { row: 1, .. })
            ),
            "NaN with finite limits"
        );
    }

    #[test]
    fn missingness_is_per_cell() {
        // M7: each variable of a record is trimmed on its own.
        let d = parse_gslib(MIXED).unwrap();
        let policy = Extraction {
            nan: OnNan::Missing,
            infinite: OnInfinite::Missing,
            ..Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap())
        };
        let a = d.extract(2, &policy).unwrap();
        let b = d.extract(3, &policy).unwrap();
        assert_eq!(a, vec![Some(1.0), None, Some(3.0), None, None, None]);
        assert_eq!(
            b,
            vec![Some(5.0), Some(6.0), None, Some(7.0), Some(8.0), Some(9.0)]
        );
    }

    #[test]
    fn samples_keep_present_values_and_their_rows() {
        let d = parse_gslib(MIXED).unwrap();
        let policy = Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap());
        let present = d.samples(0, 1, 1, 3, &policy).unwrap();
        assert_eq!(present.rows(), [0, 1, 3, 4, 5]);
        assert_eq!(present.total_rows(), 6);
        let vs: Vec<f64> = present.samples().iter().map(Sample::value).collect();
        assert_eq!(vs, vec![5.0, 6.0, 7.0, 8.0, 9.0]);
        assert_eq!(present.samples()[1].location(), Point3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn coordinates_are_never_trimmed() {
        // D-022: a coordinate equal to a sentinel is still a coordinate.
        let text = "t\n2\nx\nv\n-999 1\n0 2\n";
        let policy = Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap());
        let present = parse_gslib(text)
            .unwrap()
            .samples(0, 0, 0, 1, &policy)
            .unwrap();
        assert_eq!(present.rows(), [0, 1]);
        assert_eq!(
            present.samples()[0].location(),
            Point3::new(-999.0, -999.0, -999.0)
        );
        // A non-finite coordinate is an error naming its row and column,
        // even in a row whose value is missing.
        let nan_text = "t\n4\nx\ny\nz\nv\n0 0 0 1\n0 0 0 2\n0 0 0 -999\n";
        for (column, row) in [(1, 3), (2, 2), (3, 1)] {
            let d = parse_gslib(nan_text).unwrap();
            let mut values: Vec<f64> = d.rows().flatten().copied().collect();
            values[(row - 1) * 4 + column - 1] = f64::NAN;
            let d = GslibData::new("t", d.names().to_vec(), values).unwrap();
            match d.samples(0, 1, 2, 3, &policy) {
                Err(Error::NonFinite {
                    row: r, column: c, ..
                }) => {
                    assert_eq!((r, c), (row, column), "NaN coordinate");
                }
                other => panic!("expected a NaN coordinate error, got {other:?}"),
            }
        }
    }

    #[test]
    fn missing_codes_are_finite() {
        assert!(MissingCode::new(f64::NAN).is_err(), "NaN code");
        assert!(MissingCode::new(f64::INFINITY).is_err(), "infinite code");
        let bits = |c: MissingCode| c.value().to_bits();
        assert_eq!(bits(MissingCode::new(-1.0).unwrap()), (-1.0_f64).to_bits());
        assert_eq!(bits(MissingCode::MINUS_999), (-999.0_f64).to_bits());
        assert_eq!(bits(MissingCode::MINUS_99), (-99.0_f64).to_bits());
    }

    #[test]
    fn optional_columns_write_the_code() {
        // M8: None becomes the code; present values are written unchanged.
        let d = parse_gslib(SAMPLE).unwrap();
        let out = d
            .with_optional_column(
                "Score",
                &[Some(0.5), None, Some(-1.25)],
                MissingCode::MINUS_99,
            )
            .unwrap();
        assert_eq!(out.names().last().map(String::as_str), Some("Score"));
        assert_eq!(out.column(3).unwrap(), vec![0.5, -99.0, -1.25]);
        assert_eq!(
            out.column(2).unwrap(),
            d.column(2).unwrap(),
            "other columns kept"
        );
        let text = format_gslib(&out);
        assert!(text.lines().any(|l| l == "2 3 20 -99"), "written as {text}");
    }

    #[test]
    fn optional_columns_reject_ambiguous_or_bad_results() {
        let d = parse_gslib(SAMPLE).unwrap();
        // M3: one entry per row.
        assert!(
            matches!(
                d.with_optional_column("v", &[Some(1.0)], MissingCode::MINUS_999),
                Err(Error::LengthMismatch { .. })
            ),
            "wrong length"
        );
        // A result equal to the code would read back as missing.
        assert!(
            matches!(
                d.with_optional_column(
                    "v",
                    &[None, Some(-999.0), Some(1.0)],
                    MissingCode::MINUS_999
                ),
                Err(Error::InvalidInput { .. })
            ),
            "result equals the code"
        );
        assert!(
            d.with_optional_column("v", &[None, Some(-999.0), Some(1.0)], MissingCode::MINUS_99)
                .is_ok(),
            "fine with another code"
        );
        // A non-finite result is an error naming its row and new column.
        assert!(
            matches!(
                d.with_optional_column(
                    "v",
                    &[None, Some(1.0), Some(f64::NAN)],
                    MissingCode::MINUS_999
                ),
                Err(Error::NonFinite {
                    row: 3,
                    column: 4,
                    ..
                })
            ),
            "NaN result"
        );
    }

    #[test]
    fn results_align_with_rows() {
        // M3: results computed on the present samples go back to their rows.
        let d = parse_gslib(MIXED).unwrap();
        let policy = Extraction {
            nan: OnNan::Missing,
            infinite: OnInfinite::Missing,
            ..Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap())
        };
        let present = d.samples(0, 1, 1, 2, &policy).unwrap();
        assert_eq!(present.rows(), [0, 2]);
        let doubled: Vec<f64> = present.samples().iter().map(|s| 2.0 * s.value()).collect();
        let aligned = present.align(&doubled).unwrap();
        assert_eq!(aligned, vec![Some(2.0), None, Some(6.0), None, None, None]);
        assert_eq!(aligned.len(), d.nrows());
        assert!(
            matches!(present.align(&[1.0]), Err(Error::LengthMismatch { .. })),
            "one result per sample"
        );
        let out = d
            .with_optional_column("2a", &aligned, MissingCode::MINUS_999)
            .unwrap();
        assert_eq!(
            out.column(4).unwrap(),
            vec![2.0, -999.0, 6.0, -999.0, -999.0, -999.0]
        );
    }

    proptest! {
        #[test]
        fn optional_columns_roundtrip_through_trimming(
            raw in prop::collection::vec(
                prop::option::of(prop::num::f64::NORMAL | prop::num::f64::ZERO),
                0..30,
            ),
        ) {
            // M3 + M8: write with the code, format, parse, extract with
            // limits that exclude it: the same rows are missing, and the
            // present values are bit-identical.
            let results: Vec<Option<f64>> = raw
                .into_iter()
                .map(|v| v.filter(|x| *x > -998.0 && *x < 1.0e21))
                .collect();
            let n = results.len();
            let base = GslibData::new(
                "t",
                vec!["i".to_owned()],
                (0..n).map(|i| f64::from(u32::try_from(i).unwrap())).collect(),
            )
            .unwrap();
            let written = base
                .with_optional_column("r", &results, MissingCode::MINUS_999)
                .unwrap();
            let read = parse_gslib(&format_gslib(&written)).unwrap();
            let policy = Extraction::new(Trimming::limits(-998.0, 1.0e21).unwrap());
            let back = read.extract(1, &policy).unwrap();
            prop_assert_eq!(back.len(), n);
            for (b, r) in back.iter().zip(&results) {
                prop_assert_eq!(b.map(f64::to_bits), r.map(f64::to_bits));
            }
        }
    }

    proptest! {
        #[test]
        fn format_then_parse_roundtrips(
            ncols in 1_usize..5,
            // M1: sentinel codes and huge values survive unchanged.
            raw in prop::collection::vec(
                prop_oneof![
                    4 => prop::num::f64::NORMAL | prop::num::f64::ZERO,
                    1 => prop::sample::select(vec![-999.0, -99.0, -1.0, 1.0e21, -1.0e21]),
                ],
                0..40,
            ),
        ) {
            let whole_rows = raw.len() - raw.len() % ncols;
            let values: Vec<f64> = raw.into_iter().take(whole_rows).collect();
            let names = (0..ncols).map(|i| format!("v{i}")).collect();
            let d = GslibData::new("title", names, values).unwrap();
            prop_assert_eq!(parse_gslib(&format_gslib(&d)).unwrap(), d);
        }
    }
}
