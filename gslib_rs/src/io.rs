//! GSLIB (simplified Geo-EAS) ASCII files.
//!
//! [`parse_gslib`] and [`format_gslib`] are pure; [`read_gslib`] and
//! [`write_gslib`] are thin wrappers that touch the filesystem.

use std::fmt::Write as _;
use std::path::Path;

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::sample::Sample;

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

    /// Samples from three coordinate columns and a value column (0-based).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if a column does not exist or a row has a
    /// non-finite coordinate or value.
    pub fn samples(
        &self,
        x_col: usize,
        y_col: usize,
        z_col: usize,
        value_col: usize,
    ) -> Result<Vec<Sample>> {
        crate::sample::zip_samples(&self.points(x_col, y_col, z_col)?, &self.column(value_col)?)
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
        assert!(d.samples(0, 1, 1, 3).is_err(), "no value column");
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

    proptest! {
        #[test]
        fn format_then_parse_roundtrips(
            ncols in 1_usize..5,
            raw in prop::collection::vec(prop::num::f64::NORMAL | prop::num::f64::ZERO, 0..40),
        ) {
            let whole_rows = raw.len() - raw.len() % ncols;
            let values: Vec<f64> = raw.into_iter().take(whole_rows).collect();
            let names = (0..ncols).map(|i| format!("v{i}")).collect();
            let d = GslibData::new("title", names, values).unwrap();
            prop_assert_eq!(parse_gslib(&format_gslib(&d)).unwrap(), d);
        }
    }
}
