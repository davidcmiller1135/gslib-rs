//! Parity tests of gslib-rs against Fortran GSLIB reference fixtures.
//!
//! The tests live in `tests/` and read the committed fixtures in
//! `tests/fixtures/gslib/` at the workspace root (`docs/spec/harness.md`, R2
//! and R6). They never need network access, gfortran, or GSLIB itself. This
//! library holds what they share: readers for GSLIB's output formats and the
//! comparison rules (R7, D-010).

pub mod compare;
pub mod gamv;
pub mod gslib_math;
pub mod vmodel;

use std::path::PathBuf;

/// `tests/fixtures/gslib/` at the workspace root.
#[must_use]
pub fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/gslib")
}

/// Parse whitespace-separated numeric rows (e.g. a transformation table),
/// skipping blank lines; every row must have `width` values.
///
/// # Errors
/// A message naming the first bad line.
pub fn parse_columns(text: &str, width: usize) -> Result<Vec<Vec<f64>>, String> {
    text.lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(i, line)| {
            let row = line
                .split_whitespace()
                .map(|t| {
                    t.parse::<f64>()
                        .map_err(|_| format!("line {}: {t:?} is not a number", i + 1))
                })
                .collect::<Result<Vec<f64>, String>>()?;
            if row.len() == width {
                Ok(row)
            } else {
                Err(format!(
                    "line {}: expected {width} values, got {}",
                    i + 1,
                    row.len()
                ))
            }
        })
        .collect()
}
