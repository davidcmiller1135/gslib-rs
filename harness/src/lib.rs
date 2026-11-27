//! Parity tests of gslib-rs against Fortran GSLIB reference fixtures.
//!
//! The tests live in `tests/` and read the committed fixtures in
//! `tests/fixtures/gslib/` at the workspace root (`docs/spec/harness.md`, R2
//! and R6). They never need network access, gfortran, or GSLIB itself. This
//! library holds what they share: readers for GSLIB's output formats and the
//! comparison rules (R7, D-010).

pub mod compare;
pub mod vmodel;

use std::path::PathBuf;

/// `tests/fixtures/gslib/` at the workspace root.
#[must_use]
pub fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/gslib")
}
