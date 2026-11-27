//! Parity tests of gslib-rs against Fortran GSLIB reference fixtures.
//!
//! This crate has no library code: the tests live in `tests/`, and read the
//! committed fixtures in `tests/fixtures/gslib/` at the workspace root
//! (`docs/spec/harness.md`, R2 and R6). They never need network access,
//! gfortran, or GSLIB itself.
