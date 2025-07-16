//! GSLIB-RS library

/// Core utilities: RNG, sorting, solvers, string handling, etc.
pub mod utils;

/// I/O: reading & writing GSLIB ASCII files
pub mod io;

/// Data transforms: normal-score, back-transform, rotation, etc.
pub mod transform;

/// Variogram: empirical calculation and model fitting
pub mod variogram;

/// Kriging: ordinary, simple, cokriging
pub mod kriging;

/// Simulation: SGS, SIS, SASIM stubs
pub mod simulation;

/// Plotting (requires the `plot` feature)

pub mod plot;
