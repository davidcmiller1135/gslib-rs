//! Stochastic simulation: sequential Gaussian (SGS), sequential indicator
//! (SIS), and simulated annealing (SASIM).
//!
//! Every simulation takes its random number generator as an argument; none
//! keeps hidden global state. Seed one (e.g. `rand_chacha::ChaCha8Rng`) for
//! reproducible runs, and give each thread its own.

mod grid;
mod sasim;
mod sequential;

pub use grid::{Axis, Grid};
pub use sasim::{
    Conditioning, CorrelationObjective, HistogramObjective, IndicatorVariogramObjective,
    Objectives, Schedule, VariogramObjective, sasim,
};
pub use sequential::{sgsim, sisim};
