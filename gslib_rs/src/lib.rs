//! # `gslib_rs`
//!
//! A Rust port of [GSLIB](http://www.gslib.com/), the Geostatistical Software
//! Library: variogram modelling, kriging, and stochastic simulation.
//!
//! The crate is organised as a pure numerical core plus thin adapters:
//!
//! - **Domain types:** [`geometry::Point3`], [`sample::Sample`],
//!   [`model::VariogramModel`], [`model::CovarianceModel`].
//! - **Pure algorithms:** [`variogram`], [`kriging`], [`transform`],
//!   [`simulation`]. They take their inputs, including any random number
//!   generator, as arguments, return [`Result`] instead of panicking, and
//!   touch no global state.
//! - **Adapters:** [`io`] (GSLIB files; parsing itself is pure) and `plot`
//!   (PNG/SVG output; behind the default `plot` feature).
//!
//! ```
//! use gslib_rs::geometry::Point3;
//! use gslib_rs::kriging::ordinary_kriging;
//! use gslib_rs::model::{Kind, Structure, VariogramModel};
//! use gslib_rs::sample::Sample;
//!
//! # fn main() -> gslib_rs::Result<()> {
//! let model = VariogramModel::new(0.1, vec![Structure::isotropic(Kind::Spherical, 0.9, 10.0)?])?;
//! let data = [
//!     Sample::new(Point3::new(0.0, 0.0, 0.0), 1.0)?,
//!     Sample::new(Point3::new(4.0, 0.0, 0.0), 3.0)?,
//! ];
//! let estimate = ordinary_kriging(&data, Point3::new(2.0, 0.0, 0.0), &model)?;
//! assert!((estimate.value - 2.0).abs() < 1e-9);
//! # Ok(())
//! # }
//! ```

pub mod error;
pub mod geometry;
pub mod io;
pub mod kriging;
mod linalg;
pub mod model;
mod num;
#[cfg(feature = "parallel")]
pub mod parallel;
#[cfg(feature = "plot")]
pub mod plot;
pub mod random;
pub mod sample;
pub mod search;
pub mod simulation;
pub mod transform;
pub mod variogram;

pub use error::{Error, Result};
