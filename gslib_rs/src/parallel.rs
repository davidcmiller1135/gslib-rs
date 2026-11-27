//! Parallel maps over the core's pure per-item functions
//! (`docs/spec/parallelism.md`; feature `parallel`, D-030).
//!
//! The algorithms stay sequential and pure per item (one target, one block,
//! one realization). This module only runs many items at once, and its
//! results never depend on the thread count (`AGENTS.md`, precedence 1):
//!
//! - **Order:** results come back in item order, whatever order the threads
//!   finish in.
//! - **Errors (P2):** if any item fails, the error returned is the one from
//!   the *first failing item in item order*, not whichever thread failed
//!   first. Every item is therefore evaluated even when one fails, which
//!   keeps the outcome independent of scheduling.
//!
//! The threads come from `rayon`'s global pool (sized by `RAYON_NUM_THREADS`
//! or the number of cores) unless the call runs inside [`with_threads`].
//!
//! # Example
//!
//! Ordinary kriging at many targets, each from its own search:
//!
//! ```
//! # fn main() -> gslib_rs::Result<()> {
//! use gslib_rs::geometry::Point3;
//! use gslib_rs::kriging::ordinary_kriging;
//! use gslib_rs::model::{Kind, Structure, VariogramModel};
//! use gslib_rs::parallel::try_map;
//! use gslib_rs::sample::Sample;
//! use gslib_rs::search::{DataSearch, SearchEllipsoid, SuperBlock};
//!
//! let data = vec![
//!     Sample::new(Point3::new(0.0, 0.0, 0.0), 1.0)?,
//!     Sample::new(Point3::new(10.0, 0.0, 0.0), 2.0)?,
//!     Sample::new(Point3::new(0.0, 10.0, 0.0), 3.0)?,
//! ];
//! let model = VariogramModel::new(0.0, vec![Structure::isotropic(Kind::Spherical, 1.0, 20.0)?])?;
//! let search = DataSearch::new(SearchEllipsoid::isotropic(50.0)?, 1, 8)?;
//! let index = SuperBlock::new(&data);
//! let targets = [Point3::new(2.0, 2.0, 0.0), Point3::new(8.0, 1.0, 0.0)];
//! let estimates = try_map(&targets, |&t| {
//!     ordinary_kriging(&search.select_samples(&index, t)?, t, &model)
//! })?;
//! assert_eq!(estimates.len(), 2);
//! # Ok(())
//! # }
//! ```

use std::num::NonZeroUsize;

use rayon::prelude::*;

use crate::error::{Error, Result};

/// Apply `f` to every item in parallel and return the results in item
/// order, or the first error in item order (P2).
///
/// # Errors
/// The error of the first item, in item order, for which `f` fails.
pub fn try_map<I, T, F>(items: &[I], f: F) -> Result<Vec<T>>
where
    I: Sync,
    T: Send,
    F: Fn(&I) -> Result<T> + Sync,
{
    // `collect` into a `Vec` keeps item order; the sequential pass then
    // picks the first error by position, not by time.
    let results: Vec<Result<T>> = items.par_iter().map(&f).collect();
    results.into_iter().collect()
}

/// Run `work` in a dedicated pool of exactly `threads` threads.
///
/// Parallel maps inside `work` use that many threads. Results do not depend
/// on the choice; this exists to cap resource use and to test that they do
/// not.
///
/// # Errors
/// [`Error::InvalidInput`] if the operating system cannot provide the pool.
pub fn with_threads<R, W>(threads: NonZeroUsize, work: W) -> Result<R>
where
    R: Send,
    W: FnOnce() -> R + Send,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .build()
        .map_err(|e| Error::invalid("thread pool", e.to_string()))?;
    Ok(pool.install(work))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::kriging::{Estimate, ordinary_kriging};
    use crate::model::{Kind, Structure, VariogramModel};
    use crate::sample::Sample;
    use crate::search::{DataSearch, SearchEllipsoid, SuperBlock};

    fn pools() -> [NonZeroUsize; 3] {
        [1, 2, 8].map(|n| NonZeroUsize::new(n).unwrap())
    }

    /// Scattered data and a model for kriging tests.
    fn field() -> (Vec<Sample>, VariogramModel) {
        let data = (0..300_u32)
            .map(|i| {
                let f = f64::from(i);
                let (x, y) = (
                    (f * 0.618_034).fract() * 100.0,
                    (f * 0.414_214).fract() * 100.0,
                );
                Sample::new(Point3::new(x, y, 0.0), y.mul_add(0.01, (x * 0.05).sin())).unwrap()
            })
            .collect();
        let model = VariogramModel::new(
            0.05,
            vec![Structure::isotropic(Kind::Spherical, 1.0, 40.0).unwrap()],
        )
        .unwrap();
        (data, model)
    }

    fn bits(estimates: &[Estimate]) -> Vec<(u64, u64)> {
        estimates
            .iter()
            .map(|e| (e.value.to_bits(), e.variance.to_bits()))
            .collect()
    }

    #[test]
    fn kriging_is_identical_for_any_thread_count() {
        // P1, P2: 400 targets with ordinary kriging from 16 neighbours each.
        let (data, model) = field();
        let index = SuperBlock::new(&data);
        let search = DataSearch::new(SearchEllipsoid::isotropic(30.0).unwrap(), 1, 16).unwrap();
        let targets: Vec<Point3> = (0..400_u32)
            .map(|i| {
                let f = f64::from(i);
                Point3::new(
                    (f * 0.732_051).fract() * 100.0,
                    (f * 0.236_068).fract() * 100.0,
                    0.0,
                )
            })
            .collect();
        let krige = |&t: &Point3| ordinary_kriging(&search.select_samples(&index, t)?, t, &model);
        let sequential: Vec<Estimate> = targets.iter().map(krige).collect::<Result<_>>().unwrap();
        for threads in pools() {
            let parallel = with_threads(threads, || try_map(&targets, krige))
                .unwrap()
                .unwrap();
            assert_eq!(bits(&parallel), bits(&sequential), "{threads} threads");
        }
    }

    #[test]
    fn the_first_error_in_item_order_wins() {
        // P2: items 50 and 150 fail; item 50 is made slow, so with several
        // threads item 150 usually fails first in time. The reported error
        // must still be item 50's.
        let items: Vec<u32> = (0..200).collect();
        let f = |&i: &u32| -> Result<u32> {
            if i == 50 {
                let started = std::time::Instant::now();
                while started.elapsed() < std::time::Duration::from_millis(20) {
                    std::hint::spin_loop();
                }
            }
            if i == 50 || i == 150 {
                return Err(Error::invalid("item", format!("item {i} failed")));
            }
            Ok(i * 2)
        };
        for threads in pools() {
            match with_threads(threads, || try_map(&items, f)).unwrap() {
                Err(Error::InvalidInput { reason, .. }) => {
                    assert_eq!(reason, "item 50 failed", "{threads} threads");
                }
                other => panic!("expected item 50's error, got {other:?}"),
            }
        }
        // And without errors, results in item order.
        let ok = try_map(&items, |&i| Ok(i * 2)).unwrap();
        assert_eq!(ok, items.iter().map(|i| i * 2).collect::<Vec<_>>());
    }

    #[test]
    fn kriging_errors_name_the_first_failing_target() {
        // P2 with a real failure: targets beyond the data's reach find too
        // few data; the error names the first such target in order.
        let (data, model) = field();
        let index = SuperBlock::new(&data);
        let search = DataSearch::new(SearchEllipsoid::isotropic(10.0).unwrap(), 3, 8).unwrap();
        let mut targets: Vec<Point3> = (0..100_u32)
            .map(|i| Point3::new(f64::from(i), 50.0, 0.0))
            .collect();
        targets.insert(30, Point3::new(500.0, 500.0, 0.0));
        targets.insert(80, Point3::new(-400.0, 0.0, 0.0));
        let krige = |&t: &Point3| ordinary_kriging(&search.select_samples(&index, t)?, t, &model);
        for threads in pools() {
            match with_threads(threads, || try_map(&targets, krige)).unwrap() {
                Err(Error::TooFewData { target, .. }) => {
                    assert_eq!(target, Point3::new(500.0, 500.0, 0.0), "{threads} threads");
                }
                other => panic!("expected a too-few-data error, got {other:?}"),
            }
        }
    }
}
