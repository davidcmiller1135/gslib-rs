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
//! Entry points: [`try_map`] for any per-item function (e.g. kriging at many
//! targets), [`experimental_variogram`] for the pair loop's blocks, and
//! [`realizations`] for stochastic simulation, where realization `k` draws
//! from its own stream `Seed::realization(k)` (D-037), so the result does not
//! depend on scheduling.
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
use crate::geometry::Point3;
use crate::random::{Seed, Stream};
use crate::sample::Sample;
use crate::variogram::{Lag, LagScheme, LagSums, PairPlan};

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

/// The experimental semivariogram of `samples`, with its row blocks summed
/// in parallel: bit-identical to [`crate::variogram::experimental_variogram`]
/// for any thread count (D-031, P3).
///
/// Each block of rows is summed sequentially, as in the sequential path,
/// and the block sums are combined in block order.
#[must_use]
pub fn experimental_variogram<L>(samples: &[Sample], lags: &L) -> Vec<Lag>
where
    L: LagScheme + Sync + ?Sized,
{
    let points: Vec<Point3> = samples.iter().map(Sample::location).collect();
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    let plan = PairPlan::new(&points, &values, lags);
    let blocks: Vec<usize> = (0..plan.block_count()).collect();
    // `collect` into a `Vec` keeps block order.
    let sums: Vec<Vec<LagSums>> = blocks.par_iter().map(|&b| plan.block_sums(b)).collect();
    plan.finish(sums)
}

/// Run `count` realizations in parallel, each on its own stream.
///
/// Calls `simulate(k, stream)` with `stream = seed.realization(k)` for `k`
/// in `0..count`. Results come back in realization order, or the first
/// error in realization order (P2, P5).
///
/// Each realization depends only on `(seed, k)`, so the results are the
/// same for any thread count, and realization `k` equals a run of
/// realization `k` alone.
///
/// # Errors
/// The error of the first realization, in order, for which `simulate`
/// fails.
pub fn realizations<T, F>(seed: Seed, count: usize, simulate: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(usize, &mut Stream) -> Result<T> + Sync,
{
    let indices: Vec<usize> = (0..count).collect();
    try_map(&indices, |&k| {
        let stream = u64::try_from(k)
            .map_err(|_| Error::invalid("realization", format!("index {k} exceeds u64")))?;
        simulate(k, &mut seed.realization(stream))
    })
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
    use crate::model::{CovarianceModel, Kind, Structure, VariogramModel};
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

    #[test]
    fn variogram_is_identical_for_any_thread_count() {
        // P1, P3: several row blocks (1 500 data, 6 blocks), both lag
        // schemes; parallel equals sequential bit for bit.
        use crate::variogram::{CentredLags, LagClasses};
        let samples: Vec<Sample> = (0..1_500_u32)
            .map(|i| {
                let f = f64::from(i);
                let (x, y) = (
                    (f * 0.618_034).fract() * 200.0,
                    (f * 0.414_214).fract() * 200.0,
                );
                Sample::new(Point3::new(x, y, (f * 0.1).fract()), (x * 0.03).sin() * y).unwrap()
            })
            .collect();
        let lag_bits = |lags: &[Lag]| -> Vec<[u64; 5]> {
            lags.iter()
                .map(|l| {
                    let b = |v: Option<f64>| v.map_or(u64::MAX, f64::to_bits);
                    [
                        l.center.to_bits(),
                        b(l.mean_distance),
                        u64::try_from(l.pairs).unwrap(),
                        b(l.gamma),
                        b(l.pair_mean),
                    ]
                })
                .collect()
        };
        let centred = CentredLags::new(12, 4.0, 3.0).unwrap();
        let classes = LagClasses::new(60.0, 15).unwrap();
        let sequential = (
            lag_bits(&crate::variogram::experimental_variogram(
                &samples, &centred,
            )),
            lag_bits(&crate::variogram::experimental_variogram(
                &samples, &classes,
            )),
        );
        assert!(sequential.0.iter().any(|l| l[2] > 1000), "a real workload");
        for threads in pools() {
            let parallel = with_threads(threads, || {
                (
                    lag_bits(&experimental_variogram(&samples, &centred)),
                    lag_bits(&experimental_variogram(&samples, &classes)),
                )
            })
            .unwrap();
            assert_eq!(parallel, sequential, "{threads} threads");
        }
        // Degenerate inputs.
        assert_eq!(
            lag_bits(&experimental_variogram(&[], &centred)),
            lag_bits(&crate::variogram::experimental_variogram(&[], &centred))
        );
    }

    #[test]
    #[ignore = "timing evidence; run with --release --features parallel --ignored --nocapture"]
    fn variogram_timing() {
        use crate::variogram::CentredLags;
        let samples: Vec<Sample> = (0..20_000_u32)
            .map(|i| {
                let f = f64::from(i);
                let (x, y) = (
                    (f * 0.618_034).fract() * 1000.0,
                    (f * 0.414_214).fract() * 1000.0,
                );
                Sample::new(Point3::new(x, y, 0.0), (x * 0.01).sin()).unwrap()
            })
            .collect();
        for lags in [
            CentredLags::new(10, 5.0, 2.5).unwrap(),
            CentredLags::new(10, 20.0, 10.0).unwrap(),
        ] {
            let started = std::time::Instant::now();
            let sequential = crate::variogram::experimental_variogram(&samples, &lags);
            let one = started.elapsed();
            for threads in pools() {
                let timer = std::time::Instant::now();
                let parallel =
                    with_threads(threads, || experimental_variogram(&samples, &lags)).unwrap();
                let took = timer.elapsed();
                assert_eq!(parallel, sequential);
                println!("{lags:?}: sequential {one:?}, {threads} threads {took:?}");
            }
        }
    }

    #[test]
    fn realizations_are_independent_of_threads_and_batch() {
        // P1, P5: 12 sgsim realizations on a 15 x 15 grid with a datum.
        use crate::search::{DataSearch, NodeSearch, SearchEllipsoid, TableExtent};
        use crate::simulation::{Axis, DataConditioning, Grid, Searches, sgsim};
        let axis = |count| Axis {
            count,
            first_center: 0.5,
            spacing: 1.0,
        };
        let grid = Grid::new(axis(15), axis(15), axis(1)).unwrap();
        let (_, variogram) = field();
        let model = CovarianceModel::try_from(variogram).unwrap();
        let data = [Sample::new(Point3::new(7.5, 7.5, 0.0), 1.0).unwrap()];
        let search = DataSearch::new(SearchEllipsoid::isotropic(10.0).unwrap(), 0, 8).unwrap();
        let nodes = NodeSearch::new(
            &grid,
            search.ellipsoid(),
            &model,
            12,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap();
        let searches = Searches {
            data: DataConditioning::Search(&search),
            nodes: &nodes,
        };
        let simulate = |_: usize, rng: &mut Stream| sgsim(&data, &grid, &model, 0.0, searches, rng);
        let bits = |runs: &[Vec<f64>]| -> Vec<Vec<u64>> {
            runs.iter()
                .map(|r| r.iter().map(|v| v.to_bits()).collect())
                .collect()
        };
        let seed = Seed(69069);
        let sequential: Vec<Vec<f64>> = (0..12_u64)
            .map(|k| {
                sgsim(
                    &data,
                    &grid,
                    &model,
                    0.0,
                    searches,
                    &mut seed.realization(k),
                )
                .unwrap()
            })
            .collect();
        for threads in pools() {
            let parallel = with_threads(threads, || realizations(seed, 12, simulate))
                .unwrap()
                .unwrap();
            assert_eq!(bits(&parallel), bits(&sequential), "{threads} threads");
        }
        // Realization 7 alone.
        let alone = sgsim(
            &data,
            &grid,
            &model,
            0.0,
            searches,
            &mut seed.realization(7),
        )
        .unwrap();
        assert_eq!(bits(&[alone]), bits(&sequential[7..8]));
        // The first failing realization, in order, is reported.
        let failing = |k: usize, _: &mut Stream| -> Result<usize> {
            if k == 3 || k == 9 {
                Err(Error::invalid(
                    "realization",
                    format!("realization {k} failed"),
                ))
            } else {
                Ok(k)
            }
        };
        match realizations(seed, 12, failing) {
            Err(Error::InvalidInput { reason, .. }) => assert_eq!(reason, "realization 3 failed"),
            other => panic!("expected realization 3's error, got {other:?}"),
        }
    }
}
