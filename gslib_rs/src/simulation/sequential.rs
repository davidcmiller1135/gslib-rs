//! Sequential Gaussian and sequential indicator simulation on a [`Grid`].
//!
//! Each node, visited in random order, is kriged from two searches, as in
//! GSLIB `sgsim` (`docs/spec/search.md`, D-026): the conditioning data, found
//! by a [`DataSearch`] over a [`SuperBlock`] of the data, then the previously
//! simulated nodes, found by a [`NodeSearch`]. Data come first in the kriging
//! system, nearest first; nodes follow in the node search's order.
//!
//! A node within [`COINCIDENT`] of a datum takes the datum's value and is not
//! simulated. As in GSLIB, such a node is not an informed node for later
//! node searches, since the datum itself is found by the data search: one
//! location never enters a kriging system twice.

use rand::seq::SliceRandom;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, StandardNormal};

use super::Grid;
use crate::error::{Error, Result};
use crate::geometry::{COINCIDENT, Point3};
use crate::kriging::simple_kriging;
use crate::model::CovarianceModel;
use crate::num::mean;
use crate::sample::Sample;
use crate::search::{DataSearch, NodeSearch, SuperBlock};

/// The two searches of sequential simulation.
#[derive(Debug, Clone, Copy)]
pub struct Searches<'a> {
    /// The conditioning data.
    pub data: &'a DataSearch,
    /// The previously simulated nodes; built for the same grid.
    pub nodes: &'a NodeSearch,
}

/// Visit the grid's nodes in random order, computing each value from the
/// data and the informed nodes the searches select, and return values in
/// node order.
fn sequential<R, F>(
    data: &[Sample],
    grid: &Grid,
    searches: Searches<'_>,
    rng: &mut R,
    mut draw: F,
) -> Result<Vec<f64>>
where
    R: Rng + ?Sized,
    F: FnMut(&[Sample], Point3, &mut R) -> Result<f64>,
{
    if searches.nodes.counts() != grid.counts() {
        return Err(Error::invalid(
            "node search",
            format!(
                "built for a {:?} grid, used with {:?}",
                searches.nodes.counts(),
                grid.counts()
            ),
        ));
    }
    // Nodes holding a datum take its value and are never simulated.
    let mut fixed: Vec<Option<f64>> = vec![None; grid.len()];
    for d in data {
        if let Some(i) = grid.node_index(d.location())
            && grid.node(i).distance(d.location()) < COINCIDENT
            && let Some(slot) = fixed.get_mut(i)
        {
            *slot = Some(d.value());
        }
    }
    let data_index = SuperBlock::new(data);
    let mut path: Vec<usize> = (0..grid.len()).collect();
    path.shuffle(rng);

    let mut values = vec![0.0; grid.len()];
    let mut informed = vec![false; grid.len()];
    for index in path {
        let node = grid.node(index);
        let value = if let Some(v) = fixed.get(index).copied().flatten() {
            v
        } else {
            let mut conditioning = searches.data.select_samples(&data_index, node)?;
            for n in searches
                .nodes
                .select(grid, index, |i| informed.get(i).copied().unwrap_or(false))?
            {
                if let Some(&v) = values.get(n) {
                    conditioning.push(Sample::new(grid.node(n), v)?);
                }
            }
            let v = draw(&conditioning, node, rng)?;
            if let Some(flag) = informed.get_mut(index) {
                *flag = true;
            }
            v
        };
        if let Some(slot) = values.get_mut(index) {
            *slot = value;
        }
    }
    Ok(values)
}

/// One realization of **sequential Gaussian simulation**.
///
/// `data` must already be in normal-score units with the given stationary
/// `mean` (0 after [`crate::transform::normal_scores`]), and `model` should
/// have a sill of 1. Each node draws from N(SK estimate, SK variance).
///
/// Values are returned in node order (x fastest).
///
/// # Errors
/// - [`Error::SingularSystem`] if a kriging system fails (e.g. duplicate data
///   locations);
/// - [`Error::TooFewData`] if the data search requires more data than it
///   finds;
/// - [`Error::InvalidInput`] if the node search was built for another grid.
pub fn sgsim<R>(
    data: &[Sample],
    grid: &Grid,
    model: &CovarianceModel,
    mean: f64,
    searches: Searches<'_>,
    rng: &mut R,
) -> Result<Vec<f64>>
where
    R: Rng + ?Sized,
{
    sequential(data, grid, searches, rng, |conditioning, node, rng| {
        let estimate = simple_kriging(conditioning, node, model, mean)?;
        let deviate: f64 = StandardNormal.sample(rng);
        Ok(estimate.variance.sqrt().mul_add(deviate, estimate.value))
    })
}

/// One realization of **sequential indicator simulation** for a single
/// threshold. Returns `true` where the simulated value is below `threshold`.
///
/// Indicators are kriged by simple kriging around the global proportion p₀
/// of data below the threshold: p* = p₀ + Σ wᵢ(iᵢ − p₀), clamped to [0, 1].
/// `model` is the indicator covariance.
///
/// # Errors
/// - [`Error::InvalidInput`] if `data` is empty (p₀ is undefined) or the
///   node search was built for another grid;
/// - [`Error::SingularSystem`] if a kriging system fails;
/// - [`Error::TooFewData`] if the data search requires more data than it
///   finds.
pub fn sisim<R>(
    data: &[Sample],
    grid: &Grid,
    threshold: f64,
    model: &CovarianceModel,
    searches: Searches<'_>,
    rng: &mut R,
) -> Result<Vec<bool>>
where
    R: Rng + ?Sized,
{
    let indicators = data
        .iter()
        .map(|s| Sample::new(s.location(), indicator(s.value() < threshold)))
        .collect::<Result<Vec<_>>>()?;
    let values: Vec<f64> = indicators.iter().map(Sample::value).collect();
    let p0 =
        mean(&values).ok_or_else(|| Error::invalid("data", "sisim needs at least one datum"))?;

    let simulated = sequential(
        &indicators,
        grid,
        searches,
        rng,
        |conditioning, node, rng| {
            let p = simple_kriging(conditioning, node, model, p0)?
                .value
                .clamp(0.0, 1.0);
            Ok(indicator(rng.random::<f64>() < p))
        },
    )?;
    Ok(simulated.into_iter().map(|v| v > 0.5).collect())
}

const fn indicator(below: bool) -> f64 {
    if below { 1.0 } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Kind, Structure, VariogramModel};
    use crate::search::{SearchEllipsoid, TableExtent};
    use crate::simulation::Axis;
    use approx::assert_abs_diff_eq;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }
    /// A 1-D grid: `count` nodes from `first` with `spacing`.
    fn line(count: usize, first: f64, spacing: f64) -> Grid {
        let one = Axis {
            count: 1,
            first_center: 0.0,
            spacing: 1.0,
        };
        Grid::new(
            Axis {
                count,
                first_center: first,
                spacing,
            },
            one,
            one,
        )
        .unwrap()
    }
    fn exp_cov(range: f64) -> CovarianceModel {
        CovarianceModel::try_from(
            VariogramModel::new(
                0.0,
                vec![Structure::isotropic(Kind::Exponential, 1.0, range).unwrap()],
            )
            .unwrap(),
        )
        .unwrap()
    }
    fn near() -> DataSearch {
        DataSearch::new(SearchEllipsoid::unlimited(), 0, 12).unwrap()
    }
    fn nodes(grid: &Grid, model: &CovarianceModel, max: usize) -> NodeSearch {
        NodeSearch::new(
            grid,
            &SearchEllipsoid::unlimited(),
            model,
            max,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap()
    }
    fn rng(seed: u64) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(seed)
    }

    #[test]
    fn sgsim_reproduces_data_at_data_locations() {
        let data = [s(0.0, 1.5), s(5.0, -0.5)];
        let g = line(6, 0.0, 1.0);
        let m = exp_cov(3.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &DataSearch::unlimited(),
            nodes: &ns,
        };
        let r = sgsim(&data, &g, &m, 0.0, searches, &mut rng(1)).unwrap();
        assert_abs_diff_eq!(r[0], 1.5);
        assert_abs_diff_eq!(r[5], -0.5);
        assert_eq!(r.len(), 6);
    }

    #[test]
    fn sgsim_same_seed_same_realization() {
        let data = [s(0.0, 1.0)];
        let g = line(20, 1.0, 1.0);
        let m = exp_cov(5.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &near(),
            nodes: &ns,
        };
        let run = |seed| sgsim(&data, &g, &m, 0.0, searches, &mut rng(seed)).unwrap();
        assert_eq!(run(7), run(7));
        assert_ne!(run(7), run(8));
    }

    #[test]
    fn sgsim_unconditional_moments_are_standard_normal() {
        // Nodes far apart ⇒ independent N(0,1) draws.
        let g = line(400, 0.0, 1000.0);
        let m = exp_cov(1.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &near(),
            nodes: &ns,
        };
        let r = sgsim(&[], &g, &m, 0.0, searches, &mut rng(3)).unwrap();
        let mu = mean(&r).unwrap();
        let var = r.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / 399.0;
        assert!(mu.abs() < 0.2, "mean {mu}");
        assert!((var - 1.0).abs() < 0.2, "variance {var}");
    }

    #[test]
    fn bounded_searches_match_unlimited_when_they_cover_everything() {
        // 7 nodes and 2 data: 12 data and 12 nodes already cover everything.
        let data = [s(0.0, 1.0), s(4.0, -1.0)];
        let g = line(7, 0.5, 0.5);
        let m = exp_cov(3.0);
        let (all, twelve) = (nodes(&g, &m, 100), nodes(&g, &m, 12));
        let a = sgsim(
            &data,
            &g,
            &m,
            0.0,
            Searches {
                data: &DataSearch::unlimited(),
                nodes: &all,
            },
            &mut rng(5),
        )
        .unwrap();
        let b = sgsim(
            &data,
            &g,
            &m,
            0.0,
            Searches {
                data: &near(),
                nodes: &twelve,
            },
            &mut rng(5),
        )
        .unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert_abs_diff_eq!(x, y, epsilon = 1e-12);
        }
    }

    #[test]
    fn the_node_search_shapes_the_realization() {
        // With no data, max_nodes = 1 and max_nodes = 12 krige from
        // different neighbourhoods, so the same seed gives different fields;
        // and the node search must belong to the grid.
        let g = line(30, 0.0, 1.0);
        let m = exp_cov(10.0);
        let (one, twelve) = (nodes(&g, &m, 1), nodes(&g, &m, 12));
        let run = |ns: &NodeSearch| {
            sgsim(
                &[],
                &g,
                &m,
                0.0,
                Searches {
                    data: &near(),
                    nodes: ns,
                },
                &mut rng(9),
            )
        };
        assert_ne!(run(&one).unwrap(), run(&twelve).unwrap());
        let other = line(31, 0.0, 1.0);
        let wrong = nodes(&other, &m, 12);
        assert!(run(&wrong).is_err(), "node search for another grid");
    }

    #[test]
    fn coincident_nodes_take_the_datum_and_are_not_reused() {
        // A datum on node 2: the node holds it exactly, and kriging its
        // neighbours does not see that location twice (which would make
        // the system singular).
        let data = [s(2.0, 0.7)];
        let g = line(5, 0.0, 1.0);
        let m = exp_cov(4.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &DataSearch::unlimited(),
            nodes: &ns,
        };
        for seed in 0..20 {
            let r = sgsim(&data, &g, &m, 0.0, searches, &mut rng(seed)).unwrap();
            assert_abs_diff_eq!(r[2], 0.7);
        }
    }

    #[test]
    fn sisim_constant_data_gives_constant_field() {
        let data = [s(0.0, 0.5)];
        let g = line(2, 1.0, 1.0);
        let m = exp_cov(3.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &near(),
            nodes: &ns,
        };
        for seed in 0..5 {
            assert_eq!(
                sisim(&data, &g, 1.0, &m, searches, &mut rng(seed)).unwrap(),
                vec![true, true]
            );
            assert_eq!(
                sisim(&data, &g, 0.0, &m, searches, &mut rng(seed)).unwrap(),
                vec![false, false]
            );
        }
    }

    #[test]
    fn sisim_far_node_uses_global_proportion() {
        // Half the data below threshold; a node far from everything is a fair coin.
        let data = [s(0.0, 0.0), s(1.0, 2.0)];
        let g = line(1, 1e6, 1.0);
        let m = exp_cov(1.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &near(),
            nodes: &ns,
        };
        let mut r = rng(11);
        let n = 2000_u32;
        let hits: u32 = (0..n)
            .map(|_| u32::from(sisim(&data, &g, 1.0, &m, searches, &mut r).unwrap()[0]))
            .sum();
        let freq = f64::from(hits) / f64::from(n);
        assert!((freq - 0.5).abs() < 0.05, "freq {freq}");
    }

    #[test]
    fn sisim_needs_data() {
        let g = line(1, 0.0, 1.0);
        let m = exp_cov(1.0);
        let ns = nodes(&g, &m, 12);
        let searches = Searches {
            data: &near(),
            nodes: &ns,
        };
        assert!(
            sisim(&[], &g, 1.0, &m, searches, &mut rng(0)).is_err(),
            "no data"
        );
    }
}
