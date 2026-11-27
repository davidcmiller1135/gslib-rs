//! Sequential Gaussian and sequential indicator simulation.
//!
//! Each node is kriged from the data and previously simulated nodes chosen by
//! a [`Search`] neighbourhood. With [`Search::unlimited`], cost grows as
//! O(N⁴) in the number of nodes; a bounded search keeps it near O(N²).

use rand::seq::SliceRandom;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, StandardNormal};

use crate::error::{Error, Result};
use crate::geometry::{COINCIDENT, Point3};
use crate::kriging::simple_kriging;
use crate::model::CovarianceModel;
use crate::num::mean;
use crate::sample::Sample;
use crate::search::Search;

/// Visit `nodes` in random order, computing each value from the conditioning
/// set so far (data plus earlier nodes), and return values in node order.
///
/// A node that coincides with a conditioning location takes that location's
/// value directly: re-kriging there would only duplicate a row of the system.
fn sequential<R, F>(
    data: &[Sample],
    nodes: &[Point3],
    search: Search,
    rng: &mut R,
    mut draw: F,
) -> Result<Vec<f64>>
where
    R: Rng + ?Sized,
    F: FnMut(&[Sample], Point3, &mut R) -> Result<f64>,
{
    let mut path: Vec<(usize, Point3)> = nodes.iter().copied().enumerate().collect();
    path.shuffle(rng);

    let mut conditioning = data.to_vec();
    conditioning.reserve(nodes.len());
    let mut simulated = Vec::with_capacity(nodes.len());
    for (index, node) in path {
        let existing = conditioning
            .iter()
            .find(|s| s.location().distance(node) < COINCIDENT)
            .map(Sample::value);
        let value = if let Some(v) = existing {
            v
        } else {
            let v = draw(&search.select(&conditioning, node), node, rng)?;
            // Only new locations join the conditioning set: a duplicate
            // location would make every later kriging system singular.
            conditioning.push(Sample::new(node, v)?);
            v
        };
        simulated.push((index, value));
    }
    simulated.sort_by_key(|&(i, _)| i);
    Ok(simulated.into_iter().map(|(_, v)| v).collect())
}

/// One realization of **sequential Gaussian simulation**.
///
/// `data` must already be in normal-score units with the given stationary
/// `mean` (0 after [`crate::transform::normal_scores`]), and `model` should
/// have a sill of 1. Each node draws from N(SK estimate, SK variance).
///
/// # Errors
/// [`Error::SingularSystem`] if a kriging system fails (e.g. duplicate data
/// locations).
pub fn sgsim<R>(
    data: &[Sample],
    nodes: &[Point3],
    model: &CovarianceModel,
    mean: f64,
    search: Search,
    rng: &mut R,
) -> Result<Vec<f64>>
where
    R: Rng + ?Sized,
{
    sequential(data, nodes, search, rng, |conditioning, node, rng| {
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
/// [`Error::InvalidInput`] if `data` is empty (p₀ is undefined), or
/// [`Error::SingularSystem`] if a kriging system fails.
pub fn sisim<R>(
    data: &[Sample],
    nodes: &[Point3],
    threshold: f64,
    model: &CovarianceModel,
    search: Search,
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
        nodes,
        search,
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
    use approx::assert_abs_diff_eq;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }
    fn at(x: f64) -> Point3 {
        Point3::new(x, 0.0, 0.0)
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
    fn near() -> Search {
        Search::new(12, f64::INFINITY).unwrap()
    }
    fn rng(seed: u64) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(seed)
    }

    #[test]
    fn sgsim_reproduces_data_at_data_locations() {
        let data = [s(0.0, 1.5), s(5.0, -0.5)];
        let nodes = [at(0.0), at(2.0), at(5.0), at(0.0)];
        let r = sgsim(
            &data,
            &nodes,
            &exp_cov(3.0),
            0.0,
            Search::unlimited(),
            &mut rng(1),
        )
        .unwrap();
        assert_abs_diff_eq!(r[0], 1.5);
        assert_abs_diff_eq!(r[2], -0.5);
        assert_abs_diff_eq!(r[3], 1.5);
    }

    #[test]
    fn sgsim_same_seed_same_realization() {
        let data = [s(0.0, 1.0)];
        let nodes: Vec<Point3> = (1..20).map(|i| at(f64::from(i))).collect();
        let a = sgsim(&data, &nodes, &exp_cov(5.0), 0.0, near(), &mut rng(7)).unwrap();
        let b = sgsim(&data, &nodes, &exp_cov(5.0), 0.0, near(), &mut rng(7)).unwrap();
        let c = sgsim(&data, &nodes, &exp_cov(5.0), 0.0, near(), &mut rng(8)).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn sgsim_unconditional_moments_are_standard_normal() {
        // Nodes far apart ⇒ independent N(0,1) draws.
        let nodes: Vec<Point3> = (0..400).map(|i| at(f64::from(i) * 1000.0)).collect();
        let r = sgsim(&[], &nodes, &exp_cov(1.0), 0.0, near(), &mut rng(3)).unwrap();
        let m = mean(&r).unwrap();
        let var = r.iter().map(|v| (v - m).powi(2)).sum::<f64>() / 399.0;
        assert!(m.abs() < 0.2, "mean {m}");
        assert!((var - 1.0).abs() < 0.2, "variance {var}");
    }

    #[test]
    fn bounded_search_matches_unlimited_when_it_covers_everything() {
        let data = [s(0.0, 1.0), s(4.0, -1.0)];
        let nodes: Vec<Point3> = (1..8).map(|i| at(f64::from(i) * 0.5)).collect();
        let a = sgsim(
            &data,
            &nodes,
            &exp_cov(3.0),
            0.0,
            Search::unlimited(),
            &mut rng(5),
        )
        .unwrap();
        let b = sgsim(&data, &nodes, &exp_cov(3.0), 0.0, near(), &mut rng(5)).unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert_abs_diff_eq!(x, y, epsilon = 1e-12);
        }
    }

    #[test]
    fn sisim_constant_data_gives_constant_field() {
        let data = [s(0.0, 0.5)];
        let nodes = [at(1.0), at(2.0)];
        for seed in 0..5 {
            assert_eq!(
                sisim(&data, &nodes, 1.0, &exp_cov(3.0), near(), &mut rng(seed)).unwrap(),
                vec![true, true]
            );
            assert_eq!(
                sisim(&data, &nodes, 0.0, &exp_cov(3.0), near(), &mut rng(seed)).unwrap(),
                vec![false, false]
            );
        }
    }

    #[test]
    fn sisim_far_node_uses_global_proportion() {
        // Half the data below threshold; a node far from everything is a fair coin.
        let data = [s(0.0, 0.0), s(1.0, 2.0)];
        let mut r = rng(11);
        let n = 2000_u32;
        let hits: u32 = (0..n)
            .map(|_| {
                u32::from(sisim(&data, &[at(1e6)], 1.0, &exp_cov(1.0), near(), &mut r).unwrap()[0])
            })
            .sum();
        let freq = f64::from(hits) / f64::from(n);
        assert!((freq - 0.5).abs() < 0.05, "freq {freq}");
    }

    #[test]
    fn sisim_needs_data() {
        assert!(
            sisim(&[], &[at(0.0)], 1.0, &exp_cov(1.0), near(), &mut rng(0)).is_err(),
            "no data"
        );
    }
}
