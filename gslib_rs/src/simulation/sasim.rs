//! Simulated annealing simulation (GSLIB `sasim`).
//!
//! The field starts as random draws from the data distribution and is
//! perturbed one node at a time, accepting changes by the Metropolis rule,
//! until it matches the chosen objectives (histogram, variogram, indicator
//! variogram, correlation with secondary data).
//!
//! The objective is recomputed in full after every perturbation (O(N²) per
//! trial for variogram objectives). GSLIB updates it incrementally; that is
//! on the roadmap and will not change results.

use rand::seq::IndexedRandom;
use rand::{Rng, RngExt};

use super::grid::Grid;
use crate::error::{Error, Result, check_lengths};
use crate::geometry::Point3;
use crate::num::{floor_index, mean, to_f64};
use crate::sample::Sample;
use crate::variogram::{Lag, LagClasses, semivariogram};

/// Match the histogram of the data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HistogramObjective {
    /// Relative weight of this component.
    pub weight: f64,
    /// Number of histogram classes.
    pub bins: usize,
    /// Class range; `None` uses the data minimum and maximum.
    pub range: Option<(f64, f64)>,
}

/// Match the experimental variogram of the data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariogramObjective {
    /// Relative weight of this component.
    pub weight: f64,
    /// Lag classes for both the data and the simulated field.
    pub lags: LagClasses,
}

/// Match the experimental indicator variogram of the data at a threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndicatorVariogramObjective {
    /// Relative weight of this component.
    pub weight: f64,
    /// Lag classes for both the data and the simulated field.
    pub lags: LagClasses,
    /// Values below this are indicator 1.
    pub threshold: f64,
}

/// Match the data's correlation with a secondary variable known at every node.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationObjective {
    /// Relative weight of this component.
    pub weight: f64,
    /// Secondary values, one per grid node in grid order.
    pub secondary: Vec<f64>,
}

/// The objective function components. At least one must be set.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Objectives {
    /// Histogram reproduction.
    pub histogram: Option<HistogramObjective>,
    /// Variogram reproduction.
    pub variogram: Option<VariogramObjective>,
    /// Indicator variogram reproduction.
    pub indicator_variogram: Option<IndicatorVariogramObjective>,
    /// Correlation with secondary variables.
    pub correlations: Vec<CorrelationObjective>,
}

/// The annealing schedule (GSLIB `t0, redfac, ka, k, num, Omin`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Schedule {
    /// Starting temperature, > 0.
    pub initial_temperature: f64,
    /// Temperature multiplier per level, in (0, 1).
    pub reduction_factor: f64,
    /// Move to the next level after this many accepted perturbations.
    pub accepted_per_level: usize,
    /// Or after this many attempted perturbations.
    pub trials_per_level: usize,
    /// Maximum number of temperature levels.
    pub max_levels: usize,
    /// Stop once the temperature drops below this.
    pub min_temperature: f64,
}

/// Whether data falling inside grid cells are frozen at their values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Conditioning {
    /// Data define the targets only; every node is perturbed.
    Unconditional,
    /// Each datum is assigned to the node whose cell contains it and is
    /// never perturbed. Data outside the grid are ignored.
    Conditional,
}

/// One component's target and weight, resolved from the data.
#[derive(Debug, Clone)]
enum Component {
    Histogram {
        weight: f64,
        bins: usize,
        min: f64,
        max: f64,
        target: Vec<f64>,
    },
    Variogram {
        weight: f64,
        lags: LagClasses,
        threshold: Option<f64>,
        target: Vec<Lag>,
    },
    Correlation {
        weight: f64,
        secondary: Vec<f64>,
        target: f64,
    },
}

impl Component {
    const fn weight(&self) -> f64 {
        match *self {
            Self::Histogram { weight, .. }
            | Self::Variogram { weight, .. }
            | Self::Correlation { weight, .. } => weight,
        }
    }

    /// Unweighted mismatch between `field` and this component's target.
    fn mismatch(&self, field: &[f64], nodes: &[Point3]) -> f64 {
        match self {
            Self::Histogram {
                bins,
                min,
                max,
                target,
                ..
            } => histogram(field, *bins, *min, *max)
                .iter()
                .zip(target)
                .map(|(p, t)| (p - t).powi(2))
                .sum(),
            Self::Variogram {
                lags,
                threshold,
                target,
                ..
            } => {
                let values = threshold.map_or_else(|| field.to_vec(), |t| indicators(field, t));
                semivariogram(nodes, &values, lags)
                    .iter()
                    .zip(target)
                    .filter_map(|(sim, data)| Some((sim.gamma? - data.gamma?).powi(2)))
                    .sum()
            }
            Self::Correlation {
                secondary, target, ..
            } => {
                let pairs: Vec<(f64, f64)> = field
                    .iter()
                    .copied()
                    .zip(secondary.iter().copied())
                    .collect();
                // A constant field has no correlation; treat it as zero.
                (correlation(&pairs).unwrap_or(0.0) - target).powi(2)
            }
        }
    }
}

/// One realization of **simulated annealing** on `grid`.
///
/// Each component's mismatch is divided by its value for the initial field,
/// so weights are relative (as in GSLIB).
///
/// # Errors
/// [`Error::InvalidInput`] for empty data, no objectives, an invalid schedule
/// or weight, or a correlation target that cannot be computed;
/// [`Error::LengthMismatch`] if a secondary variable does not match the grid.
pub fn sasim<R>(
    data: &[Sample],
    grid: &Grid,
    objectives: &Objectives,
    schedule: &Schedule,
    conditioning: Conditioning,
    rng: &mut R,
) -> Result<Vec<f64>>
where
    R: Rng + ?Sized,
{
    validate_schedule(schedule)?;
    let values: Vec<f64> = data.iter().map(Sample::value).collect();
    if values.is_empty() {
        return Err(Error::invalid("data", "sasim needs at least one datum"));
    }
    let nodes = grid.nodes();
    let on_grid: Vec<(usize, f64)> = data
        .iter()
        .filter_map(|s| grid.node_index(s.location()).map(|i| (i, s.value())))
        .collect();
    let components = resolve(data, &values, &on_grid, objectives, nodes.len())?;

    // Initial field: random draws from the data, with frozen data in place.
    let mut field = (0..nodes.len())
        .map(|_| values.choose(rng).copied())
        .collect::<Option<Vec<f64>>>()
        .ok_or_else(|| Error::invalid("data", "sasim needs at least one datum"))?;
    let mut frozen = vec![false; nodes.len()];
    if conditioning == Conditioning::Conditional {
        for &(i, v) in &on_grid {
            if let (Some(slot), Some(f)) = (field.get_mut(i), frozen.get_mut(i)) {
                *slot = v;
                *f = true;
            }
        }
    }
    let free: Vec<usize> = frozen
        .iter()
        .enumerate()
        .filter_map(|(i, &f)| (!f).then_some(i))
        .collect();

    // Normalise each component by its initial mismatch (1 if already perfect).
    let scales: Vec<f64> = components
        .iter()
        .map(|c| {
            let m = c.mismatch(&field, &nodes);
            if m > 0.0 { c.weight() / m } else { c.weight() }
        })
        .collect();
    let energy = |state: &[f64]| -> f64 {
        components
            .iter()
            .zip(&scales)
            .map(|(c, s)| s * c.mismatch(state, &nodes))
            .sum()
    };

    let mut current = energy(&field);
    let mut temperature = schedule.initial_temperature;
    for _ in 0..schedule.max_levels {
        let mut accepted = 0;
        let mut trials = 0;
        while accepted < schedule.accepted_per_level && trials < schedule.trials_per_level {
            trials += 1;
            let (Some(&i), Some(&proposal)) = (free.choose(rng), values.choose(rng)) else {
                return Ok(field); // nothing can be perturbed
            };
            let Some(old) = field.get(i).copied() else {
                continue;
            };
            set(&mut field, i, proposal);
            let candidate = energy(&field);
            let delta = candidate - current;
            if delta <= 0.0 || rng.random::<f64>() < (-delta / temperature).exp() {
                current = candidate;
                accepted += 1;
            } else {
                set(&mut field, i, old);
            }
        }
        temperature *= schedule.reduction_factor;
        if temperature < schedule.min_temperature || current <= 0.0 {
            break;
        }
    }
    Ok(field)
}

fn set(field: &mut [f64], i: usize, value: f64) {
    if let Some(slot) = field.get_mut(i) {
        *slot = value;
    }
}

fn validate_schedule(s: &Schedule) -> Result<()> {
    if !(s.initial_temperature.is_finite() && s.initial_temperature > 0.0) {
        return Err(Error::invalid(
            "schedule",
            "initial temperature must be finite and positive",
        ));
    }
    if !(s.reduction_factor > 0.0 && s.reduction_factor < 1.0) {
        return Err(Error::invalid(
            "schedule",
            "reduction factor must be in (0, 1)",
        ));
    }
    if !(s.min_temperature.is_finite() && s.min_temperature >= 0.0) {
        return Err(Error::invalid(
            "schedule",
            "minimum temperature must be finite and >= 0",
        ));
    }
    Ok(())
}

fn check_weight(weight: f64) -> Result<f64> {
    if weight.is_finite() && weight > 0.0 {
        Ok(weight)
    } else {
        Err(Error::invalid(
            "objective weight",
            format!("must be finite and positive, got {weight}"),
        ))
    }
}

/// Compute every requested component's target from the data.
fn resolve(
    data: &[Sample],
    values: &[f64],
    on_grid: &[(usize, f64)],
    objectives: &Objectives,
    grid_len: usize,
) -> Result<Vec<Component>> {
    let points: Vec<Point3> = data.iter().map(Sample::location).collect();
    let mut out = Vec::new();

    if let Some(h) = objectives.histogram {
        if h.bins == 0 {
            return Err(Error::invalid(
                "histogram objective",
                "needs at least one bin",
            ));
        }
        let (min, max) = h.range.unwrap_or_else(|| {
            values
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                    (lo.min(v), hi.max(v))
                })
        });
        if !(min.is_finite() && max.is_finite() && min <= max) {
            return Err(Error::invalid(
                "histogram range",
                "must be finite with min <= max",
            ));
        }
        out.push(Component::Histogram {
            weight: check_weight(h.weight)?,
            bins: h.bins,
            min,
            max,
            target: histogram(values, h.bins, min, max),
        });
    }
    if let Some(v) = objectives.variogram {
        out.push(Component::Variogram {
            weight: check_weight(v.weight)?,
            lags: v.lags,
            threshold: None,
            target: semivariogram(&points, values, &v.lags),
        });
    }
    if let Some(v) = objectives.indicator_variogram {
        out.push(Component::Variogram {
            weight: check_weight(v.weight)?,
            lags: v.lags,
            threshold: Some(v.threshold),
            target: semivariogram(&points, &indicators(values, v.threshold), &v.lags),
        });
    }
    for c in &objectives.correlations {
        check_lengths("secondary variable/grid", c.secondary.len(), grid_len)?;
        // Pair each on-grid datum with the secondary value at *its* node.
        let pairs: Vec<(f64, f64)> = on_grid
            .iter()
            .filter_map(|&(i, v)| c.secondary.get(i).map(|&y| (v, y)))
            .collect();
        let target = correlation(&pairs).ok_or_else(|| {
            Error::invalid(
                "correlation objective",
                "needs at least two on-grid data with non-zero variance in both variables",
            )
        })?;
        out.push(Component::Correlation {
            weight: check_weight(c.weight)?,
            secondary: c.secondary.clone(),
            target,
        });
    }
    if out.is_empty() {
        return Err(Error::invalid(
            "objectives",
            "at least one objective is required",
        ));
    }
    Ok(out)
}

/// Class proportions of `values` in `bins` equal classes over `[min, max]`;
/// values outside fall in the end classes.
fn histogram(values: &[f64], bins: usize, min: f64, max: f64) -> Vec<f64> {
    let mut counts = vec![0_usize; bins];
    let span = max - min;
    for &v in values {
        let class = if span > 0.0 {
            floor_index((v - min) / span * to_f64(bins)).map_or(0, |c| c.min(bins - 1))
        } else {
            0
        };
        if let Some(c) = counts.get_mut(class) {
            *c += 1;
        }
    }
    let n = to_f64(values.len().max(1));
    counts.into_iter().map(|c| to_f64(c) / n).collect()
}

fn indicators(values: &[f64], threshold: f64) -> Vec<f64> {
    values
        .iter()
        .map(|&v| if v < threshold { 1.0 } else { 0.0 })
        .collect()
}

/// Pearson correlation, or `None` with fewer than two pairs or zero variance.
fn correlation(pairs: &[(f64, f64)]) -> Option<f64> {
    if pairs.len() < 2 {
        return None;
    }
    let xs: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = pairs.iter().map(|p| p.1).collect();
    let (mx, my) = (mean(&xs)?, mean(&ys)?);
    let (sxy, sxx, syy) = pairs
        .iter()
        .fold((0.0, 0.0, 0.0), |(sxy, sxx, syy), &(x, y)| {
            let (dx, dy) = (x - mx, y - my);
            (
                dx.mul_add(dy, sxy),
                dx.mul_add(dx, sxx),
                dy.mul_add(dy, syy),
            )
        });
    let denom = (sxx * syy).sqrt();
    (denom > 0.0).then(|| sxy / denom)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::grid::Axis;
    use approx::assert_abs_diff_eq;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn line_grid(n: usize) -> Grid {
        let one = Axis {
            count: 1,
            first_center: 0.0,
            spacing: 1.0,
        };
        Grid::new(
            Axis {
                count: n,
                first_center: 0.0,
                spacing: 1.0,
            },
            one,
            one,
        )
        .unwrap()
    }
    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }
    fn schedule() -> Schedule {
        Schedule {
            initial_temperature: 1.0,
            reduction_factor: 0.5,
            accepted_per_level: 50,
            trials_per_level: 500,
            max_levels: 20,
            min_temperature: 1e-6,
        }
    }
    fn histogram_only() -> Objectives {
        Objectives {
            histogram: Some(HistogramObjective {
                weight: 1.0,
                bins: 2,
                range: None,
            }),
            ..Objectives::default()
        }
    }
    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(42)
    }

    #[test]
    fn constant_data_gives_constant_field() {
        let field = sasim(
            &[s(0.0, 2.5)],
            &line_grid(5),
            &histogram_only(),
            &schedule(),
            Conditioning::Unconditional,
            &mut rng(),
        )
        .unwrap();
        assert!(field.iter().all(|&v| (v - 2.5).abs() < 1e-12), "{field:?}");
    }

    #[test]
    fn histogram_range_defaults_to_data_range() {
        // Regression: hist_range = None used to panic.
        let data = [s(0.0, 0.0), s(1.0, 1.0)];
        let field = sasim(
            &data,
            &line_grid(20),
            &histogram_only(),
            &schedule(),
            Conditioning::Unconditional,
            &mut rng(),
        )
        .unwrap();
        let ones = field.iter().filter(|&&v| v > 0.5).count();
        assert_eq!(
            ones, 10,
            "histogram objective should reach a 50/50 split: {field:?}"
        );
    }

    #[test]
    fn conditioning_freezes_data_at_their_nodes() {
        // Regression: conditioned nodes used to keep random initial values.
        let data = [s(0.0, 7.0), s(3.0, 1.0)];
        for seed in 0..5 {
            let mut r = ChaCha8Rng::seed_from_u64(seed);
            let field = sasim(
                &data,
                &line_grid(4),
                &histogram_only(),
                &schedule(),
                Conditioning::Conditional,
                &mut r,
            )
            .unwrap();
            assert_abs_diff_eq!(field[0], 7.0);
            assert_abs_diff_eq!(field[3], 1.0);
        }
    }

    #[test]
    fn correlation_target_pairs_data_with_their_own_nodes() {
        // Regression: targets were built from a HashSet in arbitrary order.
        // Data equal the secondary at their nodes, so the target is exactly 1.
        let secondary: Vec<f64> = (0..10).map(f64::from).collect();
        let data: Vec<Sample> = [2, 5, 7, 9]
            .iter()
            .map(|&i| s(f64::from(i), f64::from(i)))
            .collect();
        let objectives = Objectives {
            correlations: vec![CorrelationObjective {
                weight: 1.0,
                secondary,
            }],
            ..Objectives::default()
        };
        let on_grid: Vec<(usize, f64)> = data
            .iter()
            .map(|d| (line_grid(10).node_index(d.location()).unwrap(), d.value()))
            .collect();
        let values: Vec<f64> = data.iter().map(Sample::value).collect();
        let comps = resolve(&data, &values, &on_grid, &objectives, 10).unwrap();
        let Component::Correlation { target, .. } = &comps[0] else {
            panic!("expected a correlation component, got {:?}", comps[0]);
        };
        assert_abs_diff_eq!(*target, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn rejects_bad_configuration() {
        let data = [s(0.0, 1.0)];
        let g = line_grid(3);
        assert!(
            sasim(
                &[],
                &g,
                &histogram_only(),
                &schedule(),
                Conditioning::Unconditional,
                &mut rng()
            )
            .is_err(),
            "no data"
        );
        assert!(
            sasim(
                &data,
                &g,
                &Objectives::default(),
                &schedule(),
                Conditioning::Unconditional,
                &mut rng()
            )
            .is_err(),
            "no objectives"
        );
        let bad = Schedule {
            reduction_factor: 1.5,
            ..schedule()
        };
        assert!(
            sasim(
                &data,
                &g,
                &histogram_only(),
                &bad,
                Conditioning::Unconditional,
                &mut rng()
            )
            .is_err(),
            "schedule"
        );
        let wrong_len = Objectives {
            correlations: vec![CorrelationObjective {
                weight: 1.0,
                secondary: vec![0.0; 2],
            }],
            ..Objectives::default()
        };
        assert!(
            sasim(
                &data,
                &g,
                &wrong_len,
                &schedule(),
                Conditioning::Unconditional,
                &mut rng()
            )
            .is_err(),
            "secondary length"
        );
    }

    #[test]
    fn histogram_proportions() {
        assert_eq!(
            histogram(&[0.0, 0.4, 0.6, 1.0, 5.0], 2, 0.0, 1.0),
            vec![0.4, 0.6]
        );
        assert_eq!(histogram(&[3.0, 3.0], 3, 3.0, 3.0), vec![1.0, 0.0, 0.0]);
    }
}
