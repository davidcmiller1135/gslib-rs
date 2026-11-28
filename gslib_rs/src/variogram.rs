//! Experimental (empirical) semivariograms.
//!
//! Two ways to group pairs by separation, both behind [`LagScheme`]:
//! - [`LagClasses`]: equal-width, non-overlapping classes over `[0, max]`;
//! - [`CentredLags`]: GSLIB `gamv`'s classes, centred on multiples of a lag
//!   spacing with a tolerance (they may overlap), plus a zero-distance class.
//!
//! # Pair enumeration (`docs/spec/parallelism.md`, D-031, D-032)
//!
//! Pairs are not formed by a full double loop. Rows are taken in fixed
//! blocks of 256 rows; each row's partners (later rows only, so
//! each unordered pair once) come from a super-block search within the lag
//! scheme's [`LagScheme::max_distance`], widened by a relative 1e-9 so no
//! pair at the cut-off is lost to rounding. The lag scheme still decides
//! every pair, so exactly the same pairs are counted as by a full loop.
//! Partners arrive in block order: grid blocks in linear order, rows
//! ascending within a block (D-034). Each row block sums sequentially, and
//! row-block sums are combined in order.
//!
//! The block layout and partner order are part of the numerical definition:
//! the sums are the same for any thread count, and differ from a plain
//! double loop only in the last bits (the `gamv` fixtures' 1e-12 holds).

use std::ops::Range;

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::num::{floor_index, to_f64};
use crate::sample::Sample;
use crate::search::Blocks;

mod directional;

pub use directional::{Direction, DirectionWindow, PairOrder, Separation};

/// A way of assigning pair separations to lag classes.
pub trait LagScheme {
    /// Number of classes.
    fn len(&self) -> usize;

    /// Whether there are no classes (never true for a validated scheme).
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The nominal distance of a class.
    fn center(&self, class: usize) -> f64;

    /// The classes a separation `h` belongs to: a possibly empty range,
    /// since classes may overlap.
    fn classes(&self, h: f64) -> Range<usize>;

    /// The largest separation any class admits: every `h` with a non-empty
    /// [`LagScheme::classes`] satisfies `h ≤ max_distance()`. Pair
    /// enumeration skips pairs farther apart (D-032).
    fn max_distance(&self) -> f64;
}

/// Equal-width distance classes covering `[0, max_distance]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LagClasses {
    width: f64,
    count: usize,
    max_distance: f64,
}

impl LagClasses {
    /// `count` classes of width `max_distance / count`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `count` is zero or `max_distance` is not
    /// finite and positive.
    pub fn new(max_distance: f64, count: usize) -> Result<Self> {
        if count == 0 {
            return Err(Error::invalid("lag count", "must be at least 1"));
        }
        if !(max_distance.is_finite() && max_distance > 0.0) {
            return Err(Error::invalid(
                "maximum lag distance",
                format!("must be finite and positive, got {max_distance}"),
            ));
        }
        Ok(Self {
            width: max_distance / to_f64(count),
            count,
            max_distance,
        })
    }

    /// Number of classes.
    #[must_use]
    pub const fn count(&self) -> usize {
        self.count
    }
}

impl LagScheme for LagClasses {
    fn len(&self) -> usize {
        self.count
    }

    fn center(&self, class: usize) -> f64 {
        self.width * (to_f64(class) + 0.5)
    }

    fn max_distance(&self) -> f64 {
        self.max_distance
    }

    fn classes(&self, h: f64) -> Range<usize> {
        if h > self.max_distance {
            return 0..0;
        }
        // h == max_distance belongs to the last class.
        floor_index(h / self.width).map_or(0..0, |i| {
            let c = i.min(self.count - 1);
            c..c + 1
        })
    }
}

/// GSLIB `gamv`'s lag classes (`vendor/gslib-2.0/src/GAMV/GAMV.FOR`, F-008).
///
/// Class 0 holds coincident pairs (separation ≤ 1e-20, GSLIB's `EPSLON`).
/// Class `k + 1`, for `k = 0..=count`, holds separations within `tolerance`
/// of `k · spacing`, inclusive at both ends, so classes overlap when
/// `tolerance > spacing / 2` and a pair then counts in several. Pairs farther
/// apart than `(count + 0.5) · spacing` are ignored.
///
/// Unlike GSLIB, a datum is never paired with itself: GSLIB's class 0 also
/// counts every datum against itself, diluting γ there towards zero
/// (ROADMAP §2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CentredLags {
    count: usize,
    spacing: f64,
    tolerance: f64,
}

/// GSLIB `gamv`'s `EPSLON`: separations up to this are "zero".
const ZERO_SEPARATION: f64 = 1e-20;

impl CentredLags {
    /// `count` lags of `spacing` (GSLIB `nlag`, `xlag`) with half-width
    /// `tolerance` (`xltol`).
    ///
    /// GSLIB replaces a non-positive tolerance with `spacing / 2`; here it is
    /// an error (explicit failure, `AGENTS.md`).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `count` is zero, or `spacing` or
    /// `tolerance` is not finite and positive.
    pub fn new(count: usize, spacing: f64, tolerance: f64) -> Result<Self> {
        if count == 0 {
            return Err(Error::invalid("lag count", "must be at least 1"));
        }
        for (what, v) in [("lag spacing", spacing), ("lag tolerance", tolerance)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(Error::invalid(
                    what,
                    format!("must be finite and positive, got {v}"),
                ));
            }
        }
        Ok(Self {
            count,
            spacing,
            tolerance,
        })
    }
}

impl LagScheme for CentredLags {
    fn len(&self) -> usize {
        self.count + 2
    }

    fn center(&self, class: usize) -> f64 {
        class
            .checked_sub(1)
            .map_or(0.0, |k| self.spacing * to_f64(k))
    }

    /// GSLIB's cut-off `(nlag + 0.5) · xlag`.
    fn max_distance(&self) -> f64 {
        (to_f64(self.count) + 0.5) * self.spacing
    }

    fn classes(&self, h: f64) -> Range<usize> {
        let limit = self.max_distance();
        if h.is_nan() || h > limit {
            return 0..0;
        }
        if h <= ZERO_SEPARATION {
            return 0..1;
        }
        let within = |k: usize| {
            let c = self.spacing * to_f64(k);
            h >= c - self.tolerance && h <= c + self.tolerance
        };
        let first = (0..=self.count).find(|&k| within(k));
        let last = (0..=self.count).rev().find(|&k| within(k));
        match (first, last) {
            (Some(a), Some(b)) => a + 1..b + 2,
            _ => 0..0,
        }
    }
}

/// One lag of an experimental variogram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lag {
    /// Nominal distance of the class.
    pub center: f64,
    /// Mean separation of the pairs in the class (`None` if there are none).
    pub mean_distance: Option<f64>,
    /// Number of pairs in the class.
    pub pairs: usize,
    /// Semivariance γ, or `None` if the class has no pairs.
    pub gamma: Option<f64>,
    /// Mean of the values in the class's pairs (both ends of every pair);
    /// GSLIB's head and tail means for an omnidirectional variogram.
    pub pair_mean: Option<f64>,
}

/// Omnidirectional experimental semivariogram of `samples`.
#[must_use]
pub fn experimental_variogram<L: LagScheme + ?Sized>(samples: &[Sample], lags: &L) -> Vec<Lag> {
    let points: Vec<Point3> = samples.iter().map(Sample::location).collect();
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    semivariogram(&points, &values, lags)
}

/// Rows per block in the pair enumeration (D-031). A constant of the
/// algorithm, never derived from the thread count: changing it changes the
/// last bits of the sums, and is checked against the `gamv` fixtures.
pub(crate) const ROWS_PER_BLOCK: usize = 256;

/// Relative widening of the pruning radius (D-032).
const PRUNE_MARGIN: f64 = 1e-9;

/// Running sums for one lag class.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct LagSums {
    gamma: f64,
    distance: f64,
    sum: f64,
    pairs: usize,
}

/// Everything one block of rows needs: the data, the lag scheme, and the
/// block grid used to find partners.
pub(crate) struct PairPlan<'a, L: LagScheme + ?Sized> {
    points: &'a [Point3],
    values: &'a [f64],
    lags: &'a L,
    blocks: Blocks,
    /// The pruning radius: the scheme's largest separation, widened.
    reach: f64,
}

impl<'a, L: LagScheme + ?Sized> PairPlan<'a, L> {
    /// A plan over `points` and `values` (zipped; any excess in the longer
    /// slice is ignored).
    pub(crate) fn new(points: &'a [Point3], values: &'a [f64], lags: &'a L) -> Self {
        let n = points.len().min(values.len());
        let (points, values) = (
            points.get(..n).unwrap_or(&[]),
            values.get(..n).unwrap_or(&[]),
        );
        // A lag scheme without a usable cut-off is searched without one; the
        // scheme still decides every pair.
        let widened = lags.max_distance() * (1.0 + PRUNE_MARGIN);
        let reach = if widened > 0.0 {
            widened
        } else {
            f64::INFINITY
        };
        Self {
            points,
            values,
            lags,
            blocks: Blocks::new(points),
            reach,
        }
    }

    /// The number of row blocks.
    pub(crate) const fn block_count(&self) -> usize {
        self.points.len().div_ceil(ROWS_PER_BLOCK)
    }

    /// The sums contributed by the pairs whose first row is in block
    /// `block`. Partners come in block order (D-034).
    pub(crate) fn block_sums(&self, block: usize) -> Vec<LagSums> {
        let mut acc = vec![LagSums::default(); self.lags.len()];
        let start = block.saturating_mul(ROWS_PER_BLOCK);
        let end = start.saturating_add(ROWS_PER_BLOCK).min(self.points.len());
        for i in start..end {
            let (Some(&head), Some(&head_value)) = (self.points.get(i), self.values.get(i)) else {
                continue;
            };
            self.blocks.visit_within(head, self.reach, |row, tail| {
                if row <= i {
                    return;
                }
                let Some(&tail_value) = self.values.get(row) else {
                    return;
                };
                let h = head.distance(tail);
                let d = head_value - tail_value;
                for class in self.lags.classes(h) {
                    if let Some(a) = acc.get_mut(class) {
                        a.gamma = (0.5 * d).mul_add(d, a.gamma);
                        a.distance += h;
                        a.sum += head_value + tail_value;
                        a.pairs += 1;
                    }
                }
            });
        }
        acc
    }

    /// Combine block sums, in block order, into the experimental variogram.
    pub(crate) fn finish(&self, blocks: impl IntoIterator<Item = Vec<LagSums>>) -> Vec<Lag> {
        let mut total = vec![LagSums::default(); self.lags.len()];
        for sums in blocks {
            for (t, b) in total.iter_mut().zip(sums) {
                t.gamma += b.gamma;
                t.distance += b.distance;
                t.sum += b.sum;
                t.pairs += b.pairs;
            }
        }
        lags_from(self.lags, total)
    }
}

fn lags_from<L: LagScheme + ?Sized>(lags: &L, sums: Vec<LagSums>) -> Vec<Lag> {
    sums.into_iter()
        .enumerate()
        .map(|(class, a)| {
            let n = (a.pairs > 0).then(|| to_f64(a.pairs));
            Lag {
                center: lags.center(class),
                mean_distance: n.map(|n| a.distance / n),
                pairs: a.pairs,
                gamma: n.map(|n| a.gamma / n),
                pair_mean: n.map(|n| a.sum / (2.0 * n)),
            }
        })
        .collect()
}

/// Experimental semivariogram of `values` at `points` (zipped; any excess in
/// the longer slice is ignored). Each unordered pair of distinct data counts
/// once in every class its separation belongs to. Pairs are enumerated as
/// described in the module documentation.
pub(crate) fn semivariogram<L: LagScheme + ?Sized>(
    points: &[Point3],
    values: &[f64],
    lags: &L,
) -> Vec<Lag> {
    let plan = PairPlan::new(points, values, lags);
    plan.finish((0..plan.block_count()).map(|b| plan.block_sums(b)))
}

/// The plain double loop over every pair: the reference the pruned
/// enumeration must match (P4).
#[cfg(test)]
pub(crate) fn semivariogram_all_pairs<L: LagScheme + ?Sized>(
    points: &[Point3],
    values: &[f64],
    lags: &L,
) -> Vec<Lag> {
    let data: Vec<(Point3, f64)> = points.iter().copied().zip(values.iter().copied()).collect();
    let mut acc = vec![LagSums::default(); lags.len()];
    for (i, &(p, v)) in data.iter().enumerate() {
        for &(q, w) in data.iter().skip(i + 1) {
            let h = p.distance(q);
            let d = v - w;
            for class in lags.classes(h) {
                if let Some(a) = acc.get_mut(class) {
                    a.gamma = (0.5 * d).mul_add(d, a.gamma);
                    a.distance += h;
                    a.sum += v + w;
                    a.pairs += 1;
                }
            }
        }
    }
    lags_from(lags, acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn s(x: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, 0.0, 0.0), v).unwrap()
    }

    #[test]
    fn two_points_one_lag() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(1.0, 1.0)],
            &LagClasses::new(1.0, 1).unwrap(),
        );
        assert_eq!(lags.len(), 1);
        assert_abs_diff_eq!(lags[0].center, 0.5);
        assert_eq!(lags[0].pairs, 1);
        assert_abs_diff_eq!(lags[0].gamma.unwrap(), 0.5);
        assert_abs_diff_eq!(lags[0].mean_distance.unwrap(), 1.0);
    }

    #[test]
    fn empty_lags_are_none_not_zero() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(3.5, 1.0)],
            &LagClasses::new(4.0, 4).unwrap(),
        );
        assert_eq!(
            lags.iter().map(|l| l.pairs).collect::<Vec<_>>(),
            vec![0, 0, 0, 1]
        );
        assert_eq!(lags[0].gamma, None);
        assert!(lags[3].gamma.is_some(), "the populated lag has a value");
    }

    #[test]
    fn pairs_beyond_max_distance_are_excluded() {
        let lags = experimental_variogram(
            &[s(0.0, 0.0), s(5.0, 1.0)],
            &LagClasses::new(4.0, 2).unwrap(),
        );
        assert!(lags.iter().all(|l| l.pairs == 0), "{lags:?}");
    }

    #[test]
    fn linear_trend_gives_quadratic_gamma() {
        // z = x on a unit grid ⇒ γ(h) = h²/2 exactly.
        let samples: Vec<Sample> = (0..10).map(|i| s(f64::from(i), f64::from(i))).collect();
        let lags = experimental_variogram(&samples, &LagClasses::new(3.5, 7).unwrap());
        for lag in lags.iter().filter(|l| l.pairs > 0) {
            let h = lag.mean_distance.unwrap();
            assert_abs_diff_eq!(lag.gamma.unwrap(), h * h / 2.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn invalid_classes_rejected() {
        assert!(LagClasses::new(1.0, 0).is_err(), "zero lags");
        assert!(LagClasses::new(0.0, 3).is_err(), "zero distance");
        assert!(LagClasses::new(f64::NAN, 3).is_err(), "NaN distance");
    }

    fn centred(count: usize, spacing: f64, tolerance: f64) -> CentredLags {
        CentredLags::new(count, spacing, tolerance).unwrap()
    }

    #[test]
    fn centred_lags_follow_gamv() {
        // nlag = 3, xlag = 10, xltol = 5: classes 0 (zero), then centres 0,
        // 10, 20, 30; GAMV.FOR's bounds are inclusive at both ends.
        let lags = centred(3, 10.0, 5.0);
        assert_eq!(lags.len(), 5);
        assert_eq!(lags.classes(0.0), 0..1, "coincident pairs");
        assert_eq!(lags.classes(3.0), 1..2, "centre 0");
        assert_eq!(lags.classes(5.0), 1..3, "on a shared bound: both classes");
        assert_eq!(lags.classes(12.0), 2..3);
        assert_eq!(lags.classes(35.0), 4..5, "up to (nlag + 0.5) * xlag");
        assert_eq!(lags.classes(35.5), 0..0, "beyond the limit");
        assert_eq!(lags.classes(f64::NAN), 0..0, "NaN is never classified");
        assert_eq!(
            (0..5).map(|c| lags.center(c)).collect::<Vec<_>>(),
            vec![0.0, 0.0, 10.0, 20.0, 30.0]
        );
    }

    #[test]
    fn centred_lags_overlap_with_wide_tolerance() {
        // xltol > xlag / 2: 10 lies within 8 of centre 10 only; 14 lies
        // within 8 of both 10 and 20.
        let lags = centred(4, 10.0, 8.0);
        assert_eq!(lags.classes(10.0), 2..3);
        assert_eq!(lags.classes(14.0), 2..4, "counted in two classes");
    }

    #[test]
    fn centred_variogram_counts_each_distinct_pair_once_per_class() {
        let data = [s(0.0, 1.0), s(0.0, 3.0), s(10.0, 2.0), s(14.0, 6.0)];
        let lags = experimental_variogram(&data, &centred(2, 10.0, 8.0));
        // Zero class: one coincident pair (1, 3), not the four self-pairs.
        assert_eq!(lags[0].pairs, 1);
        assert_abs_diff_eq!(lags[0].gamma.unwrap(), 2.0);
        assert_abs_diff_eq!(lags[0].pair_mean.unwrap(), 2.0);
        // Separation 4 (10-14) is within 8 of centre 0 and of 10.
        assert_eq!(lags[1].pairs, 1);
        // Class 2 is centred on 10 with tolerance 8, i.e. [2, 18].
        assert_eq!(lags[2].pairs, 5, "pairs at 10, 10, 14, 14 and 4");
        // Class 3 is centred on 20, i.e. [12, 28].
        assert_eq!(lags[3].pairs, 2, "pairs at 14 and 14");
    }

    #[test]
    fn centred_lags_reject_invalid() {
        assert!(CentredLags::new(0, 1.0, 0.5).is_err(), "zero lags");
        assert!(CentredLags::new(3, 0.0, 0.5).is_err(), "zero spacing");
        assert!(
            CentredLags::new(3, 1.0, 0.0).is_err(),
            "GSLIB's default tolerance must be explicit"
        );
    }

    fn at(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    /// Pair counts must match exactly; sums only to rounding. Recursive
    /// summation of `n` terms has error at most about `(n − 1)·ε·Σ|xᵢ|`, so
    /// summing in another order can change a sum of like-signed terms by up
    /// to about `2·n·ε` relative. The tolerance is that, at least 1e-12.
    fn assert_same_pairs(pruned: &[Lag], full: &[Lag]) {
        assert_eq!(pruned.len(), full.len());
        for (k, (a, b)) in pruned.iter().zip(full).enumerate() {
            assert_eq!(a.pairs, b.pairs, "lag {k}");
            for (x, y) in [
                (a.gamma, b.gamma),
                (a.mean_distance, b.mean_distance),
                (a.pair_mean, b.pair_mean),
            ] {
                match (x, y) {
                    (Some(x), Some(y)) => {
                        let tolerance = (2.0 * to_f64(b.pairs) * f64::EPSILON).max(1e-12);
                        assert!(
                            (x - y).abs() <= tolerance * y.abs().max(1.0),
                            "lag {k}: {x} vs {y} (tolerance {tolerance:e})"
                        );
                    }
                    (None, None) => {}
                    other => panic!("lag {k}: {other:?}"),
                }
            }
        }
    }

    #[test]
    fn pairs_exactly_at_the_cut_off_are_kept() {
        // P4: GSLIB's cut-off (nlag + 0.5) xlag = 2.5, and a pair exactly
        // there (lag 2's tolerance reaches it). Also one ulp beyond.
        let lags = centred(2, 1.0, 0.5);
        let beyond = f64::from_bits(2.5_f64.to_bits() + 1);
        let points = [at(0.0, 0.0, 0.0), at(2.5, 0.0, 0.0), at(0.0, beyond, 0.0)];
        let values = [0.0, 1.0, 5.0];
        let pruned = semivariogram(&points, &values, &lags);
        assert_same_pairs(&pruned, &semivariogram_all_pairs(&points, &values, &lags));
        assert_eq!(
            pruned.iter().map(|l| l.pairs).sum::<usize>(),
            1,
            "only the pair at 2.5"
        );
        // The same for equal-width classes at their maximum.
        let classes = LagClasses::new(2.5, 5).unwrap();
        let by_class = semivariogram(&points, &values, &classes);
        assert_same_pairs(
            &by_class,
            &semivariogram_all_pairs(&points, &values, &classes),
        );
        assert_eq!(by_class.last().map(|l| l.pairs), Some(1));
    }

    #[test]
    fn pruning_spans_many_blocks() {
        // More rows than one block, on a lattice (many ties and exact
        // separations), with overlapping centred lags.
        let points: Vec<Point3> = (0..900_u32)
            .map(|i| at(f64::from(i % 30), f64::from(i.div_euclid(30)), 0.0))
            .collect();
        let values: Vec<f64> = (0..900_u32).map(|i| f64::from(i % 7) * 0.3).collect();
        assert!(points.len() > 3 * ROWS_PER_BLOCK);
        let lags = centred(6, 1.5, 1.0);
        assert_same_pairs(
            &semivariogram(&points, &values, &lags),
            &semivariogram_all_pairs(&points, &values, &lags),
        );
    }

    #[test]
    fn no_data_and_one_datum_give_empty_lags() {
        let lags = centred(3, 1.0, 0.5);
        for n in [0, 1] {
            let points = vec![at(0.0, 0.0, 0.0); n];
            let values = vec![1.0; n];
            let got = semivariogram(&points, &values, &lags);
            assert!(got.iter().all(|l| l.pairs == 0), "{n} data");
            assert_eq!(got.len(), lags.len());
        }
    }

    #[test]
    #[ignore = "timing evidence for pruning; run with --release --ignored --nocapture"]
    fn pruning_timing() {
        // 20 000 scattered 2-D data, GSLIB lags reaching 5% of the domain.
        let points: Vec<Point3> = (0..20_000_u32)
            .map(|i| {
                let f = f64::from(i);
                at(
                    (f * 0.618_034).fract() * 1000.0,
                    (f * 0.414_214).fract() * 1000.0,
                    0.0,
                )
            })
            .collect();
        let values: Vec<f64> = points.iter().map(|p| (p.x * 0.01).sin()).collect();
        let lags = centred(10, 5.0, 2.5);
        let started = std::time::Instant::now();
        let pruned = semivariogram(&points, &values, &lags);
        let fast = started.elapsed();
        let restarted = std::time::Instant::now();
        let full = semivariogram_all_pairs(&points, &values, &lags);
        let slow = restarted.elapsed();
        assert_same_pairs(&pruned, &full);
        println!("all pairs {slow:?}, pruned {fast:?}");
    }

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn pruning_counts_exactly_the_same_pairs(
            raw in prop::collection::vec(
                prop_oneof![
                    (-30.0_f64..30.0, -30.0_f64..30.0, -3.0_f64..3.0),
                    (-6_i8..=6, -6_i8..=6, -1_i8..=1).prop_map(|(x, y, z)| (f64::from(x), f64::from(y), f64::from(z))),
                ],
                0..120,
            ),
            count in 1_usize..8,
            spacing in 0.5_f64..6.0,
            tolerance_ratio in 0.1_f64..1.5,
            max in 1.0_f64..40.0,
        ) {
            // P4: both lag schemes, the zero-distance class (duplicates from
            // the lattice), overlapping lags (tolerance above half the
            // spacing), and every pair compared with the full double loop.
            let points: Vec<Point3> = raw.iter().map(|&(x, y, z)| at(x, y, z)).collect();
            let values: Vec<f64> = raw.iter().map(|&(x, y, z)| x.mul_add(0.1, y.sin()) + z).collect();
            let centred = centred(count, spacing, spacing * tolerance_ratio);
            assert_same_pairs(
                &semivariogram(&points, &values, &centred),
                &semivariogram_all_pairs(&points, &values, &centred),
            );
            let classes = LagClasses::new(max, count).unwrap();
            assert_same_pairs(
                &semivariogram(&points, &values, &classes),
                &semivariogram_all_pairs(&points, &values, &classes),
            );
        }
    }
}
