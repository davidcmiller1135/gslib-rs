//! Search neighbourhoods: which conditioning data to use at each location
//! (`docs/spec/search.md`).
//!
//! A search has two parts:
//! - a [`Neighbourhood`] finds the **candidates**, every datum inside the
//!   [`SearchEllipsoid`], in a fixed order. Implementations differ only in
//!   speed: [`BruteForce`] is the definition, and any accelerator must return
//!   exactly its result (D-025).
//! - [`DataSearch`] applies GSLIB's rules to the candidates (`KT3D.FOR`,
//!   `SRCHSUPR.FOR`; F-012): the octant limit, then the `ndmax` cut, then the
//!   `ndmin` check.
//!
//! # Formulation
//!
//! - **Membership (S1):** a datum at `x` is a candidate for target `t` iff
//!   `q(x − t) ≤ radius²`, where `q` is the squared anisotropic distance
//!   ([`Anisotropy::squared_scaled_distance`], GSLIB's `sqdist`). GSLIB
//!   compares the same squared quantities, so no square root can move a
//!   datum across the boundary.
//! - **Order (S2, D-024):** candidates are sorted by `(q, row)`, a total
//!   order. **Divergence from GSLIB** (ROADMAP §2): GSLIB orders equidistant
//!   data by an unstable quicksort over its super-block order, so which tied
//!   datum survives a cut can differ.
//! - **Octants (S3):** with `per_octant = n`, candidates are scanned in
//!   order and kept while their octant holds fewer than `n`. GSLIB stops
//!   the scan once `8n` are kept, which changes nothing (every octant is
//!   then full) but saves time. The octant of a datum follows the signs of `datum −
//!   target` with GSLIB's boundary rules (see [`octant`]).
//! - **Truncation (S4):** the first `max_data` of what remains.
//! - **Too few data (S5):** fewer than `min_data` is
//!   [`Error::TooFewData`]. GSLIB's behaviour (`kt3d` writes `UNEST`) is the
//!   opt-in [`OnTooFewData::Unestimated`].

use std::num::NonZeroUsize;

use crate::error::{Error, Result};
use crate::geometry::{Angles, Anisotropy, Point3, Ranges};
use crate::sample::Sample;

/// The search ellipsoid: GSLIB's `radius, radius1, radius2` and
/// `sang1, sang2, sang3`, independent of any variogram model's anisotropy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SearchEllipsoid {
    radius: f64,
    anisotropy: Anisotropy,
}

impl SearchEllipsoid {
    /// An ellipsoid with radii `radii` (`major` is GSLIB's `radius`, `minor`
    /// its `radius1`, `vertical` its `radius2`) oriented by `angles`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if an angle is not finite or a radius is not
    /// finite and positive.
    pub fn new(angles: Angles, radii: Ranges) -> Result<Self> {
        Ok(Self {
            radius: radii.major,
            anisotropy: Anisotropy::new(angles, radii)?,
        })
    }

    /// A sphere. `f64::INFINITY` means no distance limit.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `radius` is NaN or not positive.
    pub fn isotropic(radius: f64) -> Result<Self> {
        if radius.is_nan() || radius <= 0.0 {
            return Err(Error::invalid(
                "search radius",
                format!("must be positive, got {radius}"),
            ));
        }
        Ok(Self {
            radius,
            anisotropy: Anisotropy::isotropic(),
        })
    }

    /// No distance limit.
    #[must_use]
    pub const fn unlimited() -> Self {
        Self {
            radius: f64::INFINITY,
            anisotropy: Anisotropy::isotropic(),
        }
    }

    /// The radius along the major axis.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// The squared anisotropic distance `q` between two points, in units
    /// where the ellipsoid is the sphere of radius [`SearchEllipsoid::radius`].
    #[must_use]
    pub fn squared_distance(&self, a: Point3, b: Point3) -> f64 {
        self.anisotropy.squared_scaled_distance(a, b)
    }

    /// Whether a squared distance lies inside the ellipsoid (`q ≤ radius²`).
    #[must_use]
    pub fn admits(&self, squared_distance: f64) -> bool {
        squared_distance <= self.radius * self.radius
    }
}

/// A datum inside the search ellipsoid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Candidate {
    /// Its 0-based index in the data the neighbourhood was built from.
    pub row: usize,
    /// Its squared anisotropic distance `q` from the target.
    pub squared_distance: f64,
    /// `datum − target`, for octants.
    pub offset: [f64; 3],
}

/// Finds the candidates for a target: every datum inside the ellipsoid,
/// sorted by `(squared_distance, row)` (S1, S2).
///
/// Implementations are interchangeable: for the same data and arguments
/// they must return the same sequence as [`BruteForce`] (D-025).
pub trait Neighbourhood {
    /// The candidates for `target` under `ellipsoid`.
    fn candidates(&self, ellipsoid: &SearchEllipsoid, target: Point3) -> Vec<Candidate>;
}

/// Every datum is examined: the definition of a search, O(N) per target.
#[derive(Debug, Clone, Copy)]
pub struct BruteForce<'a> {
    samples: &'a [Sample],
}

impl<'a> BruteForce<'a> {
    /// A neighbourhood over `samples`; candidate rows index this slice.
    #[must_use]
    pub const fn new(samples: &'a [Sample]) -> Self {
        Self { samples }
    }
}

impl Neighbourhood for BruteForce<'_> {
    fn candidates(&self, ellipsoid: &SearchEllipsoid, target: Point3) -> Vec<Candidate> {
        let mut found: Vec<Candidate> = self
            .samples
            .iter()
            .enumerate()
            .filter_map(|(row, s)| {
                let q = ellipsoid.squared_distance(s.location(), target);
                ellipsoid.admits(q).then(|| Candidate {
                    row,
                    squared_distance: q,
                    offset: s.location().delta(target),
                })
            })
            .collect();
        sort_candidates(&mut found);
        found
    }
}

/// Sort by `(squared_distance, row)` (D-024).
fn sort_candidates(candidates: &mut [Candidate]) {
    candidates.sort_unstable_by(|a, b| {
        a.squared_distance
            .total_cmp(&b.squared_distance)
            .then(a.row.cmp(&b.row))
    });
}

/// GSLIB's octant (1–8) of an offset `datum − target` (`SRCHSUPR.FOR`).
///
/// 1–4 for `dz ≥ 0` and 5–8 for `dz < 0`, then, by the signs of `dx, dy`:
/// 1 if `dx ≤ 0, dy > 0`; 2 if `dx > 0, dy ≥ 0`; 3 if `dx < 0, dy ≤ 0`; 4
/// otherwise (including `dx = dy = 0`).
#[must_use]
pub fn octant([dx, dy, dz]: [f64; 3]) -> u8 {
    let quadrant = if dx <= 0.0 && dy > 0.0 {
        1
    } else if dx > 0.0 && dy >= 0.0 {
        2
    } else if dx < 0.0 && dy <= 0.0 {
        3
    } else {
        4
    };
    if dz < 0.0 { quadrant + 4 } else { quadrant }
}

/// The rows of `candidates` (in order) kept by GSLIB's octant rule: at most
/// `n` per octant (S3). Stopping at `8n` is GSLIB's shortcut: every octant
/// is full by then.
fn octant_limited(candidates: &[Candidate], n: NonZeroUsize) -> Vec<usize> {
    let limit = n.get();
    let mut counts = [0_usize; 8];
    let mut rows = Vec::new();
    for c in candidates {
        if rows.len() == limit.saturating_mul(8) {
            break;
        }
        if let Some(count) = counts.get_mut(usize::from(octant(c.offset) - 1)) {
            *count += 1;
            if *count <= limit {
                rows.push(c.row);
            }
        }
    }
    rows
}

/// What a search does when it finds fewer than `min_data` data (S5).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OnTooFewData {
    /// Fail with [`Error::TooFewData`].
    #[default]
    Error,
    /// Report the location as unestimated (GSLIB `kt3d`'s `UNEST`).
    Unestimated,
}

impl OnTooFewData {
    /// Apply the policy to a search result: under
    /// [`OnTooFewData::Unestimated`], [`Error::TooFewData`] becomes `None`.
    ///
    /// # Errors
    /// Any other error, and [`Error::TooFewData`] under
    /// [`OnTooFewData::Error`].
    pub fn resolve<T>(self, result: Result<T>) -> Result<Option<T>> {
        match (self, result) {
            (_, Ok(v)) => Ok(Some(v)),
            (Self::Unestimated, Err(Error::TooFewData { .. })) => Ok(None),
            (_, Err(e)) => Err(e),
        }
    }
}

/// GSLIB's data search: an ellipsoid, `ndmin`, `ndmax`, and optionally a
/// limit per octant (`noct`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DataSearch {
    ellipsoid: SearchEllipsoid,
    min_data: usize,
    max_data: usize,
    per_octant: Option<NonZeroUsize>,
}

impl DataSearch {
    /// Keep at most `max_data` candidates and require at least `min_data`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `max_data` is zero or `min_data` exceeds
    /// it (every search would fail).
    pub fn new(ellipsoid: SearchEllipsoid, min_data: usize, max_data: usize) -> Result<Self> {
        if max_data == 0 {
            return Err(Error::invalid("search", "max_data must be at least 1"));
        }
        if min_data > max_data {
            return Err(Error::invalid(
                "search",
                format!("min_data {min_data} exceeds max_data {max_data}"),
            ));
        }
        Ok(Self {
            ellipsoid,
            min_data,
            max_data,
            per_octant: None,
        })
    }

    /// Every datum, however many and however far (cost grows quickly).
    #[must_use]
    pub const fn unlimited() -> Self {
        Self {
            ellipsoid: SearchEllipsoid::unlimited(),
            min_data: 0,
            max_data: usize::MAX,
            per_octant: None,
        }
    }

    /// Keep at most `per_octant` data in each octant (GSLIB's `noct`).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `min_data` exceeds `8 · per_octant`
    /// (every search would fail).
    pub fn with_octants(self, per_octant: NonZeroUsize) -> Result<Self> {
        if self.min_data > per_octant.get().saturating_mul(8) {
            return Err(Error::invalid(
                "search",
                format!(
                    "min_data {} exceeds 8 octants x {per_octant}",
                    self.min_data
                ),
            ));
        }
        Ok(Self {
            per_octant: Some(per_octant),
            ..self
        })
    }

    /// The search ellipsoid.
    #[must_use]
    pub const fn ellipsoid(&self) -> &SearchEllipsoid {
        &self.ellipsoid
    }

    /// The rows of the data to use at `target`, nearest first (S1–S5).
    ///
    /// # Errors
    /// [`Error::TooFewData`] if fewer than `min_data` remain; see
    /// [`OnTooFewData::resolve`] for GSLIB's alternative.
    pub fn select(&self, neighbourhood: &impl Neighbourhood, target: Point3) -> Result<Vec<usize>> {
        let candidates = neighbourhood.candidates(&self.ellipsoid, target);
        let kept: Vec<usize> = self.per_octant.map_or_else(
            || {
                candidates
                    .iter()
                    .map(|c| c.row)
                    .take(self.max_data)
                    .collect()
            },
            |n| {
                let mut rows = octant_limited(&candidates, n);
                rows.truncate(self.max_data);
                rows
            },
        );
        if kept.len() < self.min_data {
            return Err(Error::TooFewData {
                target,
                found: kept.len(),
                required: self.min_data,
            });
        }
        Ok(kept)
    }

    /// The samples to use at `target`, nearest first: [`DataSearch::select`]
    /// over a [`BruteForce`] neighbourhood of `samples`.
    ///
    /// # Errors
    /// As [`DataSearch::select`].
    pub fn select_samples(&self, samples: &[Sample], target: Point3) -> Result<Vec<Sample>> {
        Ok(self
            .select(&BruteForce::new(samples), target)?
            .into_iter()
            .filter_map(|row| samples.get(row).copied())
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn at(x: f64, y: f64, z: f64) -> Sample {
        Sample::new(Point3::new(x, y, z), x).unwrap()
    }

    fn rows(search: &DataSearch, samples: &[Sample], target: Point3) -> Vec<usize> {
        search.select(&BruteForce::new(samples), target).unwrap()
    }

    fn sphere(radius: f64, max_data: usize) -> DataSearch {
        DataSearch::new(SearchEllipsoid::isotropic(radius).unwrap(), 0, max_data).unwrap()
    }

    #[test]
    fn membership_is_inclusive_along_each_rotated_axis() {
        // S1: azimuth 30, dip -20 (downward), plunge 0; radii 10, 4, 2. The
        // major axis points along azimuth 30 from north, dipping 20 degrees
        // down; points at exactly the radius along each axis are inside.
        let (az, dip) = (30.0_f64.to_radians(), (-20.0_f64).to_radians());
        let major = [az.sin() * dip.cos(), az.cos() * dip.cos(), dip.sin()];
        let minor = [az.cos(), -az.sin(), 0.0];
        // The third axis is perpendicular to both (its sign is irrelevant).
        let vertical = [
            major[2].mul_add(-minor[1], major[1] * minor[2]),
            major[0].mul_add(-minor[2], major[2] * minor[0]),
            major[1].mul_add(-minor[0], major[0] * minor[1]),
        ];
        let ellipsoid = SearchEllipsoid::new(
            Angles {
                azimuth: 30.0,
                dip: -20.0,
                plunge: 0.0,
            },
            Ranges {
                major: 10.0,
                minor: 4.0,
                vertical: 2.0,
            },
        )
        .unwrap();
        let o = Point3::default();
        for (axis, r) in [(major, 10.0), (minor, 4.0), (vertical, 2.0)] {
            let p = |f: f64| Point3::new(axis[0] * r * f, axis[1] * r * f, axis[2] * r * f);
            let q = |f: f64| ellipsoid.squared_distance(p(f), o);
            assert!(ellipsoid.admits(q(0.999)), "just inside, axis {axis:?}");
            assert!((q(1.0) - 100.0).abs() < 1e-9, "on the surface: {}", q(1.0));
            assert!(!ellipsoid.admits(q(1.001)), "just outside, axis {axis:?}");
        }
        // Negative dip is downward: the major axis tip is below the origin.
        assert!(major[2] < 0.0, "dip -20 points down");
    }

    #[test]
    fn membership_boundary_is_exact() {
        // A datum at exactly the radius is kept (q <= radius²), one ulp
        // farther is not.
        let samples = [
            at(3.0, 4.0, 0.0),
            at(f64::from_bits(5.0_f64.to_bits() + 1), 0.0, 0.0),
        ];
        assert_eq!(rows(&sphere(5.0, 10), &samples, Point3::default()), [0]);
    }

    #[test]
    fn candidates_are_sorted_by_distance_then_row() {
        // S2: four data at distance 1 and two at 2, in a scrambled order.
        let samples = [
            at(0.0, 2.0, 0.0),
            at(0.0, -1.0, 0.0),
            at(1.0, 0.0, 0.0),
            at(-2.0, 0.0, 0.0),
            at(0.0, 1.0, 0.0),
            at(-1.0, 0.0, 0.0),
        ];
        assert_eq!(
            rows(&sphere(10.0, 10), &samples, Point3::default()),
            [1, 2, 4, 5, 0, 3]
        );
    }

    #[test]
    fn octants_follow_gslib_boundary_rules() {
        // F-012's rules, including points on the boundaries.
        let cases = [
            ([0.0, 1.0, 0.0], 1),
            ([-1.0, 1.0, 0.0], 1),
            ([1.0, 0.0, 0.0], 2),
            ([1.0, 1.0, 0.0], 2),
            ([-1.0, 0.0, 0.0], 3),
            ([-1.0, -1.0, 0.0], 3),
            ([0.0, -1.0, 0.0], 4),
            ([1.0, -1.0, 0.0], 4),
            ([0.0, 0.0, 0.0], 4),
            ([0.0, 0.0, 1.0], 4),
            ([0.0, 1.0, -1.0], 5),
            ([1.0, 0.0, -1.0], 6),
            ([-1.0, 0.0, -1.0], 7),
            ([0.0, 0.0, -1.0], 8),
        ];
        for (offset, expected) in cases {
            assert_eq!(octant(offset), expected, "{offset:?}");
        }
    }

    #[test]
    fn octant_limit_applies_before_the_cut() {
        // S3/S4: three data in octant 2 (nearest first), one in octant 3.
        let samples = [
            at(1.0, 0.0, 0.0),
            at(2.0, 0.0, 0.0),
            at(3.0, 0.0, 0.0),
            at(-4.0, 0.0, 0.0),
        ];
        let o = Point3::default();
        let two = |max| {
            sphere(10.0, max)
                .with_octants(NonZeroUsize::new(2).unwrap())
                .unwrap()
        };
        assert_eq!(
            rows(&two(10), &samples, o),
            [0, 1, 3],
            "third in octant 2 skipped"
        );
        assert_eq!(rows(&two(2), &samples, o), [0, 1], "then the cut");
        assert_eq!(
            rows(&two(3), &samples, o),
            [0, 1, 3],
            "the cut counts only data the octants kept"
        );
        // Without octants, the cut alone.
        assert_eq!(rows(&sphere(10.0, 3), &samples, o), [0, 1, 2]);
    }

    #[test]
    fn octant_limit_caps_the_total_at_eight_per_octant() {
        // Three data per octant, two allowed: 16 kept of 24, the nearest
        // two in each octant. (GSLIB's scan also stops at 8n; that is only
        // a shortcut, since every later candidate would be rejected.)
        let mut samples = Vec::new();
        for r in [1.0, 2.0, 3.0] {
            for (x, y) in [(-1.0, 1.0), (1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)] {
                for z in [1.0, -1.0] {
                    samples.push(at(x * r, y * r, z * r));
                }
            }
        }
        let search = sphere(100.0, 100)
            .with_octants(NonZeroUsize::new(2).unwrap())
            .unwrap();
        let kept = rows(&search, &samples, Point3::default());
        assert_eq!(kept.len(), 16);
        assert!(
            kept.iter().all(|&r| r < 16),
            "the farthest ring is never kept"
        );
    }

    #[test]
    fn truncation_keeps_the_nearest() {
        // S4: max_data smaller, equal and larger than the candidate count.
        let samples = [
            at(0.0, 0.0, 0.0),
            at(10.0, 0.0, 0.0),
            at(4.0, 0.0, 0.0),
            at(6.0, 0.0, 0.0),
            at(-3.0, 0.0, 0.0),
        ];
        let t = Point3::new(5.0, 0.0, 0.0);
        assert_eq!(rows(&sphere(f64::INFINITY, 2), &samples, t), [2, 3]);
        assert_eq!(
            rows(&sphere(f64::INFINITY, 5), &samples, t),
            [2, 3, 0, 1, 4]
        );
        assert_eq!(rows(&sphere(f64::INFINITY, 9), &samples, t).len(), 5);
        assert_eq!(rows(&sphere(2.0, 9), &samples, t), [2, 3], "radius");
    }

    #[test]
    fn too_few_data_is_an_error_unless_unestimated() {
        // S5.
        let samples = [at(0.0, 0.0, 0.0), at(10.0, 0.0, 0.0)];
        let search = DataSearch::new(SearchEllipsoid::isotropic(5.0).unwrap(), 2, 4).unwrap();
        let t = Point3::new(1.0, 0.0, 0.0);
        let result = search.select(&BruteForce::new(&samples), t);
        assert!(
            matches!(
                result,
                Err(Error::TooFewData {
                    found: 1,
                    required: 2,
                    ..
                })
            ),
            "{result:?}"
        );
        let resolved =
            OnTooFewData::Unestimated.resolve(search.select(&BruteForce::new(&samples), t));
        assert_eq!(resolved.unwrap(), None);
        assert!(
            OnTooFewData::Error
                .resolve(search.select(&BruteForce::new(&samples), t))
                .is_err(),
            "default policy"
        );
        // Exactly min_data is enough.
        let enough = DataSearch::new(SearchEllipsoid::isotropic(20.0).unwrap(), 2, 4).unwrap();
        assert_eq!(
            enough.select(&BruteForce::new(&samples), t).unwrap(),
            [0, 1]
        );
        // min_data = 0 never fails, even with nothing in range.
        let none = sphere(0.5, 3).select(&BruteForce::new(&samples), Point3::new(5.0, 0.0, 0.0));
        assert_eq!(none.unwrap(), Vec::<usize>::new());
        // Other errors pass through either policy.
        let other: Result<()> = Err(Error::NoCovariance);
        assert!(
            OnTooFewData::Unestimated.resolve(other).is_err(),
            "not a search error"
        );
    }

    #[test]
    fn rejects_invalid() {
        assert!(SearchEllipsoid::isotropic(0.0).is_err(), "zero radius");
        assert!(SearchEllipsoid::isotropic(f64::NAN).is_err(), "NaN radius");
        assert!(
            SearchEllipsoid::isotropic(f64::INFINITY).is_ok(),
            "unlimited"
        );
        let e = SearchEllipsoid::unlimited();
        assert!(DataSearch::new(e, 0, 0).is_err(), "max_data 0");
        assert!(DataSearch::new(e, 3, 2).is_err(), "min above max");
        let one = NonZeroUsize::new(1).unwrap();
        assert!(
            DataSearch::new(e, 9, 20)
                .unwrap()
                .with_octants(one)
                .is_err(),
            "min above 8 per-octant slots"
        );
        assert!(DataSearch::new(e, 8, 20).unwrap().with_octants(one).is_ok());
    }

    #[test]
    fn select_samples_returns_the_selected_data_in_order() {
        let samples = [at(5.0, 0.0, 0.0), at(1.0, 0.0, 0.0), at(3.0, 0.0, 0.0)];
        let got = sphere(10.0, 2)
            .select_samples(&samples, Point3::default())
            .unwrap();
        assert_eq!(got, vec![samples[1], samples[2]]);
        assert_eq!(
            DataSearch::unlimited()
                .select_samples(&samples, Point3::default())
                .unwrap()
                .len(),
            3
        );
    }

    #[test]
    fn searching_in_parallel_matches_sequential() {
        // S7: no hidden state; results do not depend on thread count.
        let samples: Vec<Sample> = (0..200_u32)
            .map(|i| {
                let f = f64::from(i);
                at((f * 7.3) % 50.0, (f * 3.1) % 40.0, (f * 1.7) % 5.0)
            })
            .collect();
        let search = sphere(15.0, 12)
            .with_octants(NonZeroUsize::new(3).unwrap())
            .unwrap();
        let targets: Vec<Point3> = (0..40_u32)
            .map(|i| Point3::new(f64::from(i), f64::from(i) * 0.8, 2.0))
            .collect();
        let sequential: Vec<Vec<usize>> = targets
            .iter()
            .map(|&t| rows(&search, &samples, t))
            .collect();
        let parallel: Vec<Vec<usize>> = std::thread::scope(|scope| {
            #[expect(
                clippy::needless_collect,
                reason = "every thread must start before any is joined"
            )]
            let handles: Vec<_> = targets
                .chunks(7)
                .map(|chunk| {
                    scope.spawn(|| {
                        chunk
                            .iter()
                            .map(|&t| rows(&search, &samples, t))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect()
        });
        assert_eq!(parallel, sequential);
    }

    proptest! {
        #[test]
        fn order_does_not_depend_on_input_permutation(
            points in prop::collection::vec((-3_i8..=3, -3_i8..=3, -1_i8..=1), 1..30),
            shift in 0_usize..30,
        ) {
            // S2: integer coordinates give many exact ties; rotating the
            // input permutes rows, and the selected locations, in order,
            // must be the same up to that relabelling only where distances
            // differ. Here: sorted (q, location) sequences agree.
            let samples: Vec<Sample> = points
                .iter()
                .map(|&(x, y, z)| at(f64::from(x), f64::from(y), f64::from(z)))
                .collect();
            let mut rotated = samples.clone();
            rotated.rotate_left(shift % samples.len());
            let e = SearchEllipsoid::isotropic(2.5).unwrap();
            let t = Point3::new(0.25, -0.5, 0.0);
            let key = |s: &[Sample]| -> Vec<(u64, [u64; 3])> {
                BruteForce::new(s)
                    .candidates(&e, t)
                    .iter()
                    .map(|c| {
                        let p = s[c.row].location().delta(Point3::default());
                        (c.squared_distance.to_bits(), p.map(f64::to_bits))
                    })
                    .collect()
            };
            let mut a = key(&samples);
            let mut b = key(&rotated);
            // Within a tie the rows differ, so compare as multisets per q.
            a.sort_unstable();
            b.sort_unstable();
            prop_assert_eq!(a, b);
            // And each sequence is non-decreasing in q, with rows ascending
            // within ties.
            let c = BruteForce::new(&rotated).candidates(&e, t);
            for w in c.windows(2) {
                prop_assert!(
                    w[0].squared_distance < w[1].squared_distance
                        || (w[0].squared_distance.to_bits() == w[1].squared_distance.to_bits()
                            && w[0].row < w[1].row)
                );
            }
        }
    }
}
