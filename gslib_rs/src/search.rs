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
use crate::num::{floor_index, to_f64};
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

    /// The smallest factor by which the search metric stretches a direction:
    /// `q(h) ≥ smallest_scale()² · |h|²` ([`Anisotropy::smallest_scale`]).
    #[must_use]
    pub fn smallest_scale(&self) -> f64 {
        self.anisotropy.smallest_scale()
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
    /// The data the neighbourhood was built from; candidate rows index it.
    fn samples(&self) -> &[Sample];

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
    fn samples(&self) -> &[Sample] {
        self.samples
    }

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

/// Data binned into a regular grid of blocks over their bounding box (the
/// idea of GSLIB's `setsupr`/`srchsupr`), so a search examines only the
/// blocks that can hold candidates.
///
/// **Exact** (D-025, S6): it returns exactly the [`BruteForce`] sequence. A
/// block is searched unless a *lower bound* on the squared distance from
/// the target to any of its data exceeds `radius²`. The bound is `s² · d²`,
/// where `s` is [`SearchEllipsoid::smallest_scale`] and `d` the Euclidean
/// distance from the target to the tight bounding box of the block's own
/// data, reduced by a relative margin of 1e-9 so rounding cannot exclude a
/// datum on the boundary. **Divergence from GSLIB** (`picksup`, F-012):
/// GSLIB tests 27 lattice points, which may exclude blocks holding data
/// inside a rotated anisotropic ellipsoid (hypothesis, ROADMAP §6).
#[derive(Debug, Clone)]
pub struct SuperBlock {
    origin: [f64; 3],
    size: [f64; 3],
    counts: [usize; 3],
    /// The data, by row.
    samples: Vec<Sample>,
    /// Rows sorted by block; block `b` holds `rows[starts[b]..starts[b + 1]]`.
    rows: Vec<usize>,
    starts: Vec<usize>,
    /// Tight bounding box `(min, max)` of each block's data.
    boxes: Vec<([f64; 3], [f64; 3])>,
}

/// Relative margin on the block lower bound (see [`SuperBlock`]).
const BOUND_MARGIN: f64 = 1e-9;
/// About this many data per block when [`SuperBlock::new`] chooses counts.
const DATA_PER_BLOCK: usize = 4;
/// At most this many blocks along an axis.
const MAX_BLOCKS_PER_AXIS: usize = 1024;

impl SuperBlock {
    /// Blocks chosen for `samples`: about four data per block, shaped to the
    /// data's extent (one block along an axis where all data share a
    /// coordinate).
    #[must_use]
    pub fn new(samples: &[Sample]) -> Self {
        let (lo, hi) = extent(samples);
        let extents = zip3(hi, lo, |h, l| h - l);
        let active: Vec<f64> = extents.iter().copied().filter(|&e| e > 0.0).collect();
        let blocks = samples.len().div_euclid(DATA_PER_BLOCK).max(1);
        let counts = if active.is_empty() {
            [1; 3]
        } else {
            // Cubes of side `side` would give about `blocks` blocks.
            let volume: f64 = active.iter().product();
            let side = (volume / to_f64(blocks)).powf(1.0 / to_f64(active.len()));
            extents.map(|e| {
                if e > 0.0 && side > 0.0 {
                    floor_index(e / side)
                        .map_or(MAX_BLOCKS_PER_AXIS, |n| n + 1)
                        .min(MAX_BLOCKS_PER_AXIS)
                } else {
                    1
                }
            })
        };
        Self::build(samples, lo, extents, counts)
    }

    /// Exactly `counts` blocks along x, y and z (each at least 1).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if a count is zero or above 1024.
    pub fn with_counts(samples: &[Sample], counts: [usize; 3]) -> Result<Self> {
        if counts.iter().any(|&n| n == 0 || n > MAX_BLOCKS_PER_AXIS) {
            return Err(Error::invalid(
                "super-block counts",
                format!("each must be 1 to {MAX_BLOCKS_PER_AXIS}, got {counts:?}"),
            ));
        }
        let (lo, hi) = extent(samples);
        Ok(Self::build(samples, lo, zip3(hi, lo, |h, l| h - l), counts))
    }

    /// The number of blocks along x, y and z.
    #[must_use]
    pub const fn counts(&self) -> [usize; 3] {
        self.counts
    }

    fn build(samples: &[Sample], origin: [f64; 3], extents: [f64; 3], counts: [usize; 3]) -> Self {
        let [e0, e1, e2] = extents;
        let [n0, n1, n2] = counts;
        let width = |e: f64, n: usize| if e > 0.0 { e / to_f64(n) } else { 1.0 };
        let size = [width(e0, n0), width(e1, n1), width(e2, n2)];
        let mut grid = Self {
            origin,
            size,
            counts,
            samples: samples.to_vec(),
            rows: Vec::new(),
            starts: Vec::new(),
            boxes: Vec::new(),
        };
        let nblocks = counts.iter().product::<usize>();
        let mut keyed: Vec<(usize, usize)> = grid
            .samples
            .iter()
            .enumerate()
            .map(|(row, s)| (grid.block_of(coords(s.location())), row))
            .collect();
        keyed.sort_unstable();
        grid.starts = vec![0; nblocks + 1];
        for &(b, _) in &keyed {
            if let Some(next) = grid.starts.get_mut(b + 1) {
                *next += 1;
            }
        }
        for b in 0..nblocks {
            let before = grid.starts.get(b).copied().unwrap_or(0);
            if let Some(next) = grid.starts.get_mut(b + 1) {
                *next += before;
            }
        }
        grid.rows = keyed.iter().map(|&(_, row)| row).collect();
        grid.boxes = (0..nblocks)
            .map(|b| {
                grid.block_rows(b).iter().fold(
                    ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
                    |(lo, hi), &row| {
                        let c = grid
                            .samples
                            .get(row)
                            .map_or([0.0; 3], |s| coords(s.location()));
                        (zip3(lo, c, f64::min), zip3(hi, c, f64::max))
                    },
                )
            })
            .collect();
        grid
    }

    /// The block index along each axis of a point (clamped to the grid).
    fn cell(&self, c: [f64; 3]) -> [usize; 3] {
        let [x, y, z] = zip3(zip3(c, self.origin, |c, o| c - o), self.size, |d, s| d / s);
        let [nx, ny, nz] = self.counts;
        [
            clamped_index(x, nx),
            clamped_index(y, ny),
            clamped_index(z, nz),
        ]
    }

    const fn linear([i, j, k]: [usize; 3], [nx, ny, _]: [usize; 3]) -> usize {
        i + nx * (j + ny * k)
    }

    fn block_of(&self, c: [f64; 3]) -> usize {
        Self::linear(self.cell(c), self.counts)
    }

    fn block_rows(&self, b: usize) -> &[usize] {
        let (start, end) = (
            self.starts.get(b).copied().unwrap_or(0),
            self.starts.get(b + 1).copied().unwrap_or(0),
        );
        self.rows.get(start..end).unwrap_or(&[])
    }
}

impl Neighbourhood for SuperBlock {
    fn samples(&self) -> &[Sample] {
        &self.samples
    }

    fn candidates(&self, ellipsoid: &SearchEllipsoid, target: Point3) -> Vec<Candidate> {
        let t = coords(target);
        let scale = ellipsoid.smallest_scale();
        let limit = ellipsoid.radius() * ellipsoid.radius();
        // Every candidate lies within Euclidean distance radius / scale.
        let reach = ellipsoid.radius() / scale * (1.0 + BOUND_MARGIN);
        let [i0, j0, k0] = self.cell(t.map(|x| x - reach));
        let [i1, j1, k1] = self.cell(t.map(|x| x + reach));
        let mut found = Vec::new();
        for k in k0..=k1 {
            for j in j0..=j1 {
                for i in i0..=i1 {
                    let b = Self::linear([i, j, k], self.counts);
                    let Some(&(bmin, bmax)) = self.boxes.get(b) else {
                        continue;
                    };
                    // Squared Euclidean distance from the target to the box.
                    let below = zip3(bmin, t, |m, x| (m - x).max(0.0));
                    let above = zip3(t, bmax, |x, m| (x - m).max(0.0));
                    let [g0, g1, g2] = zip3(below, above, f64::max);
                    let gap = g0.mul_add(g0, g1.mul_add(g1, g2 * g2));
                    if scale * scale * gap * (1.0 - BOUND_MARGIN) > limit {
                        continue;
                    }
                    for &row in self.block_rows(b) {
                        let Some(p) = self.samples.get(row).map(Sample::location) else {
                            continue;
                        };
                        let q = ellipsoid.squared_distance(p, target);
                        if ellipsoid.admits(q) {
                            found.push(Candidate {
                                row,
                                squared_distance: q,
                                offset: p.delta(target),
                            });
                        }
                    }
                }
            }
        }
        sort_candidates(&mut found);
        found
    }
}

/// Combine two triples component by component.
fn zip3(a: [f64; 3], b: [f64; 3], f: impl Fn(f64, f64) -> f64) -> [f64; 3] {
    let ([a0, a1, a2], [b0, b1, b2]) = (a, b);
    [f(a0, b0), f(a1, b1), f(a2, b2)]
}

fn coords(p: Point3) -> [f64; 3] {
    p.delta(Point3::default())
}

/// The bounding box of the samples' locations (empty: all zeros).
fn extent(samples: &[Sample]) -> ([f64; 3], [f64; 3]) {
    if samples.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    samples.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(lo, hi), s| {
            let c = coords(s.location());
            (zip3(lo, c, f64::min), zip3(hi, c, f64::max))
        },
    )
}

/// `floor(x)` clamped to `0..count` (NaN and negatives give 0).
fn clamped_index(x: f64, count: usize) -> usize {
    let last = count.saturating_sub(1);
    if x.is_nan() || x <= 0.0 {
        0
    } else {
        floor_index(x).map_or(last, |i| i.min(last))
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
    /// mapped back to the neighbourhood's data.
    ///
    /// # Errors
    /// As [`DataSearch::select`].
    pub fn select_samples(
        &self,
        neighbourhood: &impl Neighbourhood,
        target: Point3,
    ) -> Result<Vec<Sample>> {
        let samples = neighbourhood.samples();
        Ok(self
            .select(neighbourhood, target)?
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
            .select_samples(&BruteForce::new(&samples), Point3::default())
            .unwrap();
        assert_eq!(got, vec![samples[1], samples[2]]);
        assert_eq!(
            DataSearch::unlimited()
                .select_samples(&BruteForce::new(&samples), Point3::default())
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

    fn ellipsoid(az: f64, dip: f64, plunge: f64, radii: [f64; 3]) -> SearchEllipsoid {
        let [major, minor, vertical] = radii;
        SearchEllipsoid::new(
            Angles {
                azimuth: az,
                dip,
                plunge,
            },
            Ranges {
                major,
                minor,
                vertical,
            },
        )
        .unwrap()
    }

    #[test]
    fn super_block_handles_degenerate_data() {
        // S6: no data, one datum, duplicates, all data at one point.
        let e = SearchEllipsoid::isotropic(5.0).unwrap();
        let t = Point3::new(1.0, 1.0, 0.0);
        let cases: [Vec<Sample>; 4] = [
            vec![],
            vec![at(0.0, 0.0, 0.0)],
            vec![at(2.0, 1.0, 0.0), at(2.0, 1.0, 0.0), at(0.0, 0.0, 0.0)],
            vec![at(3.0, 3.0, 3.0); 5],
        ];
        for samples in cases {
            let index = SuperBlock::new(&samples);
            assert_eq!(
                index.candidates(&e, t),
                BruteForce::new(&samples).candidates(&e, t),
                "{samples:?}"
            );
        }
        assert_eq!(SuperBlock::new(&[]).counts(), [1, 1, 1]);
    }

    #[test]
    fn super_block_counts_follow_the_data() {
        // About four data per block, one block along a flat axis.
        let flat: Vec<Sample> = (0..400_u32)
            .map(|i| at(f64::from(i % 20), f64::from(i.div_euclid(20)), 0.0))
            .collect();
        let [nx, ny, nz] = SuperBlock::new(&flat).counts();
        assert_eq!(nz, 1, "2-D data");
        assert!((80..=150).contains(&(nx * ny)), "{nx} x {ny}");
        assert!(SuperBlock::with_counts(&flat, [0, 1, 1]).is_err(), "zero");
        assert!(
            SuperBlock::with_counts(&flat, [1, 1025, 1]).is_err(),
            "too many"
        );
    }

    #[test]
    fn super_block_finds_data_gslib_may_miss() {
        // The F-012 probe geometry: blocks of 10, azimuth 20, minor radius
        // 0.2 x 25. GSLIB's lattice test may exclude blocks holding data
        // inside this ellipsoid; the exact bound must not.
        let samples: Vec<Sample> = (0..=100_u32)
            .flat_map(|i| (0..=100_u32).map(move |j| at(f64::from(i), f64::from(j), 0.0)))
            .collect();
        let index = SuperBlock::with_counts(&samples, [10, 10, 1]).unwrap();
        let e = ellipsoid(20.0, 0.0, 0.0, [25.0, 5.0, 5.0]);
        for t in [Point3::new(45.0, 45.0, 0.0), Point3::new(55.5, 41.25, 0.0)] {
            let got = index.candidates(&e, t);
            assert!(!got.is_empty());
            assert_eq!(got, BruteForce::new(&samples).candidates(&e, t), "{t:?}");
        }
    }

    #[test]
    #[ignore = "timing evidence for the super-block; run with --release --ignored --nocapture"]
    fn super_block_timing() {
        // 20 000 scattered 2-D data, 2 000 targets, a search radius of 5% of
        // the domain: the case GSLIB's super-blocks exist for.
        let samples: Vec<Sample> = (0..20_000_u32)
            .map(|i| {
                let f = f64::from(i);
                at(
                    (f * 0.618_034).fract() * 1000.0,
                    (f * 0.414_214).fract() * 1000.0,
                    0.0,
                )
            })
            .collect();
        let targets: Vec<Point3> = (0..2_000_u32)
            .map(|i| {
                let f = f64::from(i);
                Point3::new(
                    (f * 0.732_051).fract() * 1000.0,
                    (f * 0.236_068).fract() * 1000.0,
                    0.0,
                )
            })
            .collect();
        let search = sphere(50.0, 24);
        let building = std::time::Instant::now();
        let index = SuperBlock::new(&samples);
        let build = building.elapsed();
        let time = |n: &dyn Fn(Point3) -> Vec<usize>| {
            let started = std::time::Instant::now();
            let total: usize = targets.iter().map(|&t| n(t).len()).sum();
            (started.elapsed(), total)
        };
        let (brute, nb) = time(&|t| search.select(&BruteForce::new(&samples), t).unwrap());
        let (fast, nf) = time(&|t| search.select(&index, t).unwrap());
        assert_eq!(nb, nf);
        println!(
            "brute force {brute:?}, super-block {fast:?} (+ {build:?} to build, {:?} blocks)",
            index.counts()
        );
    }

    proptest! {
        #[test]
        fn super_block_equals_brute_force(
            points in prop::collection::vec(
                prop_oneof![
                    (-20.0_f64..20.0, -20.0_f64..20.0, -5.0_f64..5.0),
                    (-4_i8..=4, -4_i8..=4, -1_i8..=1)
                        .prop_map(|(x, y, z)| (f64::from(x) * 5.0, f64::from(y) * 5.0, f64::from(z))),
                ],
                0..80,
            ),
            (az, dip, plunge) in (0.0_f64..360.0, -90.0_f64..90.0, -90.0_f64..90.0),
            (radius, r1, r2) in (0.5_f64..40.0, 0.05_f64..2.0, 0.05_f64..2.0),
            target in (-60.0_f64..60.0, -60.0_f64..60.0, -10.0_f64..10.0),
            counts in prop::option::of((1_usize..7, 1_usize..7, 1_usize..4)),
        ) {
            // S6, D-025: ratios from 0.05 to 2 (so the smallest scale can be
            // below 1), targets inside and outside the data, any block grid.
            let samples: Vec<Sample> = points.iter().map(|&(x, y, z)| at(x, y, z)).collect();
            let index = counts.map_or_else(
                || SuperBlock::new(&samples),
                |c| SuperBlock::with_counts(&samples, c.into()).unwrap(),
            );
            let e = ellipsoid(az, dip, plunge, [radius, radius * r1, radius * r2]);
            let t = Point3::new(target.0, target.1, target.2);
            prop_assert_eq!(index.candidates(&e, t), BruteForce::new(&samples).candidates(&e, t));
            let unlimited = SearchEllipsoid::unlimited();
            prop_assert_eq!(
                index.candidates(&unlimited, t),
                BruteForce::new(&samples).candidates(&unlimited, t)
            );
        }
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
