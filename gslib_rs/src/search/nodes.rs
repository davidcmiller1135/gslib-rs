//! The search for previously simulated grid nodes in sequential simulation
//! (GSLIB `sgsim`'s `ctable` and `srchnd`; `docs/spec/search.md` S8–S10,
//! D-026, F-012, F-013).
//!
//! # Formulation
//!
//! - **Offset table (S8):** every grid offset `δ ≠ 0` (in cells) whose
//!   separation `h` satisfies `q(h) ≤ radius²` under the search ellipsoid,
//!   ordered by decreasing covariance `C(h)`, then increasing `q(h)`, then
//!   the offset's linear index (x fastest over the table's box). GSLIB sorts
//!   by `−C(h) + 10⁻¹⁰·q(h)` with an unstable sort; the explicit
//!   lexicographic key makes the order total (D-024, D-026).
//! - **Selection (S9):** from node `n`, offsets are scanned in table order;
//!   offsets leaving the grid are skipped; a node is kept if it is informed
//!   (simulated) and its octant holds fewer than `per_octant` kept nodes; the
//!   scan stops at `max_nodes`. Octants follow `srchnd`'s convention, which
//!   uses `target − node` and splits on `dz > 0` (F-012), unlike the data
//!   search.
//! - **Extent (S10, F-013):** by default the table covers the whole
//!   ellipsoid (within the grid). **Divergence from GSLIB** (ROADMAP §2):
//!   GSLIB's table is capped at `min((MAXCT − 1)/2, n − 1)` offsets per axis,
//!   ±35/±35/0 in the stock build, so it never searches vertically and
//!   silently ignores nodes inside the radius beyond the cap. That cap is the
//!   opt-in [`TableExtent::Gslib`].

use std::cmp::Ordering;
use std::num::NonZeroUsize;

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::model::CovarianceModel;
use crate::num::{floor_index, to_f64};
use crate::sample::Sample;
use crate::simulation::Grid;

use super::SearchEllipsoid;

/// How far the offset table reaches (S10).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TableExtent {
    /// The whole search ellipsoid, within the grid.
    #[default]
    Ellipsoid,
    /// GSLIB's compile-time cap: at most `min((max − 1)/2, n − 1)` offsets
    /// along each axis, for table sizes `max` (`MAXCTX`, `MAXCTY`, `MAXCTZ`;
    /// 71, 71, 1 in the stock `sgsim`).
    Gslib {
        /// Table sizes along x, y and z.
        max: [usize; 3],
    },
}

/// One entry of the offset table.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Offset {
    cells: [i64; 3],
    covariance: f64,
    squared_distance: f64,
    linear: usize,
}

/// The search for previously simulated nodes of one grid (S8–S10).
#[derive(Debug, Clone, PartialEq)]
pub struct NodeSearch {
    counts: [usize; 3],
    table: Vec<Offset>,
    max_nodes: usize,
    per_octant: Option<NonZeroUsize>,
}

/// The largest offset along one axis: within the grid, the ellipsoid's
/// reach, and the extent's cap.
fn reach(count: usize, spacing: f64, radius: f64, scale: f64, cap: Option<usize>) -> usize {
    let grid = count.saturating_sub(1);
    // Every offset inside the ellipsoid is within Euclidean radius / scale.
    let ellipsoid = floor_index(radius / scale / spacing).unwrap_or(usize::MAX);
    cap.map_or(grid, |c| grid.min(c)).min(ellipsoid)
}

fn signed(n: usize) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

fn to_f64_signed(n: i64) -> f64 {
    let magnitude = to_f64(usize::try_from(n.unsigned_abs()).unwrap_or(usize::MAX));
    if n < 0 { -magnitude } else { magnitude }
}

impl NodeSearch {
    /// The offset table for `grid`, ranked by `model`'s covariance within
    /// `ellipsoid`, keeping at most `max_nodes` nodes (GSLIB `nodmax`) and
    /// optionally at most `per_octant` per octant (`noct`).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `max_nodes` is zero.
    pub fn new(
        grid: &Grid,
        ellipsoid: &SearchEllipsoid,
        model: &CovarianceModel,
        max_nodes: usize,
        per_octant: Option<NonZeroUsize>,
        extent: TableExtent,
    ) -> Result<Self> {
        if max_nodes == 0 {
            return Err(Error::invalid(
                "node search",
                "max_nodes must be at least 1",
            ));
        }
        let counts = grid.counts();
        let [dx, dy, dz] = grid.spacing();
        let caps = match extent {
            TableExtent::Ellipsoid => [None; 3],
            TableExtent::Gslib { max } => max.map(|m| Some(m.saturating_sub(1).div_euclid(2))),
        };
        let scale = ellipsoid.smallest_scale();
        let [cx, cy, cz] = caps;
        let [nx, ny, nz] = counts;
        let r = ellipsoid.radius();
        let (ex, ey, ez) = (
            reach(nx, dx, r, scale, cx),
            reach(ny, dy, r, scale, cy),
            reach(nz, dz, r, scale, cz),
        );
        let (wx, wy) = (2 * ex + 1, 2 * ey + 1);
        let mut table = Vec::new();
        for k in 0..=2 * ez {
            for j in 0..=2 * ey {
                for i in 0..=2 * ex {
                    let cells = [
                        signed(i) - signed(ex),
                        signed(j) - signed(ey),
                        signed(k) - signed(ez),
                    ];
                    if cells == [0; 3] {
                        continue;
                    }
                    let [ci, cj, ck] = cells.map(to_f64_signed);
                    let h = Point3::new(ci * dx, cj * dy, ck * dz);
                    let q = ellipsoid.squared_distance(Point3::default(), h);
                    if ellipsoid.admits(q) {
                        table.push(Offset {
                            cells,
                            covariance: model.covariance(Point3::default(), h),
                            squared_distance: q,
                            linear: i + wx * (j + wy * k),
                        });
                    }
                }
            }
        }
        table.sort_unstable_by(order);
        Ok(Self {
            counts,
            table,
            max_nodes,
            per_octant,
        })
    }

    /// The ordered offset table, in cells (S8).
    #[must_use]
    pub fn offsets(&self) -> Vec<[i64; 3]> {
        self.table.iter().map(|o| o.cells).collect()
    }

    /// The grid shape this search was built for.
    #[must_use]
    pub const fn counts(&self) -> [usize; 3] {
        self.counts
    }

    /// The informed nodes to use when simulating node `index` of `grid`, in
    /// table order (S9). `informed(i)` says whether node `i` already holds a
    /// simulated value.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `grid` does not have the shape this search
    /// was built for, or `index` is outside it.
    pub fn select(
        &self,
        grid: &Grid,
        index: usize,
        informed: impl Fn(usize) -> bool,
    ) -> Result<Vec<usize>> {
        if grid.counts() != self.counts {
            return Err(Error::invalid(
                "node search",
                format!(
                    "built for a {:?} grid, used with {:?}",
                    self.counts,
                    grid.counts()
                ),
            ));
        }
        if index >= grid.len() {
            return Err(Error::invalid(
                "node search",
                format!("node {index} is outside the grid"),
            ));
        }
        let target = grid.cell_of(index).map(signed);
        let mut per_octant = [0_usize; 8];
        let mut kept = Vec::new();
        for o in &self.table {
            if kept.len() == self.max_nodes {
                break;
            }
            let [tx, ty, tz] = target;
            let [ox, oy, oz] = o.cells;
            let cell = [tx + ox, ty + oy, tz + oz];
            let Some(node) = cell
                .iter()
                .map(|&c| usize::try_from(c).ok())
                .collect::<Option<Vec<usize>>>()
                .and_then(|c| <[usize; 3]>::try_from(c).ok())
                .and_then(|c| grid.index_of(c))
            else {
                continue;
            };
            if !informed(node) {
                continue;
            }
            if let Some(limit) = self.per_octant {
                // srchnd: octant of target − node.
                let slot = usize::from(node_octant([-ox, -oy, -oz]) - 1);
                if let Some(count) = per_octant.get_mut(slot) {
                    *count += 1;
                    if *count > limit.get() {
                        continue;
                    }
                }
            }
            kept.push(node);
        }
        Ok(kept)
    }
}

/// The table order: decreasing covariance, increasing squared distance,
/// increasing linear index.
fn order(a: &Offset, b: &Offset) -> Ordering {
    b.covariance
        .total_cmp(&a.covariance)
        .then(a.squared_distance.total_cmp(&b.squared_distance))
        .then(a.linear.cmp(&b.linear))
}

/// `srchnd`'s octant (1–8) of `target − node`, in cells.
///
/// 1–4 when `dz > 0`, 5–8 otherwise, with the quadrant rules of the data
/// search: 1 if `dx ≤ 0, dy > 0`; 2 if `dx > 0, dy ≥ 0`; 3 if
/// `dx < 0, dy ≤ 0`; else 4.
#[must_use]
pub const fn node_octant([dx, dy, dz]: [i64; 3]) -> u8 {
    let quadrant = if dx <= 0 && dy > 0 {
        1
    } else if dx > 0 && dy >= 0 {
        2
    } else if dx < 0 && dy <= 0 {
        3
    } else {
        4
    };
    if dz > 0 { quadrant } else { quadrant + 4 }
}

/// What assigning data to nodes does when two data share a node (S11,
/// D-028).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OnSharedNode {
    /// Fail, naming both data and the node.
    #[default]
    Error,
    /// GSLIB's rule (`SGSIM.FOR` lines 869–881): keep the datum with the
    /// smaller L1 distance to the node; on an exact tie the later datum
    /// wins. The other datum is dropped from the simulation.
    KeepClosest,
}

/// What assigning data to nodes does with a datum outside the grid (S12,
/// F-013).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutsideGrid {
    /// Fail, naming the datum.
    #[default]
    Error,
    /// GSLIB's behaviour: `getindx` clamps each index to the grid, so the
    /// datum goes to the nearest edge node, however far away it lies.
    ClampToEdge,
}

/// The two-part search's first step (GSLIB `sgsim` with `sstrat = 1`):
/// move each datum to its nearest node. Returns, per node (x fastest), the
/// row of the datum it holds.
///
/// A datum's node is the one whose cell contains it (GSLIB's `getindx` per
/// axis). Data are taken in row order.
///
/// # Errors
/// [`Error::InvalidInput`], naming the rows and node, when two data share a
/// node under [`OnSharedNode::Error`], or a datum lies outside the grid
/// under [`OutsideGrid::Error`].
pub fn assign_to_nodes(
    data: &[Sample],
    grid: &Grid,
    shared: OnSharedNode,
    outside: OutsideGrid,
) -> Result<Vec<Option<usize>>> {
    let mut held: Vec<Option<usize>> = vec![None; grid.len()];
    for (row, d) in data.iter().enumerate() {
        let (cell, inside) = grid.clamped_cell(d.location());
        if !inside && outside == OutsideGrid::Error {
            return Err(Error::invalid(
                "data on the grid",
                format!(
                    "datum {} at {:?} lies outside the simulation grid",
                    row + 1,
                    d.location()
                ),
            ));
        }
        let Some(node) = grid.index_of(cell) else {
            continue;
        };
        let Some(slot) = held.get_mut(node) else {
            continue;
        };
        match (*slot, shared) {
            (None, _) => *slot = Some(row),
            (Some(earlier), OnSharedNode::Error) => {
                return Err(Error::invalid(
                    "data on the grid",
                    format!(
                        "data {} and {} both fall on node {} at {:?}",
                        earlier + 1,
                        row + 1,
                        node + 1,
                        grid.node(node)
                    ),
                ));
            }
            (Some(earlier), OnSharedNode::KeepClosest) => {
                let at = grid.node(node);
                let l1 = |r: usize| {
                    data.get(r).map_or(f64::INFINITY, |s| {
                        let [dx, dy, dz] = s.location().delta(at);
                        dx.abs() + dy.abs() + dz.abs()
                    })
                };
                // GSLIB: `if(test.le.test2) sim(ind) = real(id)`.
                if l1(row) <= l1(earlier) {
                    *slot = Some(row);
                }
            }
        }
    }
    Ok(held)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Angles, Ranges};
    use crate::model::{Kind, Structure, VariogramModel};
    use crate::simulation::Axis;
    use proptest::prelude::*;

    fn grid(n: [usize; 3], spacing: [f64; 3]) -> Grid {
        let axis = |count, width| Axis {
            count,
            first_center: 0.5,
            spacing: width,
        };
        Grid::new(
            axis(n[0], spacing[0]),
            axis(n[1], spacing[1]),
            axis(n[2], spacing[2]),
        )
        .unwrap()
    }

    fn model(range: f64, minor: f64, azimuth: f64) -> CovarianceModel {
        let s = Structure::new(
            Kind::Spherical,
            1.0,
            Angles {
                azimuth,
                dip: 0.0,
                plunge: 0.0,
            },
            Ranges {
                major: range,
                minor,
                vertical: range,
            },
        )
        .unwrap();
        CovarianceModel::try_from(VariogramModel::new(0.0, vec![s]).unwrap()).unwrap()
    }

    fn sphere(r: f64) -> SearchEllipsoid {
        SearchEllipsoid::isotropic(r).unwrap()
    }

    /// The definition of S8, from every pair of grid nodes rather than the
    /// table: offsets of all nodes from a corner node, filtered and sorted.
    fn table_by_definition(
        g: &Grid,
        ellipsoid: &SearchEllipsoid,
        covariance: &CovarianceModel,
    ) -> Vec<[i64; 3]> {
        let [nx, ny, nz] = g.counts().map(signed);
        let [dx, dy, dz] = g.spacing();
        let mut all: Vec<([i64; 3], f64, f64, i64)> = Vec::new();
        for k in -(nz - 1)..nz {
            for j in -(ny - 1)..ny {
                for i in -(nx - 1)..nx {
                    if [i, j, k] == [0; 3] {
                        continue;
                    }
                    let sep = Point3::new(
                        to_f64_signed(i) * dx,
                        to_f64_signed(j) * dy,
                        to_f64_signed(k) * dz,
                    );
                    let q = ellipsoid.squared_distance(Point3::default(), sep);
                    if ellipsoid.admits(q) {
                        // Linear order: k slowest, then j, then i, as the table's.
                        all.push((
                            [i, j, k],
                            covariance.covariance(Point3::default(), sep),
                            q,
                            (k * 10_000 + j) * 10_000 + i,
                        ));
                    }
                }
            }
        }
        all.sort_by(|a, b| {
            b.1.total_cmp(&a.1)
                .then(a.2.total_cmp(&b.2))
                .then(a.3.cmp(&b.3))
        });
        all.into_iter().map(|t| t.0).collect()
    }

    #[test]
    fn table_is_ranked_by_covariance_then_distance() {
        // S8 on a small anisotropic case: the major axis (azimuth 90, along
        // x) has the larger covariance, so (±1, 0) come before (0, ±1).
        let g = grid([5, 5, 1], [1.0, 1.0, 1.0]);
        let ns = NodeSearch::new(
            &g,
            &sphere(3.0),
            &model(10.0, 3.0, 90.0),
            4,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap();
        let t = ns.offsets();
        assert_eq!(&t[..2], &[[-1, 0, 0], [1, 0, 0]], "{t:?}");
        assert!(!t.contains(&[0, 0, 0]), "never the node itself");
        assert!(
            t.iter().all(|&[i, j, _]| i * i + j * j <= 9),
            "within the radius"
        );
        assert_eq!(
            t,
            table_by_definition(&g, &sphere(3.0), &model(10.0, 3.0, 90.0))
        );
    }

    #[test]
    fn node_octants_follow_srchnd() {
        // target − node; upper octants only for dz > 0.
        assert_eq!(node_octant([0, 1, 1]), 1);
        assert_eq!(node_octant([1, 0, 1]), 2);
        assert_eq!(node_octant([-1, 0, 1]), 3);
        assert_eq!(node_octant([0, -1, 1]), 4);
        assert_eq!(node_octant([0, 0, 1]), 4);
        assert_eq!(node_octant([0, 1, 0]), 5, "dz = 0 is a lower octant here");
        assert_eq!(node_octant([1, 1, -1]), 6);
        assert_eq!(node_octant([-1, -1, 0]), 7);
        assert_eq!(node_octant([1, -1, 0]), 8);
    }

    /// The definition of S9: every informed node within the ellipsoid, by
    /// table order, with the octant cap and `max_nodes`.
    fn select_by_definition(
        ns: &NodeSearch,
        g: &Grid,
        index: usize,
        informed: &[bool],
        per_octant: Option<usize>,
        max: usize,
    ) -> Vec<usize> {
        let t = g.cell_of(index).map(signed);
        let mut ranked: Vec<(usize, usize, [i64; 3])> = (0..g.len())
            .filter(|&m| m != index && informed[m])
            .filter_map(|m| {
                let c = g.cell_of(m).map(signed);
                let d = [c[0] - t[0], c[1] - t[1], c[2] - t[2]];
                ns.offsets()
                    .iter()
                    .position(|&o| o == d)
                    .map(|rank| (rank, m, d))
            })
            .collect();
        ranked.sort_unstable();
        let mut counts = [0_usize; 8];
        let mut kept = Vec::new();
        for (_, m, d) in ranked {
            if kept.len() == max {
                break;
            }
            if let Some(limit) = per_octant {
                let slot = usize::from(node_octant([-d[0], -d[1], -d[2]]) - 1);
                counts[slot] += 1;
                if counts[slot] > limit {
                    continue;
                }
            }
            kept.push(m);
        }
        kept
    }

    proptest! {
        #[test]
        fn table_and_selection_match_their_definitions(
            n in (1_usize..7, 1_usize..7, 1_usize..4),
            spacing in (0.5_f64..3.0, 0.5_f64..3.0, 0.5_f64..3.0),
            radius in 0.5_f64..12.0,
            (range, ratio, azimuth) in (1.0_f64..15.0, 0.2_f64..1.0, 0.0_f64..180.0),
            informed in prop::collection::vec(any::<bool>(), 150),
            index_seed in 0_usize..1000,
            max in 1_usize..12,
            per_octant in prop::option::of(1_usize..3),
        ) {
            // S8 and S9 against independent definitions.
            let g = grid([n.0, n.1, n.2], [spacing.0, spacing.1, spacing.2]);
            let e = sphere(radius);
            let m = model(range, range * ratio, azimuth);
            let ns = NodeSearch::new(&g, &e, &m, max, per_octant.and_then(NonZeroUsize::new), TableExtent::Ellipsoid)
                .unwrap();
            prop_assert_eq!(ns.offsets(), table_by_definition(&g, &e, &m));
            let index = index_seed % g.len();
            let flags: Vec<bool> = (0..g.len()).map(|i| informed[i % informed.len()]).collect();
            let got = ns.select(&g, index, |i| flags[i]).unwrap();
            prop_assert_eq!(got, select_by_definition(&ns, &g, index, &flags, per_octant, max));
        }
    }

    #[test]
    fn gslib_extent_caps_the_table() {
        // S10: a 1-D grid of 100 nodes and an unbounded search. By default
        // a node 40 cells away is found; GSLIB's 71/71/1 table stops at 35.
        let g = grid([100, 1, 1], [1.0, 1.0, 1.0]);
        let e = SearchEllipsoid::unlimited();
        let m = model(200.0, 200.0, 0.0);
        let only_40 = |i: usize| i == 40;
        let full = NodeSearch::new(&g, &e, &m, 5, None, TableExtent::Ellipsoid).unwrap();
        assert_eq!(full.select(&g, 0, only_40).unwrap(), [40]);
        let gslib =
            NodeSearch::new(&g, &e, &m, 5, None, TableExtent::Gslib { max: [71, 71, 1] }).unwrap();
        assert_eq!(gslib.select(&g, 0, only_40).unwrap(), Vec::<usize>::new());
        assert_eq!(gslib.select(&g, 0, |i| i == 35).unwrap(), [35]);
        // And no vertical search at all under GSLIB's MAXCTZ = 1.
        let column = grid([1, 1, 10], [1.0, 1.0, 1.0]);
        let up = NodeSearch::new(
            &column,
            &e,
            &m,
            5,
            None,
            TableExtent::Gslib { max: [71, 71, 1] },
        )
        .unwrap();
        assert!(up.offsets().is_empty(), "nctz = 0");
        let ours = NodeSearch::new(&column, &e, &m, 5, None, TableExtent::Ellipsoid).unwrap();
        assert_eq!(ours.select(&column, 0, |i| i == 9).unwrap(), [9]);
    }

    #[test]
    fn selection_validates_its_grid() {
        let g = grid([3, 3, 1], [1.0, 1.0, 1.0]);
        let ns = NodeSearch::new(
            &g,
            &sphere(2.0),
            &model(5.0, 5.0, 0.0),
            4,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap();
        let other = grid([4, 3, 1], [1.0, 1.0, 1.0]);
        assert!(ns.select(&other, 0, |_| true).is_err(), "different shape");
        assert!(ns.select(&g, 9, |_| true).is_err(), "outside the grid");
        assert!(
            NodeSearch::new(
                &g,
                &sphere(2.0),
                &model(5.0, 5.0, 0.0),
                0,
                None,
                TableExtent::Ellipsoid
            )
            .is_err(),
            "max_nodes 0"
        );
        // max_nodes caps the result; nearest by covariance first.
        let two = NodeSearch::new(
            &g,
            &sphere(2.0),
            &model(5.0, 5.0, 0.0),
            2,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap();
        let got = two.select(&g, 4, |_| true).unwrap();
        assert_eq!(got.len(), 2);
        assert!(
            got.iter().all(|&i| [1, 3, 5, 7].contains(&i)),
            "edge neighbours first: {got:?}"
        );
    }

    fn datum(x: f64, y: f64, v: f64) -> Sample {
        Sample::new(Point3::new(x, y, 0.0), v).unwrap()
    }

    #[test]
    fn data_go_to_their_nearest_nodes() {
        // Nodes at 0.5, 1.5, 2.5 on x and y.
        let g = grid([3, 3, 1], [1.0, 1.0, 1.0]);
        let data = [datum(0.6, 0.4, 1.0), datum(2.4, 1.6, 2.0)];
        let held = assign_to_nodes(&data, &g, OnSharedNode::Error, OutsideGrid::Error).unwrap();
        assert_eq!(held[0], Some(0));
        assert_eq!(held[2 + 3], Some(1), "node (2, 1)");
        assert_eq!(held.iter().flatten().count(), 2);
    }

    #[test]
    fn shared_nodes_are_errors_or_keep_the_closest() {
        // S11: two data in node (1, 1)'s cell, the second closer.
        let g = grid([3, 3, 1], [1.0, 1.0, 1.0]);
        let data = [datum(1.2, 1.9, 1.0), datum(1.4, 1.5, 2.0)];
        let err = assign_to_nodes(&data, &g, OnSharedNode::Error, OutsideGrid::Error).unwrap_err();
        let message = err.to_string();
        assert!(
            message.contains("data 1 and 2") && message.contains("node 5"),
            "{message}"
        );
        let held =
            assign_to_nodes(&data, &g, OnSharedNode::KeepClosest, OutsideGrid::Error).unwrap();
        assert_eq!(held[4], Some(1), "the closer, later datum");
        // The closer first: kept.
        let swapped = [data[1], data[0]];
        let first_closer =
            assign_to_nodes(&swapped, &g, OnSharedNode::KeepClosest, OutsideGrid::Error).unwrap();
        assert_eq!(first_closer[4], Some(0), "the closer, earlier datum");
        // An exact tie: the later wins (test.le.test2).
        let tie = [datum(1.25, 1.5, 1.0), datum(1.75, 1.5, 2.0)];
        let tied =
            assign_to_nodes(&tie, &g, OnSharedNode::KeepClosest, OutsideGrid::Error).unwrap();
        assert_eq!(tied[4], Some(1), "later datum on a tie");
    }

    #[test]
    fn data_outside_the_grid_are_errors_or_clamped() {
        // S12: one datum beyond each face, cells (0..3) x (0..3) spanning
        // [0, 3] on each axis.
        let g = grid([3, 3, 1], [1.0, 1.0, 1.0]);
        let cases = [
            (datum(-0.2, 1.5, 1.0), 3),  // x below: node (0, 1)
            (datum(3.7, 1.5, 1.0), 5),   // x above: node (2, 1)
            (datum(1.5, -40.0, 1.0), 1), // y below: node (1, 0)
            (datum(1.5, 3.1, 1.0), 7),   // y above: node (1, 2)
        ];
        for (d, node) in cases {
            let err =
                assign_to_nodes(&[d], &g, OnSharedNode::Error, OutsideGrid::Error).unwrap_err();
            assert!(err.to_string().contains("datum 1"), "{err}");
            let held =
                assign_to_nodes(&[d], &g, OnSharedNode::Error, OutsideGrid::ClampToEdge).unwrap();
            assert_eq!(held[node], Some(0), "{d:?} clamps to node {node}");
        }
        // z faces, on a 1 x 1 x 3 column.
        let column = grid([1, 1, 3], [1.0, 1.0, 1.0]);
        for (z, node) in [(-1.0, 0), (9.0, 2)] {
            let d = Sample::new(Point3::new(0.5, 0.5, z), 1.0).unwrap();
            assert!(
                assign_to_nodes(&[d], &column, OnSharedNode::Error, OutsideGrid::Error).is_err()
            );
            let held =
                assign_to_nodes(&[d], &column, OnSharedNode::Error, OutsideGrid::ClampToEdge)
                    .unwrap();
            assert_eq!(held[node], Some(0));
        }
        // Exactly on the outer edge of the first and last cells: inside, as
        // getindx's int(v + 1.5) gives.
        let edges = [datum(0.0, 0.0, 1.0), datum(2.9999, 2.9999, 2.0)];
        assert!(assign_to_nodes(&edges, &g, OnSharedNode::Error, OutsideGrid::Error).is_ok());
    }
}
