//! Regular grids in GSLIB's convention.

use crate::error::{Error, Result};
use crate::geometry::Point3;
use crate::num::{floor_index, to_f64};

/// One axis of a regular grid: `count` cells of width `spacing`, the first
/// centred at `first_center` (GSLIB `nx`, `xmn`, `xsiz`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis {
    /// Number of cells.
    pub count: usize,
    /// Coordinate of the first cell's centre.
    pub first_center: f64,
    /// Cell width.
    pub spacing: f64,
}

impl Axis {
    fn center(&self, i: usize) -> f64 {
        self.spacing.mul_add(to_f64(i), self.first_center)
    }

    /// Index of the cell containing `x`, if inside the axis.
    fn cell(&self, x: f64) -> Option<usize> {
        floor_index((x - self.first_center) / self.spacing + 0.5).filter(|&i| i < self.count)
    }
}

/// A regular 3-D grid. Nodes are ordered x fastest, then y, then z, as in
/// GSLIB grid files.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid {
    x: Axis,
    y: Axis,
    z: Axis,
}

impl Grid {
    /// # Errors
    /// [`Error::InvalidInput`] if any axis has no cells, a non-finite origin,
    /// or a spacing that is not finite and positive.
    pub fn new(x: Axis, y: Axis, z: Axis) -> Result<Self> {
        for axis in [x, y, z] {
            if axis.count == 0 {
                return Err(Error::invalid("grid", "every axis needs at least one cell"));
            }
            if !(axis.first_center.is_finite() && axis.spacing.is_finite() && axis.spacing > 0.0) {
                return Err(Error::invalid(
                    "grid",
                    "origins must be finite and spacings positive",
                ));
            }
        }
        if x.count
            .checked_mul(y.count)
            .and_then(|n| n.checked_mul(z.count))
            .is_none()
        {
            return Err(Error::invalid("grid", "too many nodes"));
        }
        Ok(Self { x, y, z })
    }

    /// Number of nodes.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.x.count * self.y.count * self.z.count
    }

    /// Whether the grid has no nodes (never true for a validated grid).
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Node coordinates, x fastest.
    #[must_use]
    pub fn nodes(&self) -> Vec<Point3> {
        (0..self.z.count)
            .flat_map(|iz| {
                (0..self.y.count).flat_map(move |iy| {
                    (0..self.x.count).map(move |ix| {
                        Point3::new(self.x.center(ix), self.y.center(iy), self.z.center(iz))
                    })
                })
            })
            .collect()
    }

    /// Index of the node whose cell contains `p`, or `None` outside the grid.
    #[must_use]
    pub fn node_index(&self, p: Point3) -> Option<usize> {
        let ix = self.x.cell(p.x)?;
        let iy = self.y.cell(p.y)?;
        let iz = self.z.cell(p.z)?;
        Some(ix + self.x.count * (iy + self.y.count * iz))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axis(count: usize) -> Axis {
        Axis {
            count,
            first_center: 0.5,
            spacing: 1.0,
        }
    }

    #[test]
    fn nodes_are_x_fastest() {
        let g = Grid::new(axis(2), axis(2), axis(1)).unwrap();
        let n = g.nodes();
        assert_eq!(n.len(), 4);
        assert_eq!(n[1], Point3::new(1.5, 0.5, 0.5));
        assert_eq!(n[2], Point3::new(0.5, 1.5, 0.5));
    }

    #[test]
    fn node_index_matches_nodes() {
        let g = Grid::new(axis(3), axis(2), axis(2)).unwrap();
        for (i, p) in g.nodes().into_iter().enumerate() {
            assert_eq!(g.node_index(p), Some(i));
        }
        assert_eq!(g.node_index(Point3::new(-0.1, 0.5, 0.5)), None);
        assert_eq!(g.node_index(Point3::new(3.1, 0.5, 0.5)), None);
    }

    #[test]
    fn rejects_empty_axes() {
        assert!(Grid::new(axis(0), axis(1), axis(1)).is_err(), "zero cells");
    }
}
