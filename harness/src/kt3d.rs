//! Reading `kt3d` output, and the parts of GSLIB the comparison must
//! reproduce rather than share (F-012, F-014; `docs/spec/search.md` step 3).

use gslib_rs::geometry::{Angles, Point3, Ranges};
use gslib_rs::model::VariogramModel;
use gslib_rs::sample::Sample;
use gslib_rs::search::{BruteForce, Candidate, Neighbourhood, SearchEllipsoid};
use oracle_case::kt3d::Params;
use oracle_case::vmodel::Structure;

/// `DEG2RAD` as declared in `GSLIB/SETROT.FOR` line 28 (truncated π).
#[expect(
    clippy::approx_constant,
    reason = "deliberately GSLIB's truncated pi from SETROT.FOR, not std's PI"
)]
pub const SETROT_DEG2RAD: f64 = 3.141_592_654 / 180.0;

/// `kt3d`'s super-block limits (`KT3D.INC`: `MAXSBX = MAXSBY = MAXSBZ = 21`).
pub const MAX_SUPER_BLOCKS: u32 = 21;

/// GSLIB's angles, converted so gslib-rs reproduces `setrot`'s radians.
///
/// `setrot` uses `(90 − az)·D`, `−dip·D` and `plunge·D` with `D` =
/// [`SETROT_DEG2RAD`]. Converting first removes the documented π divergence
/// (ROADMAP §2) from the comparison, so what remains is arithmetic.
///
/// # Errors
/// If the azimuth is 270° or more: `setrot` then uses `(450 − az)·D`, which
/// with a truncated π is not an exact turn away, so no single conversion
/// reproduces both branches. Fixtures stay below 270°.
pub fn gslib_angles([azimuth, dip, plunge]: [f64; 3]) -> Result<Angles, String> {
    if azimuth >= 270.0 {
        return Err(format!(
            "azimuth {azimuth} >= 270 takes setrot's other branch; use a smaller azimuth"
        ));
    }
    // gslib-rs converts x degrees to x·π/180; GSLIB converts to x·D.
    let factor = SETROT_DEG2RAD / 1.0_f64.to_radians();
    Ok(Angles {
        azimuth: (90.0 - azimuth).mul_add(-factor, 90.0),
        dip: dip * factor,
        plunge: plunge * factor,
    })
}

/// The gslib-rs model for a case's nugget and structures, with GSLIB's
/// angle conversion ([`gslib_angles`]).
///
/// # Errors
/// A message if gslib-rs rejects the model or an azimuth is 270° or more.
pub fn model(nugget: f64, structures: &[Structure]) -> Result<VariogramModel, String> {
    let rescaled: Vec<Structure> = structures
        .iter()
        .map(|s| {
            let a = gslib_angles(s.angles)?;
            Ok(Structure {
                angles: [a.azimuth, a.dip, a.plunge],
                ..*s
            })
        })
        .collect::<Result<_, String>>()?;
    crate::vmodel::model_from(nugget, &rescaled)
}

/// The case's search ellipsoid, with GSLIB's angle conversion.
///
/// # Errors
/// A message if gslib-rs rejects it or the azimuth is 270° or more.
pub fn ellipsoid(p: &Params) -> Result<SearchEllipsoid, String> {
    let [major, minor, vertical] = p.search_radii;
    SearchEllipsoid::new(
        gslib_angles(p.search_angles)?,
        Ranges {
            major,
            minor,
            vertical,
        },
    )
    .map_err(|e| e.to_string())
}

/// Grid node `index` (0-based, x fastest), as `kt3d` locates it: 1-based
/// indices and `xmn + (ix − 1)·xsiz` etc.
#[must_use]
#[expect(
    clippy::suboptimal_flops,
    reason = "GSLIB is built with -ffp-contract=off: reproduce its unfused arithmetic"
)]
pub fn node(p: &Params, index: usize) -> ([u32; 3], Point3) {
    let [nx, ny, _] = p
        .grid_counts
        .map(|n| usize::try_from(n).unwrap_or(1).max(1));
    let (ix, rest) = (index % nx, index.div_euclid(nx));
    let (iy, iz) = (rest % ny, rest.div_euclid(ny));
    let ijk = [ix, iy, iz].map(|i| u32::try_from(i).unwrap_or(u32::MAX));
    let [i0, i1, i2] = ijk;
    let [x0, y0, z0] = p.grid_origin;
    let [dx, dy, dz] = p.grid_spacing;
    let location = Point3::new(
        x0 + f64::from(i0) * dx,
        y0 + f64::from(i1) * dy,
        z0 + f64::from(i2) * dz,
    );
    (ijk.map(|i| i.saturating_add(1)), location)
}

/// `kt3d.out`: `(estimate, variance)` per node, `−999` where unestimated.
///
/// # Errors
/// A message naming the first malformed line.
pub fn parse_output(text: &str) -> Result<Vec<(f64, f64)>, String> {
    let mut lines = text.lines();
    lines.next().ok_or("no title line")?;
    let nvar = lines.next().ok_or("no variable count")?.trim();
    if nvar != "2" {
        return Err(format!("expected 2 variables, got {nvar:?}"));
    }
    for want in ["Estimate", "EstimationVariance"] {
        let got = lines.next().ok_or("missing a variable name")?.trim();
        if got != want {
            return Err(format!("expected {want:?}, got {got:?}"));
        }
    }
    lines
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(i, l)| {
            let v: Vec<f64> = l
                .split_whitespace()
                .map(str::parse)
                .collect::<Result<_, _>>()
                .map_err(|e| format!("data line {}: {e}", i + 1))?;
            match v.as_slice() {
                [e, s] => Ok((*e, *s)),
                _ => Err(format!("data line {}: expected 2 values", i + 1)),
            }
        })
        .collect()
}

/// One solved system in `kt3d.dbg` (`idbg = 2`).
#[derive(Debug, Clone, PartialEq)]
pub struct Block {
    /// 1-based node indices.
    pub node: [u32; 3],
    /// The node's coordinates.
    pub location: [f64; 3],
    /// GSLIB's printed Lagrange multiplier (OK only): `s(na+1)·unbias`,
    /// which equals −μ of gslib-rs's variogram form.
    pub lagrange: Option<f64>,
    /// `x y z value weight` per datum, in GSLIB's order.
    pub data: Vec<[f64; 5]>,
    /// The estimate and variance.
    pub estimate: f64,
    /// See `estimate`.
    pub variance: f64,
}

fn numbers(text: &str) -> Result<Vec<f64>, String> {
    text.split_whitespace()
        .map(|t| t.parse().map_err(|e| format!("{t:?}: {e}")))
        .collect()
}

/// The `BLOCK:` sections of `kt3d.dbg`, in file order.
///
/// # Errors
/// A message naming the first malformed line.
pub fn parse_debug(text: &str) -> Result<Vec<Block>, String> {
    let mut blocks: Vec<Block> = Vec::new();
    let mut lines = text.lines().enumerate();
    while let Some((n, line)) = lines.next() {
        let Some(rest) = line.trim_start().strip_prefix("BLOCK:") else {
            continue;
        };
        let at = |e: String| format!("kt3d.dbg line {}: {e}", n + 1);
        let (ijk, xyz) = rest
            .split_once(" at ")
            .ok_or_else(|| at("no ' at '".into()))?;
        let ijk: Vec<u32> = ijk
            .split_whitespace()
            .map(str::parse)
            .collect::<Result<_, _>>()
            .map_err(|e| at(format!("node index: {e}")))?;
        let xyz = numbers(xyz).map_err(at)?;
        let (Ok(node), Ok(location)) = (
            <[u32; 3]>::try_from(ijk.as_slice()),
            <[f64; 3]>::try_from(xyz.as_slice()),
        ) else {
            return Err(at("expected 3 indices and 3 coordinates".into()));
        };
        let mut block = Block {
            node,
            location,
            lagrange: None,
            data: Vec::new(),
            estimate: f64::NAN,
            variance: f64::NAN,
        };
        for (m, content) in lines.by_ref() {
            let here = |e: String| format!("kt3d.dbg line {}: {e}", m + 1);
            let entry = content.trim();
            if let Some(v) = entry.strip_prefix("Lagrange :") {
                block.lagrange = Some(v.trim().parse().map_err(|e| here(format!("{e}")))?);
            } else if let Some(v) = entry.strip_prefix("estimate, variance") {
                let pair = numbers(v).map_err(here)?;
                let &[estimate, variance] = pair.as_slice() else {
                    return Err(here("expected estimate and variance".into()));
                };
                (block.estimate, block.variance) = (estimate, variance);
                break;
            } else if !(entry.starts_with("BLOCK EST") || entry.is_empty()) {
                let row = numbers(entry).map_err(here)?;
                let datum = <[f64; 5]>::try_from(row.as_slice())
                    .map_err(|_| here(format!("expected 5 values, got {entry:?}")))?;
                block.data.push(datum);
            }
        }
        if block.estimate.is_nan() {
            return Err(at("block has no estimate line".into()));
        }
        blocks.push(block);
    }
    Ok(blocks)
}

/// GSLIB's `getindx`: `int((loc − min)/siz + 1.5)` clamped to `1..=n`.
fn getindx(n: u32, min: f64, size: f64, loc: f64) -> u32 {
    let v = (loc - min) / size + 1.5;
    // The largest k in 1..=n with k <= v, else 1 (int() truncates; for
    // v >= 1 that is the floor).
    (1..=n).rev().find(|&k| f64::from(k) <= v).unwrap_or(1)
}

/// `kt3d`'s super-block search (`setsupr`, `picksup`, `srchsupr`).
///
/// A datum is examined only if its super-block's offset from the target's
/// is one `picksup` selected. Candidates are otherwise exactly gslib-rs's
/// (same metric, same order), so `DataSearch::select` over this
/// neighbourhood predicts GSLIB's selection, up to D-024 ties.
#[derive(Debug)]
pub struct GslibSuperBlocks<'a> {
    samples: &'a [Sample],
    axes: [Axis; 3],
    blocks: Vec<[u32; 3]>,
    offsets: Vec<[i64; 3]>,
}

/// One axis of the super-block grid.
#[derive(Debug, Clone, Copy)]
struct Axis {
    count: u32,
    origin: f64,
    size: f64,
}

impl Axis {
    /// `setsupr`: `nsup = min(n, 21)`, `sizsup = n·siz/nsup`,
    /// `mnsup = (mn − siz/2) + sizsup/2`.
    #[expect(
        clippy::suboptimal_flops,
        reason = "GSLIB is built with -ffp-contract=off: reproduce its unfused arithmetic"
    )]
    fn new(count: u32, origin: f64, spacing: f64) -> Self {
        let n = count.clamp(1, MAX_SUPER_BLOCKS);
        let size = f64::from(count) * spacing / f64::from(n);
        Self {
            count: n,
            origin: (origin - 0.5 * spacing) + 0.5 * size,
            size,
        }
    }

    fn index(self, loc: f64) -> u32 {
        getindx(self.count, self.origin, self.size, loc)
    }

    /// Offsets `-(n − 1)..=(n − 1)`.
    fn offsets(self) -> std::ops::RangeInclusive<i64> {
        let n = i64::from(self.count) - 1;
        -n..=n
    }
}

/// `picksup`'s smallest `sqdist` from the origin to the 3×3×3 lattice
/// `(a − b)·size/2 + offset`, with `a, b ∈ {−1, 1}` per axis.
#[expect(
    clippy::suboptimal_flops,
    reason = "GSLIB is built with -ffp-contract=off: reproduce its unfused arithmetic"
)]
fn lattice_min(ellipsoid: &SearchEllipsoid, axes: [Axis; 3], offset: [i64; 3]) -> f64 {
    let to_f = |v: i64| f64::from(i32::try_from(v).unwrap_or(i32::MAX));
    let [ax, ay, az] = axes;
    let [ox, oy, oz] = offset;
    let (xo, yo, zo) = (to_f(ox) * ax.size, to_f(oy) * ay.size, to_f(oz) * az.size);
    let ends = [-1.0, 1.0];
    let mut shortest = f64::INFINITY;
    for a1 in ends {
        for b1 in ends {
            for c1 in ends {
                for a2 in ends {
                    for b2 in ends {
                        for c2 in ends {
                            let d = Point3::new(
                                (a1 - a2) * 0.5 * ax.size + xo,
                                (b1 - b2) * 0.5 * ay.size + yo,
                                (c1 - c2) * 0.5 * az.size + zo,
                            );
                            shortest =
                                shortest.min(ellipsoid.squared_distance(Point3::default(), d));
                        }
                    }
                }
            }
        }
    }
    shortest
}

impl<'a> GslibSuperBlocks<'a> {
    /// The super-block grid `kt3d` builds from its estimation grid, and the
    /// offsets `picksup` selects for `ellipsoid`.
    #[must_use]
    pub fn new(p: &Params, samples: &'a [Sample], ellipsoid: &SearchEllipsoid) -> Self {
        let [nx, ny, nz] = p.grid_counts;
        let [x0, y0, z0] = p.grid_origin;
        let [dx, dy, dz] = p.grid_spacing;
        let axes = [
            Axis::new(nx, x0, dx),
            Axis::new(ny, y0, dy),
            Axis::new(nz, z0, dz),
        ];
        let radsqd = ellipsoid.radius() * ellipsoid.radius();
        let [axx, axy, axz] = axes;
        let mut offsets = Vec::new();
        for i in axx.offsets() {
            for j in axy.offsets() {
                for k in axz.offsets() {
                    if lattice_min(ellipsoid, axes, [i, j, k]) <= radsqd {
                        offsets.push([i, j, k]);
                    }
                }
            }
        }
        let mut grid = Self {
            samples,
            axes,
            blocks: Vec::new(),
            offsets,
        };
        grid.blocks = samples.iter().map(|s| grid.block(s.location())).collect();
        grid
    }

    fn block(&self, p: Point3) -> [u32; 3] {
        let [ax, ay, az] = self.axes;
        [ax.index(p.x), ay.index(p.y), az.index(p.z)]
    }

    /// Whether GSLIB examines datum `row` when estimating at `target`.
    #[must_use]
    pub fn examines(&self, row: usize, target: Point3) -> bool {
        let [ti, tj, tk] = self.block(target);
        self.blocks.get(row).is_some_and(|&[bi, bj, bk]| {
            let offset = [
                i64::from(bi) - i64::from(ti),
                i64::from(bj) - i64::from(tj),
                i64::from(bk) - i64::from(tk),
            ];
            self.offsets.contains(&offset)
        })
    }
}

impl Neighbourhood for GslibSuperBlocks<'_> {
    fn samples(&self) -> &[Sample] {
        self.samples
    }

    fn candidates(&self, ellipsoid: &SearchEllipsoid, target: Point3) -> Vec<Candidate> {
        BruteForce::new(self.samples)
            .candidates(ellipsoid, target)
            .into_iter()
            .filter(|c| self.examines(c.row, target))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getindx_matches_gslib() {
        // GETINDX.FOR: int((loc - min)/siz + 1.5), clamped.
        assert_eq!(getindx(5, 0.0, 10.0, 0.0), 1);
        assert_eq!(getindx(5, 0.0, 10.0, 4.9), 1);
        assert_eq!(getindx(5, 0.0, 10.0, 5.0), 2);
        assert_eq!(getindx(5, 0.0, 10.0, 1e9), 5);
        assert_eq!(getindx(5, 0.0, 10.0, -1e9), 1);
    }

    #[test]
    fn gslib_angles_reproduce_setrot_radians() {
        let a = gslib_angles([30.0, -20.0, 10.0]).unwrap();
        let alpha = (90.0 - a.azimuth).to_radians();
        assert!(60.0_f64.mul_add(-SETROT_DEG2RAD, alpha).abs() < 1e-15);
        assert!(
            20.0_f64
                .mul_add(-SETROT_DEG2RAD, (-a.dip).to_radians())
                .abs()
                < 1e-15
        );
        assert!(
            10.0_f64
                .mul_add(-SETROT_DEG2RAD, a.plunge.to_radians())
                .abs()
                < 1e-15
        );
        assert!(gslib_angles([270.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn parses_output_and_debug() {
        let out = "KT3D ESTIMATES WITH: t\n2\nEstimate\nEstimationVariance\n  1.5E+000  2.0E-001\n -9.99E+002 -9.99E+002\n";
        assert_eq!(
            parse_output(out).unwrap(),
            vec![(1.5, 0.2), (-999.0, -999.0)]
        );
        let dbg = "  \n Block Covariance:    1.0\n  \n        \n BLOCK:            2           1           1  at    40.0        20.0        0.0     \n        \n   Lagrange :   -2.5E-002\n   BLOCK EST: x,y,z,vr,wt \n  1.0E+001  2.0E+001  0.0E+000  3.0E+000  6.0E-001\n  4.0E+001  5.0E+000  0.0E+000  1.0E+000  4.0E-001\n   estimate, variance     2.2       0.41     \n";
        let blocks = parse_debug(dbg).unwrap();
        assert_eq!(blocks.len(), 1);
        let b = &blocks[0];
        assert_eq!(b.node, [2, 1, 1]);
        assert_eq!(b.lagrange, Some(-0.025));
        assert_eq!(b.data.len(), 2);
        assert_eq!((b.estimate, b.variance), (2.2, 0.41));
    }
}
