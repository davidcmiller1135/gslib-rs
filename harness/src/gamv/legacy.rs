//! A `gamv` case as GSLIB reads it, for gslib-rs (`docs/spec/gamv.md`, plan
//! step 6; D-048, D-050).
//!
//! This adapter reproduces what `GAMV.FOR`'s `readparm` does before the pair
//! loop, so that gslib-rs's kernel sees the same data:
//!
//! - **Rows** are kept iff some listed variable has `tmin ≤ v < tmax`
//!   (l. 257).
//! - **Pairing** accepts a value iff `tmin ≤ v ≤ tmax` (l. 623): a datum at
//!   `tmax` is used.
//! - **Sills** are population variances over `tmin ≤ v < tmax` (ll. 293–303);
//!   a variable with no such value has none (GSLIB's −999).
//! - **Indicators** (types 9, 10) are built from the tail variable
//!   (ll. 308–334): a value with `v ≤ tmin` or `v > tmax` becomes
//!   `tmin − 1e-20`, which then still passes the pair test whenever that
//!   subtraction rounds back to `tmin` (F-025); otherwise `I = 0` if `v < c`
//!   else 1 (type 9), or `int(v + 0.5) = int(c + 0.5)` (type 10). The sill is
//!   `p0(1 − p0)` for type 9 and 0 for type 10.
//! - **Standardisation** (`isill = 1`) applies to types 1, 9 and 10 with
//!   tail = head and a positive sill (ll. 814–820).
//! - **Tolerances** `atol, dtol ≤ 0` become 45° (D-044).
//!
//! Self-pairs are not reproduced (D-017); [`Legacy::self_pairs`] gives the
//! count GSLIB adds to the zero-distance class.

use std::path::Path;

use gslib_rs::geometry::Point3;
use gslib_rs::io::parse_gslib;
use gslib_rs::variogram::{
    CentredLags, Data, Direction, GamvPolicies, GamvRequest, VariogramKind, VariogramSpec,
};
use oracle_case::gamv::Params;

/// `EPSLON` in `GAMV.INC`.
const EPSLON: f64 = 1.0e-20;

/// One variogram as GSLIB lists it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variogram {
    /// GSLIB's `ivtype`.
    pub gslib_type: i32,
    /// The names in GSLIB's title (12 characters each).
    pub tail_name: String,
    /// See `tail_name`.
    pub head_name: String,
}

/// A case, ready for gslib-rs.
#[derive(Debug, Clone)]
pub struct Legacy {
    /// The data, trimmed as GSLIB trims it.
    pub data: Data,
    /// The request, with every legacy policy.
    pub request: GamvRequest,
    /// Each variogram's type and names.
    pub variograms: Vec<Variogram>,
    /// For each variogram: the data with both its values present and in the
    /// measure's domain, which GSLIB pairs with themselves.
    pub self_pair_data: Vec<usize>,
    /// The sills, as GSLIB computed them (`None`: GSLIB's −999).
    pub sills: Vec<Option<f64>>,
    /// The data locations (a copy of what `data` holds).
    pub points: Vec<Point3>,
    /// The kernel's columns (a copy of what `data` holds).
    pub columns: Vec<Vec<Option<f64>>>,
    /// The upper trimming limit.
    pub tmax: f64,
}

impl Legacy {
    /// GSLIB's self-pairs in the zero-distance class for variogram `v` and
    /// a direction: one per datum, two if omnidirectional (except type 2).
    #[must_use]
    pub fn self_pairs(&self, v: usize, omnidirectional: bool) -> usize {
        let n = self.self_pair_data.get(v).copied().unwrap_or(0);
        let doubled = omnidirectional && self.variograms.get(v).is_some_and(|g| g.gslib_type != 2);
        if doubled { 2 * n } else { n }
    }
}

/// GSLIB's sill: `Σv²/n − (Σv/n)²`, unfused as in `GAMV.FOR` ll. 268–297.
#[expect(
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    reason = "GSLIB's unfused formula: mean of squares minus squared mean"
)]
fn population_variance(values: impl Iterator<Item = f64>) -> Option<f64> {
    let (n, sum, sq) = values.fold((0_u32, 0.0, 0.0), |(n, s, q), v| (n + 1, s + v, q + v * v));
    (n > 0).then(|| {
        let n = f64::from(n);
        let avg = sum / n;
        sq / n - avg * avg
    })
}

/// GSLIB's indicator column and sill from the raw values of its source
/// variable (`GAMV.FOR` ll. 308–334). `pairable` is the pair test.
fn indicator_column(
    raw: &[f64],
    gslib_type: i32,
    cut: f64,
    [tmin, tmax]: [f64; 2],
    pairable: impl Fn(f64) -> Option<f64>,
) -> (Vec<Option<f64>>, f64) {
    let marker = tmin - EPSLON;
    let (mut zeros, mut ones) = (0.0_f64, 0.0_f64);
    let column = raw
        .iter()
        .map(|&x| {
            if x <= tmin || x > tmax {
                // GSLIB stores tmin − EPSLON; the pair test rejects it only
                // if that is really below tmin (F-025).
                pairable(marker)
            } else if gslib_type == 9 {
                if x < cut {
                    zeros += 1.0;
                    Some(0.0)
                } else {
                    ones += 1.0;
                    Some(1.0)
                }
            } else {
                let same = (x + 0.5).trunc().to_bits() == (cut + 0.5).trunc().to_bits();
                Some(if same { 1.0 } else { 0.0 })
            }
        })
        .collect();
    // Type 10 never counts p0 or p1, so its sill is 0.
    let p0 = zeros / (zeros + ones).max(1.0);
    (column, p0 * (1.0 - p0))
}

/// The data GSLIB pairs with themselves for each variogram: both values
/// present and, for types 6 and 7, in the measure's domain.
fn self_pair_data(
    specs: &[VariogramSpec],
    variograms: &[Variogram],
    columns: &[Vec<Option<f64>>],
    n: usize,
) -> Vec<usize> {
    specs
        .iter()
        .zip(variograms)
        .map(|(s, g)| {
            let at = |c: usize, i: usize| columns.get(c).and_then(|c| c.get(i)).copied().flatten();
            (0..n)
                .filter(|&i| match (at(s.tail, i), at(s.head, i)) {
                    (Some(t), Some(h)) => match g.gslib_type {
                        6 => (t + h).abs() > EPSLON,
                        7 => t > EPSLON && h > EPSLON,
                        _ => true,
                    },
                    _ => false,
                })
                .count()
        })
        .collect()
}

fn locations(xs: &[f64], ys: &[f64], zs: &[f64]) -> Vec<Point3> {
    xs.iter()
        .zip(ys)
        .zip(zs)
        .map(|((&x, &y), &z)| Point3::new(x, y, z))
        .collect()
}

/// The directions, with GSLIB's 45° for a non-positive tolerance (D-044).
fn directions(p: &Params) -> Result<Vec<Direction>, String> {
    p.directions
        .clone()
        .unwrap_or_else(|| vec![[0.0, 90.0, 1.0e21, 0.0, 90.0, 1.0e21]])
        .iter()
        .map(|&[azm, atol, bh, dip, dtol, bv]| {
            let legacy = |t: f64| if t <= 0.0 { 45.0 } else { t };
            Direction::new(azm, legacy(atol), bh, dip, legacy(dtol), bv).map_err(|e| e.to_string())
        })
        .collect()
}

fn twelve(name: &str) -> String {
    name.chars().take(12).collect()
}

fn kind(t: i32) -> Result<VariogramKind, String> {
    Ok(match t {
        1 | 9 | 10 => VariogramKind::Semivariogram,
        2 => VariogramKind::CrossSemivariogram,
        3 => VariogramKind::Covariance,
        -3 => VariogramKind::CovarianceAsVariogram,
        4 => VariogramKind::Correlogram,
        5 => VariogramKind::GeneralRelative,
        6 => VariogramKind::PairwiseRelative,
        7 => VariogramKind::Logarithmic,
        8 => VariogramKind::Madogram,
        other => return Err(format!("variogram type {other}")),
    })
}

/// Build gslib-rs's inputs for a case in `dir`.
///
/// # Errors
/// A message if the data file cannot be read or does not fit the case.
pub fn build(dir: &Path, p: &Params) -> Result<Legacy, String> {
    let text = std::fs::read_to_string(dir.join(&p.data)).map_err(|e| e.to_string())?;
    let file = parse_gslib(&text).map_err(|e| e.to_string())?;
    let [tmin, tmax] = p.trimming_limits.unwrap_or([-1.0e21, 1.0e21]);
    let raw = |c: u32| -> Result<Vec<f64>, String> {
        if c == 0 {
            return Ok(vec![0.0; file.nrows()]);
        }
        let index = usize::try_from(c - 1).map_err(|e| e.to_string())?;
        file.column(index).map_err(|e| e.to_string())
    };
    let columns = p.columns()?;
    let vars: Vec<Vec<f64>> = columns.iter().map(|&c| raw(c)).collect::<Result<_, _>>()?;
    let in_range = |v: f64| v >= tmin && v < tmax;
    let keep: Vec<bool> = (0..file.nrows())
        .map(|row| {
            vars.iter()
                .any(|col| col.get(row).is_some_and(|&v| in_range(v)))
        })
        .collect();
    let kept = |col: &[f64]| -> Vec<f64> {
        col.iter()
            .zip(&keep)
            .filter_map(|(&v, &k)| k.then_some(v))
            .collect()
    };
    let [cx, cy, cz] = p.coordinate_columns;
    let points = locations(&kept(&raw(cx)?), &kept(&raw(cy)?), &kept(&raw(cz)?));
    let pairable = |v: f64| (v >= tmin && v <= tmax).then_some(v);
    let names = file.names();
    let mut kernel_columns: Vec<Vec<Option<f64>>> = Vec::new();
    let mut sills: Vec<Option<f64>> = Vec::new();
    let mut var_names: Vec<String> = Vec::new();
    for (&c, col) in columns.iter().zip(&vars) {
        let col = kept(col);
        sills.push(population_variance(
            col.iter().copied().filter(|&v| in_range(v)),
        ));
        kernel_columns.push(col.iter().map(|&v| pairable(v)).collect());
        let index = usize::try_from(c - 1).map_err(|e| e.to_string())?;
        var_names.push(names.get(index).map_or_else(String::new, |n| twelve(n)));
    }

    let mut specs = Vec::new();
    let mut variograms = Vec::new();
    let mut cuts = 0;
    for v in p.variogram_list() {
        let index = |i: u32| usize::try_from(i - 1).map_err(|e| e.to_string());
        let (tail, head) = if let Some(cut) = v.cutoff {
            cuts += 1;
            let raw_source = kept(vars.get(index(v.tail)?).ok_or("indicator source")?);
            let (column, sill) = indicator_column(&raw_source, v.kind, cut, [tmin, tmax], pairable);
            sills.push(Some(sill));
            kernel_columns.push(column);
            var_names.push(format!("Indicator {cuts:>2}"));
            (kernel_columns.len() - 1, kernel_columns.len() - 1)
        } else {
            (index(v.tail)?, index(v.head)?)
        };
        let gslib_kind = kind(v.kind)?;
        let standardize = p.standardize
            && tail == head
            && matches!(v.kind, 1 | 9 | 10)
            && sills.get(tail).copied().flatten().is_some_and(|s| s > 0.0);
        specs.push(VariogramSpec {
            tail,
            head,
            kind: gslib_kind,
            standardize,
        });
        variograms.push(Variogram {
            gslib_type: v.kind,
            tail_name: var_names.get(tail).cloned().unwrap_or_default(),
            head_name: var_names.get(head).cloned().unwrap_or_default(),
        });
    }

    let self_pair_data = self_pair_data(&specs, &variograms, &kernel_columns, points.len());
    let directions = directions(p)?;
    let nlag = usize::try_from(p.nlag).map_err(|e| e.to_string())?;
    let request = GamvRequest {
        lags: CentredLags::new(nlag, p.lag_spacing, p.lag_tolerance).map_err(|e| e.to_string())?,
        directions,
        variograms: specs,
        sills: sills.clone(),
        policies: GamvPolicies::gslib(),
    };
    let data = Data::new(points.clone(), kernel_columns.clone()).map_err(|e| e.to_string())?;
    Ok(Legacy {
        data,
        request,
        variograms,
        self_pair_data,
        sills,
        points,
        columns: kernel_columns,
        tmax,
    })
}
