//! Synthetic input data for reference cases (`docs/spec/harness.md`, R4).
//!
//! Each suite is a pure function of a fixed seed, so `cargo xtask gen-data`
//! always writes byte-identical files; `--check` verifies the committed copies
//! in `tests/fixtures/data/`. Every suite says what it targets, and the tests
//! below check that it really contains those features.
//!
//! Files are GSLIB (Geo-EAS) format with columns `X Y Z V` (plus `W` where
//! noted). Coordinates are written with 4 decimals and values with 6, so the
//! text, not the generator's `f64`s, is the data: GSLIB and gslib-rs both read
//! the same decimal strings.

use std::fmt::Write as _;

use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rand_distr::{Distribution, StandardNormal};

/// One generated data set.
pub(crate) struct Suite {
    /// File stem: the data is written to `<name>.dat`.
    pub(crate) name: &'static str,
    /// Seed for [`ChaCha8Rng::seed_from_u64`].
    pub(crate) seed: u64,
    /// What the suite is for (copied into the data README).
    pub(crate) note: &'static str,
    /// Produces the rows `[x, y, z, v, (w)]`.
    build: fn(&mut ChaCha8Rng) -> Rows,
}

/// Column names and rows of a data set.
pub(crate) struct Rows {
    names: &'static [&'static str],
    rows: Vec<Vec<f64>>,
}

impl Suite {
    /// The rows this suite generates.
    pub(crate) fn rows(&self) -> Rows {
        (self.build)(&mut ChaCha8Rng::seed_from_u64(self.seed))
    }

    /// The suite as GSLIB file text.
    pub(crate) fn text(&self) -> String {
        let data = self.rows();
        let mut out = String::new();
        // Writing to a String cannot fail.
        let _ = writeln!(
            out,
            "gslib-rs synthetic suite '{}' (seed {})",
            self.name, self.seed
        );
        let _ = writeln!(out, "{}", data.names.len());
        for name in data.names {
            let _ = writeln!(out, "{name}");
        }
        for row in &data.rows {
            let fields: Vec<String> = row
                .iter()
                .enumerate()
                .map(|(j, v)| {
                    if j < 3 {
                        format!("{v:.4}")
                    } else {
                        format!("{v:.6}")
                    }
                })
                .collect();
            let _ = writeln!(out, "{}", fields.join(" "));
        }
        out
    }
}

const XYZV: &[&str] = &["X", "Y", "Z", "V"];
const XYZVW: &[&str] = &["X", "Y", "Z", "V", "W"];

fn uniform(rng: &mut ChaCha8Rng, lo: f64, hi: f64) -> f64 {
    rng.random::<f64>().mul_add(hi - lo, lo)
}

fn normal(rng: &mut ChaCha8Rng) -> f64 {
    StandardNormal.sample(rng)
}

/// Round to `decimals` places, as the file will store it.
fn rounded(v: f64, decimals: i32) -> f64 {
    let scale = 10_f64.powi(decimals);
    (v * scale).round() / scale
}

/// Rows of `[x, y, 0, v]` at random locations in a square.
fn scattered(
    rng: &mut ChaCha8Rng,
    n: u32,
    side: f64,
    value: impl Fn(&mut ChaCha8Rng, f64, f64) -> f64,
) -> Vec<Vec<f64>> {
    (0..n)
        .map(|_| {
            let x = rounded(uniform(rng, 0.0, side), 4);
            let y = rounded(uniform(rng, 0.0, side), 4);
            let v = rounded(value(rng, x, y), 6);
            vec![x, y, 0.0, v]
        })
        .collect()
}

fn smooth(rng: &mut ChaCha8Rng) -> Rows {
    let rows = scattered(rng, 60, 100.0, |r, x, y| {
        0.2_f64.mul_add(normal(r), (x / 20.0).sin() + (y / 25.0).cos() + 2.0)
    });
    Rows { names: XYZV, rows }
}

fn ties(rng: &mut ChaCha8Rng) -> Rows {
    // Values from a small set, so many are exactly equal.
    let levels = [1.0, 2.0, 2.5, 4.0];
    let rows = scattered(rng, 40, 50.0, |r, _, _| {
        levels
            .get(r.random_range(0..levels.len()))
            .copied()
            .unwrap_or(4.0)
    });
    Rows { names: XYZV, rows }
}

fn skewed(rng: &mut ChaCha8Rng) -> Rows {
    // Lognormal: exp(N(0, 1)), strongly right-skewed with a long upper tail.
    let rows = scattered(rng, 80, 100.0, |r, _, _| normal(r).exp());
    Rows { names: XYZV, rows }
}

fn weighted(rng: &mut ChaCha8Rng) -> Rows {
    let rows = scattered(rng, 50, 100.0, |r, x, _| x / 10.0 + normal(r))
        .into_iter()
        .map(|mut row| {
            row.push(rounded(uniform(rng, 0.2, 3.0), 6));
            row
        })
        .collect();
    Rows { names: XYZVW, rows }
}

fn duplicates(rng: &mut ChaCha8Rng) -> Rows {
    // 30 scattered points, then 5 of them repeated with different values.
    let mut rows = scattered(rng, 30, 40.0, |r, _, _| 10.0 + normal(r));
    let repeats: Vec<Vec<f64>> = rows
        .iter()
        .step_by(6)
        .map(|row| {
            let mut copy = row.clone();
            if let Some(v) = copy.get_mut(3) {
                *v = rounded(*v + 1.0 + normal(rng).abs(), 6);
            }
            copy
        })
        .collect();
    rows.extend(repeats);
    Rows { names: XYZV, rows }
}

fn grid_3d(rng: &mut ChaCha8Rng) -> Rows {
    // 6 x 5 x 4 regular grid (spacing 2, 2, 1), x fastest, with values that
    // vary faster vertically than horizontally.
    let mut rows = Vec::new();
    for iz in 0..4_u32 {
        for iy in 0..5_u32 {
            for ix in 0..6_u32 {
                let (x, y, z) = (2.0 * f64::from(ix), 2.0 * f64::from(iy), f64::from(iz));
                let v = 0.1_f64.mul_add(normal(rng), (x / 6.0).sin() + (z / 1.5).cos());
                rows.push(vec![x, y, z, rounded(v, 6)]);
            }
        }
    }
    Rows { names: XYZV, rows }
}

fn clustered(rng: &mut ChaCha8Rng) -> Rows {
    // 40 points in a dense cluster with high values, 20 spread out with low.
    let mut rows = Vec::new();
    for k in 0..60_u32 {
        let in_cluster = k < 40;
        let (lo, hi) = if in_cluster {
            (10.0, 25.0)
        } else {
            (0.0, 100.0)
        };
        let x = rounded(uniform(rng, lo, hi), 4);
        let y = rounded(uniform(rng, lo, hi), 4);
        let base = if in_cluster { 8.0 } else { 2.0 };
        rows.push(vec![x, y, 0.0, rounded(base + normal(rng), 6)]);
    }
    Rows { names: XYZV, rows }
}

fn constant(rng: &mut ChaCha8Rng) -> Rows {
    let rows = scattered(rng, 20, 30.0, |_, _, _| 3.5);
    Rows { names: XYZV, rows }
}

fn tiny(_: &mut ChaCha8Rng) -> Rows {
    Rows {
        names: XYZV,
        rows: vec![vec![0.0, 0.0, 0.0, 1.0], vec![3.0, 4.0, 0.0, 2.0]],
    }
}

fn with_sentinels(rng: &mut ChaCha8Rng) -> Rows {
    // Every fifth value is GSLIB's customary missing-value code, -999.
    let mut rows = scattered(rng, 40, 60.0, |r, _, _| 5.0 + normal(r));
    for row in rows.iter_mut().step_by(5) {
        if let Some(v) = row.get_mut(3) {
            *v = -999.0;
        }
    }
    Rows { names: XYZV, rows }
}

/// Every suite, in the order the README lists them.
pub(crate) const SUITES: [Suite; 10] = [
    Suite {
        name: "smooth",
        seed: 1,
        note: "60 scattered 2-D points with a smooth trend plus noise and no tied values: general-purpose input for nscore, gamv and kt3d.",
        build: smooth,
    },
    Suite {
        name: "ties",
        seed: 2,
        note: "40 points whose values take only 4 levels, so most are tied: nscore/backtr tie handling (GSLIB's random despiking, F-005).",
        build: ties,
    },
    Suite {
        name: "skewed",
        seed: 3,
        note: "80 points with lognormal values (strong positive skew, long upper tail): normal-score transform and back-transform tails.",
        build: skewed,
    },
    Suite {
        name: "weighted",
        seed: 4,
        note: "50 points with a positive weight column W: weighted nscore and declustering weights.",
        build: weighted,
    },
    Suite {
        name: "duplicates",
        seed: 5,
        note: "30 points plus 5 repeats at identical locations with different values: zero-distance pairs in gamv and singular systems in kt3d.",
        build: duplicates,
    },
    Suite {
        name: "grid_3d",
        seed: 6,
        note: "A 6 x 5 x 4 regular 3-D grid (x fastest), varying faster vertically: directional and vertical gamv lags.",
        build: grid_3d,
    },
    Suite {
        name: "clustered",
        seed: 7,
        note: "40 high-valued points in a dense cluster plus 20 sparse low-valued points: preferential sampling for declus.",
        build: clustered,
    },
    Suite {
        name: "constant",
        seed: 8,
        note: "20 points with one identical value (zero variance): degenerate input for every program.",
        build: constant,
    },
    Suite {
        name: "tiny",
        seed: 9,
        note: "Two points 5 apart: the smallest non-trivial input.",
        build: tiny,
    },
    Suite {
        name: "with_sentinels",
        seed: 10,
        note: "40 points where every fifth value is -999: trimming limits and missing values (semantics pending the missing-values spec).",
        build: with_sentinels,
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    fn suite(name: &str) -> &'static Suite {
        SUITES.iter().find(|s| s.name == name).unwrap()
    }

    fn column(name: &str, j: usize) -> Vec<f64> {
        suite(name).rows().rows.iter().map(|r| r[j]).collect()
    }

    fn distinct(values: &[f64]) -> usize {
        let mut v = values.to_vec();
        v.sort_by(f64::total_cmp);
        v.dedup();
        v.len()
    }

    #[test]
    fn generation_is_deterministic() {
        for s in &SUITES {
            assert_eq!(s.text(), s.text(), "{} differs between runs", s.name);
        }
    }

    #[test]
    fn names_and_seeds_are_unique() {
        let names: std::collections::BTreeSet<_> = SUITES.iter().map(|s| s.name).collect();
        let seeds: std::collections::BTreeSet<_> = SUITES.iter().map(|s| s.seed).collect();
        assert_eq!(names.len(), SUITES.len(), "duplicate suite name");
        assert_eq!(seeds.len(), SUITES.len(), "duplicate seed");
    }

    #[test]
    fn written_text_parses_back_to_the_same_rows() {
        for s in &SUITES {
            let parsed = gslib_rs::io::parse_gslib(&s.text()).unwrap();
            let rows = s.rows();
            assert_eq!(parsed.nrows(), rows.rows.len(), "{}", s.name);
            for (a, b) in parsed.rows().zip(&rows.rows) {
                assert_eq!(
                    a,
                    b.as_slice(),
                    "{}: rounding must make text and rows agree",
                    s.name
                );
            }
        }
    }

    #[test]
    fn smooth_has_no_ties() {
        let v = column("smooth", 3);
        assert_eq!(distinct(&v), v.len(), "smooth must have no tied values");
    }

    #[test]
    fn ties_are_mostly_tied() {
        let v = column("ties", 3);
        assert!(distinct(&v) <= 4, "ties has at most 4 distinct values");
        assert!(distinct(&v) >= 3, "but enough levels to rank");
    }

    #[test]
    fn skewed_is_positively_skewed() {
        let v = column("skewed", 3);
        let n = f64::from(u32::try_from(v.len()).unwrap());
        let mean = v.iter().sum::<f64>() / n;
        let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let skew = v.iter().map(|x| (x - mean).powi(3)).sum::<f64>() / n / var.powf(1.5);
        assert!(skew > 1.0, "skewness {skew}");
        assert!(v.iter().all(|&x| x > 0.0), "lognormal values are positive");
    }

    #[test]
    fn weights_are_positive() {
        assert!(
            column("weighted", 4).iter().all(|&w| w > 0.0),
            "weights must be positive"
        );
    }

    #[test]
    fn duplicates_share_locations_with_different_values() {
        let rows = suite("duplicates").rows().rows;
        let mut pairs = 0;
        for (i, a) in rows.iter().enumerate() {
            for b in rows.iter().skip(i + 1) {
                if a[..3] == b[..3] {
                    assert!(
                        (a[3] - b[3]).abs() > 0.5,
                        "coincident points must differ in value"
                    );
                    pairs += 1;
                }
            }
        }
        assert_eq!(pairs, 5);
    }

    #[test]
    fn grid_is_regular_and_x_fastest() {
        let rows = suite("grid_3d").rows().rows;
        assert_eq!(rows.len(), 120);
        assert_eq!(rows[1][..3], [2.0, 0.0, 0.0]);
        assert_eq!(rows[6][..3], [0.0, 2.0, 0.0]);
        assert_eq!(rows[30][..3], [0.0, 0.0, 1.0]);
    }

    #[test]
    fn clustered_is_preferentially_sampled() {
        let rows = suite("clustered").rows().rows;
        let in_box = |r: &Vec<f64>| (10.0..=25.0).contains(&r[0]) && (10.0..=25.0).contains(&r[1]);
        assert!(
            rows.iter().filter(|r| in_box(r)).count() >= 40,
            "dense cluster present"
        );
        let mean = |f: &dyn Fn(&&Vec<f64>) -> bool| {
            let v: Vec<f64> = rows.iter().filter(f).map(|r| r[3]).collect();
            v.iter().sum::<f64>() / f64::from(u32::try_from(v.len()).unwrap())
        };
        assert!(
            mean(&|r| in_box(r)) > mean(&|r| !in_box(r)) + 2.0,
            "cluster has high values"
        );
    }

    #[test]
    fn constant_has_zero_variance() {
        assert_eq!(distinct(&column("constant", 3)), 1);
    }

    #[test]
    fn sentinels_appear_every_fifth_value() {
        let v = column("with_sentinels", 3);
        let is_sentinel = |x: f64| (x + 999.0).abs() < 1e-9;
        assert_eq!(v.iter().filter(|&&x| is_sentinel(x)).count(), 8);
        assert!(
            v.iter().filter(|&&x| !is_sentinel(x)).all(|&x| x > -100.0),
            "no accidental sentinels"
        );
    }
}
