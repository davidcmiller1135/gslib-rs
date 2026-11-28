//! gslib-rs's search and point kriging against GSLIB `kt3d` fixtures
//! (`docs/spec/search.md` step 3, harness step 7; F-012, F-014, D-027).
//!
//! Per grid node:
//! 1. **Search.** GSLIB's selection, read from `kt3d.dbg`, must equal
//!    gslib-rs's `DataSearch` over [`GslibSuperBlocks`] (gslib-rs's metric
//!    and rules, restricted to the super-blocks `picksup` lets GSLIB
//!    examine). Where that differs from gslib-rs's own (brute-force)
//!    selection, GSLIB skipped data inside its search ellipsoid; the case
//!    must declare the `picksup` deviation. Tied candidates at a cut are
//!    compared as a valid choice (D-024) when the case declares
//!    `search_ties`.
//! 2. **Kriging.** gslib-rs kriges GSLIB's own selection, in GSLIB's order,
//!    and must reproduce GSLIB's estimate, variance, weights and Lagrange
//!    multiplier within `K · cond(A) · ε` (D-010), with `cond(A)` the larger
//!    of the two systems' condition numbers (GSLIB's covariance form,
//!    gslib-rs's).
//! 3. **Special cases.** Too few data: GSLIB's UNEST where gslib-rs's
//!    search fails. One datum: GSLIB's shortcut (F-014), checked as an exact
//!    relation.

#[cfg(test)]
mod kt3d {
    use std::num::NonZeroUsize;

    use gslib_rs::geometry::Point3;
    use gslib_rs::io::{Extraction, Trimming, parse_gslib};
    use gslib_rs::kriging::{Solution, ordinary_kriging_solution, simple_kriging_solution};
    use gslib_rs::model::{CovarianceModel, VariogramModel};
    use gslib_rs::sample::Sample;
    use gslib_rs::search::{BruteForce, DataSearch, SearchEllipsoid};
    use harness::compare::{Comparison, DIRECT, Tolerance};
    use harness::kt3d::{
        Block, GslibSuperBlocks, ellipsoid, model, node, parse_debug, parse_output,
    };
    use nalgebra::DMatrix;
    use oracle_case::Fixture;
    use oracle_case::kt3d::{DEBUG_FILE, Kriging, OUTPUT_FILE, Params};

    /// GSLIB's unestimated code (`KT3D.INC`).
    const UNEST: f64 = -999.0;
    /// `k` in the solver bound `k · cond(A) · ε` (D-010), calibrated on these
    /// fixtures: the worst observed `|Δ| / (cond·ε·scale)` is 10.4, so 100
    /// leaves a factor of about 10 (D-035).
    const K: f64 = 100.0;

    /// What one fixture's comparison needs.
    struct Setup {
        p: Params,
        samples: Vec<Sample>,
        variogram: VariogramModel,
        covariance: CovarianceModel,
        search: DataSearch,
        ellipsoid: SearchEllipsoid,
    }

    fn setup(f: &Fixture) -> Setup {
        let p = Params::from_table(&f.case.parameters).unwrap();
        let data = parse_gslib(&std::fs::read_to_string(f.dir.join(&p.data)).unwrap()).unwrap();
        let trimming = p
            .trimming_limits
            .map_or(Trimming::NONE, |[lo, hi]| Trimming::limits(lo, hi).unwrap());
        let values = data
            .extract(
                usize::try_from(p.variable_column - 1).unwrap(),
                &Extraction::new(trimming),
            )
            .unwrap();
        // A coordinate column of 0 means the grid origin in kt3d (F-014).
        let coordinate = |axis: usize| -> Vec<f64> {
            match p.coordinate_columns[axis] {
                0 => vec![p.grid_origin[axis]; data.nrows()],
                c => data.column(usize::try_from(c - 1).unwrap()).unwrap(),
            }
        };
        let (xs, ys, zs) = (coordinate(0), coordinate(1), coordinate(2));
        let samples = (0..data.nrows())
            .filter_map(|r| {
                values[r].map(|v| Sample::new(Point3::new(xs[r], ys[r], zs[r]), v).unwrap())
            })
            .collect();
        let variogram = model(p.nugget, &p.structures).unwrap();
        let covariance = CovarianceModel::try_from(variogram.clone()).unwrap();
        let ellipsoid = ellipsoid(&p).unwrap();
        let mut search = DataSearch::new(
            ellipsoid,
            usize::try_from(p.min_data).unwrap(),
            usize::try_from(p.max_data).unwrap(),
        )
        .unwrap();
        if let Some(n) = NonZeroUsize::new(usize::try_from(p.per_octant).unwrap()) {
            search = search.with_octants(n).unwrap();
        }
        Setup {
            p,
            samples,
            variogram,
            covariance,
            search,
            ellipsoid,
        }
    }

    /// The rows of GSLIB's selected data, matched by value (exact) and
    /// location (GSLIB's debug coordinates are `xa + xloc − xsiz/2`, so
    /// within rounding of the data file's).
    fn gslib_rows(block: &Block, samples: &[Sample]) -> Result<Vec<usize>, String> {
        let mut used = vec![false; samples.len()];
        block
            .data
            .iter()
            .map(|&[x, y, z, v, _]| {
                let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * (1.0 + b.abs());
                let row = (0..samples.len()).find(|&r| {
                    let s = samples[r];
                    let l = s.location();
                    !used[r]
                        && s.value().to_bits() == v.to_bits()
                        && close(x, l.x)
                        && close(y, l.y)
                        && close(z, l.z)
                });
                row.inspect(|&r| used[r] = true)
                    .ok_or_else(|| format!("no datum matches GSLIB's ({x}, {y}, {z}) = {v}"))
            })
            .collect()
    }

    /// Whether GSLIB's selection `g` is a valid choice for gslib-rs's `ours`
    /// under D-024: the same number of data at each distance, the same data
    /// strictly inside the farthest selected distance, and in order.
    fn valid_tie_choice(g: &[usize], ours: &[usize], q: &dyn Fn(usize) -> f64) -> bool {
        let bits = |rows: &[usize]| {
            let mut b: Vec<u64> = rows.iter().map(|&r| q(r).to_bits()).collect();
            b.sort_unstable();
            b
        };
        let ordered = g.windows(2).all(|w| q(w[0]) <= q(w[1]));
        let cut = ours.iter().map(|&r| q(r)).fold(0.0, f64::max);
        let inner = |rows: &[usize]| {
            let mut v: Vec<usize> = rows.iter().copied().filter(|&r| q(r) < cut).collect();
            v.sort_unstable();
            v
        };
        g.len() == ours.len() && ordered && bits(g) == bits(ours) && inner(g) == inner(ours)
    }

    /// The 2-norm condition number of a square matrix.
    fn condition(a: &DMatrix<f64>) -> f64 {
        let s = a.singular_values();
        let max = s.iter().copied().fold(0.0, f64::max);
        let min = s.iter().copied().fold(f64::INFINITY, f64::min);
        max / min
    }

    /// cond(A) for GSLIB's system (covariances, OK bordered by `unbias =
    /// C(0)`) and for gslib-rs's (OK in variograms bordered by 1; SK in
    /// covariances), the larger of the two.
    fn conditioning(s: &Setup, data: &[Sample]) -> f64 {
        let n = data.len();
        let ok = s.p.kriging == Kriging::Ordinary;
        let size = if ok { n + 1 } else { n };
        let c0 = s.covariance.sill();
        let border = |i: usize, j: usize, unit: f64| match (data.get(i), data.get(j)) {
            (Some(_), Some(_)) => None,
            (None, None) => Some(0.0),
            _ => Some(unit),
        };
        let gslib = DMatrix::from_fn(size, size, |i, j| {
            border(i, j, c0).unwrap_or_else(|| {
                s.covariance
                    .covariance(data[i].location(), data[j].location())
            })
        });
        let ours = DMatrix::from_fn(size, size, |i, j| {
            border(i, j, 1.0).unwrap_or_else(|| {
                if ok {
                    s.variogram.gamma(data[i].location(), data[j].location())
                } else {
                    s.covariance
                        .covariance(data[i].location(), data[j].location())
                }
            })
        });
        condition(&gslib).max(condition(&ours))
    }

    /// The largest `|Δ| / (cond · ε · scale)` seen per quantity (calibration).
    #[derive(Default)]
    struct Calibration(Vec<(&'static str, f64)>);

    impl Calibration {
        fn note(&mut self, q: &'static str, ratio: f64) {
            match self.0.iter_mut().find(|(n, _)| *n == q) {
                Some((_, r)) => *r = r.max(ratio),
                None => self.0.push((q, ratio)),
            }
        }
    }

    fn solve(s: &Setup, data: &[Sample], t: Point3) -> Solution {
        match s.p.kriging {
            Kriging::Ordinary => ordinary_kriging_solution(data, t, &s.variogram),
            Kriging::Simple => simple_kriging_solution(data, t, &s.covariance, s.p.sk_mean),
        }
        .unwrap()
    }

    /// Compare a solved system: GSLIB's block against gslib-rs on GSLIB's
    /// selection.
    fn compare_block(
        cmp: &mut Comparison<'_>,
        cal: &mut Calibration,
        s: &Setup,
        block: &Block,
        rows: &[usize],
        t: Point3,
        at: &dyn Fn() -> String,
    ) {
        let data: Vec<Sample> = rows.iter().map(|&r| s.samples[r]).collect();
        let ours = solve(s, &data, t);
        let cond = conditioning(s, &data);
        cal.note("cond(A)", cond);
        let bound = cond * f64::EPSILON;
        let tolerance = Tolerance {
            relative: K * bound,
            absolute: 0.0,
        };
        let value_scale = data
            .iter()
            .map(|d| d.value().abs())
            .fold(0.0, f64::max)
            .max(1e-300);
        let sill = s.covariance.sill();
        let mut check = |q: &'static str, o: f64, g: f64, scale: f64| {
            cal.note(q, (o - g).abs() / (bound * scale));
            cmp.value(q, o, g, scale, tolerance, at);
        };
        check("estimate", ours.estimate.value, block.estimate, value_scale);
        // gslib-rs clamps a negative variance to 0; GSLIB does not.
        check(
            "variance",
            ours.estimate.variance,
            block.variance.max(0.0),
            sill,
        );
        for (w, d) in ours.weights.iter().zip(&block.data) {
            check("weight", *w, d[4], 1.0);
        }
        if let (Some(mu), Some(l)) = (ours.lagrange, block.lagrange) {
            // GSLIB prints s(na+1)·unbias, which is −μ of the variogram form.
            check("lagrange", -mu, l, sill);
        }
        cmp.structure(ours.lagrange.is_some() == block.lagrange.is_some(), || {
            format!("{}: Lagrange multiplier present in one only", at())
        });
    }

    /// Compare one fixture; returns its mismatches.
    #[expect(
        clippy::too_many_lines,
        reason = "one pass over the nodes, case by case"
    )]
    fn check(fixture: &Fixture) -> Vec<String> {
        let name = fixture
            .dir
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let s = setup(fixture);
        let read = |file: &str| std::fs::read_to_string(fixture.dir.join(file)).unwrap();
        let out = parse_output(&read(OUTPUT_FILE)).unwrap();
        let blocks = parse_debug(&read(DEBUG_FILE)).unwrap();
        let gslib_blocks = GslibSuperBlocks::new(&s.p, &s.samples, &s.ellipsoid);
        let brute = BruteForce::new(&s.samples);
        let mut cmp = Comparison::new(&fixture.case);
        let mut cal = Calibration::default();
        let nodes =
            s.p.grid_counts
                .iter()
                .map(|&n| usize::try_from(n).unwrap())
                .product::<usize>();
        cmp.structure(out.len() == nodes, || {
            format!("{} output rows, expected {nodes}", out.len())
        });
        let (mut misses, mut solved, mut unestimated, mut single, mut ties) = (0, 0, 0, 0, 0);
        for (index, &(est, var)) in out.iter().enumerate().take(nodes) {
            let (ijk, t) = node(&s.p, index);
            let at = || format!("node {ijk:?}");
            let block = blocks.iter().find(|b| b.node == ijk);
            let predicted = s.search.select(&gslib_blocks, t);
            let ours = s.search.select(&brute, t);
            if predicted.as_ref().ok() != ours.as_ref().ok() {
                if misses == 0 {
                    let skipped: Vec<String> = ours
                        .iter()
                        .flatten()
                        .filter(|r| !predicted.as_ref().is_ok_and(|p| p.contains(r)))
                        .map(|&r| {
                            let d = s.samples[r];
                            format!(
                                "row {r} at {:?}, {:.3} of the radius",
                                d.location(),
                                s.ellipsoid.squared_distance(d.location(), t).sqrt()
                                    / s.ellipsoid.radius()
                            )
                        })
                        .collect();
                    println!(
                        "{name}: first skip at node {ijk:?} {t:?}: GSLIB wrote ({est}, {var}); \
                         skipped {skipped:?}"
                    );
                }
                misses += 1;
                cmp.structure(cmp.is_expected_deviation("picksup"), || {
                    format!(
                        "{}: GSLIB's super-blocks skip data gslib-rs selects, undeclared",
                        at()
                    )
                });
            }
            match predicted {
                Err(_) => {
                    unestimated += 1;
                    cmp.structure(
                        est.to_bits() == UNEST.to_bits() && var.to_bits() == UNEST.to_bits(),
                        || format!("{}: too few data, but GSLIB wrote {est}, {var}", at()),
                    );
                    cmp.structure(block.is_none(), || {
                        format!("{}: unexpected debug block", at())
                    });
                }
                Ok(rows) if rows.len() == 1 && s.p.kriging == Kriging::Ordinary => {
                    // F-014: kt3d counts OK's unbiasedness as a drift term
                    // and leaves 1 <= na <= mdt = 1 unestimated.
                    single += 1;
                    cmp.structure(cmp.is_expected_deviation("ok_single_datum"), || {
                        format!("{}: OK with one datum, undeclared", at())
                    });
                    cmp.structure(
                        est.to_bits() == UNEST.to_bits() && var.to_bits() == UNEST.to_bits(),
                        || format!("{}: OK with one datum, but GSLIB wrote {est}, {var}", at()),
                    );
                    cmp.structure(block.is_none(), || {
                        format!("{}: unexpected debug block", at())
                    });
                }
                Ok(rows) if rows.len() == 1 => {
                    // F-014: kt3d returns the datum and 2(C(0) − C(h)).
                    single += 1;
                    let d = s.samples[rows[0]];
                    let c0 = s.covariance.sill();
                    let ch = s.covariance.covariance(d.location(), t);
                    cmp.value(
                        "single_datum_value",
                        d.value(),
                        est,
                        d.value().abs().max(1.0),
                        DIRECT,
                        at,
                    );
                    cmp.value(
                        "single_datum_variance",
                        2.0 * (c0 - ch),
                        var,
                        c0,
                        DIRECT,
                        at,
                    );
                    cmp.structure(block.is_none(), || {
                        format!("{}: unexpected debug block", at())
                    });
                    let o = solve(&s, &[d], t).estimate;
                    let quantity = if s.p.kriging == Kriging::Simple {
                        "sk_single_datum"
                    } else {
                        "estimate"
                    };
                    cmp.value(quantity, o.value, est, d.value().abs().max(1.0), DIRECT, at);
                }
                Ok(rows) => {
                    solved += 1;
                    let Some(block) = block else {
                        cmp.structure(false, || {
                            format!("{}: no debug block for a solved system", at())
                        });
                        continue;
                    };
                    cmp.value(
                        "node",
                        t.x,
                        block.location[0],
                        s.p.grid_spacing[0],
                        DIRECT,
                        at,
                    );
                    cmp.value(
                        "node",
                        t.y,
                        block.location[1],
                        s.p.grid_spacing[1],
                        DIRECT,
                        at,
                    );
                    cmp.value(
                        "node",
                        t.z,
                        block.location[2],
                        s.p.grid_spacing[2],
                        DIRECT,
                        at,
                    );
                    cmp.structure(est.to_bits() == block.estimate.to_bits(), || {
                        format!("{}: kt3d.out and kt3d.dbg disagree", at())
                    });
                    let selected = match gslib_rows(block, &s.samples) {
                        Ok(selected) => selected,
                        Err(e) => {
                            cmp.structure(false, || format!("{}: {e}", at()));
                            continue;
                        }
                    };
                    let distance =
                        |r: usize| s.ellipsoid.squared_distance(s.samples[r].location(), t);
                    if selected != rows {
                        ties += 1;
                        let ok = cmp.is_expected_deviation("search_ties")
                            && valid_tie_choice(&selected, &rows, &distance);
                        cmp.structure(ok, || {
                            format!(
                                "{}: GSLIB selected rows {selected:?}, predicted {rows:?}",
                                at()
                            )
                        });
                    }
                    compare_block(&mut cmp, &mut cal, &s, block, &selected, t, &at);
                }
            }
        }
        cmp.structure(blocks.len() == solved, || {
            format!("{} debug blocks, {solved} solved systems", blocks.len())
        });
        println!(
            "{name}: {solved} solved, {single} single-datum, {unestimated} unestimated, \
             {misses} nodes where GSLIB's super-blocks skip in-range data, \
             {ties} where GSLIB chose other tied data"
        );
        for (quantity, r) in &cal.0 {
            if *quantity == "cond(A)" {
                println!("{name}: largest cond(A) = {r:.3e}");
            } else {
                println!("{name}: worst {quantity} |Δ| / (cond·ε·scale) = {r:.3e}");
            }
        }
        cmp.failures()
            .iter()
            .map(|m| format!("{name}: {m}"))
            .collect()
    }

    #[test]
    fn point_kriging_matches_gslib() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "kt3d")
            .collect();
        assert!(!fixtures.is_empty(), "no kt3d fixtures found");
        let failures: Vec<String> = fixtures.iter().flat_map(check).collect();
        assert!(
            failures.is_empty(),
            "{} mismatches:\n{}",
            failures.len(),
            failures
                .iter()
                .take(40)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
