//! `sgsim` statistical parity (`docs/spec/harness.md` step 8, R8, D-011).
//!
//! For each statistical case, GSLIB's N realizations (from the fixture) and
//! N realizations from gslib-rs (stream k of the case's seed, D-037) are
//! reduced to the same per-realization statistics:
//! - the experimental semivariogram along x and along y at lags 1–6 (pairs
//!   of grid nodes h apart in that direction);
//! - the proportion of nodes below each decile of N(0, 1).
//!
//! For a stationary Gaussian field both are unbiased for the model: `E[γ̂(h)]
//! = γ(h)` and `E[proportion below z] = Φ(z)`. (Sample quantiles would be
//! biased by each realization's own mean, which is why "per quantile" is
//! taken at the model's deciles.)
//!
//! Checks, each within 5 standard errors (D-011):
//! - GSLIB's ensemble mean against gslib-rs's, with the combined standard
//!   error `sqrt(s²_gslib/N + s²_ours/N)`;
//! - unconditional cases only: gslib-rs's ensemble mean against the model
//!   (R8's independent check; GSLIB's distance from the model is reported);
//! - conditional cases: every realization of both reproduces every datum
//!   at its node, exactly.

#[cfg(test)]
mod sgsim_stats {
    use gslib_rs::geometry::Point3;
    use gslib_rs::io::parse_gslib;
    use gslib_rs::model::CovarianceModel;
    use gslib_rs::random::Seed;
    use gslib_rs::sample::Sample;
    use gslib_rs::search::{DataSearch, NodeSearch, SearchEllipsoid, TableExtent};
    use gslib_rs::simulation::{Axis, DataConditioning, Grid, Searches, sgsim};
    use gslib_rs::transform::standard_normal_quantile;
    use harness::compare::{Comparison, DIRECT};
    use oracle_case::Fixture;
    use oracle_case::sgsim::{OUTPUT_FILE, Params};

    /// The comparison bound in standard errors (D-011).
    const SIGMAS: f64 = 5.0;
    /// Lags compared along each axis.
    const LAGS: usize = 6;

    fn grid(p: &Params) -> Grid {
        let axis = |a: usize| Axis {
            count: usize::try_from(p.grid_counts[a]).unwrap(),
            first_center: p.grid_origin[a],
            spacing: p.grid_spacing[a],
        };
        Grid::new(axis(0), axis(1), axis(2)).unwrap()
    }

    fn samples(f: &Fixture, p: &Params) -> Vec<Sample> {
        let Some(file) = &p.data else {
            return Vec::new();
        };
        let data = parse_gslib(&std::fs::read_to_string(f.dir.join(file)).unwrap()).unwrap();
        let col = |c: u32| data.column(usize::try_from(c - 1).unwrap()).unwrap();
        let [cx, cy, cz] = p.coordinate_columns;
        let (xs, ys, zs, vs) = (col(cx), col(cy), col(cz), col(p.variable_column));
        (0..data.nrows())
            .map(|i| Sample::new(Point3::new(xs[i], ys[i], zs[i]), vs[i]).unwrap())
            .collect()
    }

    /// gslib-rs's realizations for the case, with GSLIB's settings: SK with
    /// mean 0, an isotropic search for data (up to `max_data`) and nodes (up
    /// to `max_nodes`), GSLIB's node-table limits.
    fn ours(p: &Params, grid: &Grid, data: &[Sample]) -> Vec<Vec<f64>> {
        let model = CovarianceModel::try_from(
            harness::vmodel::model_from(p.nugget, &p.structures).unwrap(),
        )
        .unwrap();
        let ellipsoid = SearchEllipsoid::isotropic(p.search_radius).unwrap();
        let search = DataSearch::new(
            ellipsoid,
            usize::try_from(p.min_data).unwrap(),
            usize::try_from(p.max_data).unwrap(),
        )
        .unwrap();
        let nodes = NodeSearch::new(
            grid,
            &ellipsoid,
            &model,
            usize::try_from(p.max_nodes).unwrap(),
            None,
            TableExtent::Gslib { max: [71, 71, 1] },
        )
        .unwrap();
        let searches = Searches {
            data: DataConditioning::Search(&search),
            nodes: &nodes,
        };
        (0..u64::from(p.realizations))
            .map(|k| {
                sgsim(
                    data,
                    grid,
                    &model,
                    0.0,
                    searches,
                    &mut Seed(u64::from(p.seed)).realization(k),
                )
                .unwrap()
            })
            .collect()
    }

    /// One realization's statistics, in a fixed order with their names.
    fn statistics(values: &[f64], counts: [usize; 2]) -> Vec<(String, f64)> {
        let [nx, ny] = counts;
        let at = |i: usize, j: usize| values[i + nx * j];
        let mut out = Vec::new();
        for h in 1..=LAGS {
            let (mut sx, mut nxp) = (0.0, 0.0);
            for j in 0..ny {
                for i in 0..nx - h {
                    sx = (at(i, j) - at(i + h, j)).powi(2).mul_add(0.5, sx);
                    nxp += 1.0;
                }
            }
            out.push((format!("gamma_x[{h}]"), sx / nxp));
            let (mut sy, mut nyp) = (0.0, 0.0);
            for j in 0..ny - h {
                for i in 0..nx {
                    sy = (at(i, j) - at(i, j + h)).powi(2).mul_add(0.5, sy);
                    nyp += 1.0;
                }
            }
            out.push((format!("gamma_y[{h}]"), sy / nyp));
        }
        for d in 1..=9_u32 {
            let q = f64::from(d) / 10.0;
            let z = standard_normal_quantile(q).unwrap();
            let below = values.iter().filter(|&&v| v < z).count();
            out.push((
                format!("below_q{d}"),
                f64::from(u32::try_from(below).unwrap())
                    / f64::from(u32::try_from(values.len()).unwrap()),
            ));
        }
        out
    }

    /// Ensemble mean and the standard error of that mean.
    fn mean_se(xs: &[f64]) -> (f64, f64) {
        let n = f64::from(u32::try_from(xs.len()).unwrap());
        let mean = xs.iter().sum::<f64>() / n;
        let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (mean, (var / n).sqrt())
    }

    fn check(f: &Fixture) -> Vec<String> {
        let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
        let p = Params::from_table(&f.case.parameters).unwrap();
        let g = grid(&p);
        let data = samples(f, &p);
        let out = parse_gslib(&std::fs::read_to_string(f.dir.join(OUTPUT_FILE)).unwrap()).unwrap();
        let all = out.column(0).unwrap();
        let gslib: Vec<Vec<f64>> = all.chunks_exact(g.len()).map(<[f64]>::to_vec).collect();
        let ours = ours(&p, &g, &data);
        let mut cmp = Comparison::new(&f.case);
        cmp.structure(
            gslib.len() == ours.len() && all.len() == gslib.len() * g.len(),
            || {
                format!(
                    "{} GSLIB values for {} realizations of {} nodes",
                    all.len(),
                    ours.len(),
                    g.len()
                )
            },
        );
        let [nx, ny, _] = g.counts();
        let per = |runs: &[Vec<f64>]| -> Vec<Vec<(String, f64)>> {
            runs.iter().map(|r| statistics(r, [nx, ny])).collect()
        };
        let (sg, so) = (per(&gslib), per(&ours));
        let model = harness::vmodel::model_from(p.nugget, &p.structures).unwrap();
        let unconditional = data.is_empty();
        let (mut worst_parity, mut worst_ours, mut worst_gslib) = (0.0_f64, 0.0_f64, 0.0_f64);
        for (k, (stat, _)) in sg[0].iter().enumerate() {
            let column =
                |s: &[Vec<(String, f64)>]| -> Vec<f64> { s.iter().map(|r| r[k].1).collect() };
            let ((mg, eg), (mo, eo)) = (mean_se(&column(&sg)), mean_se(&column(&so)));
            let z = (mg - mo).abs() / eg.hypot(eo);
            worst_parity = worst_parity.max(z);
            cmp.structure(z <= SIGMAS, || {
                format!("{stat}: GSLIB {mg:.4} vs gslib-rs {mo:.4}, {z:.2} standard errors")
            });
            if unconditional {
                let target = model_value(stat, &model);
                let (zo, zg) = ((mo - target).abs() / eo, (mg - target).abs() / eg);
                worst_ours = worst_ours.max(zo);
                worst_gslib = worst_gslib.max(zg);
                cmp.structure(zo <= SIGMAS, || {
                    format!(
                        "{stat}: gslib-rs {mo:.4} vs model {target:.4}, {zo:.2} standard errors"
                    )
                });
            }
        }
        // Data reproduction at their nodes, in every realization of both.
        for d in &data {
            let node = g.node_index(d.location()).unwrap();
            for (r, (a, b)) in gslib.iter().zip(&ours).enumerate() {
                let at = || format!("datum at {:?}, realization {}", d.location(), r + 1);
                cmp.value("data_gslib", d.value(), a[node], 1.0, DIRECT, at);
                cmp.value("data_ours", d.value(), b[node], 1.0, DIRECT, at);
            }
        }
        println!(
            "{name}: {} statistics; worst GSLIB vs gslib-rs {worst_parity:.2} SE{}",
            sg[0].len(),
            if unconditional {
                format!("; vs model: gslib-rs {worst_ours:.2} SE, GSLIB {worst_gslib:.2} SE")
            } else {
                String::new()
            }
        );
        cmp.failures()
            .iter()
            .map(|m| format!("{name}: {m}"))
            .collect()
    }

    /// The model's value of a statistic: γ(h) for variograms (isotropic
    /// model, unit spacing), the decile for proportions.
    fn model_value(stat: &str, model: &gslib_rs::model::VariogramModel) -> f64 {
        if let Some(rest) = stat.strip_prefix("below_q") {
            return f64::from(rest.parse::<u32>().unwrap()) / 10.0;
        }
        let h: f64 = stat
            .split_once('[')
            .and_then(|(_, rest)| rest.strip_suffix(']'))
            .unwrap()
            .parse()
            .unwrap();
        model.gamma(Point3::default(), Point3::new(h, 0.0, 0.0))
    }

    /// The statistical cases: those with something to search.
    fn statistical(p: &Params) -> bool {
        let shortest = p.grid_spacing.iter().copied().fold(f64::INFINITY, f64::min);
        p.data.is_some() || p.search_radius >= shortest
    }

    #[test]
    fn ensembles_match_gslib_and_the_model() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "sgsim")
            .filter(|f| statistical(&Params::from_table(&f.case.parameters).unwrap()))
            .collect();
        assert!(!fixtures.is_empty(), "no statistical sgsim fixtures found");
        let failures: Vec<String> = fixtures.iter().flat_map(check).collect();
        assert!(
            failures.is_empty(),
            "{} mismatches:\n{}",
            failures.len(),
            failures
                .iter()
                .take(30)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
