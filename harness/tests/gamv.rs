//! gslib-rs's experimental variogram against GSLIB `gamv` fixtures (plan
//! step 6): omnidirectional semivariograms with GSLIB's centred lags.
//!
//! Two documented differences are checked as exact relations, not skipped
//! (F-008):
//! - **Pair counts:** GSLIB's omnidirectional mode counts each pair in both
//!   orderings, so it reports twice gslib-rs's count (γ is unaffected).
//! - **The zero-distance class:** GSLIB also pairs each datum with itself.
//!   With `n` data and `m` coincident pairs of distinct data, it reports
//!   `np = 2(n + m)`, `γ = γ_ours · m / (n + m)`, and head and tail means of
//!   `(Σv + m · mean_ours) / (n + m)`.

#[cfg(test)]
mod gamv {
    use gslib_rs::geometry::Point3;
    use gslib_rs::io::parse_gslib;
    use gslib_rs::sample::Sample;
    use gslib_rs::variogram::{CentredLags, Lag, experimental_variogram};
    use harness::compare::{Comparison, DIRECT};
    use harness::gamv::{Row, parse_output};
    use oracle_case::Fixture;
    use oracle_case::gamv::{OUTPUT_FILE, Params};

    /// Summary statistics of the data, for scales and the zero-class formulas.
    struct Stats {
        count: usize,
        n: f64,
        sum: f64,
        /// Scale for γ: the data variance.
        gamma_scale: f64,
        /// Scale for means: the larger of the range and |mean|.
        value_scale: f64,
    }

    fn stats(values: &[f64]) -> Stats {
        let n = f64::from(u32::try_from(values.len()).unwrap());
        let sum: f64 = values.iter().sum();
        let mean = sum / n;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
        let (lo, hi) = values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                (a.min(v), b.max(v))
            });
        Stats {
            count: values.len(),
            n,
            sum,
            gamma_scale: variance.max(f64::MIN_POSITIVE),
            value_scale: (hi - lo).max(mean.abs()).max(f64::MIN_POSITIVE),
        }
    }

    /// The fixture's samples (column 0 means "coordinate is 0", as in GSLIB).
    fn samples(f: &Fixture, p: &Params) -> Vec<Sample> {
        let data = parse_gslib(&std::fs::read_to_string(f.dir.join(&p.data)).unwrap()).unwrap();
        let column = |c: u32| {
            if c == 0 {
                vec![0.0; data.nrows()]
            } else {
                data.column(usize::try_from(c - 1).unwrap()).unwrap()
            }
        };
        let [cx, cy, cz] = p.coordinate_columns;
        let (xs, ys, zs, vs) = (
            column(cx),
            column(cy),
            column(cz),
            column(p.variable_column),
        );
        xs.iter()
            .zip(&ys)
            .zip(&zs)
            .zip(&vs)
            .map(|(((&x, &y), &z), &v)| Sample::new(Point3::new(x, y, z), v).unwrap())
            .collect()
    }

    /// The zero-distance class, through GSLIB's self-pair relations.
    fn zero_class(cmp: &mut Comparison<'_>, row: &Row, lag: &Lag, s: &Stats, spacing: f64) {
        let at = || "lag 1 (zero distance)".to_owned();
        let m = f64::from(u32::try_from(lag.pairs).unwrap());
        let expected_pairs = u64::try_from(2 * (s.count + lag.pairs)).unwrap();
        cmp.structure(row.pairs == expected_pairs, || {
            format!(
                "zero class: GSLIB np {} != 2(n + m) = {expected_pairs}",
                row.pairs
            )
        });
        let gamma = lag.gamma.unwrap_or(0.0) * m / (s.n + m);
        let means = m.mul_add(lag.pair_mean.unwrap_or(0.0), s.sum) / (s.n + m);
        cmp.value(
            "zero_class_gamma",
            gamma,
            row.gamma,
            s.gamma_scale,
            DIRECT,
            at,
        );
        cmp.value(
            "zero_class_mean",
            means,
            row.head_mean,
            s.value_scale,
            DIRECT,
            at,
        );
        cmp.value(
            "zero_class_mean",
            means,
            row.tail_mean,
            s.value_scale,
            DIRECT,
            at,
        );
        cmp.value("distance", 0.0, row.distance, spacing, DIRECT, at);
    }

    /// A centred lag: GSLIB counts each pair twice, otherwise values match.
    fn centred_lag(
        cmp: &mut Comparison<'_>,
        k: usize,
        row: &Row,
        lag: &Lag,
        s: &Stats,
        spacing: f64,
    ) {
        let at = || format!("lag {}", k + 1);
        cmp.structure(row.pairs == u64::try_from(2 * lag.pairs).unwrap(), || {
            format!("lag {}: GSLIB np {} != 2 x {}", k + 1, row.pairs, lag.pairs)
        });
        match (lag.mean_distance, lag.gamma, lag.pair_mean) {
            (Some(d), Some(g), Some(pm)) => {
                cmp.value("distance", d, row.distance, spacing, DIRECT, at);
                cmp.value("gamma", g, row.gamma, s.gamma_scale, DIRECT, at);
                cmp.value("pair_mean", pm, row.head_mean, s.value_scale, DIRECT, at);
                cmp.value("pair_mean", pm, row.tail_mean, s.value_scale, DIRECT, at);
            }
            _ => cmp.structure(
                row.distance == 0.0 && row.gamma == 0.0 && row.pairs == 0,
                || format!("lag {}: empty in gslib-rs but GSLIB reports {row:?}", k + 1),
            ),
        }
    }

    /// Compare one fixture; returns its mismatches.
    fn check(f: &Fixture) -> Vec<String> {
        let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
        let p = Params::from_table(&f.case.parameters).unwrap();
        let samples = samples(f, &p);
        let values: Vec<f64> = samples.iter().map(Sample::value).collect();
        let s = stats(&values);
        let lags = CentredLags::new(
            usize::try_from(p.nlag).unwrap(),
            p.lag_spacing,
            p.lag_tolerance,
        )
        .unwrap();
        let ours = experimental_variogram(&samples, &lags);
        let tables =
            parse_output(&std::fs::read_to_string(f.dir.join(OUTPUT_FILE)).unwrap()).unwrap();

        let mut cmp = Comparison::new(&f.case);
        cmp.structure(tables.len() == 1, || {
            format!("{} tables, expected 1", tables.len())
        });
        let Some(table) = tables.first() else {
            return vec![format!("{name}: no output table")];
        };
        cmp.structure(table.rows.len() == ours.len(), || {
            format!("{} rows, expected {}", table.rows.len(), ours.len())
        });
        for (k, (row, lag)) in table.rows.iter().zip(&ours).enumerate() {
            cmp.structure(row.il == u32::try_from(k + 1).unwrap(), || {
                format!("row {k}: il {}", row.il)
            });
            if k == 0 {
                zero_class(&mut cmp, row, lag, &s, p.lag_spacing);
            } else {
                centred_lag(&mut cmp, k, row, lag, &s, p.lag_spacing);
            }
        }
        for (q, w) in cmp.worst() {
            println!(
                "{name}: worst {q} |Δ|/scale = {:.3e} at {}",
                w.scaled_error, w.location
            );
        }
        cmp.failures()
            .iter()
            .map(|m| format!("{name}: {m}"))
            .collect()
    }

    #[test]
    fn omnidirectional_semivariogram_matches_gslib() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "gamv")
            .collect();
        assert!(!fixtures.is_empty(), "no gamv fixtures found");
        let failures: Vec<String> = fixtures.iter().flat_map(check).collect();
        assert!(
            failures.is_empty(),
            "{} mismatches:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
