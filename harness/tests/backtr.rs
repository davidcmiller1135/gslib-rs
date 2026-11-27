//! gslib-rs's back-transform against GSLIB `backtr` fixtures (plan step 5).
//!
//! What GSLIB does (`GSLIB/BACKTR.FOR`, F-007): inside the table, linear
//! interpolation in score space (as gslib-rs); in the tails (option 1),
//! linear interpolation in *probability* towards `zmin`/`zmax`, with the
//! probabilities from `gcum`, an approximation of Φ (D-016).

#[cfg(test)]
mod backtr {
    use gslib_rs::io::parse_gslib;
    use gslib_rs::transform::{BackTransform, NormalScoreTable, Tail, standard_normal_cdf};
    use harness::compare::{Comparison, DIRECT, Tolerance};
    use harness::gslib_math::GCUM_MAX_ERROR;
    use oracle_case::backtr::{OUTPUT_FILE, Params};

    #[test]
    fn back_transform_matches_gslib() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "backtr")
            .collect();
        assert!(!fixtures.is_empty(), "no backtr fixtures found");
        let mut failures = Vec::new();
        for f in &fixtures {
            let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
            let p = Params::from_table(&f.case.parameters).unwrap();
            let read = |file: &str| std::fs::read_to_string(f.dir.join(file)).unwrap();
            let table = harness::parse_columns(&read(&p.table), 2).unwrap();
            let pairs: Vec<(f64, f64)> = table.iter().map(|r| (r[0], r[1])).collect();
            let back = BackTransform::new(
                NormalScoreTable::from_pairs(&pairs).unwrap(),
                Tail::Linear(p.zmin),
                Tail::Linear(p.zmax),
            )
            .unwrap();
            let (first, last) = (pairs[0], pairs[pairs.len() - 1]);
            let range = p.zmax - p.zmin;

            let data = parse_gslib(&read(&p.data)).unwrap();
            let scores = data
                .column(usize::try_from(p.variable_column - 1).unwrap())
                .unwrap();
            let out = parse_gslib(&read(OUTPUT_FILE)).unwrap();
            let mut cmp = Comparison::new(&f.case);
            cmp.structure(
                out.nrows() == data.nrows() && out.ncols() == data.ncols() + 1,
                || format!("output is {}x{}", out.nrows(), out.ncols()),
            );
            let gslib = out.column(out.ncols() - 1).unwrap();

            for (i, (&y, &g)) in scores.iter().zip(&gslib).enumerate() {
                let ours = back.apply(y);
                let at = || format!("datum {} (score {y})", i + 1);
                if y > first.1 && y < last.1 {
                    cmp.value("value", ours, g, range, DIRECT, at);
                } else {
                    // Tail value = edge ± span·(ratio of probabilities); each
                    // probability is off by at most ε in gcum, so the value is
                    // off by at most span·2ε/(p_edge − ε).
                    let (span, p_edge) = if y <= first.1 {
                        (first.0 - p.zmin, standard_normal_cdf(first.1))
                    } else {
                        (p.zmax - last.0, 1.0 - standard_normal_cdf(last.1))
                    };
                    let bound = span * 2.0 * GCUM_MAX_ERROR / (p_edge - GCUM_MAX_ERROR);
                    let tol = Tolerance {
                        relative: 0.0,
                        absolute: bound,
                    };
                    cmp.value("tail_value", ours, g, 1.0, tol, at);
                }
            }
            for (q, w) in cmp.worst() {
                println!(
                    "{name}: worst {q} |Δ|/scale = {:.3e} at {}",
                    w.scaled_error, w.location
                );
            }
            failures.extend(cmp.failures().iter().map(|m| format!("{name}: {m}")));
        }
        assert!(
            failures.is_empty(),
            "{} mismatches:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
