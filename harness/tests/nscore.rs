//! gslib-rs's normal-score transform against GSLIB `nscore` fixtures (plan
//! step 5).
//!
//! What GSLIB does (`NSCORE/NSCORE.FOR`, F-005, F-007): each datum becomes
//! `v + acorni·1e-6` (ACORN, seed 69069), the perturbed values are sorted,
//! and each gets `gauinv` of the midpoint of its cumulative-weight interval.
//! For data whose values differ by more than 1e-6 this is gslib-rs's
//! transform, except that `gauinv` is an approximation (D-016).

#[cfg(test)]
mod nscore {
    use gslib_rs::io::parse_gslib;
    use gslib_rs::transform::normal_scores;
    use harness::compare::{Comparison, Tolerance};
    use harness::gslib_math::GAUINV_MAX_ERROR;
    use oracle_case::nscore::{OUTPUT_FILE, Params, TABLE_FILE};

    /// Scores: absolute, from `gauinv`'s measured error.
    const SCORE: Tolerance = Tolerance {
        relative: 0.0,
        absolute: GAUINV_MAX_ERROR,
    };
    /// Table values: GSLIB stores `v + acorni·1e-6`, with acorni in [0, 1).
    const DESPIKE: Tolerance = Tolerance {
        relative: 0.0,
        absolute: 1e-6,
    };

    #[test]
    fn normal_scores_match_gslib() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "nscore")
            .collect();
        assert!(!fixtures.is_empty(), "no nscore fixtures found");
        let mut failures = Vec::new();
        for f in &fixtures {
            let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
            let p = Params::from_table(&f.case.parameters).unwrap();
            let read = |file: &str| std::fs::read_to_string(f.dir.join(file)).unwrap();
            let data = parse_gslib(&read(&p.data)).unwrap();
            let col = |c: u32| data.column(usize::try_from(c - 1).unwrap()).unwrap();
            let values = col(p.variable_column);
            let weights = (p.weight_column > 0).then(|| col(p.weight_column));
            let ours = normal_scores(&values, weights.as_deref()).unwrap();

            let out = parse_gslib(&read(OUTPUT_FILE)).unwrap();
            let table = harness::parse_columns(&read(TABLE_FILE), 2).unwrap();
            let mut cmp = Comparison::new(&f.case);
            cmp.structure(
                out.ncols() == data.ncols() + 1 && out.nrows() == data.nrows(),
                || {
                    format!(
                        "output is {}x{}, data {}x{}",
                        out.nrows(),
                        out.ncols(),
                        data.nrows(),
                        data.ncols()
                    )
                },
            );
            cmp.structure(
                out.names().last().map(String::as_str) == Some("Normal Score value"),
                || format!("last column is {:?}", out.names().last()),
            );
            cmp.structure(table.len() == values.len(), || {
                format!("table has {} rows", table.len())
            });
            let gslib = out.column(out.ncols() - 1).unwrap();

            // Scores, datum by datum.
            for (i, (&o, &g)) in ours.scores.iter().zip(&gslib).enumerate() {
                cmp.value("score", o, g, 1.0, SCORE, || format!("datum {}", i + 1));
            }
            // Distinct values keep their order, whatever happens to ties.
            for (i, (&vi, &gi)) in values.iter().zip(&gslib).enumerate() {
                for (j, (&vj, &gj)) in values.iter().zip(&gslib).enumerate() {
                    if vi < vj {
                        cmp.structure(gi < gj, || {
                            format!("data {} < {} but scores are not ordered", i + 1, j + 1)
                        });
                    }
                }
            }
            // The table: GSLIB's perturbed values against the sorted data,
            // and its scores against ours in the same order.
            let mut sorted: Vec<(f64, f64)> = values
                .iter()
                .copied()
                .zip(ours.scores.iter().copied())
                .collect();
            sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
            for (k, (row, &(v, s))) in table.iter().zip(&sorted).enumerate() {
                let at = || format!("table row {}", k + 1);
                cmp.structure(row[0] >= v, || {
                    format!("table row {}: despiking only adds", k + 1)
                });
                cmp.value("table_value", v, row[0], 1.0, DESPIKE, at);
                cmp.value("table_score", s, row[1], 1.0, SCORE, at);
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
