//! gslib-rs's variogram model against GSLIB `vmodel` fixtures (plan step 4).

#[cfg(test)]
mod vmodel {
    use gslib_rs::geometry::Point3;
    use gslib_rs::model::{CovarianceModel, Kind};
    use harness::compare::{Comparison, DIRECT};
    use harness::vmodel::{evaluation_point, model, parse_output};
    use oracle_case::vmodel::{OUTPUT_FILE, Params};

    #[test]
    fn variogram_model_matches_gslib() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "vmodel")
            .collect();
        assert!(!fixtures.is_empty(), "no vmodel fixtures found");
        let mut failures = Vec::new();
        for f in &fixtures {
            let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
            let params = Params::from_table(&f.case.parameters).unwrap();
            let ours = model(&params).unwrap();
            let has_power = ours.structures().iter().any(|s| matches!(s.kind(), Kind::Power { .. }));
            let cov_model = CovarianceModel::try_from(ours.clone()).ok();
            // Scale for γ, C and the comparison: the total sill (c0 + Σcc).
            let sill = params.nugget + params.structures.iter().map(|s| s.cc).sum::<f64>();
            let tables = parse_output(&std::fs::read_to_string(f.dir.join(OUTPUT_FILE)).unwrap()).unwrap();

            let mut cmp = Comparison::new(&f.case);
            cmp.structure(tables.len() == params.directions.len(), || {
                format!("{} direction tables, expected {}", tables.len(), params.directions.len())
            });
            let expected_rows = usize::try_from(params.nlag).unwrap() + 2;
            for (d, (dir, table)) in params.directions.iter().zip(&tables).enumerate() {
                cmp.structure(table.rows.len() == expected_rows, || {
                    format!("direction {}: {} rows, expected {expected_rows}", d + 1, table.rows.len())
                });
                for (k, row) in table.rows.iter().enumerate() {
                    let il = if k == 0 { 1 } else { u32::try_from(k).unwrap() };
                    cmp.structure(row.il == il && row.ndir as usize == params.directions.len(), || {
                        format!("direction {} row {k}: il {} ndir {}", d + 1, row.il, row.ndir)
                    });
                    let p = evaluation_point(dir, k, il);
                    let at = || format!("direction {} il {il} (row {k})", d + 1);
                    // vmodel prints the first row's h on the repeated row.
                    let printed = if k == 1 { evaluation_point(dir, 0, 1) } else { p };
                    let h = printed.distance(Point3::default());
                    cmp.value("h", h, row.h, dir.lag, DIRECT, at);
                    let gamma = ours.gamma(Point3::default(), p);
                    cmp.value("gamma", gamma, row.gamma, sill, DIRECT, at);
                    if let Some(c) = &cov_model {
                        let cov = c.covariance(Point3::default(), p);
                        cmp.value("cov", cov, row.cov, sill, DIRECT, at);
                        cmp.value("ro", cov / c.sill(), row.ro, 1.0, DIRECT, at);
                    } else {
                        assert!(has_power, "only power models lack a covariance");
                        assert!(
                            cmp.is_expected_deviation("cov") && cmp.is_expected_deviation("ro"),
                            "{name}: power model cov/ro must be declared expected deviations"
                        );
                    }
                }
            }
            for (q, w) in cmp.worst() {
                println!("{name}: worst {q} |Δ|/scale = {:.3e} at {}", w.scaled_error, w.location);
            }
            failures.extend(cmp.failures().iter().map(|m| format!("{name}: {m}")));
        }
        assert!(failures.is_empty(), "{} mismatches:\n{}", failures.len(), failures.join("\n"));
    }
}
