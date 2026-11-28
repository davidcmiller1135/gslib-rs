//! gslib-rs's GSLIB random-stream building blocks against `sgsim` fixtures
//! (`docs/spec/rng.md` RN3, RN5, RN6, RN7).
//!
//! The fixtures are unconditional runs with nothing to search, so each node
//! is `gmean + 1·gauinv(p) = gauinv(p)` (`SGSIM.FOR` lines 939–958). Per
//! realization, from one ACORN stream: one draw per node builds the random
//! path, then one draw per node in path order gives its deviate. GSLIB's
//! output must equal that, bit for bit.

#[cfg(test)]
mod sgsim {
    use gslib_rs::io::parse_gslib;
    use gslib_rs::random::{Acorn, gslib_path};
    use gslib_rs::simulation::{Axis, Grid};
    use gslib_rs::transform::gslib_gauinv;
    use harness::compare::{Comparison, DIRECT};
    use oracle_case::Fixture;
    use oracle_case::sgsim::{OUTPUT_FILE, Params};

    /// The RN7 cases: no data and a search radius below the node spacing,
    /// so every node is `gauinv(p)`. The other `sgsim` cases are statistical
    /// (`sgsim_stats.rs`).
    fn nothing_to_search(p: &Params) -> bool {
        let shortest = p.grid_spacing.iter().copied().fold(f64::INFINITY, f64::min);
        p.data.is_none() && p.search_radius < shortest
    }

    fn grid(p: &Params) -> Grid {
        let axis = |a: usize| Axis {
            count: usize::try_from(p.grid_counts[a]).unwrap(),
            first_center: p.grid_origin[a],
            spacing: p.grid_spacing[a],
        };
        Grid::new(axis(0), axis(1), axis(2)).unwrap()
    }

    /// The realizations gslib-rs predicts from GSLIB's stream.
    fn predicted(p: &Params, grid: &Grid) -> Vec<Vec<f64>> {
        let mut acorn = Acorn::gslib(u64::from(p.seed)).unwrap();
        let levels = (p.multiple_grid > 0).then(|| usize::try_from(p.multiple_grid).unwrap());
        (0..p.realizations)
            .map(|_| {
                let path = gslib_path(&mut acorn, grid, levels);
                let mut values = vec![f64::NAN; grid.len()];
                for node in path {
                    values[node] = gslib_gauinv(acorn.next_acorni());
                }
                values
            })
            .collect()
    }

    fn check(fixture: &Fixture) -> Vec<String> {
        let name = fixture
            .dir
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let p = Params::from_table(&fixture.case.parameters).unwrap();
        let g = grid(&p);
        let out =
            parse_gslib(&std::fs::read_to_string(fixture.dir.join(OUTPUT_FILE)).unwrap()).unwrap();
        let gslib = out.column(0).unwrap();
        let ours: Vec<f64> = predicted(&p, &g).concat();
        let mut cmp = Comparison::new(&fixture.case);
        cmp.structure(out.names() == ["value"], || {
            format!("columns {:?}", out.names())
        });
        cmp.structure(gslib.len() == ours.len(), || {
            format!("{} values, expected {}", gslib.len(), ours.len())
        });
        let mut exact = 0;
        for (i, (&o, &v)) in ours.iter().zip(&gslib).enumerate() {
            let (r, node) = (i.div_euclid(g.len()), i % g.len());
            let at = || format!("realization {}, node {}", r + 1, node + 1);
            cmp.value("deviate", o, v, 1.0, DIRECT, at);
            cmp.structure(o.to_bits() == v.to_bits(), || {
                format!("{}: not bit-identical", at())
            });
            exact += usize::from(o.to_bits() == v.to_bits());
        }
        println!("{name}: {exact} of {} values bit-identical", gslib.len());
        cmp.failures()
            .iter()
            .map(|m| format!("{name}: {m}"))
            .collect()
    }

    #[test]
    fn gslib_streams_are_reproduced() {
        let fixtures: Vec<_> = oracle_case::discover(&harness::fixture_root())
            .unwrap()
            .into_iter()
            .filter(|f| f.case.program == "sgsim")
            .filter(|f| nothing_to_search(&Params::from_table(&f.case.parameters).unwrap()))
            .collect();
        assert!(!fixtures.is_empty(), "no sgsim stream fixtures found");
        let failures: Vec<String> = fixtures.iter().flat_map(check).collect();
        assert!(
            failures.is_empty(),
            "{} mismatches:\n{}",
            failures.len(),
            failures
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
