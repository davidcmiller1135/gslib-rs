//! `variogram::gamv` against every GSLIB `gamv` fixture (`docs/spec/gamv.md`,
//! P1, P2): directions, every variogram type, several variables, indicators
//! and standardisation, with GSLIB's legacy reading (`harness::gamv::legacy`).
//!
//! Differences that D-017 keeps are checked as exact relations:
//! - **Pair counts.** GSLIB's `np` is gslib-rs's `orderings` (both
//!   orderings of an omnidirectional pair, D-043); in the zero-distance class
//!   GSLIB adds one self-pair per datum with both values present, two when
//!   omnidirectional (except type 2).
//! - **The zero-distance class's values** are compared only where GSLIB's
//!   class holds nothing but self-pairs and gslib-rs's is empty: there is no
//!   gslib-rs value to compare then, and the class's separation must be 0.
//!   The omnidirectional semivariogram cases (gamv.rs) check the self-pair
//!   relations of F-008 in full.
//!
//! - **A datum at exactly `tmax` in an omnidirectional direction** (F-026),
//!   only where a case declares the expected deviation
//!   `omnidirectional_at_tmax`: GSLIB tests the swapped ordering against
//!   `[tmin, tmax)`, so it loses one ordering per pair involving such a datum.
//!   Its count is then checked exactly, and values only in lags without such
//!   pairs.
//!
//! Values use D-015's 1e-12, relative to a scale per table: the largest of
//! |γ|, `|m_t·m_h|` (covariances are differences of such products) and the
//! sills for type −3; the largest |mean|; `xlag` for distances.

#[cfg(test)]
mod gamv_directional {
    use gslib_rs::variogram::{DirectionalVariogram, GamvLag, LagScheme, Separation, gamv};
    use harness::compare::{Comparison, DIRECT};
    use harness::gamv::legacy::{Legacy, build};
    use harness::gamv::{Row, Table, parse_output};
    use oracle_case::Fixture;
    use oracle_case::gamv::{OUTPUT_FILE, Params};

    /// GSLIB's title for a variogram type (`writeout`).
    fn title(gslib_type: i32) -> &'static str {
        match gslib_type.abs() {
            1 => "Semivariogram",
            2 => "Cross Semivariogram",
            3 => "Covariance",
            4 => "Correlogram",
            5 => "General Relative",
            6 => "Pairwise Relative",
            7 => "Variogram of Logarithms",
            8 => "Semimadogram",
            _ => "Indicator 1/2 Variogram",
        }
    }

    /// GSLIB's value for an undefined lag (D-049).
    const fn undefined(gslib_type: i32) -> f64 {
        if gslib_type == -3 { -999.0 } else { 0.0 }
    }

    struct Scales {
        gamma: f64,
        mean: f64,
        variance: f64,
    }

    fn scales(table: &Table, legacy: &Legacy, v: usize) -> Scales {
        let max = |f: &dyn Fn(&Row) -> f64| {
            table
                .rows
                .iter()
                .skip(1)
                .map(f)
                .fold(f64::MIN_POSITIVE, f64::max)
        };
        let spec = legacy.request.variograms.get(v);
        let sill = |i: Option<usize>| {
            i.and_then(|i| legacy.sills.get(i).copied().flatten())
                .unwrap_or(0.0)
                .abs()
        };
        let sills = sill(spec.map(|s| s.tail)).max(sill(spec.map(|s| s.head)));
        Scales {
            gamma: max(&|r: &Row| r.gamma.abs().max((r.head_mean * r.tail_mean).abs())).max(sills),
            mean: max(&|r: &Row| r.head_mean.abs().max(r.tail_mean.abs())),
            variance: max(&|r: &Row| r.variances.map_or(0.0, |(a, b)| a.abs().max(b.abs()))),
        }
    }

    struct Check<'a, 'c> {
        cmp: &'a mut Comparison<'c>,
        name: &'a str,
        spacing: f64,
    }

    impl Check<'_, '_> {
        fn lag(&mut self, at: &str, row: &Row, lag: &GamvLag, ty: i32, s: &Scales) {
            let at = || at.to_owned();
            let gslib_np = row.pairs;
            self.cmp
                .structure(gslib_np == u64::try_from(lag.orderings).unwrap(), || {
                    format!(
                        "{}: {}: GSLIB np {gslib_np} != orderings {}",
                        self.name,
                        at(),
                        lag.orderings
                    )
                });
            if lag.orderings == 0 {
                self.cmp
                    .structure(row.distance == 0.0 && row.gamma == 0.0, || {
                        format!("{}: {}: empty here, GSLIB {row:?}", self.name, at())
                    });
                return;
            }
            let ours = |x: Option<f64>| x.unwrap_or(f64::NAN);
            self.cmp.value(
                "distance",
                ours(lag.mean_distance),
                row.distance,
                self.spacing,
                DIRECT,
                at,
            );
            match lag.gamma {
                Some(g) => self.cmp.value("gamma", g, row.gamma, s.gamma, DIRECT, at),
                None => self
                    .cmp
                    .structure(row.gamma.to_bits() == undefined(ty).to_bits(), || {
                        format!(
                            "{}: {}: undefined here, GSLIB γ {}",
                            self.name,
                            at(),
                            row.gamma
                        )
                    }),
            }
            // GSLIB's hm is the tail mean, tm the head mean (`writeout`).
            self.cmp.value(
                "tail_mean",
                ours(lag.tail_mean),
                row.head_mean,
                s.mean,
                DIRECT,
                at,
            );
            self.cmp.value(
                "head_mean",
                ours(lag.head_mean),
                row.tail_mean,
                s.mean,
                DIRECT,
                at,
            );
            if let Some((hv, tv)) = row.variances {
                self.cmp.value(
                    "tail_variance",
                    ours(lag.tail_variance),
                    hv,
                    s.variance,
                    DIRECT,
                    at,
                );
                self.cmp.value(
                    "head_variance",
                    ours(lag.head_variance),
                    tv,
                    s.variance,
                    DIRECT,
                    at,
                );
            }
        }
    }

    /// The expected deviation F-026 introduces.
    const AT_TMAX: &str = "omnidirectional_at_tmax";

    /// Per lag: the orderings GSLIB loses in omnidirectional direction `d`
    /// because a value of the variogram's pair equals `tmax` (F-026), counted
    /// independently with gslib-rs's public pair primitives. Only for
    /// variograms with tail = head.
    fn lost_at_tmax(legacy: &Legacy, variogram: usize, direction: usize) -> Vec<usize> {
        let spec = legacy.request.variograms[variogram];
        assert_eq!(spec.tail, spec.head, "F-026 relation needs tail = head");
        let column = &legacy.columns[spec.tail];
        let window = legacy.request.directions[direction].window(legacy.request.policies.pi);
        let mut lost = vec![0; legacy.request.lags.len()];
        let count = legacy.points.len();
        for i in 0..count {
            for j in i..count {
                let (Some(a), Some(b)) = (column[i], column[j]) else {
                    continue;
                };
                if a.to_bits() != legacy.tmax.to_bits() && b.to_bits() != legacy.tmax.to_bits() {
                    continue;
                }
                let sep = Separation::new(legacy.points[i], legacy.points[j]);
                if window.accept(&sep).is_none() {
                    continue;
                }
                for k in legacy.request.lags.classes(sep.distance()) {
                    lost[k] += 1;
                }
            }
        }
        lost
    }

    fn compare(
        f: &Fixture,
        name: &str,
        p: &Params,
        legacy: &Legacy,
        ours: &[DirectionalVariogram],
        tables: &[Table],
    ) -> Vec<String> {
        let mut cmp = Comparison::new(&f.case);
        cmp.structure(tables.len() == ours.len(), || {
            format!("{name}: {} tables, expected {}", tables.len(), ours.len())
        });
        let directions = p
            .directions
            .clone()
            .unwrap_or_else(|| vec![[0.0, 90.0, 0.0, 0.0, 0.0, 0.0]]);
        for (dv, table) in ours.iter().zip(tables) {
            let g = &legacy.variograms[dv.variogram];
            let expected = format!(
                "{:<24}tail:{:<12} head:{:<12}  direction {:>2}",
                title(g.gslib_type),
                g.tail_name,
                g.head_name,
                dv.direction + 1
            );
            cmp.structure(table.title == expected.trim(), || {
                format!(
                    "{name}: title {:?}, expected {:?}",
                    table.title,
                    expected.trim()
                )
            });
            cmp.structure(table.rows.len() == dv.lags.len(), || {
                format!(
                    "{name}: {} rows, expected {}",
                    table.rows.len(),
                    dv.lags.len()
                )
            });
            let omni = directions[dv.direction][1] >= 90.0;
            let lost = if omni && cmp.is_expected_deviation(AT_TMAX) {
                lost_at_tmax(legacy, dv.variogram, dv.direction)
            } else {
                vec![0; dv.lags.len()]
            };
            let s = scales(table, legacy, dv.variogram);
            let mut check = Check {
                cmp: &mut cmp,
                name,
                spacing: p.lag_spacing,
            };
            for (k, (row, lag)) in table.rows.iter().zip(&dv.lags).enumerate() {
                let at = format!(
                    "variogram {} direction {} lag {}",
                    dv.variogram + 1,
                    dv.direction + 1,
                    k + 1
                );
                if k == 0 {
                    let selfs = legacy.self_pairs(dv.variogram, omni);
                    let with_selfs = u64::try_from(lag.orderings + selfs - lost[0]).unwrap();
                    check.cmp.structure(row.pairs == with_selfs, || {
                        format!(
                            "{name}: {at}: GSLIB np {} != orderings + self-pairs {with_selfs}",
                            row.pairs
                        )
                    });
                    check.cmp.structure(row.distance.abs() <= 1e-20, || {
                        format!("{name}: {at}: zero class at distance {}", row.distance)
                    });
                } else if lost[k] > 0 {
                    // F-026: GSLIB summed fewer orderings; only the count can
                    // be compared.
                    let remaining = u64::try_from(lag.orderings - lost[k]).unwrap();
                    check.cmp.structure(row.pairs == remaining, || {
                        format!(
                            "{name}: {at}: GSLIB np {} != orderings − lost at tmax {remaining}",
                            row.pairs
                        )
                    });
                } else {
                    check.lag(&at, row, lag, g.gslib_type, &s);
                }
            }
        }
        for (q, w) in cmp.worst() {
            println!(
                "{name}: worst {q} |Δ|/scale = {:.3e} at {}",
                w.scaled_error, w.location
            );
        }
        cmp.failures().to_vec()
    }

    /// Compare one fixture; returns its mismatches.
    fn check(f: &Fixture) -> Vec<String> {
        let name = f.dir.file_name().unwrap().to_string_lossy().into_owned();
        let p = Params::from_table(&f.case.parameters).unwrap();
        let legacy = build(&f.dir, &p).unwrap();
        let text = std::fs::read_to_string(f.dir.join(OUTPUT_FILE)).unwrap();
        let tables = parse_output(&text).unwrap();
        let non_finite = tables
            .iter()
            .flat_map(|t| &t.rows)
            .any(|r| !(r.gamma.is_finite() && r.head_mean.is_finite() && r.tail_mean.is_finite()));
        match (gamv(&legacy.data, &legacy.request), non_finite) {
            // P2: where GSLIB writes NaN or Infinity, gslib-rs must fail.
            (Err(_), true) => {
                println!("{name}: GSLIB wrote a non-finite value; gslib-rs fails, as required");
                Vec::new()
            }
            (Ok(_), true) => vec![format!(
                "{name}: GSLIB wrote a non-finite value, gslib-rs did not fail"
            )],
            (Err(e), false) => vec![format!("{name}: gslib-rs failed: {e}")],
            (Ok(ours), false) => compare(f, &name, &p, &legacy, &ours, &tables),
        }
    }

    #[test]
    fn gamv_matches_gslib() {
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
