//! End-to-end tests of the `gslib-rs` binary.

#[cfg(test)]
mod cli {
    use std::path::{Path, PathBuf};

    use assert_cmd::Command;
    use predicates::str::contains;
    use tempfile::TempDir;

    /// A GSLIB file with columns X, Y, Z, V.
    const DATA: &str = "Test Data\n4\nX\nY\nZ\nV\n0 0 0 1\n10 0 0 2\n0 10 0 3\n";

    /// Model flags: exponential, sill 1, practical range 30.
    const MODEL: [&str; 2] = ["--structure", "exp:1:30"];

    fn write(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let path = dir.path().join(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    fn bin() -> Command {
        Command::cargo_bin("gslib-rs").unwrap()
    }

    /// Parse a GSLIB file into (names, rows).
    fn read(path: &Path) -> (Vec<String>, Vec<Vec<f64>>) {
        let text = std::fs::read_to_string(path).unwrap();
        let mut lines = text.lines().skip(1);
        let nvar: usize = lines.next().unwrap().trim().parse().unwrap();
        let names = lines.by_ref().take(nvar).map(str::to_owned).collect();
        let rows = lines
            .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
            .collect();
        (names, rows)
    }

    fn column(rows: &[Vec<f64>], j: usize) -> Vec<f64> {
        rows.iter().map(|r| r[j]).collect()
    }

    fn realization(prefix: &Path, i: usize) -> PathBuf {
        PathBuf::from(format!("{}_{i}.gslib", prefix.display()))
    }

    #[test]
    fn variogram_writes_csv_and_plot_with_model() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let csv = dir.path().join("v.csv");
        let png = dir.path().join("v.png");
        bin()
            .args(["variogram", "--max-dist", "20", "--n-lags", "2"])
            .arg("--input")
            .arg(&data)
            .arg("--out-csv")
            .arg(&csv)
            .arg("--plot")
            .arg(&png)
            .args(MODEL)
            .assert()
            .success();
        let text = std::fs::read_to_string(&csv).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "lag,mean_distance,pairs,gamma");
        assert_eq!(lines.len(), 3, "{text}");
        // Two pairs at distance 10 and one at 14.1: all in the second lag.
        assert_eq!(lines[1], "5,,0,", "first lag is empty, not zero");
        assert!(lines[2].starts_with("15,"), "{}", lines[2]);
        assert!(lines[2].contains(",3,"), "three pairs: {}", lines[2]);
        assert!(std::fs::metadata(&png).unwrap().len() > 0, "plot is empty");
    }

    #[test]
    fn variogram_defaults_to_stdout() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        bin()
            .args(["variogram", "-d", "20", "-n", "1", "-i"])
            .arg(&data)
            .assert()
            .success()
            .stdout(contains("lag,mean_distance,pairs,gamma"))
            .stdout(contains(",3,"));
    }

    #[test]
    fn ordinary_kriging_honours_data_and_reports_variance() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let targets = write(&dir, "t.dat", "targets\n3\nX\nY\nZ\n0 0 0\n5 5 0\n");
        let out = dir.path().join("k.dat");
        bin()
            .args(["krige", "--kind", "ordinary"])
            .arg("--data")
            .arg(&data)
            .arg("--targets")
            .arg(&targets)
            .arg("--out")
            .arg(&out)
            .args(MODEL)
            .assert()
            .success();
        let (names, rows) = read(&out);
        assert_eq!(names, ["X", "Y", "Z", "estimate", "variance"]);
        assert!((rows[0][3] - 1.0).abs() < 1e-9, "exact at datum: {rows:?}");
        assert!(rows[0][4].abs() < 1e-9, "zero variance at datum: {rows:?}");
        assert!(
            rows[1][4] > 0.0,
            "positive variance away from data: {rows:?}"
        );
    }

    #[test]
    fn simple_kriging_far_away_returns_the_mean() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let targets = write(&dir, "t.dat", "targets\n3\nX\nY\nZ\n1e6 0 0\n");
        let out = dir.path().join("k.dat");
        bin()
            .args(["krige", "--kind", "simple", "--mean", "7.5"])
            .arg("--data")
            .arg(&data)
            .arg("--targets")
            .arg(&targets)
            .arg("--out")
            .arg(&out)
            .args(MODEL)
            .assert()
            .success();
        let (_, rows) = read(&out);
        assert!((rows[0][3] - 7.5).abs() < 1e-9, "{rows:?}");
        assert!(
            (rows[0][4] - 1.0).abs() < 1e-9,
            "variance is the sill: {rows:?}"
        );
    }

    fn sgsim(dir: &TempDir, prefix: &str, seed: &str) -> PathBuf {
        let data = write(dir, "d.dat", DATA);
        let prefix = dir.path().join(prefix);
        bin()
            .args(["sgsim", "--nsim", "2", "--seed", seed])
            .args(["--grid-x", "3:0:5", "--grid-y", "3:0:5"])
            .arg("--data")
            .arg(&data)
            .arg("--out-prefix")
            .arg(&prefix)
            .args(MODEL)
            .assert()
            .success();
        prefix
    }

    #[test]
    fn sgsim_appends_back_transformed_realizations() {
        let dir = TempDir::new().unwrap();
        let prefix = sgsim(&dir, "sgs", "1");
        for i in 0..2 {
            let (names, rows) = read(&realization(&prefix, i));
            assert_eq!(names, ["X", "Y", "Z", &format!("sim_{i}")]);
            assert_eq!(rows.len(), 9, "3 x 3 nodes");
            assert_eq!(
                (column(&rows, 0)[1], column(&rows, 1)[3]),
                (5.0, 5.0),
                "x fastest, spacing 5"
            );
            let sim = column(&rows, 3);
            assert!((sim[0] - 1.0).abs() < 1e-9, "datum reproduced: {sim:?}");
            assert!(
                sim.iter().all(|&v| (1.0..=3.0).contains(&v)),
                "clamped to data range: {sim:?}"
            );
        }
    }

    #[test]
    fn sgsim_is_reproducible_from_its_seed() {
        let dir = TempDir::new().unwrap();
        let a = read(&realization(&sgsim(&dir, "a", "42"), 1));
        let b = read(&realization(&sgsim(&dir, "b", "42"), 1));
        assert_eq!(a, b);
    }

    fn sisim(threshold: &str) -> Vec<f64> {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let prefix = dir.path().join("sis");
        bin()
            .args(["sisim", "--threshold", threshold])
            .args(["--grid-x", "3:1:3"])
            .arg("--data")
            .arg(&data)
            .arg("--out-prefix")
            .arg(&prefix)
            .args(MODEL)
            .assert()
            .success();
        let (names, rows) = read(&realization(&prefix, 0));
        assert_eq!(names.last().map(String::as_str), Some("ind_0"));
        column(&rows, 3)
    }

    #[test]
    fn sisim_all_below_threshold_gives_all_ones() {
        assert_eq!(sisim("5"), vec![1.0; 3]);
    }

    #[test]
    fn sisim_all_above_threshold_gives_all_zeros() {
        assert_eq!(sisim("0"), vec![0.0; 3]);
    }

    #[test]
    fn missing_value_column_is_a_clear_error() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", "coords\n3\nX\nY\nZ\n0 0 0\n");
        bin()
            .args(["variogram", "-d", "1", "-n", "1", "-i"])
            .arg(&data)
            .assert()
            .failure()
            .stderr(contains("column 4 requested but the file has 3"));
    }

    #[test]
    fn model_is_required_for_simulation() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        bin()
            .args(["sgsim", "--out-prefix", "unused", "--grid-x", "2:0:1"])
            .arg("--data")
            .arg(&data)
            .assert()
            .failure()
            .stderr(contains("total sill must be positive"));
    }

    #[test]
    fn bad_structure_spec_is_rejected_by_the_parser() {
        bin()
            .args([
                "variogram",
                "-i",
                "x",
                "-d",
                "1",
                "-n",
                "1",
                "--structure",
                "sph:1",
            ])
            .assert()
            .failure()
            .stderr(contains("expected TYPE:CC:A_HMAX"));
    }

    #[test]
    fn malformed_file_reports_line_number() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", "t\n4\nX\nY\nZ\nV\n0 0 0 1\n0 0 oops 1\n");
        bin()
            .args(["variogram", "-d", "1", "-n", "1", "-i"])
            .arg(&data)
            .assert()
            .failure()
            .stderr(contains("line 8"));
    }

    /// The `with_sentinels` suite: 40 points, every fifth value −999.
    fn sentinels() -> String {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/data/with_sentinels.dat");
        std::fs::read_to_string(path).unwrap()
    }

    /// The same file with the rows whose value is −999 removed.
    fn without_sentinel_rows(text: &str) -> String {
        let mut lines: Vec<&str> = text.lines().collect();
        lines.retain(|l| !l.trim_end().ends_with("-999.000000"));
        lines.join("\n") + "\n"
    }

    fn variogram_csv(data: &Path, extra: &[&str]) -> assert_cmd::assert::Assert {
        bin()
            .args(["variogram", "-d", "40", "-n", "4", "-i"])
            .arg(data)
            .args(extra)
            .assert()
    }

    fn stdout(a: &assert_cmd::assert::Assert) -> String {
        String::from_utf8(a.get_output().stdout.clone()).unwrap()
    }

    #[test]
    fn sentinels_need_explicit_trimming_in_every_command() {
        // M6: each command that reads data fails before computing, naming
        // the column and the first offending row.
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "s.dat", &sentinels());
        let nodes = write(&dir, "n.dat", "nodes\n3\nX\nY\nZ\n1 1 0\n");
        let out = dir.path().join("out");
        let file = |flag: &str, p: &Path| [flag.to_owned(), p.display().to_string()];
        let commands: Vec<Vec<String>> = vec![
            ["variogram", "-d", "40", "-n", "4"]
                .map(str::to_owned)
                .to_vec(),
            [
                ["krige".to_owned()].as_slice(),
                &file("--targets", &nodes),
                &file("--out", &out),
            ]
            .concat(),
            [
                ["sgsim".to_owned()].as_slice(),
                &["--grid-x".to_owned(), "3:1:1".to_owned()],
                &file("--out-prefix", &out),
            ]
            .concat(),
            [
                ["sisim".to_owned(), "--threshold".to_owned(), "1".to_owned()].as_slice(),
                &["--grid-x".to_owned(), "3:1:1".to_owned()],
                &file("--out-prefix", &out),
            ]
            .concat(),
        ];
        for args in commands {
            let data_flag = if args[0] == "variogram" {
                "--input"
            } else {
                "--data"
            };
            bin()
                .args(&args)
                .args(MODEL)
                .arg(data_flag)
                .arg(&data)
                .assert()
                .failure()
                .stderr(contains(
                    "column 4 (V), row 1: -999 looks like a missing-value code",
                ))
                .stderr(contains("--tmin"));
            assert!(!out.exists(), "{}: wrote output", args[0]);
        }
    }

    #[test]
    fn tmin_trims_sentinel_rows() {
        // Trimming with --tmin gives exactly the result on the file with
        // those rows deleted.
        let dir = TempDir::new().unwrap();
        let text = sentinels();
        let data = write(&dir, "s.dat", &text);
        let clean = write(&dir, "c.dat", &without_sentinel_rows(&text));
        assert_eq!(
            without_sentinel_rows(&text).lines().count(),
            text.lines().count() - 8,
            "8 sentinel rows"
        );
        let trimmed = stdout(&variogram_csv(&data, &["--tmin", "-998"]).success());
        let reference = stdout(&variogram_csv(&clean, &[]).success());
        assert_eq!(trimmed, reference);
        let equals_form = stdout(&variogram_csv(&data, &["--tmin=-998", "--tmax=1e21"]).success());
        assert_eq!(
            equals_form, reference,
            "--flag=value form and an explicit tmax"
        );
        // --no-trimming accepts the codes as data, which changes the result.
        let as_data = stdout(&variogram_csv(&data, &["--no-trimming"]).success());
        assert_ne!(as_data, reference, "-999 used as data");
    }

    #[test]
    fn trimming_flag_errors() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "s.dat", &sentinels());
        variogram_csv(&data, &["--no-trimming", "--tmin", "-998"])
            .failure()
            .stderr(contains("cannot be used with"));
        variogram_csv(&data, &["--tmin", "5", "--tmax", "1"])
            .failure()
            .stderr(contains("tmin 5 must be less than tmax 1"));
        variogram_csv(&data, &["--tmin", "1e6"])
            .failure()
            .stderr(contains("all 40 rows are missing after trimming"));
    }

    #[test]
    fn non_finite_values_follow_their_policies() {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", &format!("{DATA}5 5 0 NaN\n7 7 0 inf\n"));
        variogram_csv(&data, &[])
            .failure()
            .stderr(contains("row 4, column 4: NaN is not finite"))
            .stderr(contains("--nan missing"));
        variogram_csv(&data, &["--nan", "missing"])
            .failure()
            .stderr(contains("row 5, column 4: inf is not finite"));
        let clean = write(&dir, "c.dat", DATA);
        let reference = stdout(&variogram_csv(&clean, &[]).success());
        let dropped =
            stdout(&variogram_csv(&data, &["--nan", "missing", "--inf", "missing"]).success());
        assert_eq!(dropped, reference);
        let clamped =
            stdout(&variogram_csv(&data, &["--nan", "missing", "--inf", "clamp"]).success());
        let with_max = write(&dir, "m.dat", &format!("{DATA}7 7 0 3\n"));
        assert_eq!(
            clamped,
            stdout(&variogram_csv(&with_max, &[]).success()),
            "+Inf clamps to the largest value, 3"
        );
    }

    /// A command's outcome: its parsed output, or its error message.
    enum Run<T> {
        Done(T),
        Failed(String),
    }

    impl<T: std::fmt::Debug> Run<T> {
        fn unwrap(self) -> T {
            match self {
                Self::Done(v) => v,
                Self::Failed(e) => panic!("command failed: {e}"),
            }
        }

        fn unwrap_err(self) -> String {
            match self {
                Self::Done(v) => panic!("command succeeded with {v:?}"),
                Self::Failed(e) => e,
            }
        }
    }

    /// Krige one target at the origin from `data`, with extra flags; the
    /// estimate and variance written.
    fn krige_origin(data: &str, extra: &[&str]) -> Run<(f64, f64)> {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", data);
        let targets = write(&dir, "t.dat", "targets\n3\nX\nY\nZ\n0 0 0\n");
        let out = dir.path().join("k.dat");
        let result = bin()
            .args(["krige", "--structure", "sph:1:50"])
            .arg("--data")
            .arg(&data)
            .arg("--targets")
            .arg(&targets)
            .arg("--out")
            .arg(&out)
            .args(extra)
            .output()
            .unwrap();
        if !result.status.success() {
            return Run::Failed(String::from_utf8_lossy(&result.stderr).into_owned());
        }
        let (names, rows) = read(&out);
        assert_eq!(names[3..], ["estimate", "variance"]);
        Run::Done((rows[0][3], rows[0][4]))
    }

    const CROSS: &str = "cross\n4\nX\nY\nZ\nV\n5 0 0 1\n-5 0 0 1\n0 5 0 10\n0 -5 0 10\n";

    #[test]
    fn an_anisotropic_search_selects_along_its_axis() {
        // Isotropic: all four data, equal weights by symmetry: 5.5. A
        // search elongated east-west (azimuth 90, minor radius 1) keeps only
        // the two data on the x axis, both 1.
        let (iso, _) = krige_origin(CROSS, &[]).unwrap();
        assert!((iso - 5.5).abs() < 1e-9, "{iso}");
        let east_west = [
            "--radius",
            "10",
            "--radius-minor",
            "1",
            "--search-angles",
            "90:0:0",
        ];
        let (along_x, _) = krige_origin(CROSS, &east_west).unwrap();
        assert!((along_x - 1.0).abs() < 1e-9, "{along_x}");
        // North-south instead: the two 10s.
        let north_south = [
            "--radius",
            "10",
            "--radius-minor",
            "1",
            "--search-angles",
            "0:0:0",
        ];
        let (along_y, _) = krige_origin(CROSS, &north_south).unwrap();
        assert!((along_y - 10.0).abs() < 1e-9, "{along_y}");
        // An anisotropic search needs a finite radius.
        let err = krige_origin(CROSS, &["--radius-minor", "1"]).unwrap_err();
        assert!(err.contains("finite --radius"), "{err}");
    }

    #[test]
    fn octants_limit_data_per_direction() {
        // Three data east (1, 2, 3) and one west (10). The nearest two are
        // both east; one per octant takes the nearest east and the west.
        let data = "line\n4\nX\nY\nZ\nV\n1 0 0 1\n2 0 0 2\n3 0 0 3\n-4 0 0 10\n";
        let (nearest, _) = krige_origin(data, &["--max-samples", "2"]).unwrap();
        assert!((1.0..=2.0).contains(&nearest), "{nearest}");
        let (octants, _) =
            krige_origin(data, &["--max-samples", "2", "--per-octant", "1"]).unwrap();
        assert!(octants > 2.0, "the west datum enters: {octants}");
    }

    #[test]
    fn too_few_data_fail_or_are_written_as_unest() {
        // Two data required within radius 3 of the origin; none are.
        let flags = ["--min-samples", "2", "--radius", "3"];
        let err = krige_origin(CROSS, &flags).unwrap_err();
        assert!(err.contains("only 0 data within the search"), "{err}");
        let unest = [&flags[..], &["--on-too-few", "unestimated"]].concat();
        assert_eq!(krige_origin(CROSS, &unest).unwrap(), (-999.0, -999.0));
        let custom = [&unest[..], &["--unest", "-99"]].concat();
        assert_eq!(krige_origin(CROSS, &custom).unwrap(), (-99.0, -99.0));
    }

    /// Run sgsim on a 3 x 3 grid of spacing 5 with `data` and extra flags.
    fn sgsim_with(data: &str, extra: &[&str]) -> Run<Vec<f64>> {
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", data);
        let prefix = dir.path().join("s");
        let result = bin()
            .args(["sgsim", "--grid-x", "3:0:5", "--grid-y", "3:0:5"])
            .args(MODEL)
            .arg("--data")
            .arg(&data)
            .arg("--out-prefix")
            .arg(&prefix)
            .args(extra)
            .output()
            .unwrap();
        if !result.status.success() {
            return Run::Failed(String::from_utf8_lossy(&result.stderr).into_owned());
        }
        Run::Done(column(&read(&realization(&prefix, 0)).1, 3))
    }

    #[test]
    fn the_two_part_search_and_its_policies() {
        // Data near nodes 0 (0, 0) and 4 (5, 5): moved there and held.
        let near = "d\n4\nX\nY\nZ\nV\n0.4 0.3 0 1\n5.2 4.9 0 3\n2 2 0 2\n";
        let shared = sgsim_with(near, &["--assign-data"]).unwrap_err();
        assert!(
            shared.contains("both fall on node 1"),
            "(0.4, 0.3) and (2, 2) share node 1: {shared}"
        );
        let held =
            sgsim_with(near, &["--assign-data", "--on-shared-node", "keep-closest"]).unwrap();
        assert!(
            (held[0] - 1.0).abs() < 1e-9 && (held[4] - 3.0).abs() < 1e-9,
            "{held:?}"
        );
        // A datum off the grid.
        let far = "d\n4\nX\nY\nZ\nV\n0 0 0 1\n40 0 0 3\n";
        let outside = sgsim_with(far, &["--assign-data"]).unwrap_err();
        assert!(outside.contains("outside the simulation grid"), "{outside}");
        let clamped = sgsim_with(far, &["--assign-data", "--outside-grid", "clamp"]).unwrap();
        assert!(
            (clamped[2] - 3.0).abs() < 1e-9,
            "clamped to node (10, 0): {clamped:?}"
        );
        // The policies require --assign-data.
        let needs_flag = sgsim_with(far, &["--outside-grid", "clamp"]).unwrap_err();
        assert!(needs_flag.contains("--assign-data"), "{needs_flag}");
        // --per-octant also limits the node search, as GSLIB's noct does:
        // with one datum the data search is unaffected, so a different
        // realization from the same seed can only come from the nodes.
        // Linear tails keep the back-transform from mapping every score to
        // the single datum.
        let one = "d\n4\nX\nY\nZ\nV\n0 0 0 1\n";
        let tails = ["--zmin", "0", "--zmax", "10"];
        let all_nodes = sgsim_with(one, &tails).unwrap();
        let by_octant = sgsim_with(one, &[&tails[..], &["--per-octant", "1"]].concat()).unwrap();
        assert_ne!(all_nodes, by_octant, "octants change the node search");
        // GSLIB's node table runs too.
        assert_eq!(
            sgsim_with(far, &["--node-table", "gslib"]).unwrap().len(),
            9
        );
    }

    #[test]
    fn a_realization_does_not_depend_on_the_others() {
        // RN2 (D-037): realization 2 of a CLI batch equals realization 2
        // computed alone from Seed(7).realization(2), with the CLI's default
        // searches (24 data, unlimited radius, 12 nodes). With one stream
        // shared through the batch, it would depend on realizations 0 and 1.
        use gslib_rs::geometry::Point3;
        use gslib_rs::model::{CovarianceModel, Kind, Structure, VariogramModel};
        use gslib_rs::random::Seed;
        use gslib_rs::sample::Sample;
        use gslib_rs::search::{DataSearch, NodeSearch, SearchEllipsoid, TableExtent};
        use gslib_rs::simulation::{Axis, DataConditioning, Grid, Searches, sisim};

        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let prefix = dir.path().join("sis");
        bin()
            .args(["sisim", "--threshold", "2.5", "--nsim", "3", "--seed", "7"])
            .args(["--grid-x", "4:0:5", "--grid-y", "3:0:5"])
            .args(MODEL)
            .arg("--data")
            .arg(&data)
            .arg("--out-prefix")
            .arg(&prefix)
            .assert()
            .success();
        let batch = column(&read(&realization(&prefix, 2)).1, 3);

        let samples = [(0.0, 0.0, 1.0), (10.0, 0.0, 2.0), (0.0, 10.0, 3.0)]
            .map(|(x, y, v)| Sample::new(Point3::new(x, y, 0.0), v).unwrap());
        let axis = |count, spacing| Axis {
            count,
            first_center: 0.0,
            spacing,
        };
        let grid = Grid::new(axis(4, 5.0), axis(3, 5.0), axis(1, 1.0)).unwrap();
        let model = CovarianceModel::try_from(
            VariogramModel::new(
                0.0,
                vec![Structure::isotropic(Kind::Exponential, 1.0, 30.0).unwrap()],
            )
            .unwrap(),
        )
        .unwrap();
        let search = DataSearch::new(SearchEllipsoid::unlimited(), 0, 24).unwrap();
        let nodes = NodeSearch::new(
            &grid,
            &SearchEllipsoid::unlimited(),
            &model,
            12,
            None,
            TableExtent::Ellipsoid,
        )
        .unwrap();
        let searches = Searches {
            data: DataConditioning::Search(&search),
            nodes: &nodes,
        };
        let alone = sisim(
            &samples,
            &grid,
            2.5,
            &model,
            searches,
            &mut Seed(7).realization(2),
        )
        .unwrap()
        .into_iter()
        .map(|b| if b { 1.0 } else { 0.0 })
        .collect::<Vec<f64>>();
        assert_eq!(batch, alone);
    }

    #[test]
    fn realization_files_do_not_depend_on_thread_count() {
        // P1, P5 end to end: the same batch with 1 and 8 threads writes
        // byte-identical files.
        let dir = TempDir::new().unwrap();
        let data = write(&dir, "d.dat", DATA);
        let run = |threads: &str, name: &str| -> Vec<String> {
            let prefix = dir.path().join(name);
            bin()
                .env("RAYON_NUM_THREADS", threads)
                .args(["sgsim", "--nsim", "6", "--seed", "11"])
                .args(["--grid-x", "6:0:2", "--grid-y", "6:0:2"])
                .args(MODEL)
                .arg("--data")
                .arg(&data)
                .arg("--out-prefix")
                .arg(&prefix)
                .assert()
                .success();
            (0..6)
                .map(|i| std::fs::read_to_string(realization(&prefix, i)).unwrap())
                .collect()
        };
        let one = run("1", "one");
        let eight = run("8", "eight");
        assert_eq!(one, eight);
        assert_ne!(one[0], one[1], "realizations differ from each other");
    }
}
