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
}
