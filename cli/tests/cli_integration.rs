use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

/// Create a minimal 4-column GSLIB file (X, Y, Z, V).
fn make_sample_gslib() -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, "Test Data").unwrap();    // title
    writeln!(f, "4").unwrap();            // nvar
    writeln!(f, "X").unwrap();
    writeln!(f, "Y").unwrap();
    writeln!(f, "Z").unwrap();
    writeln!(f, "V").unwrap();
    writeln!(f, "0.0 0.0 0.0 1.0").unwrap();
    writeln!(f, "1.0 0.0 0.0 2.0").unwrap();
    writeln!(f, "0.0 1.0 0.0 3.0").unwrap();
    f
}

#[test]
fn variogram_cli_produces_csv_and_plot() {
    let data_file = make_sample_gslib();
    let out_csv   = NamedTempFile::new().unwrap();
    let out_png   = NamedTempFile::new().unwrap().path().with_extension("png");

    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("variogram")
       .arg("--input").arg(data_file.path())
       .arg("--max-dist").arg("2.0")
       .arg("--n-lags").arg("2")
       .arg("--out-csv").arg(out_csv.path())
       .arg("--plot").arg(&out_png);
    cmd.assert().success();

    let csv_text = std::fs::read_to_string(out_csv.path()).unwrap();
    let mut lines = csv_text.lines();
    assert_eq!(lines.next().unwrap(), "lag,semivar");
    assert!(lines.next().unwrap().contains(','));

    let meta = std::fs::metadata(&out_png).unwrap();
    assert!(meta.len() > 0);
}

#[test]
fn sgsim_cli_writes_gslib_files() {
    let data_file = make_sample_gslib();
    let sim_file  = make_sample_gslib();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("sgsim")
       .arg("--data").arg(data_file.path())
       .arg("--sim").arg(sim_file.path())
       .arg("--nsim").arg("2")
       .arg("--seed").arg("123")
       .arg("--out-prefix").arg(&prefix);
    cmd.assert().success();

    for i in 0..2 {
        let path = prefix.with_file_name(format!("out_{}.gslib", i));
        let meta = std::fs::metadata(&path).expect("missing output gslib");
        assert!(meta.len() > 0, "GSLIB file {} empty", path.display());
    }
}

#[test]
fn sisim_cli_all_ones() {
    let data_file = make_sample_gslib();
    let sim_file  = make_sample_gslib();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("sisim")
       .arg("--data").arg(data_file.path())
       .arg("--sim").arg(sim_file.path())
       .arg("--threshold").arg("5.0")
       .arg("--nsim").arg("2")
       .arg("--out-prefix").arg(&prefix);
    cmd.assert().success();

    for i in 0..2 {
        let path = prefix.with_file_name(format!("out_{}.gslib", i));
        assert_eq!(last_column(&path, &format!("ind_{i}")), vec![1.0; 3]);
    }
}

#[test]
fn sisim_cli_all_zeros() {
    let data_file = make_sample_gslib();
    let sim_file  = make_sample_gslib();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("sisim")
       .arg("--data").arg(data_file.path())
       .arg("--sim").arg(sim_file.path())
       .arg("--threshold").arg("0.0")
       .arg("--nsim").arg("1")
       .arg("--out-prefix").arg(&prefix);
    cmd.assert().success();

    let path = prefix.with_file_name("out_0.gslib");
    assert_eq!(last_column(&path, "ind_0"), vec![0.0; 3]);
}

#[test]
fn sgsim_cli_keeps_sim_columns_and_appends_result() {
    let data_file = make_sample_gslib();
    let sim_file = make_sample_gslib();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("sgs_cols");
    Command::cargo_bin("cli")
        .unwrap()
        .arg("sgsim")
        .arg("--data").arg(data_file.path())
        .arg("--sim").arg(sim_file.path())
        .arg("--seed").arg("1")
        .arg("--out-prefix").arg(&prefix)
        .assert()
        .success();

    let path = prefix.with_file_name("sgs_cols_0.gslib");
    let text = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[1], "5");
    assert_eq!(&lines[2..7], ["X", "Y", "Z", "V", "sim_0"]);
    // The sim file's own V column (1, 2, 3) must be carried through unchanged.
    let v: Vec<f64> = lines[7..]
        .iter()
        .map(|l| l.split_whitespace().nth(3).unwrap().parse().unwrap())
        .collect();
    assert_eq!(v, vec![1.0, 2.0, 3.0]);
}

#[test]
fn sgsim_cli_rejects_data_without_value_column() {
    let mut data_file = NamedTempFile::new().unwrap();
    writeln!(data_file, "Coords only\n3\nX\nY\nZ\n0 0 0").unwrap();
    let sim_file = make_sample_gslib();
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("bad");

    Command::cargo_bin("cli")
        .unwrap()
        .arg("sgsim")
        .arg("--data").arg(data_file.path())
        .arg("--sim").arg(sim_file.path())
        .arg("--out-prefix").arg(&prefix)
        .assert()
        .failure()
        .stderr(predicates::str::contains("expected at least 4"));
}

/// Read a GSLIB output file, check its last variable is named `name`, and
/// return that column's values.
fn last_column(path: &std::path::Path, name: &str) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    lines.next(); // title
    let nvar: usize = lines.next().unwrap().trim().parse().unwrap();
    let names: Vec<&str> = lines.by_ref().take(nvar).collect();
    assert_eq!(names.last().copied(), Some(name), "last variable name");
    lines
        .map(|l| {
            let fields: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(fields.len(), nvar);
            fields[nvar - 1].parse().unwrap()
        })
        .collect()
}
