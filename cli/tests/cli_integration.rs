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

    let prefix = NamedTempFile::new().unwrap().path().with_file_name("out");
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

    let prefix = NamedTempFile::new().unwrap().path().with_file_name("out");
    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("sisim")
       .arg("--data").arg(data_file.path())
       .arg("--sim").arg(sim_file.path())
       .arg("--threshold").arg("5.0")
       .arg("--nsim").arg("2")
       .arg("--out-prefix").arg(&prefix);
    cmd.assert().success();

    // Header = 1 title + 1 nvar + 5 var_names = 7 lines
    for i in 0..2 {
        let path = prefix.with_file_name(format!("out_{}.gslib", i));
        let text = std::fs::read_to_string(&path).unwrap();
        for line in text.lines().skip(7) {
            let fields: Vec<_> = line.split_whitespace().collect();
            let ind_idx = fields.len() - 1; // indicator is second-to-last
            let ind: f64 = fields[ind_idx].parse().unwrap();
            assert_eq!(ind, 1.0);
        }
    }
}

#[test]
fn sisim_cli_all_zeros() {
    let data_file = make_sample_gslib();
    let sim_file  = make_sample_gslib();

    let prefix = NamedTempFile::new().unwrap().path().with_file_name("out");
    let mut cmd = Command::cargo_bin("cli").unwrap();
    cmd.arg("sisim")
       .arg("--data").arg(data_file.path())
       .arg("--sim").arg(sim_file.path())
       .arg("--threshold").arg("0.0")
       .arg("--nsim").arg("1")
       .arg("--out-prefix").arg(&prefix);
    cmd.assert().success();

    let path = prefix.with_file_name("out_0.gslib");
    let text = std::fs::read_to_string(&path).unwrap();
    for line in text.lines().skip(7) {
        let fields: Vec<_> = line.split_whitespace().collect();
        let ind_idx = fields.len() - 1;
        let ind: f64 = fields[ind_idx].parse().unwrap();
        assert_eq!(ind, 0.0);
    }
}
