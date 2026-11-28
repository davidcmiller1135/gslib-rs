//! `cargo xtask oracle generate`: run GSLIB on every reference case and
//! write its fixture (`docs/spec/harness.md`, R1 step 4, R6).

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

use oracle_case::{backtr, gamv, kt3d, nscore, parse_case, sha256_hex, vmodel};

use super::{AppResult, BuildOutcome, FLAGS, rewrite};

/// Paths the recipe depends on: generation refuses to record a commit if
/// any of them has uncommitted changes.
const RECIPE_PATHS: [&str; 5] = ["xtask", "oracle-case", "vendor", "Cargo.toml", "Cargo.lock"];

/// What one program needs to run a case: its parameter file and outputs.
struct Run {
    par: String,
    outputs: Vec<&'static str>,
}

/// Fail unless every file the parameters name is a declared case input.
fn require_inputs(files: &[&str], inputs: &[String]) -> AppResult {
    if let Some(f) = files.iter().find(|f| !inputs.iter().any(|i| i == *f)) {
        return Err(format!("{f} is used by the parameters but not listed in inputs").into());
    }
    Ok(())
}

fn prepare(program: &str, parameters: &toml::Table, inputs: &[String]) -> AppResult<Run> {
    match program {
        "vmodel" => Ok(Run {
            par: vmodel::Params::from_table(parameters)?.to_par(),
            outputs: vec![vmodel::OUTPUT_FILE],
        }),
        "nscore" => {
            let p = nscore::Params::from_table(parameters)?;
            require_inputs(&[&p.data], inputs)?;
            Ok(Run {
                par: p.to_par(),
                outputs: vec![nscore::OUTPUT_FILE, nscore::TABLE_FILE],
            })
        }
        "backtr" => {
            let p = backtr::Params::from_table(parameters)?;
            require_inputs(&[&p.data, &p.table], inputs)?;
            Ok(Run {
                par: p.to_par(),
                outputs: vec![backtr::OUTPUT_FILE],
            })
        }
        "gamv" => {
            let p = gamv::Params::from_table(parameters)?;
            require_inputs(&[&p.data], inputs)?;
            Ok(Run {
                par: p.to_par(),
                outputs: vec![gamv::OUTPUT_FILE],
            })
        }
        "kt3d" => {
            let p = kt3d::Params::from_table(parameters)?;
            require_inputs(&[&p.data], inputs)?;
            Ok(Run {
                par: p.to_par(),
                outputs: vec![kt3d::OUTPUT_FILE, kt3d::DEBUG_FILE],
            })
        }
        other => Err(format!("no fixture generator for {other} yet").into()),
    }
}

fn git(root: &Path, args: &[&str]) -> AppResult<String> {
    let out = Command::new("git").current_dir(root).args(args).output()?;
    if !out.status.success() {
        return Err(format!(
            "git {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr)
        )
        .into());
    }
    Ok(String::from_utf8(out.stdout)?.trim().to_owned())
}

/// The commit to record, refusing if the recipe has uncommitted changes.
fn recipe_commit(root: &Path) -> AppResult<String> {
    let mut args = vec!["status", "--porcelain", "--"];
    args.extend(RECIPE_PATHS);
    let dirty = git(root, &args)?;
    if !dirty.is_empty() {
        return Err(
            format!("commit the recipe before generating fixtures; uncommitted:\n{dirty}").into(),
        );
    }
    git(root, &["rev-parse", "HEAD"])
}

/// `YYYY-MM-DDTHH:MM:SSZ` for seconds since the Unix epoch.
pub(crate) fn utc_timestamp(secs: i64) -> String {
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    // Civil-from-days (H. Hinnant's algorithm), proleptic Gregorian.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe.div_euclid(1_460) + doe.div_euclid(36_524) - doe.div_euclid(146_096))
        .div_euclid(365);
    let doy = doe - (365 * yoe + yoe.div_euclid(4) - yoe.div_euclid(100));
    let mp = (5 * doy + 2).div_euclid(153);
    let day = doy - (153 * mp + 2).div_euclid(5) + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        rem.div_euclid(3_600),
        rem.rem_euclid(3_600).div_euclid(60),
        rem.rem_euclid(60)
    )
}

fn now_utc() -> AppResult<String> {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)?
        .as_secs();
    Ok(utc_timestamp(i64::try_from(secs)?))
}

fn quote(s: &str) -> String {
    super::toml_string(s)
}

fn provenance(build: &BuildOutcome, case_sha256: &str, recipe: &Recipe) -> AppResult<String> {
    let mut p = String::new();
    writeln!(
        p,
        "# Written by `cargo xtask oracle generate`. Do not edit."
    )?;
    writeln!(p, "schema = {}", oracle_case::SCHEMA_VERSION)?;
    writeln!(
        p,
        "source_manifest_sha256 = {}",
        quote(&build.manifest_sha256)
    )?;
    writeln!(p, "compiler = {}", quote(&build.compiler))?;
    let flags: Vec<String> = FLAGS.iter().map(|f| quote(f)).collect();
    writeln!(p, "compiler_flags = [{}]", flags.join(", "))?;
    writeln!(p, "rewrite_rules_version = {}", rewrite::RULES_VERSION)?;
    writeln!(p, "case_sha256 = {}", quote(case_sha256))?;
    writeln!(p, "recipe_sha256 = {}", quote(&recipe.sha256))?;
    writeln!(p, "recipe_commit = {}", quote(&recipe.commit))?;
    writeln!(p, "generated_utc = {}", quote(&now_utc()?))?;
    let args: Vec<String> = std::env::args().map(|a| quote(&a)).collect();
    writeln!(p, "command_line = [{}]", args.join(", "))?;
    Ok(p)
}

/// Run one case and write its fixture files.
/// The recipe's identity: its content hash and the commit it ran from.
struct Recipe {
    sha256: String,
    commit: String,
}

fn generate_case(
    case_dir: &Path,
    program: &str,
    out: &Path,
    build: &BuildOutcome,
    recipe: &Recipe,
) -> AppResult {
    let case_path = case_dir.join("case.toml");
    let case_bytes = std::fs::read(&case_path)?;
    let case = parse_case(std::str::from_utf8(&case_bytes)?, &case_path)?;
    if case.program != program {
        return Err(format!(
            "{}: case is for {} but lives under {program}",
            case_dir.display(),
            case.program
        )
        .into());
    }
    let run = prepare(program, &case.parameters, &case.inputs)
        .map_err(|e| format!("{}: {e}", case_path.display()))?;

    let name = case_dir
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or("bad case directory name")?;
    let work = out.join("run").join(program).join(name);
    if work.exists() {
        std::fs::remove_dir_all(&work)?;
    }
    std::fs::create_dir_all(&work)?;
    for input in &case.inputs {
        // An input missing from the case is taken from the synthetic suites
        // and copied into the case, so the fixture holds the exact bytes used.
        let local = case_dir.join(input);
        if !local.is_file() {
            let suite = case_dir.join("../../../data").join(input);
            if !suite.is_file() {
                return Err(format!(
                    "{}: input {input} is neither in the case nor a data suite",
                    case_dir.display()
                )
                .into());
            }
            std::fs::copy(&suite, &local)?;
        }
        std::fs::copy(&local, work.join(input))?;
    }
    std::fs::write(work.join("params.par"), &run.par)?;

    let mut child = Command::new(out.join("bin").join(program))
        .current_dir(&work)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    child
        .stdin
        .take()
        .ok_or("no stdin")?
        .write_all(b"params.par\n")?;
    let result = child.wait_with_output()?;
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::write(work.join("run.log"), &log)?;
    // GSLIB reports parameter errors with STOP, which exits 0.
    if !result.status.success() || log.contains("ERROR") || log.contains("Error") {
        return Err(format!(
            "{program} failed on {}; see {}",
            case_dir.display(),
            work.join("run.log").display()
        )
        .into());
    }

    std::fs::write(case_dir.join("params.par"), &run.par)?;
    for output in &run.outputs {
        let produced = work.join(output);
        if !produced.is_file() {
            return Err(format!(
                "{program} did not write {output} for {}",
                case_dir.display()
            )
            .into());
        }
        std::fs::copy(&produced, case_dir.join(output))?;
    }
    let provenance_path = case_dir.join("provenance.toml");
    std::fs::write(
        &provenance_path,
        provenance(build, &sha256_hex(&case_bytes), recipe)?,
    )?;
    // The fixture must now pass the same checks the harness applies.
    oracle_case::load_fixture(case_dir)?;
    Ok(())
}

fn subdirs(dir: &Path) -> AppResult<Vec<std::path::PathBuf>> {
    let mut dirs = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_dir() {
            dirs.push(path);
        }
    }
    dirs.sort();
    Ok(dirs)
}

/// Generate fixtures for every case under `tests/fixtures/gslib/`, or only
/// for `program`.
pub(crate) fn generate(root: &Path, out: &Path, program: Option<&str>) -> AppResult {
    let recipe = Recipe {
        commit: recipe_commit(root)?,
        sha256: oracle_case::recipe_sha256(root)?,
    };
    let build = super::build(root, out)?;
    let fixtures = root.join("tests/fixtures/gslib");
    let mut count = 0;
    for program_dir in subdirs(&fixtures)? {
        let name = program_dir
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or("bad program directory")?;
        if program.is_some_and(|p| p != name) {
            continue;
        }
        for case_dir in subdirs(&program_dir)? {
            generate_case(&case_dir, name, out, &build, &recipe)?;
            count += 1;
            writeln!(std::io::stdout(), "generated {}", case_dir.display())?;
        }
    }
    writeln!(
        std::io::stdout(),
        "{count} fixture(s) generated; recipe {} (commit {})",
        recipe.sha256,
        recipe.commit
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timestamps() {
        assert_eq!(utc_timestamp(0), "1970-01-01T00:00:00Z");
        // 2000-02-29T12:34:56Z (leap day), checked with `date -u -d @951827696`.
        assert_eq!(utc_timestamp(951_827_696), "2000-02-29T12:34:56Z");
        // 2026-11-27T16:00:00Z, checked with `date -u -d @1795795200`.
        assert_eq!(utc_timestamp(1_795_795_200), "2026-11-27T16:00:00Z");
    }

    #[test]
    fn only_known_programs_have_generators() {
        assert!(
            prepare("sgsim", &toml::Table::new(), &[]).is_err(),
            "sgsim has no generator yet"
        );
        let t: toml::Table = "data = \"a.dat\"\nvariable_column = 4".parse().unwrap();
        assert!(
            prepare("nscore", &t, &[]).is_err(),
            "data must be a declared input"
        );
        assert!(
            prepare("nscore", &t, &["a.dat".to_owned()]).is_ok(),
            "declared input"
        );
    }
}
