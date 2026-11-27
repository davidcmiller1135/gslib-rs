//! Developer tasks for gslib-rs: the Fortran GSLIB reference recipe
//! (`docs/spec/harness.md`, R1).
//!
//! - `cargo xtask oracle verify-source [--archive PATH]` checks the vendored
//!   GSLIB 2.0 source against `vendor/gslib-2.0/MANIFEST.sha256` (and,
//!   optionally, an original archive against its pinned checksum).
//! - `cargo xtask oracle build [--out DIR]` copies the source, widens its
//!   output formats ([`rewrite`]), and builds every program with double
//!   promotion into `DIR` (default `target/oracle`).
//! - `cargo xtask oracle generate [--program NAME]` builds, then runs GSLIB on
//!   every case in `tests/fixtures/gslib/` and writes its fixture
//!   ([`generate`]).

mod generate;
mod rewrite;

use std::error::Error;
use std::fmt::Write as _;
use std::io::Write as IoWrite;
use std::path::{Path, PathBuf};
use std::process::{Command as Process, ExitCode};

use clap::{Parser, Subcommand};
use oracle_case::sha256_hex;

/// Compiler flags for the reference build.
///
/// - `-fdefault-real-8 -fdefault-double-8`: double promotion (`AGENTS.md`).
/// - `-std=legacy`: accept Fortran 77 idioms without errors.
/// - `-fno-automatic`: static local storage, as F77 compilers of the era
///   provided and as old code may rely on.
/// - `-O0 -ffp-contract=off`: no reordering and no fused multiply-add, so
///   results depend only on the source.
/// - `-ffixed-line-length-none`: widened formats may pass column 72 (no
///   vendored code line does).
const FLAGS: [&str; 7] = [
    "-fdefault-real-8",
    "-fdefault-double-8",
    "-std=legacy",
    "-fno-automatic",
    "-O0",
    "-ffp-contract=off",
    "-ffixed-line-length-none",
];

/// SHA-256 of the GSLIB 2.0 archive the source was vendored from (F-003).
const ARCHIVE_SHA256: &str = "f8a4bd8fe2b90b26f6b3f33ccf795f9d73cb2560c524128c7865b050a5f763d1";

type AppResult<T = ()> = Result<T, Box<dyn Error>>;

/// gslib-rs developer tasks.
#[derive(Debug, Parser)]
#[command(name = "xtask")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// The Fortran GSLIB reference recipe.
    #[command(subcommand)]
    Oracle(Oracle),
}

#[derive(Debug, Subcommand)]
enum Oracle {
    /// Verify the vendored GSLIB 2.0 source against its manifest.
    VerifySource {
        /// Also check an original GSLIB 2.0 archive against its pinned SHA-256.
        #[arg(long)]
        archive: Option<PathBuf>,
    },
    /// Build every GSLIB program with double promotion and widened formats.
    Build {
        /// Output directory (default: `target/oracle` in the workspace).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Build, then regenerate the fixtures under `tests/fixtures/gslib/`.
    Generate {
        /// Only this program's cases.
        #[arg(long)]
        program: Option<String>,
        /// Build directory (default: `target/oracle` in the workspace).
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            writeln!(std::io::stderr(), "error: {e}").ok();
            ExitCode::FAILURE
        }
    }
}

fn workspace_root() -> AppResult<PathBuf> {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .map(Path::to_path_buf)
        .ok_or_else(|| "cannot locate the workspace root".into())
}

fn run(cli: Cli) -> AppResult {
    match cli.command {
        Command::Oracle(Oracle::VerifySource { archive }) => {
            let vendor = workspace_root()?.join("vendor/gslib-2.0");
            let manifest = std::fs::read_to_string(vendor.join("MANIFEST.sha256"))?;
            let count = verify_tree(&vendor.join("src"), &manifest)?;
            writeln!(
                std::io::stdout(),
                "vendored GSLIB 2.0: {count} files match MANIFEST.sha256"
            )?;
            if let Some(path) = archive {
                verify_archive(&std::fs::read(&path)?)?;
                writeln!(
                    std::io::stdout(),
                    "{}: matches the pinned SHA-256",
                    path.display()
                )?;
            }
            Ok(())
        }
        Command::Oracle(Oracle::Build { out }) => {
            let root = workspace_root()?;
            build(&root, &out.unwrap_or_else(|| root.join("target/oracle"))).map(drop)
        }
        Command::Oracle(Oracle::Generate { program, out }) => {
            let root = workspace_root()?;
            let out = out.unwrap_or_else(|| root.join("target/oracle"));
            generate::generate(&root, &out, program.as_deref())
        }
    }
}

/// Run `gfortran` (or `ar`) and return its combined output on failure.
fn run_tool(program: &str, args: &[String]) -> Result<(), String> {
    let output = Process::new(program)
        .args(args)
        .output()
        .map_err(|e| format!("cannot run {program}: {e}"))?;
    if output.status.success() {
        Ok(())
    } else {
        let mut text = String::from_utf8_lossy(&output.stderr).into_owned();
        text.push_str(&String::from_utf8_lossy(&output.stdout));
        Err(text)
    }
}

fn path_arg(p: &Path) -> AppResult<String> {
    Ok(p.to_str().ok_or("non-UTF-8 path")?.to_owned())
}

/// The first `error`-looking line of compiler output, for summaries.
fn first_error(text: &str) -> String {
    text.lines()
        .find(|l| l.contains("Error") || l.contains("error"))
        .unwrap_or_else(|| text.lines().next().unwrap_or("unknown failure"))
        .trim()
        .to_owned()
}

fn toml_string(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

/// What a build produced, for `build-info.toml` and the summary.
struct BuildOutcome {
    compiler: String,
    manifest_sha256: String,
    widened: usize,
    internal_skipped: usize,
    built: Vec<String>,
    failed: Vec<(String, String)>,
}

/// The reference build (R1 steps 2 and 3).
fn build(root: &Path, out: &Path) -> AppResult<BuildOutcome> {
    let vendor = root.join("vendor/gslib-2.0");
    let manifest = std::fs::read(vendor.join("MANIFEST.sha256"))?;
    let src = vendor.join("src");
    verify_tree(&src, std::str::from_utf8(&manifest)?)?;
    let compiler = compiler_version()?;

    if out.exists() {
        std::fs::remove_dir_all(out)?;
    }
    let work = out.join("src");
    let (widened, internal_skipped) = prepare_sources(&src, &work)?;
    let library = build_library(&work, out)?;
    let (built, failed) = build_programs(&work, out, &library)?;

    let outcome = BuildOutcome {
        compiler,
        manifest_sha256: sha256_hex(&manifest),
        widened,
        internal_skipped,
        built,
        failed,
    };
    std::fs::write(out.join("build-info.toml"), build_info(&outcome)?)?;
    report(&outcome, out)?;
    Ok(outcome)
}

fn compiler_version() -> AppResult<String> {
    let output = Process::new("gfortran")
        .arg("--version")
        .output()
        .map_err(|e| {
            format!(
                "gfortran not found ({e}); install it (apt-get update && apt-get install gfortran)"
            )
        })?;
    Ok(String::from_utf8_lossy(&output.stdout)
        .lines()
        .next()
        .unwrap_or("")
        .trim()
        .to_owned())
}

/// Copy the source with lower-case names (programs `include 'kt3d.inc'`),
/// widening output formats in `.for` files. Returns the rewrite totals.
fn prepare_sources(src: &Path, work: &Path) -> AppResult<(usize, usize)> {
    let mut widened = 0;
    let mut internal_skipped = 0;
    for rel in list_files(src)? {
        let lower = rel.to_ascii_lowercase();
        let dest = work.join(&lower);
        if let Some(parent) = dest.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let bytes = std::fs::read(src.join(&rel))?;
        let is_fortran = Path::new(&lower)
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("for"));
        if is_fortran {
            let (text, r) = rewrite::widen_formats(&String::from_utf8_lossy(&bytes));
            widened += r.widened;
            internal_skipped += r.internal_skipped;
            std::fs::write(&dest, text)?;
        } else {
            std::fs::write(&dest, bytes)?;
        }
    }
    Ok((widened, internal_skipped))
}

fn flags() -> Vec<String> {
    FLAGS.iter().map(|f| (*f).to_owned()).collect()
}

/// Compile `gslib/*.for` into `libgslib.a`.
fn build_library(work: &Path, out: &Path) -> AppResult<PathBuf> {
    let obj = out.join("obj");
    std::fs::create_dir_all(&obj)?;
    let mut objects = Vec::new();
    for rel in list_files(&work.join("gslib"))? {
        let Some(stem) = rel.strip_suffix(".for") else {
            continue;
        };
        let object = obj.join(format!("{stem}.o"));
        let mut args = flags();
        args.extend([
            "-c".to_owned(),
            path_arg(&work.join("gslib").join(&rel))?,
            "-o".to_owned(),
            path_arg(&object)?,
        ]);
        run_tool("gfortran", &args).map_err(|e| format!("gslib/{rel}: {}", first_error(&e)))?;
        objects.push(path_arg(&object)?);
    }
    let library = out.join("libgslib.a");
    let mut ar_args = vec!["rcs".to_owned(), path_arg(&library)?];
    ar_args.extend(objects);
    run_tool("ar", &ar_args).map_err(|e| format!("ar: {}", first_error(&e)))?;
    Ok(library)
}

/// Programs built, and programs that failed with their first error line.
type ProgramResults = (Vec<String>, Vec<(String, String)>);

/// Build every `<name>/<name>.for` against the library.
fn build_programs(work: &Path, out: &Path, library: &Path) -> AppResult<ProgramResults> {
    let bin = out.join("bin");
    std::fs::create_dir_all(&bin)?;
    let mut built = Vec::new();
    let mut failed = Vec::new();
    for entry in std::fs::read_dir(work)? {
        let dir = entry?.path();
        let Some(name) = dir.file_name().and_then(|n| n.to_str()).map(str::to_owned) else {
            continue;
        };
        let source = dir.join(format!("{name}.for"));
        if !source.is_file() {
            continue;
        }
        let mut args = flags();
        args.extend([
            format!("-I{}", path_arg(&dir)?),
            format!("-J{}", path_arg(&out.join("obj"))?),
            path_arg(&source)?,
            path_arg(library)?,
            "-o".to_owned(),
            path_arg(&bin.join(&name))?,
        ]);
        match run_tool("gfortran", &args) {
            Ok(()) => built.push(name),
            Err(e) => failed.push((name, first_error(&e))),
        }
    }
    built.sort();
    failed.sort();
    Ok((built, failed))
}

fn build_info(o: &BuildOutcome) -> AppResult<String> {
    let mut info = String::new();
    writeln!(info, "# Written by `cargo xtask oracle build`.")?;
    writeln!(
        info,
        "source_manifest_sha256 = {}",
        toml_string(&o.manifest_sha256)
    )?;
    writeln!(info, "compiler = {}", toml_string(&o.compiler))?;
    let quoted: Vec<String> = FLAGS.iter().map(|f| toml_string(f)).collect();
    writeln!(info, "compiler_flags = [{}]", quoted.join(", "))?;
    writeln!(info, "rewrite_rules_version = {}", rewrite::RULES_VERSION)?;
    writeln!(info, "formats_widened = {}", o.widened)?;
    writeln!(info, "internal_writes_skipped = {}", o.internal_skipped)?;
    let names: Vec<String> = o.built.iter().map(|n| toml_string(n)).collect();
    writeln!(info, "built = [{}]", names.join(", "))?;
    for (name, err) in &o.failed {
        writeln!(info, "\n[failed.{name}]\nerror = {}", toml_string(err))?;
    }
    Ok(info)
}

fn report(o: &BuildOutcome, out: &Path) -> AppResult {
    let mut stdout = std::io::stdout().lock();
    writeln!(stdout, "compiler: {}", o.compiler)?;
    writeln!(
        stdout,
        "formats widened: {}; internal writes left alone: {}",
        o.widened, o.internal_skipped
    )?;
    writeln!(
        stdout,
        "built {} programs into {}",
        o.built.len(),
        out.join("bin").display()
    )?;
    for (name, err) in &o.failed {
        writeln!(stdout, "FAILED {name}: {err}")?;
    }
    if o.failed.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} program(s) failed to build; see {}",
            o.failed.len(),
            out.join("build-info.toml").display()
        )
        .into())
    }
}

fn verify_archive(bytes: &[u8]) -> AppResult {
    let actual = sha256_hex(bytes);
    if actual == ARCHIVE_SHA256 {
        Ok(())
    } else {
        Err(format!("archive SHA-256 is {actual}, expected {ARCHIVE_SHA256}").into())
    }
}

/// Every regular file under `dir`, as paths relative to it, sorted.
fn list_files(dir: &Path) -> AppResult<Vec<String>> {
    fn walk(base: &Path, dir: &Path, out: &mut Vec<String>) -> AppResult {
        for entry in std::fs::read_dir(dir)? {
            let path = entry?.path();
            if path.is_dir() {
                walk(base, &path, out)?;
            } else {
                let rel = path.strip_prefix(base)?;
                out.push(rel.to_str().ok_or("non-UTF-8 file name")?.to_owned());
            }
        }
        Ok(())
    }
    let mut out = Vec::new();
    walk(dir, dir, &mut out)?;
    out.sort();
    Ok(out)
}

/// Check that `dir` holds exactly the files in `manifest` (`sha256sum`
/// format, paths starting `./`) with matching hashes. Returns the file count.
fn verify_tree(dir: &Path, manifest: &str) -> AppResult<usize> {
    let mut expected: Vec<(String, String)> = manifest
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let (hash, path) = line
                .split_once("  ")
                .ok_or_else(|| format!("bad manifest line: {line:?}"))?;
            Ok((path.trim_start_matches("./").to_owned(), hash.to_owned()))
        })
        .collect::<AppResult<_>>()?;
    expected.sort();

    let actual = list_files(dir)?;
    let listed: Vec<&str> = expected.iter().map(|(p, _)| p.as_str()).collect();
    if actual.iter().map(String::as_str).ne(listed.iter().copied()) {
        let extra: Vec<&String> = actual
            .iter()
            .filter(|a| !listed.contains(&a.as_str()))
            .collect();
        let missing: Vec<&&str> = listed
            .iter()
            .filter(|l| !actual.iter().any(|a| a == **l))
            .collect();
        return Err(format!(
            "file set differs from manifest; extra: {extra:?}, missing: {missing:?}"
        )
        .into());
    }
    for (path, hash) in &expected {
        let actual_hash = sha256_hex(&std::fs::read(dir.join(path))?);
        if &actual_hash != hash {
            return Err(
                format!("{path}: SHA-256 {actual_hash} does not match the manifest").into(),
            );
        }
    }
    Ok(expected.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree() -> (tempfile::TempDir, String) {
        let dir = tempfile::tempdir().unwrap();
        std::fs::create_dir_all(dir.path().join("GSLIB")).unwrap();
        std::fs::write(dir.path().join("README"), b"readme\r\n").unwrap();
        std::fs::write(dir.path().join("GSLIB/COVA3.FOR"), b"code\r\n").unwrap();
        let manifest = format!(
            "{}  ./GSLIB/COVA3.FOR\n{}  ./README\n",
            sha256_hex(b"code\r\n"),
            sha256_hex(b"readme\r\n")
        );
        (dir, manifest)
    }

    #[test]
    fn matching_tree_verifies() {
        let (dir, manifest) = tree();
        assert_eq!(verify_tree(dir.path(), &manifest).unwrap(), 2);
    }

    #[test]
    fn modified_file_is_detected() {
        let (dir, manifest) = tree();
        // Converting CRLF to LF is exactly the corruption to catch.
        std::fs::write(dir.path().join("README"), b"readme\n").unwrap();
        assert!(
            verify_tree(dir.path(), &manifest)
                .unwrap_err()
                .to_string()
                .contains("README"),
            "edit"
        );
    }

    #[test]
    fn added_or_removed_files_are_detected() {
        let (dir, manifest) = tree();
        std::fs::write(dir.path().join("EXTRA.FOR"), b"x").unwrap();
        assert!(
            verify_tree(dir.path(), &manifest)
                .unwrap_err()
                .to_string()
                .contains("EXTRA.FOR"),
            "added"
        );
        let (dir2, manifest2) = tree();
        std::fs::remove_file(dir2.path().join("README")).unwrap();
        assert!(
            verify_tree(dir2.path(), &manifest2)
                .unwrap_err()
                .to_string()
                .contains("README"),
            "removed"
        );
    }

    #[test]
    fn archive_checksum_is_pinned() {
        assert!(
            verify_archive(b"not the archive").is_err(),
            "wrong archive accepted"
        );
    }

    #[test]
    fn the_real_vendored_source_verifies() {
        let vendor = workspace_root().unwrap().join("vendor/gslib-2.0");
        let manifest = std::fs::read_to_string(vendor.join("MANIFEST.sha256")).unwrap();
        assert_eq!(verify_tree(&vendor.join("src"), &manifest).unwrap(), 145);
    }
}
