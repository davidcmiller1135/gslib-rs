//! Developer tasks for gslib-rs.
//!
//! `cargo xtask oracle verify-source [--archive PATH]` checks the vendored
//! GSLIB 2.0 source against `vendor/gslib-2.0/MANIFEST.sha256` (and,
//! optionally, an original archive against its pinned checksum). This is the
//! first stage of the reference recipe (`docs/spec/harness.md`, R1).

use std::error::Error;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use oracle_case::sha256_hex;

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
