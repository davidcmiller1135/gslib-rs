//! Schema and discovery for GSLIB reference cases ("fixtures").
//!
//! Specified in `docs/spec/harness.md`: one directory per case,
//! `tests/fixtures/gslib/<program>/<case>/`, holding `case.toml` (R5), the
//! generated `params.par`, inputs, GSLIB's verbatim outputs, and
//! `provenance.toml` (R6). This crate parses and validates those files; it
//! knows nothing about any one program's parameters or outputs.

use std::fmt::Write as _;
use std::path::{Component, Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub mod backtr;
pub mod gamv;
pub mod nscore;
pub mod vmodel;

/// Version of the case and provenance schema this crate reads.
pub const SCHEMA_VERSION: u32 = 1;

/// GSLIB 2.0 programs, as named by the directories of
/// `vendor/gslib-2.0/src/`.
pub const GSLIB_PROGRAMS: [&str; 39] = [
    "addcoord", "anneal", "backtr", "bicalib", "bigaus", "bivplt", "cokb3d", "declus", "draw",
    "ellipsim", "gam", "gamv", "gtsim", "histplt", "histsmth", "ik3d", "kb2d", "kt3d", "locmap",
    "lusim", "nscore", "pfsim", "pixelplt", "postik", "postsim", "probplt", "qpplt", "rotcoord",
    "sasim", "scatplt", "scatsmth", "sgsim", "sisim", "sisim_gs", "sisim_lm", "trans", "vargplt",
    "varmap", "vmodel",
];

/// The files whose contents define the reference recipe (D-014), relative
/// to the workspace root. Directories are walked recursively.
pub const RECIPE_FILES: [&str; 5] = [
    "xtask",
    "oracle-case",
    "vendor/gslib-2.0/MANIFEST.sha256",
    "Cargo.toml",
    "Cargo.lock",
];

/// Compiler flags every reference build must use (`AGENTS.md`, *Reference
/// build*).
pub const REQUIRED_FLAGS: [&str; 2] = ["-fdefault-real-8", "-fdefault-double-8"];

/// Everything that can make a case or fixture invalid.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A file could not be read.
    #[error("{path}: {source}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A file is not valid TOML for its schema.
    #[error("{path}: {source}")]
    Toml {
        /// The file.
        path: PathBuf,
        /// The parse error.
        #[source]
        source: Box<toml::de::Error>,
    },
    /// A file parsed, but its contents break a rule of the spec.
    #[error("{path}: {message}")]
    Invalid {
        /// The file or fixture directory.
        path: PathBuf,
        /// Which rule was broken.
        message: String,
    },
}

/// `Result` specialised to [`Error`].
pub type Result<T, E = Error> = std::result::Result<T, E>;

fn invalid(path: &Path, message: impl Into<String>) -> Error {
    Error::Invalid {
        path: path.to_owned(),
        message: message.into(),
    }
}

/// A reference case (`case.toml`, R5).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// Schema version; must equal [`SCHEMA_VERSION`].
    pub schema: u32,
    /// The GSLIB program, lower case (one of [`GSLIB_PROGRAMS`]).
    pub program: String,
    /// What behaviour this case exercises.
    pub description: String,
    /// Input files, relative to the case directory.
    #[serde(default)]
    pub inputs: Vec<String>,
    /// Program parameters. Typed per program as each gains coverage.
    #[serde(default)]
    pub parameters: toml::Table,
    /// Tolerance overrides, each justified (D-010).
    #[serde(default)]
    pub tolerances: Vec<ToleranceOverride>,
    /// Known, documented deviations from GSLIB (harness.md, *Expected
    /// deviations*).
    #[serde(default)]
    pub expected_deviations: Vec<ExpectedDeviation>,
}

/// A tolerance that replaces a default for one output quantity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ToleranceOverride {
    /// The output quantity (e.g. a column name).
    pub quantity: String,
    /// Maximum relative error.
    pub relative: Option<f64>,
    /// Maximum absolute error.
    pub absolute: Option<f64>,
    /// Why the default does not apply. Required.
    pub justification: String,
}

/// A quantity expected to differ from GSLIB, and why.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExpectedDeviation {
    /// The output quantity.
    pub quantity: String,
    /// The divergence it stems from, as listed in `docs/ROADMAP.md` §2.
    pub divergence: String,
    /// Why the outputs differ.
    pub reason: String,
}

/// How a fixture was produced (`provenance.toml`, R6 / D-009 / D-013).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// Schema version; must equal [`SCHEMA_VERSION`].
    pub schema: u32,
    /// SHA-256 of `vendor/gslib-2.0/MANIFEST.sha256` at generation time.
    pub source_manifest_sha256: String,
    /// The compiler's version line (e.g. from `gfortran --version`).
    pub compiler: String,
    /// Every compiler flag used; must include [`REQUIRED_FLAGS`].
    pub compiler_flags: Vec<String>,
    /// Version of the recipe's output-format rewrite rules.
    pub rewrite_rules_version: u32,
    /// SHA-256 of this fixture's `case.toml`.
    pub case_sha256: String,
    /// Content hash of the recipe that generated the fixture
    /// ([`recipe_sha256`]); survives squash merges (D-014).
    pub recipe_sha256: String,
    /// Git commit the recipe was run from. Informational only: a squash merge
    /// can leave it unreachable from `main` (D-014).
    pub recipe_commit: String,
    /// Generation time, UTC, as `YYYY-MM-DDTHH:MM:SSZ`.
    pub generated_utc: String,
    /// The recipe's full command line.
    pub command_line: Vec<String>,
}

/// Lower-case hexadecimal SHA-256 of `bytes`.
#[must_use]
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .fold(String::with_capacity(64), |mut s, b| {
            // Writing to a String cannot fail.
            let _ = write!(s, "{b:02x}");
            s
        })
}

fn is_hex(s: &str, len: usize) -> bool {
    s.len() == len
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn is_utc_timestamp(s: &str) -> bool {
    // YYYY-MM-DDTHH:MM:SSZ
    let b = s.as_bytes();
    b.len() == 20
        && b.iter().enumerate().all(|(i, &c)| match i {
            4 | 7 => c == b'-',
            10 => c == b'T',
            13 | 16 => c == b':',
            19 => c == b'Z',
            _ => c.is_ascii_digit(),
        })
}

/// A relative path that stays inside its directory.
fn is_contained(path: &str) -> bool {
    !path.is_empty()
        && Path::new(path)
            .components()
            .all(|c| matches!(c, Component::Normal(_)))
}

fn parse_toml<T: for<'de> Deserialize<'de>>(text: &str, path: &Path) -> Result<T> {
    toml::from_str(text).map_err(|source| Error::Toml {
        path: path.to_owned(),
        source: Box::new(source),
    })
}

/// Parse and validate a `case.toml`. `path` is used only in error messages.
///
/// # Errors
/// [`Error::Toml`] for malformed TOML or unknown fields; [`Error::Invalid`]
/// for a wrong schema version, unknown program, empty description, an input
/// path leaving the case directory, or an unjustified or empty tolerance.
pub fn parse_case(text: &str, path: &Path) -> Result<Case> {
    let case: Case = parse_toml(text, path)?;
    if case.schema != SCHEMA_VERSION {
        return Err(invalid(
            path,
            format!("schema {} is not {SCHEMA_VERSION}", case.schema),
        ));
    }
    if !GSLIB_PROGRAMS.contains(&case.program.as_str()) {
        return Err(invalid(
            path,
            format!("unknown GSLIB program {:?}", case.program),
        ));
    }
    if case.description.trim().is_empty() {
        return Err(invalid(path, "description is empty"));
    }
    if let Some(bad) = case.inputs.iter().find(|p| !is_contained(p)) {
        return Err(invalid(
            path,
            format!("input {bad:?} must be a relative path inside the case"),
        ));
    }
    for t in &case.tolerances {
        if t.justification.trim().is_empty() {
            return Err(invalid(
                path,
                format!("tolerance for {:?} has no justification", t.quantity),
            ));
        }
        let bounds = [t.relative, t.absolute];
        if bounds.iter().all(Option::is_none) {
            return Err(invalid(
                path,
                format!("tolerance for {:?} sets no bound", t.quantity),
            ));
        }
        if bounds
            .iter()
            .flatten()
            .any(|v| !(v.is_finite() && *v >= 0.0))
        {
            return Err(invalid(
                path,
                format!("tolerance for {:?} must be finite and >= 0", t.quantity),
            ));
        }
    }
    if let Some(d) = case
        .expected_deviations
        .iter()
        .find(|d| d.divergence.trim().is_empty() || d.reason.trim().is_empty())
    {
        return Err(invalid(
            path,
            format!(
                "expected deviation for {:?} needs a divergence and a reason",
                d.quantity
            ),
        ));
    }
    Ok(case)
}

/// Parse and validate a `provenance.toml`. `path` is used only in error
/// messages.
///
/// # Errors
/// [`Error::Toml`] for malformed TOML, unknown or missing fields;
/// [`Error::Invalid`] for a wrong schema version, malformed hashes, commit or
/// timestamp, an empty compiler or command line, or missing required flags.
pub fn parse_provenance(text: &str, path: &Path) -> Result<Provenance> {
    let p: Provenance = parse_toml(text, path)?;
    if p.schema != SCHEMA_VERSION {
        return Err(invalid(
            path,
            format!("schema {} is not {SCHEMA_VERSION}", p.schema),
        ));
    }
    if ![&p.source_manifest_sha256, &p.case_sha256, &p.recipe_sha256]
        .iter()
        .all(|h| is_hex(h, 64))
    {
        return Err(invalid(
            path,
            "SHA-256 values must be 64 lower-case hex digits",
        ));
    }
    if !is_hex(&p.recipe_commit, 40) {
        return Err(invalid(
            path,
            "recipe_commit must be a full 40-digit commit hash",
        ));
    }
    if !is_utc_timestamp(&p.generated_utc) {
        return Err(invalid(
            path,
            "generated_utc must look like 2026-01-31T12:00:00Z",
        ));
    }
    if p.compiler.trim().is_empty() || p.command_line.is_empty() {
        return Err(invalid(path, "compiler and command_line must not be empty"));
    }
    if let Some(flag) = REQUIRED_FLAGS
        .iter()
        .find(|f| !p.compiler_flags.iter().any(|g| g == *f))
    {
        return Err(invalid(
            path,
            format!("compiler_flags lacks {flag}: the reference build must be double-promoted"),
        ));
    }
    Ok(p)
}

fn collect_files(base: &Path, path: &Path, out: &mut Vec<(String, PathBuf)>) -> Result<()> {
    if path.is_dir() {
        if path.file_name().is_some_and(|n| n == "target") {
            return Ok(());
        }
        for dir in subdirectories(path)? {
            collect_files(base, &dir, out)?;
        }
        let entries = std::fs::read_dir(path).map_err(|source| Error::Io {
            path: path.to_owned(),
            source,
        })?;
        for entry in entries {
            let file = entry
                .map_err(|source| Error::Io {
                    path: path.to_owned(),
                    source,
                })?
                .path();
            if file.is_file() {
                collect_files(base, &file, out)?;
            }
        }
    } else {
        let rel = path
            .strip_prefix(base)
            .map_err(|_| invalid(path, "outside the workspace"))?;
        let name: Vec<&str> = rel
            .components()
            .map(|c| {
                c.as_os_str()
                    .to_str()
                    .ok_or_else(|| invalid(path, "non-UTF-8 file name"))
            })
            .collect::<Result<_>>()?;
        out.push((name.join("/"), path.to_owned()));
    }
    Ok(())
}

/// Content hash of the reference recipe ([`RECIPE_FILES`] under `root`):
/// SHA-256 over each file's relative path, a zero byte, its length (u64,
/// little-endian) and its bytes, in sorted path order (D-014).
///
/// # Errors
/// [`Error::Io`] if a recipe file cannot be read.
pub fn recipe_sha256(root: &Path) -> Result<String> {
    let mut files = Vec::new();
    for entry in RECIPE_FILES {
        collect_files(root, &root.join(entry), &mut files)?;
    }
    files.sort();
    let mut hasher = Sha256::new();
    for (name, path) in files {
        let bytes = read(&path)?;
        hasher.update(name.as_bytes());
        hasher.update([0_u8]);
        hasher.update(u64::try_from(bytes.len()).unwrap_or(u64::MAX).to_le_bytes());
        hasher.update(&bytes);
    }
    Ok(hasher
        .finalize()
        .iter()
        .fold(String::with_capacity(64), |mut s, b| {
            // Writing to a String cannot fail.
            let _ = write!(s, "{b:02x}");
            s
        }))
}

/// A validated fixture directory.
#[derive(Debug, Clone, PartialEq)]
pub struct Fixture {
    /// The case directory.
    pub dir: PathBuf,
    /// Its `case.toml`.
    pub case: Case,
    /// Its `provenance.toml`.
    pub provenance: Provenance,
}

fn read(path: &Path) -> Result<Vec<u8>> {
    std::fs::read(path).map_err(|source| Error::Io {
        path: path.to_owned(),
        source,
    })
}

fn read_text(path: &Path) -> Result<String> {
    String::from_utf8(read(path)?).map_err(|_| invalid(path, "is not valid UTF-8"))
}

/// Load and cross-check one fixture directory, `<root>/<program>/<case>/`.
///
/// # Errors
/// Any error from reading or parsing its files, or [`Error::Invalid`] if
/// `params.par` or an input is missing, the program does not match the
/// directory, or `case.toml`'s hash differs from the provenance record.
pub fn load_fixture(dir: &Path) -> Result<Fixture> {
    let case_path = dir.join("case.toml");
    let case_bytes = read(&case_path)?;
    let case_text = String::from_utf8(case_bytes.clone())
        .map_err(|_| invalid(&case_path, "is not valid UTF-8"))?;
    let case = parse_case(&case_text, &case_path)?;
    let provenance_path = dir.join("provenance.toml");
    let provenance = parse_provenance(&read_text(&provenance_path)?, &provenance_path)?;

    let program_dir = dir
        .parent()
        .and_then(Path::file_name)
        .and_then(|n| n.to_str());
    if program_dir != Some(case.program.as_str()) {
        return Err(invalid(
            dir,
            format!(
                "case is for {:?} but lives under {program_dir:?}",
                case.program
            ),
        ));
    }
    if sha256_hex(&case_bytes) != provenance.case_sha256 {
        return Err(invalid(
            dir,
            "case.toml does not match provenance.case_sha256 (regenerate the fixture)",
        ));
    }
    if !dir.join("params.par").is_file() {
        return Err(invalid(dir, "params.par is missing"));
    }
    if let Some(missing) = case.inputs.iter().find(|i| !dir.join(i).is_file()) {
        return Err(invalid(dir, format!("input {missing:?} is missing")));
    }
    Ok(Fixture {
        dir: dir.to_owned(),
        case,
        provenance,
    })
}

fn subdirectories(dir: &Path) -> Result<Vec<PathBuf>> {
    let entries = std::fs::read_dir(dir).map_err(|source| Error::Io {
        path: dir.to_owned(),
        source,
    })?;
    let mut dirs = entries
        .map(|e| {
            e.map(|entry| entry.path()).map_err(|source| Error::Io {
                path: dir.to_owned(),
                source,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    dirs.retain(|p| p.is_dir());
    dirs.sort();
    Ok(dirs)
}

/// Every fixture under `root` (`<root>/<program>/<case>/`), sorted by path.
/// A missing `root` means there are no fixtures yet.
///
/// # Errors
/// The first invalid fixture, so that one bad fixture is never skipped.
pub fn discover(root: &Path) -> Result<Vec<Fixture>> {
    if !root.exists() {
        return Ok(Vec::new());
    }
    let mut fixtures = Vec::new();
    for program in subdirectories(root)? {
        for case in subdirectories(&program)? {
            fixtures.push(load_fixture(&case)?);
        }
    }
    Ok(fixtures)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CASE: &str = r#"
schema = 1
program = "vmodel"
description = "spherical structure along the azimuth"
inputs = []

[parameters]
nlag = 10

[[tolerances]]
quantity = "gamma"
relative = 1e-10
justification = "GSLIB's truncated pi (ROADMAP §2)"
"#;

    fn provenance(case_sha: &str) -> String {
        format!(
            r#"
schema = 1
source_manifest_sha256 = "{zeros}"
compiler = "GNU Fortran (Ubuntu 13.3.0) 13.3.0"
compiler_flags = ["-O2", "-fdefault-real-8", "-fdefault-double-8"]
rewrite_rules_version = 1
case_sha256 = "{case_sha}"
recipe_sha256 = "{zeros}"
recipe_commit = "{commit}"
generated_utc = "2026-11-27T12:00:00Z"
command_line = ["cargo", "xtask", "oracle", "generate"]
"#,
            zeros = "0".repeat(64),
            commit = "a".repeat(40),
        )
    }

    fn p() -> &'static Path {
        Path::new("case.toml")
    }

    #[test]
    fn sha256_of_known_input() {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn valid_case_parses() {
        let case = parse_case(CASE, p()).unwrap();
        assert_eq!(case.program, "vmodel");
        assert_eq!(case.tolerances.len(), 1);
        assert_eq!(
            case.parameters
                .get("nlag")
                .and_then(toml::Value::as_integer),
            Some(10)
        );
    }

    fn case_error(text: &str) -> String {
        parse_case(text, p()).unwrap_err().to_string()
    }

    #[test]
    fn case_rules_are_enforced() {
        assert!(case_error(&CASE.replace("schema = 1", "schema = 2")).contains("schema 2"));
        assert!(
            case_error(&CASE.replace("\"vmodel\"", "\"vmodle\"")).contains("unknown GSLIB program")
        );
        assert!(
            case_error(&CASE.replace("inputs = []", "inputs = [\"../x.dat\"]"))
                .contains("inside the case")
        );
        assert!(
            case_error(&CASE.replace(
                "justification = \"GSLIB's truncated pi (ROADMAP §2)\"",
                "justification = \" \""
            ))
            .contains("no justification")
        );
        assert!(case_error(&CASE.replace("relative = 1e-10\n", "")).contains("sets no bound"));
        assert!(case_error(&format!("{CASE}\nsurprise = 1\n")).contains("unknown field"));
    }

    #[test]
    fn provenance_rules_are_enforced() {
        let good = provenance(&"b".repeat(64));
        assert!(parse_provenance(&good, p()).is_ok());
        let err = |text: String| parse_provenance(&text, p()).unwrap_err().to_string();
        assert!(err(good.replace(", \"-fdefault-real-8\"", "")).contains("double-promoted"));
        assert!(err(good.replace(&"b".repeat(64), "B".repeat(64).as_str())).contains("SHA-256"));
        assert!(err(good.replace("2026-11-27T12:00:00Z", "yesterday")).contains("generated_utc"));
        assert!(
            err(good.replace("rewrite_rules_version = 1\n", "")).contains("rewrite_rules_version")
        );
    }

    fn write_fixture(
        root: &Path,
        program: &str,
        case_text: &str,
        recorded_sha: Option<&str>,
    ) -> PathBuf {
        let dir = root.join(program).join("case1");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("case.toml"), case_text).unwrap();
        let sha = recorded_sha.map_or_else(|| sha256_hex(case_text.as_bytes()), str::to_owned);
        std::fs::write(dir.join("provenance.toml"), provenance(&sha)).unwrap();
        std::fs::write(dir.join("params.par"), "parameters\n").unwrap();
        dir
    }

    #[test]
    fn recipe_hash_tracks_content() {
        let root = tempfile::tempdir().unwrap();
        for f in ["xtask/src", "oracle-case/src", "vendor/gslib-2.0"] {
            std::fs::create_dir_all(root.path().join(f)).unwrap();
        }
        for f in [
            "xtask/src/main.rs",
            "oracle-case/src/lib.rs",
            "vendor/gslib-2.0/MANIFEST.sha256",
            "Cargo.toml",
            "Cargo.lock",
        ] {
            std::fs::write(root.path().join(f), f).unwrap();
        }
        let before = recipe_sha256(root.path()).unwrap();
        assert_eq!(before, recipe_sha256(root.path()).unwrap(), "deterministic");
        std::fs::create_dir_all(root.path().join("xtask/target")).unwrap();
        std::fs::write(root.path().join("xtask/target/junk"), "x").unwrap();
        assert_eq!(
            before,
            recipe_sha256(root.path()).unwrap(),
            "target/ ignored"
        );
        std::fs::write(root.path().join("xtask/src/main.rs"), "changed").unwrap();
        assert_ne!(
            before,
            recipe_sha256(root.path()).unwrap(),
            "content change detected"
        );
    }

    #[test]
    fn discovers_valid_fixtures() {
        let root = tempfile::tempdir().unwrap();
        write_fixture(root.path(), "vmodel", CASE, None);
        let found = discover(root.path()).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].case.program, "vmodel");
    }

    #[test]
    fn missing_root_means_no_fixtures() {
        assert!(
            discover(Path::new("/definitely/not/here"))
                .unwrap()
                .is_empty(),
            "no fixtures"
        );
    }

    #[test]
    fn edited_case_is_rejected() {
        let root = tempfile::tempdir().unwrap();
        write_fixture(root.path(), "vmodel", CASE, Some(&"c".repeat(64)));
        assert!(
            discover(root.path())
                .unwrap_err()
                .to_string()
                .contains("case_sha256"),
            "hash mismatch"
        );
    }

    #[test]
    fn program_must_match_directory() {
        let root = tempfile::tempdir().unwrap();
        write_fixture(root.path(), "kt3d", CASE, None);
        assert!(
            discover(root.path())
                .unwrap_err()
                .to_string()
                .contains("lives under"),
            "wrong dir"
        );
    }

    #[test]
    fn missing_files_are_rejected() {
        let root = tempfile::tempdir().unwrap();
        let dir = write_fixture(root.path(), "vmodel", CASE, None);
        std::fs::remove_file(dir.join("params.par")).unwrap();
        assert!(
            discover(root.path())
                .unwrap_err()
                .to_string()
                .contains("params.par"),
            "no params"
        );
        let with_input = CASE.replace("inputs = []", "inputs = [\"data.dat\"]");
        let root2 = tempfile::tempdir().unwrap();
        write_fixture(root2.path(), "vmodel", &with_input, None);
        assert!(
            discover(root2.path())
                .unwrap_err()
                .to_string()
                .contains("data.dat"),
            "no input"
        );
    }
}
