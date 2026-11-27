//! Every committed fixture must be complete and consistent (R6).

#[cfg(test)]
mod fixtures {
    use std::path::PathBuf;

    fn fixture_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/gslib")
    }

    /// Fails on the first fixture with a missing file, an invalid schema or
    /// provenance record, or a `case.toml` edited since generation.
    #[test]
    fn all_committed_fixtures_are_valid() {
        let fixtures = oracle_case::discover(&fixture_root()).unwrap();
        for f in &fixtures {
            assert!(
                f.dir.starts_with(fixture_root()),
                "fixture outside the fixture root: {}",
                f.dir.display()
            );
        }
    }
}
