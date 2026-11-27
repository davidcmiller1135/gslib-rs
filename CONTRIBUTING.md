# Contributing

## Workflow

- `main` is the only long-lived branch. It is always releasable.
- Work on a short-lived branch, open a pull request into `main`, and merge
  once CI is green.
- **Squash-merge** pull requests so `main` has one commit per change, and
  delete the branch afterwards.
- Write PR titles (the squashed commit's subject) in the imperative, with an
  optional area prefix, e.g. `kriging: add block kriging`.

## Checks

CI runs all of these; run them locally first:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy -p gslib_rs --all-targets --no-default-features -- -D warnings
cargo test --workspace
cargo doc --workspace --no-deps
cargo deny check
```

The toolchain is pinned in `rust-toolchain.toml`. To upgrade it, bump the
version there, fix any new lints, and do it in its own PR. The minimum
supported Rust version (`rust-version` in `Cargo.toml`) is tested separately
in CI.

## Code standards

- Library code returns `Result`; it never panics on bad input.
- No `unwrap`/`expect`/`panic!`/`[]` indexing outside tests; no `as` casts
  outside `gslib_rs/src/num.rs`.
- A lint exception must be `#[expect(lint, reason = "...")]` on the smallest
  item that needs it.
- Numerical behaviour that differs from Fortran GSLIB must be documented at
  the definition and listed in `docs/ROADMAP.md` §2.
- Every bug fix comes with a test that failed before the fix.

## Dependencies

Dependabot opens grouped weekly PRs for crates and GitHub Actions. New
dependencies must pass `cargo deny check`: permissive licences only, no
known advisories, crates.io only.
