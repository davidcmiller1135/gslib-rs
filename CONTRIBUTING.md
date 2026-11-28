# Contributing

The engineering rules for this project, covering precedence, GSLIB parity,
boundaries and testing, are in [`AGENTS.md`](AGENTS.md). They apply to humans
and agents alike. This file covers the mechanics.

## Workflow

- `main` is the only long-lived branch. It is always releasable.
- Work on a short-lived branch, open a pull request into `main`, and merge
  once every check below passes locally. There is no CI (D-053); the PR
  says which checks were run.
- **Squash-merge** pull requests so `main` has one commit per change, and
  delete the branch afterwards.
- Write PR titles (the squashed commit's subject) in the imperative, with an
  optional area prefix, e.g. `kriging: add block kriging`.

## Checks

There is no CI. Run all of these locally before merging (the same list as
`AGENTS.md`):

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy -p gslib_rs --all-targets --no-default-features -- -D warnings
cargo test --workspace
cargo test -p gslib_rs --no-default-features
cargo +1.89 test --workspace       # MSRV (rust-version in Cargo.toml)
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo deny check                   # cargo install cargo-deny --locked
cargo xtask oracle verify-source   # vendored GSLIB 2.0 is untouched
cargo xtask gen-data --check       # synthetic data suites match their generator
cargo xtask oracle build           # builds the Fortran reference (needs gfortran)
```

The toolchain is pinned in `rust-toolchain.toml`. To upgrade it, bump the
version there, fix any new lints, and do it in its own PR. The minimum
supported Rust version (`rust-version` in `Cargo.toml`) is checked with
`cargo +1.89 test --workspace`, above.

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

Dependabot opens a grouped weekly PR for crate updates; run the checks
locally before merging it. New
dependencies must pass `cargo deny check`: permissive licences only, no
known advisories, crates.io only.
