# Engineering & Architectural Rules (GSLIB Rust Rewrite)

Rules for anyone, human or agent, changing this repository. The goal is to
**recreate GSLIB's behaviour, then extend it**.

## Precedence & Scope

When rules conflict, the earlier item in this list wins:

1. **Numerical correctness and determinism.** Results are reproducible from
   the seed and inputs, and do not depend on the number of threads.
2. **Explicit failure over silent coercion** (see *Legacy behaviour* below).
3. **Documented GSLIB parity.**
4. **Scoped, reviewable edits.**
5. Performance, then architectural or aesthetic preferences.

Further rules:

- **Correctness before performance, numerical stability before expedience.**
  Use `f64` unless a routine or input specification requires otherwise.
- **Scoped changes only.** Read the surrounding code and tests first. Make
  small edits; no unrequested refactors or changes in unrelated modules.

## Legacy behaviour

- **Preserve GSLIB's numerical behaviour and edge cases** unless an
  intentional fix or extension is explicitly scoped and documented.
- **Silent coercions are opt-in, off by default.** Where GSLIB coerces
  silently, the library fails explicitly by default and offers the legacy
  behaviour as a named, explicit policy (e.g. `OnSingular::SkipAsUnestimated`)
  that a GSLIB-compatibility mode can select. Examples of such coercions:
  clamping anisotropy ratios, skipping an unsolvable kriging system, or faking
  a covariance for the power model.
- **Do not treat legacy agreement as proof of correctness** where a known
  GSLIB defect is in question. Fix it, document it, and test it.

## Parity levels

- **Deterministic programs** (`gamv`, `kt3d`, `nscore`, `declus`, …):
  numerical parity with the reference build (below) within documented
  tolerances.
- **Stochastic programs** (`sgsim`, `sisim`, `sasim`, …): **statistical parity
  by default**, meaning that over many realizations the reproduced histograms,
  variograms and proportions match within justified tolerances.
- **Exact stream parity** (the same realization from the same seed) is an
  **opt-in** per program. It needs the GSLIB generator (ACORN), the same
  random path, search order and arithmetic. Because this port computes in
  `f64`, exact parity is defined against the **double-promoted reference
  build**, never against stock single-precision GSLIB.

## Reference build (oracle)

- The Fortran reference is compiled with **double promotion**
  (`-fdefault-real-8 -fdefault-double-8` or equivalent): we recreate the
  behaviour, not single-precision rounding.
- Promotion is not free. Single-precision-tuned constants and thresholds may
  behave differently. Fortran output formats (e.g. `f12.4`) truncate precision,
  so widen them in the reference build or derive tolerances from the format.
  Record the compiler, flags, and any source patches with each fixture.

## Architecture & Boundaries

- **Functional core, imperative shell.** Geostatistical algorithms
  (variograms, kriging, sequential/LU/indicator/object simulation, transforms)
  live in pure modules with no I/O, global mutable state, or CLI/framework
  dependencies.
- **Isolate boundaries.** File parsing (GSLIB/Geo-EAS ASCII, parameter files,
  VTK, …), CLI/config parsing, parallel execution, and progress reporting stay
  at the edges.
- **Missing values are `Option` in memory.** Sentinels (`-999`, `UNEST`) and
  trimming limits (`tmin`/`tmax`) exist only at the I/O boundary.
- **Domain types vs. hot loops.** Configuration, coordinates and parameters
  use explicit domain types (`Grid`, `VariogramModel`, `Search`, seeds, …).
  Inner numerical loops use flat contiguous slices (`&[f64]`, `&mut [f64]`).
- **Conventions are domain law, defined once.** Grid origin is the first cell's
  centre; nodes are x-fastest; azimuth is clockwise from north; dip is negative
  downward. Use the existing types (`Grid`, `Angles`, `Anisotropy`); never
  re-derive conventions locally.
- **No ceremony.** No Clean/Hexagonal boilerplate, speculative abstractions,
  or traits unless they express a real domain concept or enable genuine
  substitution. A search-neighbourhood trait (brute force / super-block /
  k-d tree) qualifies; a single-implementation trait does not.

## Numerical & Geostatistical Integrity

- **Explicit inputs and streams.** Pass RNGs, seeds, simulation paths, search
  parameters, covariance models and anisotropy explicitly. No implicit global
  RNG or hidden mutable state. In parallel code, derive one seed per
  realization; never share a stream between threads.
- **Explicit failures.** Return typed errors (`gslib_rs::Error`) for singular or
  non-positive-definite systems, invalid models, search exhaustion,
  out-of-bounds indices and parse failures. Never silently coerce `NaN`/`Inf`,
  clamp invalid parameters, or swallow singular-matrix errors.
- **Document assumptions and divergences.** Record formulations, tolerances,
  edge cases, and intentional departures from GSLIB in module docs, and list
  every divergence in `docs/ROADMAP.md` §2.
- **No fabricated provenance.** Never invent citations, manual references,
  test oracles or expected outputs. A claim about GSLIB not yet checked against
  its source or a reference run is marked `TODO(parity-verify)` and listed in
  `docs/ROADMAP.md` §6.

## Testing & Verification

- **Every bug fix adds a regression test** that fails without the fix.
- **Test invariants and contracts** independently of the implementation:
  - exactness at conditioning points;
  - unbiasedness and weights summing to one where required;
  - variogram/covariance properties;
  - bounds and non-negativity;
  - symmetry and round trips (property tests with `proptest`).
- **Differential and fixture tests** compare against the reference build.
  Document tolerances and their source (precision, output format), and any
  expected deviation.
- **Statistical tests** use fixed seeds and tolerances justified by sample size
  (e.g. 5σ), so they are neither flaky nor vacuous.
- **Deterministic and isolated.** Fixed seeds, no shared global state, a
  per-test temp directory, and explicit degenerate and boundary cases.

## Rust & Repo Discipline

- **Follow workspace conventions:**
  - match existing idioms, the `Error` type and crate boundaries;
  - `unsafe` is forbidden;
  - library code never panics on bad input.
- **Lints are strict on purpose.** Exceptions are
  `#[expect(lint, reason = "...")]` on the smallest item. In a hot kernel,
  indexing or casts may be allowed this way only with benchmark evidence;
  otherwise use iterators, `chunks_exact`, `get()`, or `num.rs` conversions.
- **Never weaken a lint, test, tolerance or CI step to get green.**
- **Justify new crates.** Add a dependency only when needed; it must pass
  `cargo deny check`.
- **Explicit unresolved questions.** Use `TODO(<context>): <question>` only for
  a genuine ambiguity or parity choice that blocks completion;
  `TODO(parity-verify)` marks unverified claims about GSLIB.
- **Report honestly.** Run the checks below before claiming a change is done,
  and say what was and was not verified.

## Repo Specifics

- **Core (pure):** `gslib_rs/src/{geometry, model, variogram, kriging, search,
  transform, simulation, linalg, num, sample}.rs`.
- **Boundaries:**
  - `gslib_rs/src/io.rs`: pure `parse_gslib`/`format_gslib`, plus file wrappers;
  - `gslib_rs/src/plot.rs`: behind the `plot` feature;
  - `cli/`: arguments, seeding, files.
- **Harness (not part of the library):**
  - `oracle-case/`: fixture schema (`case.toml`, `provenance.toml`),
    validation, discovery;
  - `harness/`: tests-only crate checking gslib-rs against the fixtures;
  - `xtask/`: the reference recipe, e.g. `cargo xtask oracle verify-source`.
- **Checks** (all run in CI; the toolchain is pinned in `rust-toolchain.toml`):
  ```sh
  cargo fmt --all --check
  cargo clippy --workspace --all-targets -- -D warnings
  cargo clippy -p gslib_rs --all-targets --no-default-features -- -D warnings
  cargo test --workspace
  cargo doc --workspace --no-deps
  cargo deny check
  cargo xtask oracle verify-source
  ```
- **GSLIB 2.0 source (the oracle):** `vendor/gslib-2.0/src/`, vendored
  byte-for-byte. **Never modify it**: `MANIFEST.sha256` verifies it, and
  build-time changes are applied to a copy. Answer questions about GSLIB's
  behaviour from this source, citing file and routine. Its licence notice is
  in `NOTICE`.
- **Reference fixtures:** not yet built. Planned location is
  `tests/fixtures/gslib/`, each fixture with a provenance note: compiler and
  flags, source patches, input data, and parameter file
  (`docs/spec/harness.md`).
- **Roadmap, divergences, provenance ledger:** `docs/ROADMAP.md`.
