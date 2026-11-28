# Spec: Fortran comparison harness

Status: **specified** (decisions D-001 to D-013 in `decisions.md`); no
blockers. The GSLIB 2.0 source is vendored (D-013).

## Purpose

Provide verified reference outputs ("fixtures") from Fortran GSLIB so that
gslib-rs can be checked for parity (`AGENTS.md`, *Parity levels*), and so that
items in the provenance ledger (`docs/ROADMAP.md` §6) can be ticked off.

## Fixed by AGENTS.md (not up for decision here)

- The reference build uses double promotion (`-fdefault-real-8
  -fdefault-double-8` or equivalent).
- Output precision must not be limited by Fortran formats: either widen the
  formats or derive tolerances from them.
- Every fixture records the GSLIB version, compiler and flags, source patches,
  input data, and parameter file.
- Fixtures live under `tests/fixtures/gslib/`.
- Deterministic programs are compared numerically; stochastic programs
  statistically by default.

## Decisions so far

- Reference source: official download from gslib.com (D-001).
- Reference: GSLIB 2.0 (Fortran 77) under the 1996 Stanford notice (D-002,
  F-003). It is vendored unmodified in `vendor/gslib-2.0/src/`, with a
  top-level `NOTICE` (D-013).
- Oracle runs offline via a pinned recipe; tests compare against committed
  fixtures (D-003).
- Input data: our own seeded synthetic suites plus the GSLIB book data sets
  (the latter once their licence is known) (D-004).
- The oracle is the double-promoted build from source **only**; precompiled
  binaries (e.g. upc-ghs/gslib-executables, F-001) are never a test oracle
  (D-005).
- The recipe builds the vendored source with gfortran in a session (D-006,
  D-013).
- Each reference case is a versioned TOML case file. The recipe generates
  GSLIB's `.par` from it; Rust tests read the TOML (D-007).
- The recipe is a Rust `xtask`. The schema lives in an internal `oracle-case`
  crate, and the harness tests in a tests-only `harness` crate (D-008).
- Fixtures store GSLIB's outputs verbatim alongside `provenance.toml` (D-009).
- Deterministic outputs are compared with per-quantity, justified tolerances;
  structure must match exactly (D-010).
- Stochastic programs are compared ensemble against ensemble at 5 standard
  errors, plus independent invariant checks (D-011).
- Coverage order: vmodel → nscore/backtr → gamv → kt3d → sgsim (D-012).

## Requirements

- **R1 Recipe.** `cargo xtask oracle` (D-008), given gfortran:
  1. verifies `vendor/gslib-2.0/src/` against `MANIFEST.sha256`, and fails on
     any mismatch (with `--archive PATH`, also checks the original zip's
     SHA-256);
  2. widens numeric output formats with generic rewrite rules (no GSLIB source
     lines stored in the repository);
  3. builds with double promotion, recording compiler version and flags;
  4. runs each reference case and writes its fixture and provenance record.
- **R2 Offline tests.** `cargo test` never needs network access, gfortran or
  GSLIB; it reads committed fixtures only.
- **R3 Reproducible fixtures.** Re-running the recipe with the same pinned
  inputs reproduces the committed fixtures byte-for-byte, or the difference is
  explained in the provenance record.
  *Evidence:* after the step 3 recipe change, regenerating all 8 vmodel
  fixtures reproduced every `vmodel.var` and `params.par` byte-for-byte;
  only `provenance.toml` changed (recipe hash, commit, time).
- **R4 Synthetic suites.** A committed generator script, seeded, produces each
  synthetic input and a note saying what behaviour it targets.
- **R5 Case files.** Each reference case is one versioned TOML file. The
  recipe generates the program's `.par` from it and stores that `.par` with
  the fixture.
- **R6 Fixture layout.** `tests/fixtures/gslib/<program>/<case>/` contains
  `case.toml`, the inputs, `params.par`, GSLIB's outputs byte-for-byte, and
  `provenance.toml` with every field listed in D-009. With vendoring (D-013),
  the archive hash is replaced by `source_manifest_sha256`, the SHA-256 of
  `vendor/gslib-2.0/MANIFEST.sha256`. The harness fails if any field is
  missing, if `case.toml` no longer matches its recorded hash, or if the
  compiler flags lack `-fdefault-real-8`/`-fdefault-double-8`. Provenance also
  records `recipe_sha256`, a content hash of the recipe that survives squash
  merges (D-014). Implemented in `oracle-case`.
- **R7 Deterministic comparison.** Structure (shape, names, missing positions)
  must match exactly. Values: by default, relative 1e-12 for direct formulas
  and k·cond(A)·ε for kriging results (k to be calibrated, D-010). Overrides
  need a justification in `case.toml`. Failures report the maximum
  absolute/relative error and its location.
- **R8 Statistical comparison.** For each stochastic case, N realizations
  (N = 100 to start) from each implementation. Per lag and per quantile, the
  difference of ensemble means must be within 5 combined standard errors.
  gslib-rs must also reproduce the input model and histogram within 5
  standard errors on its own. Seeds are fixed; fixtures stay under 2 MB.

## Expected deviations

Known divergences between gslib-rs and GSLIB (`docs/ROADMAP.md` §2) would
make some fixtures fail by design. They are handled explicitly, never by
loosening tolerances:

- **Lag classes (gamv).** gslib-rs uses equal classes over `[0, max]`; GSLIB
  uses `xlag`/`xltol`. GSLIB-style lag classes must be implemented before gamv
  fixtures are added (plan step 6).
- **Ties (nscore).** Cases without ties must match exactly. Cases with ties
  are either declared an expected deviation in `case.toml`, with the reason,
  or compared against an opt-in legacy tie policy once one exists.
- **Search (kt3d, sgsim).** The first kriging cases use a search that
  includes every datum, so only the model and solver are compared. Cases that
  exercise GSLIB's super-block and octant search wait for that search to be
  implemented.

## Open questions

- None. (D-002 resolved by F-003 and D-013.)

## Implementation plan

Each step is one reviewable PR. The checks in `AGENTS.md` must pass, and no
step weakens a lint or tolerance.

1. **Harness scaffolding** (done). Add the `oracle-case`, `harness` and `xtask`
   crates (`publish = false`) and the `cargo xtask` alias.
   - Schema v1 types and TOML parsing in `oracle-case`.
   - Fixture discovery and provenance validation (R6) in `harness`.
   - Vendored-source verification against `MANIFEST.sha256` in
     `xtask oracle`.
   - New dependencies (e.g. `serde`, `toml`, `sha2`) justified and passing
     `cargo deny`.
   - Tests: schema round trip; rejection of a missing provenance field;
     checksum mismatch fails; `cargo test` passes with no fixtures present.
2. **Build stage** (done: `cargo xtask oracle build`, F-005). The xtask
   locates gfortran (setup runs `apt-get update` first, D-006), applies the
   generic format-widening rules, builds with double promotion, and records
   the full compiler version and flags.
   - Tests: rewrite rules on Fortran snippets written for the test (no GSLIB
     source lines).
   - The patterns the rules must match are taken from the vendored source.
3. **Synthetic data generator** (done: `cargo xtask gen-data [--check]`, 10
   suites in `tests/fixtures/data/`). Seeded and deterministic, writing a
   note per suite (R4). A case lists a suite file in `inputs`, and `generate`
   copies it into the case.
   - Tests: the same seed gives byte-identical output; each suite contains
     the edge cases it claims.
4. **vmodel** (done: 8 fixtures, F-006; tolerance default confirmed, D-015).
   - Work: `.par` generator (layout read from the source); reader for its
     output; comparison through `VariogramModel::gamma` along GSLIB's
     directions.
   - Cases: each structure type, nugget only, nested structures, azimuth,
     dip (including its sign), plunge, and anisotropy ratios.
   - Outcome: tick off the `setrot`/`cova3` items in ROADMAP §6, or fix the
     code.
5. **nscore / backtr** (done: 6 fixtures, F-007, D-016).
   - Work: `.par` generators and a transformation-table reader.
   - Cases: tie-free (exact), weighted, ties (expected deviation), and tails
     (option 1).
   - Outcome: tick off the nscore, backtr and gauinv ledger items.
   - Then a decision entry calibrating the 1e-12 default from the evidence
     (D-010).
6. **gamv** (done for the omnidirectional semivariogram: 7 fixtures, F-008,
   D-017).
   - Work: first implement GSLIB-style lag classes (`xlag`/`xltol`) in
     `variogram`, with their own unit tests; then the `.par` generator and an
     output reader.
   - Cases: omnidirectional only at this stage.
7. **kt3d** (done for simple and ordinary point kriging: nine fixtures,
   F-015; k calibrated to 100, D-035).
   - Cases: simple and ordinary kriging at points, with a search covering
     every datum; well- and ill-conditioned configurations.
   - Then a decision entry calibrating k in k·cond(A)·ε (D-010).
8. **sgsim.** Groundwork done: an `sgsim` case generator (so far only
   unconditional runs with nothing to search, RNG spec RN7) and two
   fixtures reproduced bit for bit (F-017). The stock build is 2-D only
   (`MAXZ = 1`), so statistical-parity cases stay 2-D or the recipe patches
   `SGSIM.INC` in its copy, recorded in provenance.
   - Cases: N = 100 realizations per implementation on a small grid; R8
     statistics plus the independent model-reproduction checks.
   - Fixture under 2 MB (D-011).

Nothing is blocked: every step can proceed in order. Step 4 also ticks off
provenance-ledger items against the vendored source.
