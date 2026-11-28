# Integration plan and handoff

The single entry point for picking up this project. It sequences all remaining
work and names the specs that still have to be written before code. Rules:
`AGENTS.md`. What is ported and what diverges: `docs/ROADMAP.md`. Decisions
and findings: `docs/spec/decisions.md`.

## 1. State at handoff (2026-11-27)

- **Branch `claude/code-review-project-status-orm47m`, PR #4** (open, not
  merged). It contains harness steps 1–6 and the missing-values spec.
- **Local checks:** every check in `AGENTS.md` passes. That is 170 tests,
  including the fixture parity tests, plus `cargo xtask oracle
  verify-source`, `cargo xtask gen-data --check` and `cargo xtask oracle
  status`.
- **Fixtures:** 34 under `tests/fixtures/gslib/` (vmodel 8, nscore 7,
  backtr 1, gamv 9, kt3d 9). All were made by the current recipe, and regenerating
  them gives byte-identical GSLIB outputs.
- **Specs:**

  | Spec | State |
  |---|---|
  | `docs/spec/harness.md` | Steps 1–6 done (6 for the omnidirectional semivariogram only); 7 (kt3d) and 8 (sgsim) to do |
  | `docs/spec/missing-values.md` | Decided (D-018…D-023); plan steps 1–4 to do |

- **D-017** (gamv counts each distinct pair once, with no self-pairs) was
  confirmed by the owner on 2026-11-27.
- **GitHub CI does not run** (§2). Nothing has been verified by CI. Every
  "passes" above is from local runs.

## 2. GitHub integration

What was checked, through the GitHub API, on 2026-11-27:

- **The repository is already public**, so making it public changes
  nothing. Standard GitHub-hosted runners (`ubuntu-latest`, which is all
  `ci.yml` uses) are free for public repositories, so this is not a minutes
  quota either.
- **Actions is enabled.** Dependabot's update jobs run on GitHub's
  infrastructure and succeed.
- **Every `CI` run fails within about 4 s.** Every run checked, on PR #4
  and on Dependabot's PR #2, fails the same way:
  - each job has `runner_id: 0`, no runner name and no steps;
  - job logs return 404, and check runs have no output.

  The workflow never starts on a machine, so the cause is not in
  `ci.yml`'s steps.
- **The most likely cause is the account itself, not the repository.** For
  example, a billing lock (a failed payment or an unpaid invoice locks
  Actions even for public repositories) or an Actions spending limit or
  policy. GitHub shows the exact reason only as an annotation in the web UI,
  which the API does not return.

**Owner actions (only the owner can do these):**

1. Open any failed run, e.g.
   <https://github.com/davidcmiller1135/gslib-rs/actions/runs/46063147735>,
   and read the annotation at the top of the summary page.
2. If it mentions billing or a spending limit: go to *Settings → Billing and
   plans* (personal account) and resolve the payment or raise the limit.
   Then re-run the workflow.
3. If it mentions disallowed actions: go to *Repository Settings → Actions →
   General* and allow `actions/checkout`, `Swatinem/rust-cache` and `EmbarkStudios/cargo-deny-action` (or all
   actions).
4. Once CI is green on PR #4:
   - add branch protection on `main` (require the CI jobs, require a PR,
     no force-push);
   - in *Settings → General → Pull Requests*, allow squash merging only.
5. Delete the stale branches `master` and `docs/agent-rules`. This session's
   credentials could not delete them.
6. Dependabot PR #2 (`actions/checkout` 4 → 7): merge it after CI runs.
   Check that the `oracle` job still works with the new version.

Until CI runs, every PR description must say which checks were run locally
and which were not (`AGENTS.md`, *Report honestly*).

## 3. Setting up a new session

- **Toolchain:** pinned in `rust-toolchain.toml` (1.97.0); rustup installs it
  on first use. The MSRV is 1.89 (CI `msrv` job).
- **gfortran** is needed only to regenerate fixtures, not to run tests. The
  tests use the committed fixtures. To install it:

  ```sh
  apt-get update && apt-get install -y gfortran
  ```

  The recipe records gfortran 13.3.0. A different version is allowed, but it
  is recorded in each fixture's provenance, and its outputs must be compared
  with the committed ones before they are committed (R3).
- **cargo-deny** is needed only for `cargo deny check`: `cargo install
  cargo-deny --locked`.
- **First commands:**

  ```sh
  cargo test --workspace
  cargo xtask oracle verify-source
  cargo xtask oracle status
  ```

## 4. Work queue

Steps are in dependency order. Each is one reviewable PR. **Spec first**
marks work whose decisions are not yet recorded: run the project `grill-me`
skill on it before writing code. Steps without that mark already have a spec
and acceptance tests named in it.

### Phase A: close out the current PR

| # | Work | Spec | Depends on |
|---|---|---|---|
| A2 | CI running and green on PR #4; merge (squash) | §2 above | owner actions |

### Phase B: missing values (every later program reads data through this)

| # | Work | Spec | Depends on |
|---|---|---|---|
| B1 | **Done.** Extraction API: `Trimming`, `NonFinite`, `GslibData::extract`, `samples` with row indices | `missing-values.md` step 1 (M1, M4, M5, M7, M9) | A2 |
| B2 | **Done.** Row-aligned writing: `with_optional_column` | step 2 (M3, M8) | B1 |
| B3 | **Done** (`--unest` deferred to the first row-aligned output). CLI flags and the sentinel guard | step 3 (M6) | B1, B2 |
| B4 | **Done** (F-011). Parity fixtures with trimming (`nscore`, `gamv`, a limit equal to a datum) | step 4 | B1 |

### Phase C: finish the harness plan

| # | Work | Spec | Depends on |
|---|---|---|---|
| C1 | **Spec first:** directional `gamv` (azimuth/dip, angular tolerances, bandwidths) and the other nine variogram types. Decide which types are in scope, how `ivtype` maps to a Rust type, and how the cross-variable inputs are shaped. Then fixtures. | new `docs/spec/gamv.md` | A2 (B4 for trimmed cases) |
| C2 | **Spec done** (`docs/spec/search.md`, D-024…D-028); **steps 1–3 done.** Search neighbourhood: data search with brute force and an exact super-block, `kt3d` debug fixtures, `sgsim` node search and two-part search (plan steps 1–6) | `docs/spec/search.md` | A2 |
| C3 | **Done** (F-015, D-035): `kt3d` harness step 7, simple and ordinary point kriging, with the search checked through `kt3d`'s debug output; the `picksup` miss confirmed | `harness.md` step 7 | — |
| C4 | **Spec first:** `kt3d` completeness (block discretisation, trend/KT, external drift, cross-validation, gridded output, `OnSingular` policy and its UNEST code, ledger item "unsolvable systems are skipped") | new `docs/spec/kt3d.md` | C3 |
| C5 | **Spec first:** ACORN and exact-stream parity as an opt-in. Decides the RNG trait, the `Seed` type, the random path and the per-realization seed derivation. | new `docs/spec/rng.md` | A2 |
| C6 | `sgsim` harness step 8: statistical parity, N = 100 realizations, R8 statistics, fixture < 2 MB. Also adds the covariance lookup table and spiral search. | `harness.md` step 8; needs C2, C5 | C2, C5 |
| C7 | **Spec done** (`docs/spec/parallelism.md`, D-029…D-033); **steps 1–3 done**: `rayon` in an opt-in `parallel` module; parallel kriging over targets; exact pair pruning and a fixed-block parallel variogram; parallel realizations after C5; GPU only with evidence | `docs/spec/parallelism.md` | C2 step 2 (done); realizations also need C5 |

### Phase D: remaining parity (ROADMAP §1, suggested order)

Each item is **spec first**, except where noted.

| # | Work | Notes |
|---|---|---|
| D1 | `declus` | cell-size scan, origin offsets; deterministic, fixture-friendly |
| D2 | `nscore`/`backtr` completion | reference distribution, power/hyperbolic tails, writing the transformation table, despiking as a legacy policy, CLI subcommands |
| D3 | `ik3d` and multi-threshold `sisim` | order-relation corrections, categorical data, soft data, Markov–Bayes |
| D4 | `cokb3d` | ordinary/collocated cokriging; check of the linear model of coregionalisation |
| D5 | `.par` adapter / GSLIB-compatibility mode | the one place that selects every legacy policy (D-020, D-021, D-023, `OnSingular`) |
| D6 | the rest: `gam`, `varmap`, `trans`, `addcoord`, `draw`, `histsmth`/`scatsmth`, `bicalib`, `postik`/`postsim`, `anneal`, `lusim`, `gtsim`, `pfsim`, `ellipsim`, plotting programs | one at a time, each with its own fixtures |

### Phase E: beyond GSLIB

`docs/ROADMAP.md` §5. Not started. The architecture direction in §4 applies
throughout: functional core, imperative shell, domain types.

## 5. Open provenance items (ROADMAP §6)

These must be checked against the vendored source or a fixture before code
relies on them:

- ACORN and seed 69069 in each simulation program (confirmed for `nscore`
  only): C5.
- Grid conventions (`xmn` is the first cell's centre, x-fastest order):
  `simulation/grid.rs`. Check from `sgsim`'s source or a fixture in C6.
- Search parameter names (`ndmax`, `nodmax`, `radius`): C2.
- `sasim` schedule parameters and objective normalisation.
- Unsolvable kriging systems are skipped and written as UNEST: C4.

## 6. Rules for the next session

- Read `AGENTS.md` first. It wins over this plan where they conflict.
- Specs are written with the project `grill-me` skill: one question at a
  time, a recommendation each time, and the decision recorded at once.
- GSLIB facts come from `vendor/gslib-2.0/src/`, cited by file and routine,
  or from a fixture. Never from memory.
- Never edit generated fixtures. Change `case.toml` or the recipe, commit,
  then run `cargo xtask oracle generate --program <name>`, then check
  `cargo xtask oracle status`.
- Every parity deviation is either a documented divergence (ROADMAP §2) or
  explained by a check that reproduces GSLIB's own arithmetic (D-016).
- Run mutation checks on new comparisons: plant a bug, confirm the harness
  catches it, and record the count in the finding.
- Do not merge without the owner's go-ahead.
