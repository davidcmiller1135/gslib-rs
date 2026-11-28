# Integration plan and handoff

The single entry point for picking up this project. It sequences all remaining
work and names the specs that still have to be written before code. Rules:
`AGENTS.md`. What is ported and what diverges: `docs/ROADMAP.md`. Decisions
and findings: `docs/spec/decisions.md`.

## 1. State at handoff (2026-11-28)

- **`main` holds all work to date.** PR #4 (the harness, steps 1–8; missing
  values; search; parallelism; the RNG building blocks; `sgsim` statistical
  parity) and PR #5 (directional `gamv`, F-027) were merged on 2026-11-28.
- **Checks are local** (D-053): there is no CI. Every check in `AGENTS.md`
  passes, with 46 reference fixtures, decisions D-001…D-053 and findings
  F-001…F-027.

### Milestones

| # | Milestone | State |
|---|---|---|
| M0 | Foundations | done |
| M1 | Merged on GitHub | done (PRs #4 and #5); no CI, by decision (D-053) |
| M2 | Missing values handled correctly | done (F-011) |
| M3 | Variograms match GSLIB | `gamv` done on PR #5: directions, every variogram type, several variables, indicators, 17 fixtures (F-027). `gam` and `varmap` remain (no spec yet) |
| M4 | Kriging matches GSLIB | search (steps 1–5, CLI flags) and point SK/OK done (F-015); block kriging, trend and external drift, cross-validation, `OnSingular` need a spec |
| M5 | Simulation matches statistically | RNG spec; ACORN, path and `gauinv` bit for bit (F-017); `sgsim` statistical parity (F-018); parallel realizations. Open: `sisim` parity, transformed `sgsim` cases, the covariance lookup table, 3-D `sgsim` (the stock build is 2-D only; needs a recorded patch), the optional exact-stream mode |
| M6 | Drop-in for GSLIB users | not started: `.par` files and a compatibility mode, `declus`, finishing `nscore`/`backtr` |
| M7 | Every GSLIB program | not started: `ik3d`, multi-threshold `sisim`, `cokb3d`, the rest |
| M8 | Beyond GSLIB | parallelism done; optional GPU benchmark and roadmap §5 remain |

**Specs still to write (grill sessions):** the rest of
`kt3d`; `.par` files and the compatibility mode; then `declus` and `ik3d`.

**Suggested next steps:** (A) `sisim` parity fixtures, to close M5's most
visible gap; (B) the `kt3d` spec (M4); (C) a fixture with a wide dip
tolerance, to confirm F-024.

## 2. GitHub and checks

- **There is no CI, by the owner's decision** (D-053, 2026-11-28).
  `.github/workflows/ci.yml` was removed. Before a change is merged, every
  check in `AGENTS.md` is run locally, and the PR says which were run and
  which were not (`AGENTS.md`, *Report honestly*).
- **Dependabot** still opens weekly crate-update PRs, which run on GitHub's
  side at no cost. Run the checks locally before merging one. Its GitHub
  Actions updates were removed along with the workflow, so its PR #2
  (`actions/checkout`) is obsolete and can be closed.
- **Owner housekeeping:** delete the merged or stale branches `master`,
  `docs/agent-rules` and `ccr-af700c64-3iveok`. Agent sessions cannot
  delete branches.

## 3. Setting up a new session

- **Toolchain:** pinned in `rust-toolchain.toml` (1.97.0); rustup installs it
  on first use. The MSRV is 1.89: check it with `cargo +1.89 test
  --workspace` (`rustup toolchain install 1.89 --profile minimal` first).
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
| A2 | **Done.** PR #4 merged (2026-11-28). CI dropped (D-053) | §2 above | |

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
| C1 | **Done** on PR #5 (F-027): directional `gamv` and the other nine variogram types, with fixtures. | `docs/spec/gamv.md` | A2 |
| C2 | **Done** (`docs/spec/search.md`, D-024…D-028, all six steps, with CLI flags). Search neighbourhood: data search with brute force and an exact super-block, `kt3d` debug fixtures, `sgsim` node search and two-part search (plan steps 1–6) | `docs/spec/search.md` | A2 |
| C3 | **Done** (F-015, D-035): `kt3d` harness step 7, simple and ordinary point kriging, with the search checked through `kt3d`'s debug output; the `picksup` miss confirmed | `harness.md` step 7 | — |
| C4 | **Spec first:** `kt3d` completeness (block discretisation, trend/KT, external drift, cross-validation, gridded output, `OnSingular` policy and its UNEST code, ledger item "unsolvable systems are skipped") | new `docs/spec/kt3d.md` | C3 |
| C5 | **Spec done** (`docs/spec/rng.md`, D-036…D-039); **steps 1–4 done** (F-017: ACORN, path and `gauinv` bit for bit against GSLIB): ChaCha8 by default with realization k on stream k; ACORN, GSLIB's path and `gauinv` as checked building blocks; an `sgsim` exact mode after harness step 8 | `docs/spec/rng.md` | A2 |
| C6 | **Done** (F-018, D-040): `sgsim` statistical parity (harness step 8), unconditional and conditional, within 5 SE of GSLIB and of the model | `harness.md` step 8 | — |
| C7 | **Spec done** (`docs/spec/parallelism.md`, D-029…D-033); **steps 1–4 done**: `rayon` in an opt-in `parallel` module; parallel kriging over targets; exact pair pruning and a fixed-block parallel variogram; parallel realizations after C5; GPU only with evidence | `docs/spec/parallelism.md` | C2 step 2 (done); realizations also need C5 |

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
