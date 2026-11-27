# Roadmap

Goal, in order:

1. **Parity:** everything GSLIB (2nd edition, Deutsch & Journel 1998) can do,
   gslib-rs can do, validated against the Fortran.
2. **Foundations:** clean architecture (hexagonal, DDD vocabulary, functional
   core / imperative shell, test-first) so parity work doesn't calcify.
3. **Extension:** beyond what GSLIB offers.

Status key: ✅ done · 🟡 partial · ⬜ missing.

## 1. Program parity

### Data, transforms, utilities

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| file format | 🟡 | `io::parse_gslib` / `format_gslib` | **trimming limits** (`tmin`/`tmax`) and missing values, used by every program; Fortran fixed-format output |
| `.par` files | ⬜ | (CLI flags) | read GSLIB parameter files so existing workflows run unchanged |
| `nscore` | 🟡 | `transform::normal_scores`; parity fixtures (F-007) | trimming limits; reference distribution option; rank output (negative column); random despiking as a legacy policy; writing the transformation table; CLI subcommand |
| `backtr` | 🟡 | `transform::BackTransform`; parity fixtures (F-007) | power (ω ≠ 1) and hyperbolic tails; CLI subcommand |
| `declus` | ⬜ | | cell declustering with cell-size scan and origin offsets |
| `trans` | ⬜ | | histogram identification with control on how far data may move |
| `rotcoord` | ✅ | `geometry::Rotation2` | CLI subcommand |
| `addcoord` | ⬜ | | add grid coordinates to a gridded file (`simulation::Grid::nodes` does the maths) |
| `draw` | ⬜ | | draw a random sample from a distribution |
| `histsmth`, `scatsmth` | ⬜ | | smoothed (annealed) univariate/bivariate distributions |
| `bicalib` | ⬜ | | calibration of secondary data for Markov–Bayes |
| `postik`, `postsim` | ⬜ | | post-processing of IK/simulation results (E-type, quantiles, probabilities) |

### Variograms

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `cova3` / `setrot` | ✅ | `model::VariogramModel`, `geometry::Anisotropy` | all five structure types and three-angle anisotropy; tested against hand-derived values, not yet against Fortran |
| `gamv` | 🟡 | `variogram::{experimental_variogram, CentredLags}`; omnidirectional parity fixtures (F-008) | **directions** (azimuth/dip with tolerances and bandwidths); trimming limits; sill standardization; the other 9 variogram types (cross, covariance, correlogram, general/pairwise relative, log, madogram, indicator); multiple variables; CLI on GSLIB's lags |
| `gam` | ⬜ | | variograms of gridded data along grid directions |
| `varmap` | ⬜ | | variogram maps |
| `vmodel` | 🟡 | `VariogramModel::gamma` | write model values along directions |
| `vargplt` | 🟡 | `plot::VariogramPlot` | multiple directions on one plot |

### Estimation

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `kb2d` / `kt3d` | 🟡 | `kriging::{simple_kriging, ordinary_kriging}`, `krige` CLI | **block kriging** (discretisation); **kriging with a trend** (KT, `ktsol`) and **external drift**; cross-validation / jackknife; gridded output; locally varying mean |
| search (`setsupr`/`srchsupr`) | 🟡 | `search::Search` | **anisotropic search ellipsoid**; **octant search**; **super-block search** (currently brute force, O(N) per lookup) |
| `cokb3d` | 🟡 | `kriging::simple_cokriging` | ordinary cokriging (both constraint forms); collocated cokriging; check the linear model of coregionalisation |
| `ik3d` | ⬜ | (single threshold via `sisim`'s kriging) | multiple thresholds; median IK; **order-relation corrections**; categorical variables; soft data |

### Simulation

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `sgsim` | 🟡 | `simulation::sgsim`, `sgsim` CLI | **covariance lookup table** and spiral search of simulated nodes; separate data / previously simulated search counts; multiple-grid refinement; assigning data to grid nodes; locally varying mean; collocated cokriging with a secondary variable |
| `sisim` | 🟡 | `simulation::sisim` | multiple thresholds; categorical; soft data; Markov–Bayes; order-relation corrections |
| `sasim` | 🟡 | `simulation::sasim` | conditional-distribution objective; **incremental objective updates** (currently a full recompute per perturbation); variogram objectives in several directions; CLI subcommand |
| `anneal` | ⬜ | | post-processing of an existing realization by annealing |
| `lusim` | ⬜ | | LU (Cholesky) simulation for small grids |
| `gtsim` | ⬜ | | truncated Gaussian simulation |
| `pfsim` | ⬜ | | probability-field simulation |
| `ellipsim` | ⬜ | | Boolean simulation of ellipsoids |

### Plotting

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `histplt` | 🟡 | `plot::HistogramPlot` | summary statistics on the plot; weights; log scale; CLI |
| `scatplt` | 🟡 | `plot::ScatterPlot` | statistics; log scales; CLI |
| `pixelplt` | 🟡 | `plot::PixelMap` | slices of 3-D grids; categorical palettes; colour bar; CLI |
| `probplt`, `qpplt`, `locmap`, `bivplt` | ⬜ | | |

### Suggested order

1. Trimming limits and missing values. Every GSLIB program assumes them.
2. Golden-file validation harness (§3). Do this before more algorithms.
3. Directional `gamv`, then variogram types.
4. Anisotropic and super-block search; covariance lookup table in `sgsim`.
5. `kt3d` completeness (block kriging, trend, cross-validation).
6. `declus`; `nscore`/`backtr` CLI, tails.
7. `ik3d` and multi-threshold `sisim` with order-relation corrections.
8. `.par` file adapter, so GSLIB users can switch over.
9. The remaining simulation and plotting programs.

## 2. Deliberate divergences from GSLIB

Rules (see `AGENTS.md`): explicit failure is the default, and GSLIB's silent
coercions become **opt-in legacy policies, off by default**. Each divergence
is documented at its definition and covered by tests. Descriptions of GSLIB
behaviour in this table are recalled from memory until checked off in §6.

| Topic | GSLIB | gslib-rs default | Legacy policy |
|---|---|---|---|
| Errors | stops, or writes a message and carries on | returns `Result`; never panics | (n/a) |
| Unsolvable kriging system | node skipped / written as `UNEST` | `Error::SingularSystem` | planned: `OnSingular::SkipAsUnestimated` |
| Anisotropy ratio ≤ 0 | clamped to `EPSLON` | `Error::InvalidInput` | planned: clamp policy |
| Power model covariance | `PMX − γ` with a large `PMX` | `CovarianceModel` cannot be built from it | planned, if exact parity needs it |
| Ties in `nscore` | every datum perturbed by ACORN × 1e-6 (seed 69069), so ties break randomly (F-005) | tied values share one score | planned: random despiking with ACORN |
| Normal quantile | `gauinv`: rational approximation, error ≤ 1.501e-8 (F-007) | `erfc_inv`, ~1e-16 | not planned: extension |
| Normal CDF (`backtr` tails) | `gcum`: polynomial approximation, error ≤ 7.458e-8, exactly 0 or 1 past \|x\| > 6 (F-007) | `libm::erfc`, ~1 ulp | not planned: extension |
| Back-transform tails | options 1 (linear), 2 (power), 4 (hyperbolic, upper only); no clamping option | 1 (linear) only, plus `Tail::Clamp` (an extension) | power and hyperbolic planned |
| RNG | ACORN, global state | caller-supplied `rand::Rng` (CLI: ChaCha8) | ACORN as an `Rng`, for exact-stream parity |
| Lag classes | centres at multiples of `xlag`, inclusive tolerance `xltol`, may overlap; `xltol ≤ 0` silently becomes `xlag/2` | `CentredLags` reproduces this (F-008); a non-positive tolerance is an error; equal-width `LagClasses` also available | n/a |
| Zero-distance lag (`gamv`) | also pairs each datum with itself, diluting γ there towards 0 (F-008) | distinct coincident pairs only | none (GSLIB's value is recoverable via F-008's relations); D-017 |
| Trimming interval | mostly `tmin ≤ v < tmax`; closed tests in `gamv` pair rejection, `sasim`, `ik3d`, `sisim_gs`, `bivplt` (F-009) | `tmin ≤ v < tmax` everywhere (D-019) | per program, only if a fixture shows it matters |
| NaN / ±Inf in data | NaN kept as data; ±Inf trimmed by ±1e21 limits (F-010) | explicit policies, default error, applied before trimming (D-020, M9) | `Missing` / `ClampToRange` are opt-in policies |
| `lusim` unestimated code | −1.0 (collides with real data) | −999 by default, configurable (D-023) | `--unest -1` |
| Omnidirectional pair counts (`gamv`) | each pair counted in both orderings, so `np` is doubled (γ unaffected) | each pair counted once | n/a (reporting convention); D-017 |
| π in angle conversion | `3.141592654` (`setrot`), `3.14159265` (`cova3`) | exact `π` | planned for exact-stream parity |
| Coincident points (`cova3`) | squared anisotropic distance < `1e-10` | Euclidean distance < `1e-10` | planned: GSLIB coincidence policy |
| Power model range | range slot holds ω | `Kind::Power { exponent }`; range is a unit distance | range 1 reproduces GSLIB |
| Power model nugget | **defect:** the nugget cancels out of γ (`cova3` puts `c0` in `cmax`, F-006) | the nugget applies at every h > 0 | not planned (GSLIB defect) |

## 3. Validation strategy

Parity levels (from `AGENTS.md`):

| Programs | Default parity | Opt-in |
|---|---|---|
| Deterministic (`gamv`, `kt3d`, `nscore`, `declus`, …) | numerical, within documented tolerances | (n/a) |
| Stochastic (`sgsim`, `sisim`, `sasim`, …) | statistical: histogram, variogram and proportions over many realizations | exact stream parity (ACORN, same path and search, against the promoted build) |

Planned verification work:

- **Reference build.** Fortran GSLIB compiled with double promotion
  (`-fdefault-real-8 -fdefault-double-8`), with output formats widened (or
  tolerances derived from the formats). Compiler, flags and source patches are
  recorded with every fixture. Check GSLIB's licence before vendoring its
  source or data; until then, commit generated fixtures plus the recipe.
- **Golden fixtures** in `tests/fixtures/gslib/`, from the book's data sets
  (e.g. `cluster.dat`, `true.dat`), each with a provenance note.
- **Statistical tests** for simulations: fixed seeds; tolerances justified by
  sample size.
- **Property tests** (`proptest`, already in use): symmetry, round trips,
  monotonicity, positive definiteness for valid models.
- **Performance budgets:** `criterion` benchmarks at realistic sizes, so
  regressions are caught.

## 4. Architecture direction

The guiding principle is **no ceremony**: structure has to earn its place
(`AGENTS.md`).

**Now:**
- one library crate with a pure core (`geometry`, `model`, `variogram`,
  `kriging`, `search`, `transform`, `simulation`);
- boundary modules in the same crate (`io`, and `plot` behind a feature);
- a CLI shell.

**Next, when justified:**

- **Neighbourhood trait.** Brute force, super-block and k-d tree searches are
  genuine substitutes, so they share a trait. This is the one trait currently
  planned.
- **Crate split only for a concrete reason**, such as a boundary that drags in
  heavy dependencies. The `plot` feature already keeps system libraries out of
  the core; a separate I/O crate becomes worth it if formats like VTK or
  Parquet arrive.
- **Missing values.** `Option` in memory; sentinels and trimming limits only in
  `io` (and a parameter-file adapter).
- **Domain types at the edges of algorithms, flat slices inside.** Kernels
  work on `&[f64]` / `&mut [f64]`, with grid indexing through one `Grid` type.
- **Parse, don't validate.** `Sample`, `VariogramModel`, `CovarianceModel`,
  `LagClasses`, `Search` and `Grid` can only be built valid. SASIM's
  `Schedule` is still checked at use; fix that.
- **Legacy policies as values.** For example, `OnSingular`, passed explicitly
  and selected by a GSLIB-compatibility mode in the CLI.

## 5. Beyond GSLIB

Once parity is solid:

- **Parallel realizations** (rayon), with per-realization seeds derived from
  one master seed, so results don't depend on thread count.
- **Automatic variogram fitting** (weighted least squares with constraints).
  GSLIB leaves this to the user.
- **Fast neighbourhoods** (k-d tree, approximate search) and streaming over
  large grids.
- **Non-stationarity:** locally varying anisotropy and trends.
- **Multiple-point statistics** (SNESIM-style) and Gaussian-process
  formulations of kriging.
- **Bindings:** Python (PyO3) and WebAssembly (the `plot` feature already
  keeps the core free of system dependencies).
- **Uncertainty workflows:** Bayesian model parameters and realization
  ranking.

## 6. Provenance ledger

Claims about GSLIB written from memory, not yet checked against its source or
a reference run. Each is marked `TODO(parity-verify)` in the code. Tick an
item off only after checking it against the Fortran source or a recorded
reference run, and say which.

- [x] `setrot` rotation matrix and angle conventions: matches term for term;
  dip is negative downward per `gamv`'s comment and `setrot`'s arithmetic
  (`setrot`'s own comment says otherwise). F-004.
- [x] `setrot` clamps anisotropy ratios to `EPSLON = 1.e-20`. F-004.
- [x] `cova3` formulas, types 1–5, power ω in the range slot, `PMX = 999.`:
  all match **except the Gaussian, which was wrong in gslib-rs and is now
  fixed**. F-004.
- [x] `backtr` interior interpolation (linear in score space) and the meaning
  of tail option 1 (linear in probability, via `gcum`): `BACKTR.FOR`,
  confirmed by fixtures to 1.4e-16. F-007.
- [x] `nscore` tie handling: random despiking of every datum (`EPSLON =
  1.0e-6`, ACORN seeded 69069), `NSCORE.FOR` lines 31, 142, 179. F-005.
- [x] `gauinv` algorithm and accuracy: rational approximation (Kennedy &
  Gentle 1980, per its own comment), measured error ≤ 1.501e-8. F-007.
- [~] ACORN as the generator and 69069 as its seed: confirmed in `nscore`
  (F-005); still to check in each simulation program. §2, `cli/src/main.rs`
- [ ] Grid conventions: `xmn` is the first cell's centre, x-fastest order:
  `simulation/grid.rs`
- [ ] Search parameter names (`ndmax`, `nodmax`, `radius`): `search.rs`
- [ ] `sasim` schedule parameters (`t0, redfac, ka, k, num, Omin`) and
  objective normalisation by initial values: `simulation/sasim.rs`
- [x] `gamv` lag definition: centres `(il−2)·xlag`, inclusive `±xltol`
  (overlapping), zero class `h ≤ 1e-20`, cut-off `(nlag + 0.5)·xlag`.
  `GAMV.FOR`, confirmed by fixtures including exact-boundary pairs. F-008.
- [ ] Unsolvable kriging systems are skipped and written as `UNEST`: §2
- [x] Standard-normal reference values in the `transform.rs` tests: agree
  within 4.4e-16 with Python's `statistics.NormalDist` (CPython 3.13). F-007.
- [x] GSLIB version and licence terms of the reference source: GSLIB 2.0
  (`VERSION=2.000`), 1996 Stanford notice. Verified from the archive itself
  (F-003, D-013 in `docs/spec/decisions.md`).
