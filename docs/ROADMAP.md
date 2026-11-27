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
| `nscore` | 🟡 | `transform::normal_scores` | declustering weights from file; reference distribution option; writing the transformation table; CLI subcommand |
| `backtr` | 🟡 | `transform::BackTransform` | power (ω ≠ 1) and hyperbolic tails; CLI subcommand |
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
| `gamv` | 🟡 | `variogram::experimental_variogram` | **directions** (azimuth/dip with tolerances and bandwidths); lag tolerance; the other 9 variogram types (cross, covariance, correlogram, general/pairwise relative, log, madogram, indicator); multiple variables |
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

Each of these is documented at the definition and covered by tests.

| Topic | GSLIB | gslib-rs | Why |
|---|---|---|---|
| Errors | stops or writes a message and carries on | returns `Result`; never panics | callers decide; nothing is silently wrong |
| Unsolvable kriging system | estimate skipped (or flagged) | `Error::SingularSystem` | a fallback to the mean hid duplicate data |
| Power model covariance | `PMX − γ` with `PMX = 999` | `CovarianceModel` cannot be built from it | the trick gives arbitrary results in SK/SGS |
| Ties in `nscore` | broken by sort order (despiking optional) | tied values share one score | deterministic, order-independent |
| Normal quantile | `gauinv` (Odeh–Evans, ~1e-7) | `erfc_inv` (~1e-16) + `libm::erfc` | double-precision accuracy |
| RNG | ACORN, global state | caller-supplied `rand::Rng` (CLI: ChaCha8) | reproducible, thread-safe; GSLIB realizations cannot be bit-reproduced anyway (path and search differ) |
| Lag classes | lag centres at multiples of `xlag` with `xltol` | equal classes over `[0, max]` | to be revisited with directional `gamv` |
| Power model range | range slot holds ω | `Kind::Power { exponent }`, range is a unit distance | separates the two meanings; range 1 matches GSLIB |

## 3. Validation strategy

The current tests check analytic cases and invariants. Parity needs comparison
with the Fortran.

- **Golden files.** Build Fortran GSLIB with gfortran in CI (or commit its
  outputs). Run both on the book's data sets (`cluster.dat`, `true.dat`,
  `ydata.dat`). Deterministic programs (`gamv`, `kt3d`, `nscore`, `declus`)
  should agree to ~1e-6. Check GSLIB's licence before vendoring its source or
  data.
- **Statistical tests for simulation.** Over many realizations, the mean
  histogram and variogram should reproduce the model within tolerance. This
  checks correctness rather than matching random numbers.
- **Property tests** (already started with `proptest`): symmetry, round trips,
  monotonicity, positive-definiteness of covariance matrices for valid models.
- **Performance budgets.** `criterion` benchmarks for the variogram, kriging
  and SGS on realistic sizes, run in CI so regressions are caught.

## 4. Architecture direction

Where the code is now, and where it is heading.

**Now:** one library crate with a mostly pure core (`model`, `kriging`,
`variogram`, `transform`, `simulation`, `search`), adapters in the same crate
(`io`, `plot` behind a feature), and a CLI shell.

**Next steps:**

- **Split the crates along the hexagon:**
  - `gslib-core`: domain types and pure algorithms; no filesystem, no plotting.
  - `gslib-io`: GSLIB ASCII, `.par`, CSV, later Parquet/VTK/GeoTIFF.
  - `gslib-plot`.
  - `gslib-cli`.
- **Ports as traits** in the core:
  - a `Neighbourhood`/spatial-index trait (brute force now, k-d tree or super
    block later);
  - a `CovarianceFunction` trait so non-GSLIB models can plug in;
  - sources and sinks for data and realizations.
- **Ubiquitous language (DDD).** Keep the core's names in GSLIB and
  geostatistics terms: `Sample`, `Structure`, `VariogramModel`, `Search`,
  `Grid`, `Realization`. Candidate bounded contexts: exploratory analysis,
  variography, estimation, simulation, transforms.
- **Functional core.** Algorithms take every input explicitly, including the
  random number generator, and return values. Local mutation inside a
  function is fine; shared state is not. Next: replace the remaining
  `Vec`-mutating loops in SASIM with an explicit state value.
- **Parse, don't validate.** `Sample`, `VariogramModel`, `CovarianceModel`,
  `LagClasses`, `Search` and `Grid` can only be built valid. Continue this
  for every configuration type (e.g. SASIM's `Schedule`, which is still
  checked at use).

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
