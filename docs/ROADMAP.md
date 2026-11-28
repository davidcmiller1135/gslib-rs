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
| file format | 🟡 | `io::parse_gslib` / `format_gslib`; trimming and non-finite policies on extraction (`io::Extraction`); missing results written with a code (`io::MissingCode`) | CLI trimming flags and sentinel guard; Fortran fixed-format output |
| `.par` files | ⬜ | (CLI flags) | read GSLIB parameter files so existing workflows run unchanged |
| `nscore` | 🟡 | `transform::normal_scores`; parity fixtures (F-007), with trimming (F-011) | reference distribution option; rank output (negative column); random despiking as a legacy policy; writing the transformation table; CLI subcommand |
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
| `gamv` | ✅ | `variogram::{gamv, Direction, CentredLags, indicator_continuous, indicator_categorical}`, CLI `variogram`; 17 parity fixtures (F-008, F-011, F-027) | a parallel `gamv`; flagged corrections (`docs/spec/gamv.md`) |
| `gam` | ⬜ | | variograms of gridded data along grid directions |
| `varmap` | ⬜ | | variogram maps |
| `vmodel` | 🟡 | `VariogramModel::gamma` | write model values along directions |
| `vargplt` | 🟡 | `plot::VariogramPlot` | multiple directions on one plot |

### Estimation

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `kb2d` / `kt3d` | 🟡 | `kriging::{simple_kriging, ordinary_kriging}` (with weights: `*_solution`), `krige` CLI; point SK/OK parity fixtures with search and weights (F-015) | **block kriging** (discretisation); **kriging with a trend** (KT, `ktsol`) and **external drift**; cross-validation / jackknife; gridded output; locally varying mean |
| search (`setsupr`/`srchsupr`) | 🟡 | `search::{DataSearch, SearchEllipsoid, BruteForce, SuperBlock, NodeSearch}`: anisotropic ellipsoid, `ndmin`/`ndmax`, octants, an exact super-block, `sgsim`'s node search; `kt3d` parity fixtures (F-015) and the two-part search (`assign_to_nodes`, `DataConditioning::AssignToNodes`) (`docs/spec/search.md` steps 1–5), with CLI flags | k-d tree (ROADMAP §5) |
| `cokb3d` | 🟡 | `kriging::simple_cokriging` | ordinary cokriging (both constraint forms); collocated cokriging; check the linear model of coregionalisation |
| `ik3d` | ⬜ | (single threshold via `sisim`'s kriging) | multiple thresholds; median IK; **order-relation corrections**; categorical variables; soft data |

### Simulation

| GSLIB | Status | gslib-rs | Missing for parity |
|---|---|---|---|
| `sgsim` | 🟡 | `simulation::sgsim` on a `Grid` with separate data and node searches (`search::NodeSearch`, GSLIB's covariance-ranked "spiral" table), `sgsim` CLI; statistical parity fixtures (F-018) and GSLIB's random stream reproduced bit for bit (F-017) | **covariance lookup table** for the kriging matrices; multiple-grid refinement; assigning data to grid nodes; locally varying mean; collocated cokriging with a secondary variable |
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
| RNG | ACORN, global state, one stream across all realizations (F-016) | caller-supplied `rand::Rng`; by default ChaCha8, realization k on stream k of the master seed (D-036, D-037) | ACORN as a `rand` generator with GSLIB's seeding and warm-up, for exact-stream parity (D-038; `docs/spec/rng.md`) |
| Lag classes | centres at multiples of `xlag`, inclusive tolerance `xltol`, may overlap; `xltol ≤ 0` silently becomes `xlag/2` | `CentredLags` reproduces this (F-008); a non-positive tolerance is an error; equal-width `LagClasses` also available | n/a; `sasim` should move to grid lag vectors, as GSLIB's SASIM uses (D-042) |
| Zero-distance lag (`gamv`) | also pairs each datum with itself, diluting γ there towards 0 (F-008) | distinct coincident pairs only | none (GSLIB's value is recoverable via F-008's relations); D-017 |
| Super-block selection misses data (kt3d) | `picksup` tests 27 lattice points between blocks, so a thin rotated search ellipsoid can reach into a skipped block: data inside the ellipsoid are ignored (F-012, confirmed in F-015) | exact block bound, never skips an in-range datum (D-025) | none (a defect) |
| Ordinary kriging with one datum (kt3d) | UNEST, since OK's unbiasedness counts as a drift term (`1 ≤ na ≤ mdt = 1`, F-014) | estimates with weight 1, variance 2γ(h) | possible legacy policy |
| Simple kriging with one datum (kt3d) | returns the datum and the OK variance 2(C(0) − C(h)) (F-014) | true SK: m + (C(h)/C(0))(v − m), variance C(0) − C(h)²/C(0) | none (a defect) |
| Node-search extent (`sgsim`) | the covariance table is capped at ±35 nodes in x and y and 0 in z (`MAXCTX = MAXCTY = 71`, `MAXCTZ = 1`), so nodes inside the radius beyond the cap, and every node above or below, are silently ignored (F-013) | the whole search ellipsoid (S10) | `TableExtent::Gslib { max }` |
| Node-search order (`sgsim`) | unstable sort by `−C(h) + 10⁻¹⁰·q(h)` (F-012) | decreasing `C(h)`, then `q(h)`, then offset index: a total order (D-026) | GSLIB's order: deferred to the exact-stream parity spec |
| Two data on one node (`sgsim`, `sstrat = 1`) | keeps the datum closer to the node (L1; the later on a tie) and drops the other, with only a debug-file warning (D-028) | error naming both data and the node | `OnSharedNode::KeepClosest` |
| Data outside the grid (`sgsim`, `sstrat = 1`) | `getindx` clamps the datum to the edge node, however far away (F-013) | error naming the datum | `OutsideGrid::ClampToEdge` |
| Order of equidistant search data | unstable quicksort over super-block order (F-012) | by distance, then input row (D-024) | GSLIB's order: deferred to the exact-stream parity spec |
| Trimming interval | mostly `tmin ≤ v < tmax`; `gamv` uses four intervals: reading `[tmin, tmax)`, pairing `[tmin, tmax]`, the omnidirectional swapped ordering `[tmin, tmax)`, indicators `(tmin, tmax]` (F-011, F-021, F-026); closed tests also in `sasim`, `ik3d`, `sisim_gs`, `bivplt` (F-009) | `tmin ≤ v < tmax` everywhere (D-019) | per program, only if a fixture shows it matters |
| NaN / ±Inf in data | NaN kept as data; ±Inf trimmed by ±1e21 limits (F-010) | explicit policies, default error, applied before trimming (D-020, M9) | `Missing` / `ClampToRange` are opt-in policies |
| `lusim` unestimated code | −1.0 (collides with real data) | −999 by default, configurable (D-023) | `--unest -1` |
| Omnidirectional pair counts (`gamv`) | each pair counted in both orderings, so `np` is doubled; with tail ≠ head this averages the two orderings | each pair reported once; both orderings still averaged, so values equal GSLIB's (D-043) | n/a (reporting convention); D-017 |
| Angle tolerances ≤ 0 (`gamv`) | silently 45° | `Error::InvalidInput` | `NonPositiveAngleTolerance::Gslib45` (D-044) |
| Types 6/7, omnidirectional, tail ≠ head (`gamv`) | the swapped term is computed from the unswapped values; can give `NaN`/`Inf` | uses the swapped values; non-finite → error | `OmniRelativeLog::GslibDuplicate` (D-047) |
| Undefined lag values (`gamv`) | γ = 0 (types 4, 5) or −999 (type −3); standardisation silently skipped for a sill ≤ 0 | `gamma: None`; standardising with a sill ≤ 0 is an error | GSLIB's values at output (D-049) |
| Indicators (`gamv`, types 9/10) | `I = 1` if `v ≥ c`; categories by `int(v + 0.5)`; categorical sill 0 (F-022) | `I = 1` if `v ≤ c` (IK3D's convention); exact integer categories; sill p(1−p) | gamv's conventions in the fixture adapter (D-050). Flag: `sasim` uses `v < t` |
| `gamv` edge cases | `atol = 90` drops pairs exactly perpendicular to the azimuth (F-023, confirmed by fixture); tail/head depends on row order for mixed-sign pairs (F-024) | same (reproduced) | flagged: correct later |
| Missing indicators (`gamv`) | a missing indicator is stored as `tmin − 1e-20`, which rounds to `tmin` and is paired as data (F-025, confirmed by fixture) | `None`, never paired | none (a defect; reproduced only by the fixture adapter) |
| π in angle conversion | `3.141592654` (`setrot`), `3.14159265` (`cova3`, `gamv`) | exact `π` | `PiConstant::Gslib`: planned for exact-stream parity, and needed for exact `gamv` pair counts at tolerance edges (F-020, D-046). Flagged: a symmetric boundary rule |
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
- **Parallelism at the edge.** The core stays sequential and pure per item
  (one target, one block of pairs, one realization). Parallel maps over
  those items live in a boundary module behind an opt-in feature, and their
  output never depends on the thread count (`docs/spec/parallelism.md`).

## 5. Beyond GSLIB

Once parity is solid:

- **Parallel execution** (`docs/spec/parallelism.md`, D-029–D-033), with output
  bit-identical for any thread count:
  - kriging at many targets (done: `parallel::try_map`, used by `krige`);
  - experimental variograms, with a fixed-block deterministic reduction and
    exact pruning of pairs beyond the largest lag through the super-block
    (done: pruning, D-034, gives 14× on 20 000 data with short lags and 1.6× with lags to 20% of the domain; parallel blocks add 3.9× on 4 cores);
  - realizations (done: `parallel::realizations`, realization k on stream k
    of the master seed, D-037).
- **GPU acceleration** (opt-in, evidence first): only the regular O(N²) pair
  loop is a good fit. Obstacles: slow `f64` on consumer GPUs and none in
  WebGPU, `unsafe`-based CUDA bindings, non-deterministic atomic reductions,
  and no GPU CI runners. Each kernel must be checked against the CPU result.
- **Automatic variogram fitting** (weighted least squares with constraints).
  GSLIB leaves this to the user.
- **Fast neighbourhoods:** a k-d tree as a third exact `Neighbourhood`, and
  streaming over large grids. (The exact super-block is done, about 48×
  faster than brute force on 20 000 data. An *approximate* search would
  break D-025's exactness, so it could only be a named opt-in.)
- **True-cone direction tolerance** for `gamv` (angle to the 3-D direction
  vector, cylindrical bandwidth), alongside GSLIB's window (D-045).
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
- [x] ACORN as the generator and 69069 as its seed: `nscore` (F-005), and
  every stochastic program from its source (`sgsim`, `sisim`, `sasim`,
  `lusim`, `ellipsim`, `draw`, …; 69069 in `SGSIM.PAR` and `SISIM.PAR`).
  F-016. §2, `cli/src/main.rs`
- [ ] Grid conventions: `xmn` is the first cell's centre, x-fastest order:
  `simulation/grid.rs`
- [x] Search parameter names (`ndmax`, `nodmax`, `radius`): `ndmin, ndmax`,
  `noct`, `radius radius1 radius2`, `sang1..3` in `KT3D.FOR`; `nodmax`,
  `sstrat` in `SGSIM.FOR`. F-012.
- [x] GSLIB's super-block selection (`picksup`) can miss data inside the
  search ellipsoid under a rotated anisotropic search: confirmed by the
  `kt3d/picksup_probe` fixture (150 of 1 764 nodes). F-015.
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
