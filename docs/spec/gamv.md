# Spec: gamv (directional variograms of scattered data)

Status: **specified** (decisions D-041 to D-052; findings F-019 to F-024 in
`decisions.md`). GSLIB facts were read from
`vendor/gslib-2.0/src/GAMV/GAMV.FOR`, except F-020, which is a computation
that the fixtures in plan step 6 confirm.

This spec was first written on a branch based on `main`, before PR #4 was
known. It was then ported onto PR #4. Where the two disagreed, PR #4's
earlier decisions stand (user, 2026-11-28): D-017 (no self-pairs, each pair
reported once, no legacy policy) and D-031/D-034 (the block-ordered pair
enumeration).

## Purpose

Bring `variogram` to parity with GSLIB 2.0 `gamv` for scattered 3-D data
(D-041):
- directions;
- all ten variogram types;
- multiple variables, with tail and head variables;
- indicators;
- sill standardisation;
- output.

This finishes milestone M3 (PLAN.md).

## Current state (before this work)

- `CentredLags` reproduces gamv's lag classes (F-008), behind `LagScheme`.
- `experimental_variogram`: the omnidirectional semivariogram of one
  variable, with each pair counted once and no self-pairs (D-017). It
  enumerates pairs in fixed row blocks with super-block pruning (D-031,
  D-032, D-034).
- Missing values and trimming: `GslibData::extract`, `[tmin, tmax)` (D-018
  to D-023, F-011).
- Fixtures: nine omnidirectional `gamv` cases, checked through F-008's exact
  relations.
- CLI `variogram`: `--max-dist`/`--n-lags` (equal-width `LagClasses`), plus
  the missing-value flags.

## Fixed by AGENTS.md or earlier decisions

- Deterministic program: numerical parity within D-010/D-015's tolerances,
  with exact structure.
- GSLIB's silent coercions are errors by default, with opt-in legacy
  policies (`AGENTS.md`, *Legacy behaviour*).
- Azimuth is clockwise from north; dip is negative downward (F-004).
- No self-pairs; each pair is reported once (D-017).
- Pairs are enumerated as in D-031, D-032 and D-034.
- Working principle for this spec: reproduce GSLIB first, and flag
  corrections and extensions for later (D-042).

## Decisions (summary)

| # | Topic | Default (gslib-rs) | Legacy (GSLIB) |
|---|---|---|---|
| D-041 | Scope | all of gamv | (n/a) |
| D-042 | Lags | `CentredLags`; `sasim` keeps `LagClasses` (flagged) | (n/a) |
| D-043 | Omnidirectional, tail ≠ head | both orderings averaged internally; pairs reported once | (n/a: equal to GSLIB) |
| D-044 | Tolerances | `xltol` omitted or ≤ 0 at the boundary means `xlag/2`; `atol`, `dtol ≤ 0` are errors | `NonPositiveAngleTolerance::Gslib45` |
| D-045 | Window geometry | GSLIB's four tests exactly (the true cone is a flagged extension) | (n/a) |
| D-046 | π | exact π, GSLIB's operation order (boundary rule flagged) | `PiConstant::Gslib` |
| D-047 | Types 6/7, omnidirectional | swapped values; non-finite → error | `OmniRelativeLog::GslibDuplicate` |
| D-048 | Trimming | `[tmin, tmax)` (D-019) | gamv's three intervals (fixture adapter) |
| D-049 | Undefined values | `gamma: None`; sill ≤ 0 with standardisation is an error | GSLIB's 0/−999, skipped standardisation (adapter, output) |
| D-050 | Indicators | `I = 1` if `v ≤ c`; exact integer categories; sill p(1−p) | gamv's conventions (adapter) |
| D-051 | API | `gamv(&Data, &GamvRequest)` on PR #4's enumeration | `GamvPolicies::gslib()` |
| D-052 | CLI | `--nlag --xlag --dir --vario`, CSV or gamv layout | `--gslib-compat` |

## Requirements

`EPSLON = 1e-20`. Line numbers refer to `GAMV.FOR`.

### Directions (D-044 to D-046)

- **W1.** `Direction::new(azm, atol, bandh, dip, dtol, bandv)` errors on
  non-finite input, on `atol ≤ 0` or `dtol ≤ 0`, and on `bandh < 0` or
  `bandv < 0`. The 45° substitution exists only at the boundary, under
  `NonPositiveAngleTolerance::Gslib45`.
- **W2.** A pair is accepted by a direction exactly when GSLIB's four tests
  pass, with GSLIB's comparisons (ll. 545–582):
  1. azimuth: `|dcazm| ≥ cos(atol)`, with `dcazm = 1` when the horizontal
     length is `< EPSLON`;
  2. horizontal bandwidth;
  3. dip, on the signed horizontal length, skipped for the zero-distance
     class;
  4. vertical bandwidth.

  The separation is `x_j − x_i` for data `i < j`, as in GSLIB.
- **W3.** Unit vectors and thresholds are computed as `(90 − azm)·π/180`
  and `cos(atol·π/180)`, and likewise for dip. π is exact by default, and
  `3.14159265` under `PiConstant::Gslib`.
- **W4.** A direction is omnidirectional iff `atol ≥ 90`. Each accepted pair
  then contributes both orderings to the sums, except for the
  cross-semivariogram, and is reported once (D-043).
- **W5.** Tail and head: datum i is the tail when `dcazm ≥ 0` and `dcdec ≥
  0`; otherwise j is. This is not antisymmetric (F-024).

### Measures (D-047, D-049, D-051)

- **M1.** The accumulation and finalisation of each kind match ll. 628–871:
  - the factor ½ for types 1, 2, 6, 7 and 8;
  - centring for 3 and 4;
  - `√(sill_t·sill_h) − C` for −3;
  - division by the squared lag mean, with no ½, for 5.
- **M2.** Omnidirectional types 6 and 7 use the swapped values by default.
  `OmniRelativeLog::GslibDuplicate` reproduces ll. 730 and 758.
- **M3.** A non-finite result is an error in every mode.
- **M4.** `gamma` is `None` when a lag is empty, or when:
  - a type-4 lag has `hv·tv < EPSLON`;
  - a type-5 lag has a squared mean `< EPSLON`;
  - a type −3 lag has no sill (a negative sill counts as none).

  `pairs` is still reported.
- **M5.** Standardisation (`VariogramSpec::standardize`) divides a
  semivariogram with tail = head by the variable's sill, before the factor
  ½. Requesting it for anything else, or with a sill ≤ 0 or missing, is an
  error.
- **M6.** Types 6 and 7 skip orderings outside the measure's domain and
  report them per lag in `skipped`.
- **M7.** Each lag reports `mean_distance`, `gamma`, `pairs`, `tail_mean`,
  `head_mean`, `tail_variance` and `head_variance` (correlogram only), and
  `skipped`.

### Data and indicators (D-048, D-050)

- **D1.** `Data` holds coordinates and `Option<f64>` columns of equal
  length. A pair contributes to a variogram only when the values it needs
  are all present.
- **D2.** `indicator_continuous(var, c)` gives `I = 1` if `v ≤ c`.
  `indicator_categorical(var, k)` requires integral `v` and `k` and gives
  `I = 1` if `v = k`. Each returns the variable and its sill p(1 − p), or
  `None` with no present value.

### Output and CLI (D-052)

- **O1.** Results are ordered [variogram][direction][lag], with `nlag + 2`
  lags.
- **O2.** The CLI writes long-form CSV (empty cells for `None`) or
  `--format gamv`: GSLIB's layout and titles at full precision. Empty lags
  are written as zeros, and undefined values as −999, or as GSLIB's 0/−999
  under `--gslib-compat`.
- **O3.** `--gslib-compat` selects every legacy policy listed above.

### Parity (harness)

- **P1.** Directional and multi-variable fixtures match GSLIB:
  - pair counts exactly, through D-017's relations: GSLIB's count is
    gslib-rs's, doubled for omnidirectional types other than 2. GSLIB's
    zero-distance class adds one self-pair per datum with both values
    present, also doubled when omnidirectional;
  - values within D-015's tolerances. The zero-distance class is compared
    through the self-pair relations of F-008 for semivariograms with
    tail = head, and otherwise only where it holds no coincident data.
- **P2.** Cases where the legacy formula gives `NaN`/`Inf` (D-047) are
  expected errors.

## Flagged for later (corrections and extensions)

- A true-cone window (D-045; ROADMAP §5).
- A symmetric boundary rule for edge pairs (D-046).
- `atol = 90` drops pairs exactly perpendicular to the azimuth (F-023).
- Tail/head order depends on row order when the azimuth and dip cosines
  have opposite signs (F-024).
- `sasim`: move to grid lag vectors (D-042), and adopt the `v ≤ c`
  indicator convention (D-050).
- Exact fixed-format output, with the `.par` adapter (D-052).
- A parallel `gamv` by row blocks (D-051).

## Open questions

None.

## Implementation plan

Each step is one reviewable PR (here: commits on PR #5), and every check in
`AGENTS.md` passes.

1. **Direction windows.** (**Done.**) `variogram/directional.rs`: `Direction`,
   `DirectionWindow`, `PairOrder`, `geometry::PiConstant`.
   - Tests:
     - W1 errors;
     - azimuth sense and order; vertical pairs; bandwidths at their limits;
       dip sign; dip measured in the pair's own plane; the vertical band
       when omnidirectional; the zero class skipping the dip test;
     - F-020 under both π, F-023, F-024;
     - properties: a pair and its reverse are accepted alike; tolerances
       above 90° accept everything.
2. **Kernel.** `variogram/gamv.rs`: `Data`, `GamvRequest`, `GamvPolicies`,
   `gamv` for the semivariogram on PR #4's block enumeration, with D-043.
   - Tests:
     - equal to `experimental_variogram` for one omnidirectional
       semivariogram;
     - linear trend; γ = 0 across a constant direction;
     - tail and head order;
     - output order;
     - permuting rows leaves an omnidirectional cross measure unchanged
       (D-043);
     - property: missing values never pair.
   - **Status: done.** `GamvLag::orderings` reports GSLIB's `np` (the
     orderings summed), next to `pairs`, which counts each pair once.
3. **All measures.** Types 2 to 8 and −3, standardisation, D-047, `None`
   for undefined values, errors on non-finite results.
   - Tests: hand-computed cases for each measure; the
     cross-semivariogram of a variable with itself; −3; standardisation
     errors; skip counts; both types 6/7 policies; correlogram bounds.
   - **Status: done.**
4. **Indicators.** `indicator_continuous` and `indicator_categorical`.
   - Tests: the cut-off value; the categorical sill; negative codes;
     non-integer errors; missing values; a standardised indicator
     semivariogram.
   - **Status: done.**
5. **CLI.** D-052 on PR #4's CLI, keeping its missing-value flags.
   - Tests: the CSV layout; the gamv layout under `--gslib-compat`;
     `atol = 0` with and without compat; the `--standardize` error;
     indicator conventions.
6. **Fixtures.**
   - Work:
     - extend `oracle-case`'s gamv parameters with optional directions,
       variables, variograms and `isill`. Existing cases must regenerate
       byte for byte;
     - a harness adapter with gamv's legacy trimming and indicators;
     - new cases: directions on a regular grid with edge pairs (F-020,
       F-023), dip in both signs, narrow bandwidths, every type with
       tail = head and tail ≠ head, omnidirectional cross measures,
       indicators with ties at the cut-off, missing indicator values with
       several variables, the legacy types 6/7 defect (expected error), and
       `cluster.dat` with `GAMV.PAR`.
   - Outcome: confirm or correct F-020, F-023 and F-024 against reference
     runs.
