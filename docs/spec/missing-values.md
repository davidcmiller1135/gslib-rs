# Spec: Missing values and trimming limits

Status: **specified** (decisions D-018 to D-023 in `decisions.md`; evidence
F-009, F-010).

## Purpose

Decide how gslib-rs represents, detects, and writes missing and unestimated
values, so that GSLIB data files with trimming limits and sentinel codes are
read and written correctly, and no sentinel is ever mistaken for data.

## Fixed by AGENTS.md

- Missing values are `Option` in memory. Sentinels (`-999`, `UNEST`) and
  trimming limits (`tmin`/`tmax`) exist only at the I/O boundary.
- Never silently coerce `NaN`/`Inf`.
- Divergences from GSLIB are documented in `docs/ROADMAP.md` §2.

## Evidence

F-009 (GSLIB's conventions per program), F-007, F-008.

## Decisions so far

- Missingness is decided at **extraction**: `GslibData` stays raw and
  lossless; applying trimming to a column yields `Option<f64>`; algorithms
  see only present data; row-aligned results are `Option` (D-018).
- A value is in range when `tmin ≤ v < tmax`, in every program (D-019).
- Non-finite values: separate policies for NaN (`Error` | `Missing`) and ±Inf
  (`Error` | `Missing` | `ClampToRange`), both defaulting to `Error`
  (D-020; GSLIB's behaviour in F-010).
- CLI: no limits by default, but a sentinel guard fails on −999, −99 or
  |v| ≥ 1e21 unless `--tmin`/`--tmax` or `--no-trimming` is given (D-021).
- Missingness is per cell; coordinates are never trimmed (D-022).
- Writers emit a configurable code for `None`, defaulting to GSLIB's code per
  program (−999 or −99), except `lusim`, which defaults to −999 rather than −1
  (D-023).
- **Derived rule (from D-020):** the non-finite policy applies *before*
  trimming. A ±Inf is therefore an error by default even when finite limits
  would have trimmed it, where GSLIB trims it (F-010). This is a documented
  divergence.

## Requirements

- **M1 Lossless files.** `format_gslib(parse_gslib(t))` reproduces every value
  of `t` exactly, including sentinel codes such as −999.
- **M2 One door to the algorithms.** Extracting samples or values from a
  `GslibData` requires a trimming argument (possibly "none"). No API turns a
  column into `Sample`s without deciding missingness.
- **M3 Row alignment.** A result per input row is `Vec<Option<f64>>` with
  exactly one entry per row, `None` where the input was missing.
- **M4 Half-open interval.** Given limits `(tmin, tmax)`, a value `v` is
  present iff `tmin ≤ v < tmax`: `v = tmin` is present, `v = tmax` is not.
- **M5 Non-finite values.** Extraction with the default policy fails on a NaN
  or ±Inf with an error naming the 1-based row and column. With `Missing`,
  the value becomes `None`. With `ClampToRange` (±Inf only), +Inf becomes the
  largest finite present value of the column and −Inf the smallest; if the
  column has no finite present values, extraction fails.
- **M6 Sentinel guard (CLI).** Given a data column containing −999.0, −99.0,
  or any |v| ≥ 1e21, and no `--tmin`, `--tmax` or `--no-trimming`, the
  command fails before computing anything, naming the column and the first
  offending row.
- **M7 Per-cell missingness.** In a record with several variables, each is
  trimmed independently; coordinates are never trimmed.
- **M8 Output codes.** Writing a row-aligned result replaces each `None` with
  the configured code; the default per program follows D-023. A written
  file's present values are bit-identical to the computed values.
- **M8a Unambiguous codes (derived from M8 and `AGENTS.md` precedence 2).**
  Writing fails if a present result equals the code, or is not finite: the
  first would read back as missing, the second would fail extraction. The
  error names the row, so the user can choose another code.
- **M9 Order of rules.** For each cell: the non-finite policy (M5) first,
  then trimming (M4). A cell is present only if both keep it.

## Open questions

None.

## Implementation plan

Each step is one reviewable PR; the checks in `AGENTS.md` must pass.

1. **Extraction API** (`gslib_rs::io`; done: `Trimming`, `OnNan`,
   `OnInfinite`, `Extraction`, `GslibData::extract`, `samples` returning
   `PresentSamples`, `Error::NonFinite`; 10 tests, each caught a planted
   bug in a mutation check).
   - Add `Trimming` (half-open limits, or none) and `NonFinite` policies
     (NaN: `Error`/`Missing`; ±Inf: `Error`/`Missing`/`ClampToRange`), and
     `GslibData::extract(column, &policy) -> Result<Vec<Option<f64>>>`.
   - Change `samples` to take the policy and return the present samples
     with their row indices (for M3).
   - Keep `column` as documented raw access.
   - Tests: M1 (round trip with −999, existing proptest extended); M4 (`v =
     tmin` present, `v = tmax` absent); M5 (each policy, with row and column
     in the error, clamp to the trimmed finite range, no finite values);
     M7 (two variables in one record); M9 (Inf with finite limits under
     `Error`).
2. **Row-aligned writing** (done: `GslibData::with_optional_column`,
   `MissingCode`, and `PresentSamples::align` to spread per-sample results
   back over the rows). `GslibData::with_optional_column(name,
   &[Option<f64>], code)`. Tests: M3 (one entry per row) and M8 (codes;
   present values bit-identical).
3. **CLI** (done: `cli/src/missing.rs`, on `variogram`, `krige`, `sgsim`
   and `sisim`; 4 unit and 4 end-to-end tests).
   - `--tmin`/`--tmax`/`--no-trimming`, `--nan`, `--inf`, `--unest` on
     commands that read data.
   - **Deferred: `--unest`.** No current command writes one result per data
     row (kriging and simulation write per target or node), so the flag
     would do nothing yet. It arrives with the first row-aligned output,
     e.g. an `nscore` command or a kriging `OnSingular` skip policy, using
     `MissingCode` with that program's default (D-023).
   - The guard checks the value column only: coordinates are never trimmed
     (D-022), so `--tmin` could not help with a coordinate of −999. A
     missing limit is unbounded on its side (`--tmin -998` alone keeps
     everything from −998 up).
   - The sentinel guard (M6).
   - Integration tests on the `with_sentinels` suite: fails without flags,
     names the column and row; succeeds with `--tmin -998`; writes −999 for
     trimmed rows where a row-aligned output exists.
4. **Parity fixtures with trimming** (done: four fixtures, F-011).
   - Work: `oracle-case` `nscore` and `gamv` parameters gain optional
     trimming limits (a recipe change, so regenerate).
   - Cases:
     - `nscore` on `with_sentinels` with `tmin = −998`: GSLIB writes −999
       for trimmed rows, gslib-rs returns `None`, and the scores of present
       data match;
     - `nscore` and `gamv` cases with a limit exactly equal to a datum, to
       check M4 against GSLIB's half-open tests.
   - Depends on step 1.
5. **ROADMAP §2** (done with this spec): rows for NaN/Inf handling (D-020,
   M9) and the `lusim` default code (D-023). Closed-interval programs are
   added as each is ported (D-019).
