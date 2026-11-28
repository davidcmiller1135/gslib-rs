# Spec: kt3d completeness (block kriging, drift, cross-validation)

Status: **in progress** (grill session, 2026-11-28). Decisions start at
D-054 in `decisions.md`. GSLIB facts are read from
`vendor/gslib-2.0/src/KT3D/KT3D.FOR` (`readparm`, `kt3d`),
`KT3D.INC` and `vendor/gslib-2.0/src/GSLIB/KTSOL.FOR`, cited by line.

## Purpose

Bring `kriging` to parity with GSLIB 2.0 `kt3d` beyond point SK/OK
(PLAN C4, milestone M4).

## Current state (before this work)

- `kriging::{simple_kriging, ordinary_kriging}` and their `*_solution`
  forms estimate one point from given samples. OK is written in variograms.
- The search (`docs/spec/search.md`) selects data as `kt3d` does, except
  for `picksup`'s misses (D-025, F-015). Too few data is `Error::TooFewData`,
  or `None` under `OnTooFewData::Unestimated`.
- Nine point SK/OK `kt3d` fixtures, compared within `k·cond(A)·ε`, k = 100
  (D-035).
- Linear solves (`linalg::solve_spd`, `solve_general`) fail only on a
  failed Cholesky, an exact zero pivot, or a non-finite result.

## Fixed by AGENTS.md or earlier decisions

- Deterministic program: numerical parity within documented tolerances
  (D-010, D-015, D-035).
- GSLIB's silent coercions are errors by default, with opt-in legacy
  policies (`AGENTS.md`, *Legacy behaviour*).
- Missing values are `Option` in memory; `UNEST` exists only at output
  (D-018 to D-023).
- Known `kt3d` divergences stand: SK and OK with one datum (F-014), and the
  exact super-block (D-025).

## Scope (D-054)

In scope:
- point and block kriging;
- kriging with a trend (the nine polynomial drift terms);
- external drift (`ktype = 3`);
- simple kriging with a locally varying mean (`ktype = 2`);
- trend estimation (`itrend = 1`);
- cross-validation (`koption = 1`) and jackknife (`koption = 2`);
- gridded output;
- the `OnSingular` policy and the other UNEST paths.

Deferred, but required for parity: the undocumented IK-type distribution
built from the kriging weights (negative `koption`). It belongs to the
`ik3d`/`postik` spec (PLAN D3).

## GSLIB behaviour (from the source)

- **Drift count** (`kt3d`, lines 508–519): `mdt = 1 + Σ idrif` for OK and
  KT, plus 1 for external drift; 0 for SK and LVM SK, whose `idrif` are
  forced to 0.
- **Block discretisation** (lines 528–558): `nxdis·nydis·nzdis` points
  (each clamped to ≥ 1, at most `MAXDIS = 64`), at the centres of equal
  sub-cells, as offsets from the block's lower-left corner.
- **Block covariance** (lines 566–588): `cbb = C(0)` for a point;
  otherwise the mean of `C(xᵢ, xⱼ)` over all pairs of discretisation
  points, with the nugget removed on the diagonal (`i = j`).
- **Right-hand side** (lines 792–812): the mean covariance between the datum
  and the discretisation points, with the nugget removed where a datum is
  within `EPSLON = 1e-6` (squared distance) of a point.
- **Drift rows** (lines 816–930): data coordinates relative to the block's
  lower-left corner, scaled by `resc = 4·radius²/max(covmax, 1e-4)`
  (`2·radius/…` when `radius² < 1`, lines 479–483); the right-hand side is the
  block mean of each drift function (`bv`). External drift is scaled by
  `resce = covmax/max(extest, 1e-4)`.
- **OK row** holds `unbias = C(0)` instead of 1 (line 573).
- **Solve** (`ktsol`): Gaussian elimination, partial pivoting; singular if
  `n ≤ 1` or a pivot has `|p| < 1e-11` (absolute).
- **Results** (lines 968–990): `est = Σ sᵢvᵢ` (OK, KT, external drift) or
  `m + Σ sᵢ(vᵢ − m)` (SK, LVM SK); `estv = cbb − Σ sⱼrⱼ` over all rows,
  not clamped.
- **UNEST paths:** too few data (`na < ndmin`); `1 ≤ na ≤ mdt` (lines
  707–716); a singular system (line 968); an external-drift or LVM value at
  the target outside `[tmin, tmax)` (lines 659–667, read at 657–663).
- **One datum, `mdt = 0`** (lines 722–757): `est = v₁` and
  `estv = cbb − 2·cb + C(0)` (F-014).
- **Flagged from the source, to be decided below:**
  - LVM SK (`ktype = 2`) uses the target's mean for every datum's residual:
    `skmean = extest`, then `vra(j) − skmean` (lines 975–986); the data's
    local means are never read (`readparm`, lines 304–313, reads `ve` only
    for `ktype = 3`).
  - With `ndmin = 0`, SK and no data found, the one-datum branch reads a
    stale `xa(1)` (lines 707–757).
  - Cross-validation drops every datum within L1 distance `1e-6` of the
    target, not only the datum itself (lines 685–688).
  - `itrend = 1` zeroes the data rows of the right-hand side (lines
    942–947), and the variance is computed from the same zeroed vector.

## Decisions (summary)

| # | Topic | Default (gslib-rs) | Legacy (GSLIB) |
|---|---|---|---|
| D-054 | Scope | all of `kt3d` except the IK-type distribution (deferred to D3) | (n/a) |

## Acceptance criteria

To be written as decisions are made.
