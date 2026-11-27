# Spec: Search neighbourhoods

Status: **in progress** (grill-me session started 2026-11-27; evidence F-012).

## Purpose

Decide how gslib-rs chooses the conditioning data (and, in simulation, the
previously simulated nodes) for each estimate, so that `kt3d` and `sgsim`
can reach parity with GSLIB, and the choice is deterministic, documented
and fast.

## Fixed by AGENTS.md and the source

- **Search metric:** the search has its own anisotropic ellipsoid (radius,
  two ratios, three angles), independent of the variogram model, built with
  the existing `Anisotropy` and `Angles` types. A datum is inside iff its
  squared anisotropic distance is `≤ radius²` (F-012).
- **Octants:** GSLIB's octant definition and boundary rules, applied before
  the `ndmax` cut (F-012).
- **Too few data** (`ndmin`, "search exhaustion" in `AGENTS.md`): a typed
  error by default. GSLIB's behaviour (`kt3d`: UNEST; `sgsim`: node left
  unsimulated) is an opt-in policy.
- **A neighbourhood trait:** brute force, super-block and later a k-d tree
  qualify as genuine substitution (`AGENTS.md`, *No ceremony*).
- **Inputs are explicit:** the search parameters and the data are passed
  in; there is no hidden state.

## Decisions so far

None yet.

## Open questions (in dependency order)

1. Order of equidistant data (decides what "the first `ndmax`" means).
2. Must the super-block search equal brute force exactly, or reproduce
   GSLIB's `picksup` (F-012 hypothesis)?
3. Is `sgsim`'s node search (covariance-ordered spiral, `nodmax`, two-part
   search) part of this spec?
4. Is `kt3d`'s debug output used as a direct oracle for the selected data?
