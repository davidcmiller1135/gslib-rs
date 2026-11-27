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

- **Order** (D-024): candidates are returned by ascending squared
  anisotropic distance, then ascending input row. GSLIB's tie order is an
  opt-in for exact-stream parity, deferred to that spec.
- **Accelerators are exact** (D-025): super-block (and any later
  accelerator) returns exactly the brute-force sequence. Its block
  selection may never exclude a block holding an in-range datum. GSLIB's
  possible miss (F-012) is verified with a fixture and, if confirmed,
  documented as a fixed defect.
- **Node search in scope** (D-026): `sgsim`'s search of previously
  simulated nodes is specified and built here too. Order: decreasing
  covariance, then increasing anisotropic distance, then the offset's
  linear index. It keeps at most `nodmax` nodes, with GSLIB's `srchnd`
  octant convention.
- **Oracle** (D-027): `kt3d`'s debug output is a direct oracle. Per
  estimate, the selected set (exact up to D-024 ties), its order, and the
  weights (D-010 tolerance) are compared.

## Open questions (in dependency order)

1. ~~Order of equidistant data~~ (D-024).
2. ~~Super-block vs brute force~~ (D-025).
3. ~~Node search in scope~~ (D-026).
4. ~~Debug output as oracle~~ (D-027).
5. (New, from D-026.) Two-part search (`sstrat = 1`): what happens when two
   data map to the same node? GSLIB keeps the closer one and warns, dropping
   the other.
