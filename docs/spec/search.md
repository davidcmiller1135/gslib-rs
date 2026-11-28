# Spec: Search neighbourhoods

Status: **specified** (decisions D-024 to D-028; evidence F-012, F-013).

## Purpose

Decide how gslib-rs chooses the conditioning data, and in simulation the
previously simulated nodes, for each estimate. The choice must let `kt3d`
and `sgsim` reach parity with GSLIB, and be deterministic, documented and
fast.

## Fixed by AGENTS.md and the source

- **Search metric:** the search has its own anisotropic ellipsoid (radius,
  two ratios, three angles), independent of the variogram model, built with
  the existing `Angles`, `Ranges` and `Anisotropy` types. A datum is inside
  iff its squared anisotropic distance is `≤ radius²` (F-012).
- **Octants:** GSLIB's definition and boundary rules, applied before the
  `ndmax` cut (F-012).
- **Too few data** (`ndmin`; "search exhaustion" in `AGENTS.md`): a typed
  error by default. GSLIB's behaviour (`kt3d`: UNEST; `sgsim`: node left
  unsimulated) is an opt-in policy.
- **Silent limits are opt-in** (F-013):
  - by default, the node search covers the whole search ellipsoid; GSLIB's
    capped lookup table is an opt-in;
  - in the two-part search, a datum outside the grid is an error by
    default; GSLIB clamps it to the edge node (opt-in).
- **A neighbourhood trait:** brute force, super-block and later a k-d tree
  qualify as genuine substitution (`AGENTS.md`, *No ceremony*).
- **Explicit inputs:** search parameters, data, grid and model are passed
  in, and a search holds no mutable state.

## Decisions

- **Order** (D-024): candidates are returned by ascending squared
  anisotropic distance, then ascending input row. GSLIB's tie order is an
  opt-in for exact-stream parity, deferred to that spec.
- **Accelerators are exact** (D-025): super-block (and any later
  accelerator) returns exactly the brute-force sequence, and its block
  selection may never exclude a block holding an in-range datum. GSLIB's
  possible miss (F-012) is verified with a fixture and, if confirmed,
  documented as a fixed defect.
- **Node search in scope** (D-026): `sgsim`'s search of previously
  simulated nodes is specified and built here too. Order: decreasing
  covariance, then increasing anisotropic distance, then the offset's
  linear index. It keeps at most `nodmax` nodes, using `srchnd`'s octant
  convention.
- **Oracle** (D-027): `kt3d`'s debug output is a direct oracle. Per
  estimate, the selected set (exact up to D-024 ties), its order, and the
  weights (D-010 tolerance) are compared.
- **Shared nodes** (D-028): in the two-part search, two data on one node
  are an error by default. GSLIB's keep-closest rule (L1 distance, the
  later datum wins ties) is the opt-in policy `OnSharedNode::KeepClosest`.

## Vocabulary

- **Search ellipsoid:** `radius` (major axis) with `radius1`, `radius2` and
  `sang1..3`, as an `Anisotropy` from `Angles` and `Ranges { major: radius,
  minor: radius1, vertical: radius2 }`. Write `q(h)` for the squared scaled
  distance `scaled_distance(·)²`, in units where the ellipsoid is the
  sphere of radius `radius`.
- **Data search:** an ellipsoid, `min_data` (`ndmin`), `max_data`
  (`ndmax`), and optionally `per_octant` (`noct`).
- **Node search:** an ellipsoid, a `Grid`, a `CovarianceModel`,
  `max_nodes` (`nodmax`), optionally `per_octant`, and a table extent
  (whole ellipsoid by default, GSLIB's cap as an opt-in).

## Requirements

- **S1 Membership.** Datum `i` is a candidate for target `t` iff
  `q(xᵢ − t) ≤ radius²`. Tests: points just inside, exactly on, and just
  outside the ellipsoid along each principal axis, for rotated, anisotropic
  ellipsoids; GSLIB's convention that negative dip is downward.
- **S2 Order** (D-024). Candidates are sorted by `(q, row)`, a total order.
  Tests: equidistant data come back in row order whatever their input
  permutation (property test).
- **S3 Octants.** With `per_octant = n`, the sorted candidates are scanned
  in order. Each is kept if fewer than `n` earlier ones share its octant
  (`SRCHSUPR.FOR`'s definition and boundary rules). GSLIB also stops the
  scan after `8n` are kept, which is only a shortcut: every octant is then
  full, so the result is the same. Tests: one hand-placed datum on each
  octant boundary; cap per octant; the cut counting only data the octants
  kept.
- **S4 Truncation.** The result is the first `max_data` of the sequence
  after S3. Tests: `max_data` smaller, equal and larger than the candidate
  count.
- **S5 Too few data.** If fewer than `min_data` remain, the search fails
  with `Error::TooFewData { target, found, required }`, or, under the opt-in
  policy `OnTooFewData::Unestimated`, reports the location as unestimated
  (`None` in a row-aligned result, written with the program's code,
  D-023). Tests: both policies; `min_data = 0` never fails.
- **S6 Exact accelerators** (D-025). For any data, ellipsoid, target and
  parameters, `SuperBlock` returns exactly the `BruteForce` sequence.
  Property tests cover:
  - random 3-D data;
  - rotated, strongly anisotropic ellipsoids (ratios down to 0.05);
  - targets inside and outside the data's extent;
  - duplicate locations, a single datum, and no data;
  - the F-012 probe geometry as an explicit case.
- **S7 Purity and determinism.** Searching takes `&self` and the target and
  returns a fresh result. Results do not depend on call order or thread
  count. Test: a parallel map over targets equals the sequential one.
- **S8 Node search order** (D-026). The offset table holds every grid
  offset `δ ≠ 0` with `q(δ) ≤ radius²`, ordered by decreasing `C(δ)`, then
  increasing `q(δ)`, then the offset's linear index. Test: equals a
  brute-force enumeration and sort of all offsets.
- **S9 Node selection.** At node `n`, the table is scanned in order and
  offsets landing outside the grid are skipped. A node is kept if it is
  already simulated (or holds a datum, under the two-part search) and
  passes the per-octant cap, using `srchnd`'s convention (target − node,
  F-012). The scan stops at `max_nodes`. Tests: against a brute-force
  definition on random partially simulated grids; octant boundaries.
- **S10 Table extent** (F-013). By default the table covers the whole
  ellipsoid. Under the opt-in policy `TableExtent::Gslib { max_x, max_y,
  max_z }`, it is limited per axis to `min((max − 1)/2, n − 1)` offsets.
  Tests: a node beyond ±35 offsets is found by default and missed under
  the GSLIB policy with 71/71/1.
- **S11 Two-part search, shared nodes** (D-028). Assigning data to nodes
  fails on two data mapping to one node, naming both rows and the node.
  Under `OnSharedNode::KeepClosest`, the smaller L1 distance wins, and the
  later datum on a tie. Tests: both policies, including the exact tie.
- **S12 Two-part search, outside the grid** (F-013). A datum whose nearest
  node index lies outside the grid is an error, naming its row. Under
  `OutsideGrid::ClampToEdge`, it is assigned to the edge node, as
  `getindx` does. Tests: both policies, on each face.

## Open questions

None.

## Implementation plan

Each step is one reviewable PR; the checks in `AGENTS.md` must pass.

1. **Data search types and brute force** (`gslib_rs::search`; done: 13
   tests. A mutation check caught 10 of 11 planted changes; the 11th,
   dropping GSLIB's `8n` stop, is equivalent, see S3).
   - Work:
     - `SearchEllipsoid` and `DataSearch`, plus `Error::TooFewData` and
       `OnTooFewData`;
     - the `Neighbourhood` trait (returning row indices in S2 order) with
       `BruteForce`;
     - replace today's `Search` (isotropic radius, max count) and migrate
       its callers: `kriging`, `sgsim`, `sisim` and the CLI, keeping CLI
       behaviour;
     - tick the search-parameter ledger item (names verified in F-012).
   - Tests: S1–S5, S7.
2. **Super-block accelerator** (done: `SuperBlock`, about four data per
   block over the data's extent; each block is skipped only when
   `smallest_scale² · d²` to the tight box of its own data exceeds
   `radius²`, with a 1e-9 relative margin for rounding. Equal to brute
   force over 2 000 random cases. Timing, release build, 20 000 scattered
   2-D data, 2 000 targets, radius 5% of the domain, `max_data = 24`:
   850 ms brute force against 17.6 ms, plus 1.5 ms to build, about 48×.
   All 6 planted bugs in the bound and block range were caught. The CLI's
   `krige` uses it).
   - Work:
     - a `SuperBlock` built from the data's extent, with a block-selection
       bound proven to be a lower bound (D-025);
     - a timing comparison against brute force, recorded in the PR as the
       evidence for keeping it.
   - Tests: S6.
   - Depends on 1.
3. **Search fixtures through `kt3d`** (with harness step 7; done: nine
   fixtures, F-015. GSLIB's selection equals a port of its super-block
   search at every node; the `picksup` miss is confirmed (150 of 1 764
   nodes in `picksup_probe`); k = 100 (D-035)).
   - Work:
     - a `kt3d` case generator in `oracle-case`, run with `idbg = 2`;
     - a debug-file parser in `harness/` (D-027).
   - Cases:
     - isotropic and anisotropic ellipsoids;
     - octant search;
     - an `ndmax` cut, and a tie straddling it (expected deviation, checked
       as a valid choice);
     - `ndmin` not met (GSLIB writes UNEST; compared under the opt-in
       policy);
     - the F-012 probe geometry, which confirms or withdraws the `picksup`
       hypothesis.
   - Depends on 1 and on point kriging matching `kt3d` (harness step 7);
     benefits from 2.
4. **Node search.**
   - Work: the offset table (S8, S10) and node selection (S9) for a
     `Grid`; `sgsim` moves from a node list to a `Grid`.
   - Tests: S8–S10.
   - Depends on 1.
5. **Two-part search.**
   - Work: assigning data to nodes with `OnSharedNode` and `OutsideGrid`.
   - Tests: S11, S12.
   - Depends on 4.
6. **Records.**
   - ROADMAP §2 rows: tie order (D-024); `sgsim`'s node-search extent and
     the two-part-search policies (F-013, D-028); and the `picksup` row if
     step 3 confirms it.
   - ROADMAP §1 search row updated.
   - Done as each step lands.

Blocked on `parity-verify`: only the `picksup` hypothesis (F-012), which step
3 resolves. Exact-stream tie order (D-024) waits for the RNG spec (PLAN C5).
