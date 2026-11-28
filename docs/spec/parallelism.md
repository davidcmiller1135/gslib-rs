# Spec: Parallel execution and acceleration

Status: **specified** (decisions D-029 to D-033), except realization seeds,
which wait for the RNG spec.

## Purpose

Decide how gslib-rs uses several CPU cores, and possibly GPUs, for
searches, kriging, experimental variograms and simulation, without
letting results depend on the number of threads or the hardware.

## Fixed by AGENTS.md

- **Precedence 1:** results are reproducible from the seed and inputs and
  **do not depend on the number of threads**. Performance comes last
  (precedence 5).
- **Parallel execution is a boundary concern** (*Isolate boundaries*): the
  algorithms stay pure and sequential per item; parallel maps live at the
  edges.
- **Streams:** in parallel code, derive one seed per realization and never
  share a stream between threads.
- **`unsafe` is forbidden** in our crates, and a new crate must be justified
  and pass `cargo deny check`.

## Current state (2026-11-27)

- **Searches are ready.** `DataSearch::select` takes `&self` and holds no
  mutable state; `BruteForce` and `SuperBlock` are `Send + Sync`; a test
  (search spec S7) shows that searching from several threads gives exactly
  the sequential result. The `Neighbourhood` trait does not require `Sync`,
  so generic parallel callers add the bound.
- **Kriging per target is independent.** Each target's search and solve use
  only shared, read-only inputs.
- **Experimental variograms are not ready.** `variogram::semivariogram` is a
  single O(N²) loop that adds floating-point sums in pair order. Splitting it
  naively changes that order, and so the last bits of γ, with the thread
  count. It also examines every pair, even beyond the largest lag.
- **Simulation:** a realization is sequential by construction, each node
  depending on the earlier ones. Realizations are independent apart from
  their seeds, whose derivation belongs to the RNG spec (PLAN C5).
- **No thread-pool crate** is a dependency yet.

## Decisions

- **Thread pool** (D-029): `rayon`.
- **Location** (D-030): `gslib_rs::parallel` behind an opt-in `parallel`
  feature, off by default. It maps over the core's pure per-item functions,
  and the CLI enables it.
- **Variogram reduction** (D-031): pairs are grouped into blocks fixed by
  the algorithm, not the thread count. Each block sums sequentially, and
  the partial sums are combined in block order. The sequential path uses the
  same blocks, so both agree bit for bit.
- **Pruning** (D-032, D-034): always on. Each row's partners come from the
  super-block grid within the lag scheme's largest distance, widened by a
  small margin, in block order; the lag scheme still decides each pair.
- **GPU** (D-033): not now. It is an opt-in extension, pursued only with
  benchmark evidence after the parallel CPU paths exist.

## Open questions

- **Realization seeds:** how realization `k`'s seed derives from the master
  seed. Deferred to the RNG spec (PLAN C5); it blocks plan step 5 only.

## Requirements

- **P1 Thread-count invariance.** For every parallel entry point, the output
  with 1, 2 and 8 threads is bit-identical. Tests: each entry point run under
  explicitly sized pools.
- **P2 Kriging at many targets.** The parallel map over targets returns
  bit-identically what the sequential loop returns, including errors: the
  first error in *target order*, not the first to occur. Tests: P1, plus a
  case where two targets fail.
- **P3 Deterministic variogram reduction** (D-031). The blocked pair loop
  is bit-identical for any thread count (P1), and the sequential and
  parallel paths are bit-identical to each other. Both agree with the old
  single loop within the summation bound (`2·n·ε` relative for a lag with
  `n` pairs, at least 1e-12), and all `gamv` fixtures still pass at their
  D-010 tolerances. Tests: P1; sequential equals parallel bit for bit; a
  property test against the old single loop; the harness.
- **P4 Exact pruning** (D-032). The pruned pair loop counts exactly the same
  pairs in every lag as the unpruned loop, including pairs exactly at the
  cut-off. Tests: a property test of pair counts on
  random data and lag schemes, including the zero-distance class and
  overlapping lags.
- **P5 Independent realizations.** Realization `k` depends only on the
  master seed and `k`. Running realizations 0..n in parallel gives the same
  files as running them one at a time, or running only realization `k`.
  Tests: P1; realization `k` alone equals realization `k` of a batch.
  Depends on the RNG spec.
- **P6 GPU, if ever** (D-033). Off by default behind a feature. Each GPU kernel is
  checked against the CPU result within a documented tolerance, has a CI job
  on GPU hardware, and comes with benchmark evidence that it beats the
  parallel CPU path.

## Implementation plan

Each step is one reviewable PR; the checks in `AGENTS.md` must pass, and
CI gains a `--features parallel` build where it is not already covered by
feature unification (D-030).

1. **Parallel kriging over targets** (done: `rayon` 1.12 passes `cargo deny`;
   `parallel::try_map` and `parallel::with_threads`; the `krige` command
   kriges its targets in parallel. Tests: bit-identical estimates with 1, 2
   and 8 threads; the first error in item order with a deliberately slow
   failing item; the first failing target named. Collecting straight into
   a `Result` was caught on 3 of 3 runs).
   - Work:
     - `rayon` (D-029, passing `cargo deny check`);
     - the `parallel` feature and `gslib_rs::parallel` (D-030), with a
       parallel map that krigs many targets and reports the first error in
       target order;
     - the `krige` command uses it.
   - Tests: P1 with pools of 1, 2 and 8 threads; P2.
2. **Exact pair pruning** (D-032; done).
   - Result:
     - `LagScheme::max_distance`; the search's block grid is shared as an
       internal `Blocks` over bare locations, so `sasim`'s grid nodes use it
       too;
     - `PairPlan` holds rows in fixed blocks of 256 with per-block sums,
       ready for step 3;
     - the old double loop is kept as a test reference.
   - Evidence:
     - pair counts are exact against the double loop (property test, also
       at 1 000 cases; the cut-off pair; several blocks on a lattice);
     - `gamv` fixtures unchanged (worst γ 6.3e-16);
     - the library tests caught all 5 planted bugs.
   - Timing: first committed with partners in distance order, then
     switched to block order after a benchmark (D-034). On 20 000 scattered
     data, release build: 0.98 s for all pairs against 0.07 s with lags to
     5% of the domain; 1.60 s against 0.99 s with lags to 20%.
   - Work: `LagScheme` gains its largest separation; the pair loop finds
     partners through the super-block, in D-031's fixed blocks of rows.
   - Tests: P4 as a property test over random data and both lag schemes
     (including the zero-distance class, overlapping lags and pairs exactly
     at the cut-off); the `gamv` fixtures still pass; a timing note.
3. **Deterministic parallel variogram** (D-031).
   - Work: blocks reduced in parallel, combined in block order.
   - Tests: P1; P3 (sequential equals parallel bit for bit, and agreement
     with the old single loop within 1e-12).
   - Depends on 2.
4. **Parallel realizations** (`sgsim`, `sisim`).
   - Tests: P1, P5.
   - Blocked on the RNG spec (PLAN C5).
5. **GPU evaluation** (optional, D-033).
   - Work: a benchmark of the pair loop only.
   - Tests: P6.
   - Taken further only if it beats step 3 on representative data.
