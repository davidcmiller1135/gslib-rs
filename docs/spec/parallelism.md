# Spec: Parallel execution and acceleration

Status: **draft**. The questions below carry recommendations but are not
decided; a grill-me session confirms them before implementation.

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

## Decisions so far

- **Thread pool** (D-029): `rayon`.
- **Location** (D-030): `gslib_rs::parallel` behind an opt-in `parallel`
  feature, mapping over the core's pure per-item functions; the CLI enables
  it.

## Open questions (with recommendations)

1. ~~Thread pool~~ (D-029: `rayon`). Was: `rayon` *(recommended: work stealing, `par_iter` on
   slices, mature, permissive licence, widely audited)* or only
   `std::thread::scope` (no dependency, but manual chunking everywhere).
2. ~~Where the parallel loops live~~ (D-030). Was: **Where the parallel loops live.**
   - *Recommended:* a boundary module `gslib_rs::parallel` behind an opt-in
     `parallel` feature. It holds parallel maps over the pure per-item
     functions: kriging at many targets, many realizations, the pair loop's
     blocks. The core stays sequential and dependency-free.
   - Alternative: parallel loops only in `cli/`, so library users wire their
     own.
3. **Deterministic reduction for experimental variograms.**
   - *Recommended:* split the pairs into fixed blocks of rows whose size
     does not depend on the thread count. Each block accumulates its own
     partial sums sequentially, and the partial sums are combined in block
     order. The result is bit-identical for any thread count.
   - Alternatives:
     - exact accumulation (e.g. integer or expansion arithmetic): order-free,
       but slower and more code;
     - keep the variogram sequential.
   - Either way, the parallel result may differ from today's sequential one
     in the last bits (a different but fixed order). It must still pass the
     `gamv` parity fixtures at their D-010 tolerances.
4. **Pruning pairs by distance.** Pairs beyond GSLIB's cut-off
   `(nlag + 0.5)·xlag` contribute nothing.
   - *Recommended:* enumerate pairs through the super-block, so far pairs are
     never examined, with a property test that the counted pairs are exactly
     the same.
   - This is independent of threads and is probably the larger gain for big
     data sets, since the work falls from O(N²) toward O(N·k).
5. **Parallel realizations.** Deferred to the RNG spec (PLAN C5). It must
   say how realization `k`'s seed is derived from the master seed, so that
   realization `k` depends only on `(master seed, k)`.
6. **GPU.**
   - *Recommended:* not now. Record it as an opt-in extension (ROADMAP §5),
     evaluated only with benchmark evidence. The obstacles:
     - **Double precision:** consumer GPUs run `f64` at 1/32 to 1/64 of
       their `f32` speed, and WebGPU (`wgpu`) has essentially no `f64`.
     - **Safety:** CUDA from Rust relies on `unsafe` bindings. That is
       allowed only inside a dependency, and it is heavy and vendor-specific.
     - **Determinism:** reductions with atomic additions have no fixed order.
     - **Testing:** CI has no GPU runners.
   - The one GPU-friendly workload is the regular O(N²) pair loop of the
     experimental variogram. Searches are branchy and memory-bound, and the
     kriging systems are tiny (about 24×24).

## Requirements (draft, to be confirmed with the questions)

- **P1 Thread-count invariance.** For every parallel entry point, the output
  with 1, 2 and 8 threads is bit-identical. Tests: each entry point run under
  explicitly sized pools.
- **P2 Kriging at many targets.** The parallel map over targets returns
  bit-identically what the sequential loop returns, including errors: the
  first error in *target order*, not the first to occur. Tests: P1, plus a
  case where two targets fail.
- **P3 Deterministic variogram reduction.** The blocked pair loop is
  bit-identical for any thread count (P1). It agrees with the sequential
  loop within D-010's direct-formula tolerance, and all `gamv` fixtures still
  pass. Tests: P1; a property test against the sequential loop; the harness.
- **P4 Exact pruning.** The pruned pair loop counts exactly the same pairs in
  every lag as the unpruned loop. Tests: a property test of pair counts on
  random data and lag schemes, including the zero-distance class and
  overlapping lags.
- **P5 Independent realizations.** Realization `k` depends only on the
  master seed and `k`. Running realizations 0..n in parallel gives the same
  files as running them one at a time, or running only realization `k`.
  Tests: P1; realization `k` alone equals realization `k` of a batch.
  Depends on the RNG spec.
- **P6 GPU, if ever.** Off by default behind a feature. Each GPU kernel is
  checked against the CPU result within a documented tolerance, has a CI job
  on GPU hardware, and comes with benchmark evidence that it beats the
  parallel CPU path.

## Implementation plan (draft)

1. **Grill session** on questions 1–4 and 6 (5 waits for the RNG spec).
2. **Parallel kriging over targets.**
   - Work: `rayon` behind the `parallel` feature; `parallel::krige_all`, or
     whatever question 2 decides; the `krige` command uses it.
   - Tests: P1, P2.
3. **Exact pair pruning** through the super-block.
   - Tests: P4; `gamv` fixtures unchanged.
4. **Deterministic parallel variogram.**
   - Work: the blocked pair loop.
   - Tests: P1, P3.
   - Depends on 3, which supplies the block structure.
5. **Parallel realizations** for `sgsim` and `sisim`.
   - Tests: P1, P5.
   - Depends on the RNG spec (PLAN C5).
6. **GPU evaluation** (optional).
   - Work: a benchmark of the pair loop only.
   - Tests: P6.
   - Taken further only if it beats step 4 on representative data.
