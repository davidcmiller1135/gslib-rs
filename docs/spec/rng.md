# Spec: Random numbers, seeds and streams

Status: **in progress** (grill-me session started 2026-11-28; evidence
F-016).

## Purpose

Decide which generators gslib-rs uses, how seeds become streams
(including one per realization, for parallel runs), and what the opt-in
"exact stream parity" with GSLIB requires, so that stochastic programs are
reproducible from their seed and independent of thread count.

## Fixed by AGENTS.md

- **Reproducible from the seed and inputs**, independent of thread count
  (precedence 1).
- **Explicit streams:** RNGs and seeds are passed explicitly; there is no
  global RNG. In parallel code, one seed per realization; streams are never
  shared between threads.
- **Parity levels:**
  - statistical parity is the default for stochastic programs;
  - exact stream parity is an opt-in per program, needing ACORN, the same
    random path, search order and arithmetic, against the double-promoted
    build.
- **Seeds are a domain type** (*Domain types vs. hot loops*).

## Current state

- **Library:** the algorithms take `&mut impl rand::Rng`.
- **CLI:** seeds a `ChaCha8Rng` with `seed_from_u64(--seed)`, default 69069
  (confirmed as GSLIB's example seed, F-016). One stream runs through all
  realizations.
- **`sgsim`'s path** comes from `SliceRandom::shuffle`, and its deviates
  from `rand_distr::StandardNormal`.

## Decisions

- **Generator** (D-036): ChaCha8 by default; ACORN as an opt-in
  `rand::RngCore` for exact-stream parity.
- **Streams** (D-037): realization k uses `ChaCha8Rng::seed_from_u64(master)`
  with `set_stream(k)`, so it depends only on (master, k).
- **Exact-stream parity** (D-038): now, ACORN with GSLIB's seeding and
  warm-up, GSLIB's random path, and `gauinv` deviates, each checked bit for
  bit. Later, after harness step 8, an `sgsim` exact mode: same path and
  draws, values within a solver tolerance.
- **Seeds** (D-039): a `Seed(u64)` type owning the stream rule; ACORN
  requires `0 < seed < 2³⁰`; even seeds are accepted.

## Open questions

None.

## Requirements

- **RN1 Realization streams** (D-037, D-039). `Seed(s).realization(k)`
  equals `ChaCha8Rng::seed_from_u64(s)` with `set_stream(k)`. A regression
  test pins its first outputs for two seeds and streams. The values are
  recorded from the implementation, as a guard against dependency changes,
  not as an oracle.
- **RN2 Order independence** (D-037, P5). For `sgsim` and `sisim`,
  realization k run alone equals realization k of a batch of n > k, bit for
  bit.
- **RN3 ACORN** (D-036, D-039, F-016).
  - `Acorn::new(seed)` holds `[seed, 0, …, 0]` (13 words) and fails unless
    `0 < seed < 2³⁰`. `Acorn::gslib(seed)` also discards 1 000 draws.
  - `next_acorni()` performs `ixv(i+1) = (ixv(i+1) + ixv(i)) mod 2³⁰` for
    i = 1..12 and returns `ixv(13) / 2³⁰`.
  - **Oracle:** GSLIB `nscore` seeds `Acorn::gslib(69069)` and perturbs
    each *kept* datum, in data order, to `v + acorni·10⁻⁶` (`NSCORE.FOR`
    lines 142–179; trimmed rows draw nothing). Its table is those values
    sorted. The harness checks them exactly, on every `nscore` fixture,
    including `ties`, whose random tie order this reproduces.
- **RN4 ACORN as a `rand` generator** (rand 0.10's `TryRng` with
  `Error = Infallible`, which provides `Rng`; `RngCore` is deprecated in
  0.10). `next_u32` is the 30-bit word shifted to
  the top of 32 bits; `next_u64` combines two draws. The exact mode uses
  `next_acorni` only; the trait exists so ACORN can drive any gslib-rs
  algorithm.
- **RN5 GSLIB's random path.** One `next_acorni` per node in grid order,
  then nodes visited in increasing order of their draws (ties, which need
  equal 30-bit words, by node index), with `imult` added at the
  multiple-grid nodes when that option is used.
- **RN6 `gauinv`.** GSLIB's approximation moves into the library as a
  documented function, within 1.501e-8 of the exact inverse (D-016). The
  harness uses the library's version.
- **RN7 The oracle for RN3, RN5 and RN6 together.** An unconditional `sgsim`
  fixture with a search radius below the node spacing and 2 realizations:
  - GSLIB's values equal `gauinv` of the draws along the predicted path,
    exactly;
  - the second realization continues the same stream.
- **RN8 CLI seeds.** `--seed` becomes a `Seed`. In the default mode any
  `u64` is accepted; the exact mode (when built) validates ACORN's range.

## Implementation plan

Each step is one reviewable PR; the checks in `AGENTS.md` must pass.

1. **Seeds and realization streams** (done: `gslib_rs::random::{Seed,
   Stream}`; the CLI gives realization i stream i of `--seed`, and no longer
   depends on `rand` or `rand_chacha` directly. Tests: RN1 (the definition,
   plus pinned outputs for two seeds and streams), RN2 (a CLI batch's
   realization 2 equals the library's realization 2 alone). Planted bugs:
   a shared stream and an ignored stream number were both caught).
   - Work: `gslib_rs::random` with `Seed` and `Seed::realization(k)`; the
     CLI seeds each realization from it.
   - Tests: RN1, RN2, RN8.
   - Unblocks parallel realizations (parallelism step 4).
2. **ACORN** (done: `random::Acorn`. Unit tests check the recurrence against
   its closed form `s·C(n + 11, 12) mod 2³⁰` for five seeds over 60 draws,
   including a sum of exactly 2³⁰. **GSLIB check:** the despiked values of
   all 7 `nscore` tables, `ties` and trimmed `with_sentinels` included, are
   reproduced bit for bit. Planted bugs: a 999-draw warm-up (612
   mismatches against GSLIB), order 11, and a `>` in the reduction, all
   caught once two unit tests were sharpened).
   - Work: `random::Acorn` (`new`, `gslib`, `next_acorni`, `TryRng`).
   - Tests: RN3 (unit tests from the recurrence, and the `nscore`
     despiking oracle in the harness), RN4.
3. **`gauinv` in the library.**
   - Work: move the harness port into `transform`, documented as GSLIB's
     approximation.
   - Tests: RN6 (the existing bound tests move with it).
4. **GSLIB's random path and the unconditional `sgsim` fixture.**
   - Work: `random::gslib_path`, and a minimal `sgsim` case generator in
     `oracle-case`, the start of harness step 8.
   - Tests: RN5, RN7.
5. **Later: the `sgsim` exact mode** (D-038), after harness step 8.

Nothing is blocked on `parity-verify`: F-016 is read from the source, and
RN3 and RN7 check it against GSLIB runs.
