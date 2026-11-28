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

## Decisions so far

- **Generator** (D-036): ChaCha8 by default; ACORN as an opt-in
  `rand::RngCore` for exact-stream parity.
- **Streams** (D-037): realization k uses `ChaCha8Rng::seed_from_u64(master)`
  with `set_stream(k)`, so it depends only on (master, k).

## Open questions (in dependency order)

1. ~~The default generator~~ (D-036).
2. ~~Realization streams~~ (D-037).
3. What exact stream parity covers, and when it is built.
4. Seed type and validation.
