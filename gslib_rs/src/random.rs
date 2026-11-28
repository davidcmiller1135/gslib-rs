//! Seeds and random streams (`docs/spec/rng.md`).
//!
//! The default generator is `ChaCha8` (D-036). One master [`Seed`] gives every
//! realization its own stream: realization `k` draws from `ChaCha8` keyed by
//! the seed, on stream (nonce) `k` (D-037). Realization `k` therefore depends
//! only on `(seed, k)`: it is the same run alone or in a batch, sequentially
//! or in parallel.
//!
//! **Divergence from GSLIB** (ROADMAP §2): GSLIB runs one ACORN stream
//! through all realizations (F-016), so its realization `k` depends on every
//! earlier one. Reproducing that is the opt-in exact-stream mode (D-038).

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

/// The default generator: one realization's stream (D-036).
pub type Stream = ChaCha8Rng;

/// A master seed for a run of realizations (D-039).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Seed(pub u64);

impl Seed {
    /// The stream for realization `k` (0-based): `ChaCha8` keyed by
    /// `seed_from_u64(seed)`, on stream `k` of its 2⁶⁴ (RN1).
    #[must_use]
    pub fn realization(self, k: u64) -> Stream {
        let mut rng = ChaCha8Rng::seed_from_u64(self.0);
        rng.set_stream(k);
        rng
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn first(mut rng: ChaCha8Rng) -> [u64; 3] {
        [rng.next_u64(), rng.next_u64(), rng.next_u64()]
    }

    #[test]
    fn realization_streams_are_chacha_streams() {
        // RN1: the definition.
        let mut expected = ChaCha8Rng::seed_from_u64(69069);
        expected.set_stream(3);
        assert_eq!(first(Seed(69069).realization(3)), first(expected));
        // Stream 0 is the plain seeded generator.
        assert_eq!(
            first(Seed(69069).realization(0)),
            first(ChaCha8Rng::seed_from_u64(69069))
        );
    }

    #[test]
    fn streams_differ_by_realization_and_seed() {
        let a = first(Seed(69069).realization(0));
        assert_ne!(a, first(Seed(69069).realization(1)));
        assert_ne!(a, first(Seed(69070).realization(0)));
    }

    #[test]
    fn streams_are_pinned() {
        // RN1 regression pin: values recorded from this implementation with
        // rand_chacha 0.10.0 (2026-11-28), not an independent oracle. A
        // dependency update that changes them would silently change every
        // realization; it must fail here instead (D-036).
        assert_eq!(
            first(Seed(69069).realization(0)),
            [
                8_309_600_352_294_684_511,
                16_278_930_870_170_229_085,
                16_415_774_044_364_257_802
            ]
        );
        assert_eq!(
            first(Seed(69069).realization(1)),
            [
                14_007_242_514_127_940_649,
                5_389_970_904_244_152_485,
                3_559_890_008_414_256_390
            ]
        );
        assert_eq!(
            first(Seed(1).realization(0)),
            [
                7_424_550_030_962_593_201,
                1_482_817_706_323_250_795,
                11_004_592_982_271_133_285
            ]
        );
    }
}
