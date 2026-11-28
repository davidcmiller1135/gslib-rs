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
//! earlier one. Reproducing that is the opt-in exact-stream mode (D-038),
//! built from [`Acorn`], GSLIB's generator.

use std::convert::Infallible;

use rand::{SeedableRng, TryRng};
use rand_chacha::ChaCha8Rng;

use crate::error::{Error, Result};
use crate::simulation::Grid;

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

/// ACORN's modulus, `MAXINT = 2³⁰` (`ACORNI.FOR`).
pub const ACORN_MODULUS: u32 = 1 << 30;

/// Draws GSLIB discards after seeding (`do i=1,1000`, F-016).
pub const GSLIB_WARM_UP: usize = 1000;

/// GSLIB's random number generator: ACORN of order 12, modulus 2³⁰
/// (`vendor/gslib-2.0/src/GSLIB/ACORNI.FOR`; F-016, RN3, RN4).
///
/// Each draw updates `ixv(i+1) = (ixv(i+1) + ixv(i)) mod 2³⁰` for
/// i = 1..12, in order and with already-updated values, then yields
/// `ixv(13)`. [`Acorn::next_acorni`] returns it as GSLIB's `acorni`,
/// `ixv(13) / 2³⁰` in [0, 1).
///
/// As a `rand` generator, `next_u32` is the 30-bit word shifted to the top
/// of 32 bits and `next_u64` joins two draws; only `next_acorni`
/// reproduces GSLIB's numbers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Acorn {
    ixv: [u32; 13],
}

impl Acorn {
    /// ACORN seeded as its notes require: `ixv(1) = seed`, the other 12
    /// words 0, with no draws discarded.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] unless `0 < seed < 2³⁰` (D-039). Even seeds
    /// are accepted, as GSLIB accepts them.
    pub fn new(seed: u64) -> Result<Self> {
        let word = u32::try_from(seed)
            .ok()
            .filter(|&w| w > 0 && w < ACORN_MODULUS)
            .ok_or_else(|| {
                Error::invalid(
                    "ACORN seed",
                    format!("must satisfy 0 < seed < 2^30 = {ACORN_MODULUS}, got {seed}"),
                )
            })?;
        let mut ixv = [0; 13];
        if let Some(first) = ixv.first_mut() {
            *first = word;
        }
        Ok(Self { ixv })
    }

    /// ACORN seeded as GSLIB's programs seed it: [`Acorn::new`], then
    /// [`GSLIB_WARM_UP`] draws discarded.
    ///
    /// # Errors
    /// As [`Acorn::new`].
    pub fn gslib(seed: u64) -> Result<Self> {
        let mut rng = Self::new(seed)?;
        for _ in 0..GSLIB_WARM_UP {
            rng.next_word();
        }
        Ok(rng)
    }

    /// The next 30-bit word, `ixv(13)`.
    fn next_word(&mut self) -> u32 {
        let (first, rest) = self.ixv.split_at_mut(1);
        let mut previous = first.first().copied().unwrap_or(0);
        for word in rest {
            // Both are below 2³⁰, so the sum fits in u32.
            *word += previous;
            if *word >= ACORN_MODULUS {
                *word -= ACORN_MODULUS;
            }
            previous = *word;
        }
        previous
    }

    /// The next number exactly as GSLIB's `acorni`: `ixv(13) / 2³⁰`.
    pub fn next_acorni(&mut self) -> f64 {
        f64::from(self.next_word()) / f64::from(ACORN_MODULUS)
    }
}

impl TryRng for Acorn {
    type Error = Infallible;

    fn try_next_u32(&mut self) -> std::result::Result<u32, Infallible> {
        Ok(self.next_word() << 2)
    }

    fn try_next_u64(&mut self) -> std::result::Result<u64, Infallible> {
        let high = u64::from(self.next_word() << 2);
        let low = u64::from(self.next_word() << 2);
        Ok((high << 32) | low)
    }

    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> std::result::Result<(), Infallible> {
        for chunk in dst.chunks_mut(4) {
            let bytes = (self.next_word() << 2).to_le_bytes();
            for (d, b) in chunk.iter_mut().zip(bytes) {
                *d = b;
            }
        }
        Ok(())
    }
}

/// GSLIB `sgsim`'s random path (`SGSIM.FOR` lines 814–842; RN5): the order
/// in which a realization visits the nodes of `grid`, as 0-based indices.
///
/// One [`Acorn::next_acorni`] is drawn per node in grid order (x fastest).
/// With `multiple_grid = Some(n)` (GSLIB `mults = 1, nmult = n`), `imult` is
/// added to the draw of every node on the coarse grid of level `imult`
/// (every `4·imult`-th node along each axis that has more than one), for
/// `imult = 1..=n`, exactly as GSLIB does. Nodes are then visited in
/// increasing order of their keys; equal keys, which need equal 30-bit
/// draws, go by node index (GSLIB's `sortem` leaves their order
/// unspecified).
#[must_use]
pub fn gslib_path(acorn: &mut Acorn, grid: &Grid, multiple_grid: Option<usize>) -> Vec<usize> {
    let mut keys: Vec<f64> = (0..grid.len()).map(|_| acorn.next_acorni()).collect();
    if let Some(levels) = multiple_grid {
        let [nx, ny, nz] = grid.counts();
        for imult in 1..=levels {
            let step = imult.saturating_mul(4);
            let coarse = |n: usize| n.div_euclid(step).max(1);
            let (nnx, nny, nnz) = (coarse(nx), coarse(ny), coarse(nz));
            // GSLIB's 1-based jx, jy, jz start at 1 and change only along
            // axes with more than one coarse node.
            let at = |n: usize, i: usize| if n > 1 { i * step } else { 1 };
            for iz in 1..=nnz {
                for iy in 1..=nny {
                    for ix in 1..=nnx {
                        let cell = [at(nnx, ix) - 1, at(nny, iy) - 1, at(nnz, iz) - 1];
                        if let Some(key) = grid.index_of(cell).and_then(|i| keys.get_mut(i)) {
                            *key += crate::num::to_f64(imult);
                        }
                    }
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..grid.len()).collect();
    order.sort_by(|&a, &b| {
        let key = |i: usize| keys.get(i).copied().unwrap_or(f64::INFINITY);
        key(a).total_cmp(&key(b)).then(a.cmp(&b))
    });
    order
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

    /// `C(n, k) mod 2³⁰`, exactly, via u128.
    fn binomial_mod(n: u64, k: u64) -> u64 {
        let mut c: u128 = 1;
        for i in 0..k {
            // Exact: c · (n − i) is divisible by (i + 1) at every step.
            c = (c * u128::from(n - i)).div_euclid(u128::from(i + 1));
        }
        u64::try_from(c % (1 << 30)).unwrap()
    }

    #[test]
    fn acorn_follows_its_closed_form() {
        // RN3: from [s, 0, ..., 0], each draw is a running sum, so after n
        // draws ixv(13) = s * C(n + 11, 12) mod 2^30, an independent check
        // of the recurrence and of the modular reduction.
        // 2^29 makes a sum hit 2^30 exactly on the second draw, which must
        // reduce to 0.
        for seed in [1_u64, 7, 69069, 1 << 29, (1 << 30) - 1] {
            let mut rng = Acorn::new(seed).unwrap();
            for n in 1..=60_u64 {
                let expected = (seed * binomial_mod(n + 11, 12)) % (1 << 30);
                let got = rng.next_word();
                assert_eq!(u64::from(got), expected, "seed {seed}, draw {n}");
            }
        }
        // The first draws from seed 1: 1, 13, 91 (over 2^30).
        let mut rng = Acorn::new(1).unwrap();
        let scale = f64::from(ACORN_MODULUS);
        for word in [1.0, 13.0, 91.0] {
            assert_eq!(rng.next_acorni().to_bits(), (word / scale).to_bits());
        }
    }

    #[test]
    fn acorn_seeds_are_validated_and_warmed_up() {
        assert!(Acorn::new(0).is_err(), "zero");
        assert!(Acorn::new(1 << 30).is_err(), "2^30");
        assert!(Acorn::new(u64::MAX).is_err(), "huge");
        assert!(Acorn::new(2).is_ok(), "even seeds are accepted");
        // gslib() is new() plus 1 000 discarded draws (written out, not
        // through the constant, so a wrong constant fails here).
        let mut warmed = Acorn::new(69069).unwrap();
        for _ in 0..1000 {
            warmed.next_acorni();
        }
        assert_eq!(Acorn::gslib(69069).unwrap(), warmed);
        // Values stay in [0, 1).
        let mut rng = Acorn::gslib(69069).unwrap();
        assert!((0..10_000).all(|_| (0.0..1.0).contains(&rng.next_acorni())));
    }

    #[test]
    fn acorn_drives_rand() {
        // RN4: next_u32 is the word shifted up by 2; next_u64 joins two.
        let mut a = Acorn::new(5).unwrap();
        let mut b = Acorn::new(5).unwrap();
        assert_eq!(a.next_u32(), b.next_word() << 2);
        let expected = (u64::from(b.next_word() << 2) << 32) | u64::from(b.next_word() << 2);
        assert_eq!(a.next_u64(), expected);
        let mut bytes = [0_u8; 6];
        a.fill_bytes(&mut bytes);
        let w1 = (b.next_word() << 2).to_le_bytes();
        let w2 = (b.next_word() << 2).to_le_bytes();
        assert_eq!(bytes, [w1[0], w1[1], w1[2], w1[3], w2[0], w2[1]]);
    }

    fn grid(nx: usize, ny: usize) -> Grid {
        use crate::simulation::Axis;
        let axis = |count| Axis {
            count,
            first_center: 0.5,
            spacing: 1.0,
        };
        Grid::new(axis(nx), axis(ny), axis(1)).unwrap()
    }

    #[test]
    fn the_path_sorts_nodes_by_their_draws() {
        // RN5: one draw per node in grid order; visit by increasing draw.
        let g = grid(5, 4);
        let path = gslib_path(&mut Acorn::gslib(69069).unwrap(), &g, None);
        let mut acorn = Acorn::gslib(69069).unwrap();
        let draws: Vec<f64> = (0..20).map(|_| acorn.next_acorni()).collect();
        assert!(
            path.windows(2).all(|w| draws[w[0]] <= draws[w[1]]),
            "{path:?}"
        );
        let mut sorted = path;
        sorted.sort_unstable();
        assert_eq!(sorted, (0..20).collect::<Vec<_>>(), "a permutation");
    }

    #[test]
    fn the_multiple_grid_delays_coarse_nodes() {
        // 12 x 10: level 1 adds 1 at x in {4, 8, 12}, y in {4, 8} (1-based);
        // level 2 has one coarse node per axis (12/8 = 10/8 = 1), so it adds
        // 2 at (1, 1). All other keys stay in [0, 1).
        let g = grid(12, 10);
        let path = gslib_path(&mut Acorn::gslib(5).unwrap(), &g, Some(2));
        let level_one: Vec<usize> = [4, 8, 12]
            .iter()
            .flat_map(|&x| [4, 8].map(|y| (x - 1) + 12 * (y - 1)))
            .collect();
        let n = path.len();
        assert_eq!(path[n - 1], 0, "node (1, 1) has key in [2, 3): last");
        let before_last: Vec<usize> = path[n - 7..n - 1].to_vec();
        let mut expected = level_one;
        expected.sort_unstable();
        let mut got = before_last;
        got.sort_unstable();
        assert_eq!(got, expected, "the six level-1 nodes come just before it");
        // Without the multiple grid, the same draws give another order.
        let plain = gslib_path(&mut Acorn::gslib(5).unwrap(), &g, None);
        assert_ne!(plain, path);
    }
}
