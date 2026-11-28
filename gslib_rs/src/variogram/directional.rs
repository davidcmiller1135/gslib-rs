//! GSLIB `gamv` direction windows (`docs/spec/gamv.md`).
//!
//! They reproduce `vendor/gslib-2.0/src/GAMV/GAMV.FOR`, routine `gamv`
//! (F-019), with the same comparisons and the same order of floating-point
//! operations, so a pair on a window edge is decided as in the reference
//! build:
//!
//! - sums and products are left unfused (no `mul_add`), as gfortran compiles
//!   them for a target without FMA;
//! - the reference build does not contract them either (`-ffp-contract=off`,
//!   `xtask/src/main.rs`).
//!
//! # Pairs
//!
//! A [`Separation`] is the vector from datum `i` to datum `j > i`, with
//! GSLIB's unfused squared length and distance (`GAMV.FOR` ll. 511–520).
//! Lag classes are `CentredLags` (F-008); the separation is coincident (lag
//! class 0) when its distance is at most `1e-20`.
//!
//! # Direction windows
//!
//! A [`DirectionWindow`] is GSLIB's four separate tests, **not a 3-D cone**
//! (D-045):
//!
//! 1. azimuth: the horizontal projection of the pair is within `atol` of the
//!    azimuth, in either sense; a pair with no horizontal extent passes;
//! 2. horizontal bandwidth;
//! 3. dip: measured in the vertical plane containing the pair, not the one
//!    containing the direction; skipped for lag 0;
//! 4. vertical bandwidth, applied to every direction.
//!
//! A direction with `atol ≥ 90` is omnidirectional. The variogram kernel
//! adds both orderings of its pairs to the sums and reports each pair once
//! (D-043).
//!
//! # Divergences from GSLIB
//!
//! - `atol ≤ 0` or `dtol ≤ 0` is an error here; GSLIB silently uses 45°. The
//!   legacy substitution belongs in the parameter-file adapter (D-044).
//! - π is exact unless [`PiConstant::Gslib`] is selected (D-046).
//!
//! # Known GSLIB edge cases, reproduced (F-020, F-023, F-024)
//!
//! - A pair exactly on a tolerance edge is decided by the last bits of the
//!   arithmetic, and with GSLIB's π the window is not symmetric about its
//!   azimuth.
//! - Tail and head are not antisymmetric: when the pair's azimuth and dip
//!   cosines have opposite signs, it is `Swapped` in both input orders, so
//!   the tail depends on the order of the data (F-024).
//! - `atol = 90` exactly is omnidirectional, yet the azimuth test still
//!   requires `|cos| ≥ cos(90°)`, which is slightly positive in floating
//!   point. Pairs exactly perpendicular to the azimuth can therefore be
//!   rejected (e.g. north–south pairs for `azm = 90`).

use crate::error::{Error, Result};
use crate::geometry::{PiConstant, Point3};

/// `EPSLON` in `GAMV.INC`.
pub(crate) const EPSLON: f64 = 1.0e-20;

/// `PI` in `GAMV.FOR` (l. 449).
#[expect(
    clippy::approx_constant,
    reason = "GSLIB's literal, reproduced on purpose under PiConstant::Gslib"
)]
const GAMV_PI: f64 = 3.141_592_65;

const fn pi_value(pi: PiConstant) -> f64 {
    match pi {
        PiConstant::Exact => std::f64::consts::PI,
        PiConstant::Gslib => GAMV_PI,
    }
}

/// `a·b + c·d`, unfused, as GSLIB's Fortran evaluates it.
#[expect(
    clippy::suboptimal_flops,
    reason = "unfused on purpose: matches the reference build's rounding (module docs)"
)]
fn dot2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * b + c * d
}

/// `a·b − c·d`, unfused, as GSLIB's Fortran evaluates it.
#[expect(
    clippy::suboptimal_flops,
    reason = "unfused on purpose: matches the reference build's rounding (module docs)"
)]
fn cross2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * b - c * d
}

/// The separation from datum `i` to datum `j` (GSLIB's `dx = x(j) − x(i)`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Separation {
    delta: [f64; 3],
    horizontal_squared: f64,
    distance: f64,
}

impl Separation {
    /// The separation from `tail` to `head` (GSLIB's datum `i` to datum `j`),
    /// computed as in `GAMV.FOR` ll. 511–520: unfused squares, summed in
    /// order, then the square root.
    #[must_use]
    pub fn new(tail: Point3, head: Point3) -> Self {
        let [dx, dy, dz] = head.delta(tail);
        let dxs = dx * dx;
        let dys = dy * dy;
        let dzs = dz * dz;
        let horizontal_squared = dxs + dys;
        let distance = (horizontal_squared + dzs).max(0.0).sqrt();
        Self {
            delta: [dx, dy, dz],
            horizontal_squared,
            distance,
        }
    }

    /// Separation distance `h`.
    #[must_use]
    pub const fn distance(&self) -> f64 {
        self.distance
    }

    /// True for coincident data: GSLIB's zero-distance lag, `h ≤ 1e-20`.
    #[must_use]
    pub fn is_coincident(&self) -> bool {
        self.distance <= EPSLON
    }
}

/// One `gamv` direction, as written in the parameter file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction {
    azimuth: f64,
    azimuth_tolerance: f64,
    horizontal_bandwidth: f64,
    dip: f64,
    dip_tolerance: f64,
    vertical_bandwidth: f64,
}

impl Direction {
    /// A direction from GSLIB's `azm, atol, bandh, dip, dtol, bandv`, in
    /// degrees and data units.
    ///
    /// Azimuth is clockwise from north; dip is negative downward.
    /// `atol ≥ 90` makes the direction omnidirectional.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if any value is not finite, a tolerance is not
    /// positive (GSLIB would silently use 45°, D-044), or a bandwidth is
    /// negative.
    pub fn new(
        azimuth: f64,
        azimuth_tolerance: f64,
        horizontal_bandwidth: f64,
        dip: f64,
        dip_tolerance: f64,
        vertical_bandwidth: f64,
    ) -> Result<Self> {
        let all = [
            azimuth,
            azimuth_tolerance,
            horizontal_bandwidth,
            dip,
            dip_tolerance,
            vertical_bandwidth,
        ];
        if !all.iter().all(|v| v.is_finite()) {
            return Err(Error::invalid(
                "direction",
                format!("must be finite, got {all:?}"),
            ));
        }
        for (what, t) in [
            ("azimuth tolerance", azimuth_tolerance),
            ("dip tolerance", dip_tolerance),
        ] {
            if t <= 0.0 {
                return Err(Error::invalid(
                    what,
                    format!("must be positive, got {t} (GSLIB silently uses 45°)"),
                ));
            }
        }
        for (what, b) in [
            ("horizontal bandwidth", horizontal_bandwidth),
            ("vertical bandwidth", vertical_bandwidth),
        ] {
            if b < 0.0 {
                return Err(Error::invalid(
                    what,
                    format!("must not be negative, got {b}"),
                ));
            }
        }
        Ok(Self {
            azimuth,
            azimuth_tolerance,
            horizontal_bandwidth,
            dip,
            dip_tolerance,
            vertical_bandwidth,
        })
    }

    /// True if `atol ≥ 90`.
    #[must_use]
    pub fn is_omnidirectional(&self) -> bool {
        self.azimuth_tolerance >= 90.0
    }

    /// The precomputed window, with unit vectors and cosine thresholds
    /// computed exactly as in `GAMV.FOR` ll. 460–479.
    #[must_use]
    pub fn window(&self, pi: PiConstant) -> DirectionWindow {
        let pi = pi_value(pi);
        let azimuth = (90.0 - self.azimuth) * pi / 180.0;
        let declination = (90.0 - self.dip) * pi / 180.0;
        DirectionWindow {
            ux: azimuth.cos(),
            uy: azimuth.sin(),
            cos_atol: (self.azimuth_tolerance * pi / 180.0).cos(),
            uz: declination.cos(),
            uh: declination.sin(),
            cos_dtol: (self.dip_tolerance * pi / 180.0).cos(),
            bandh: self.horizontal_bandwidth,
            bandv: self.vertical_bandwidth,
            omni: self.is_omnidirectional(),
        }
    }
}

/// Which datum of an accepted pair is the tail (GSLIB ll. 591–613).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairOrder {
    /// The pair's `tail` argument is the tail.
    AsGiven,
    /// Tail and head swap.
    Swapped,
}

/// A [`Direction`] prepared for testing pairs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionWindow {
    ux: f64,
    uy: f64,
    cos_atol: f64,
    uz: f64,
    uh: f64,
    cos_dtol: f64,
    bandh: f64,
    bandv: f64,
    omni: bool,
}

impl DirectionWindow {
    /// True if the direction is omnidirectional (`atol ≥ 90`).
    #[must_use]
    pub const fn is_omnidirectional(&self) -> bool {
        self.omni
    }

    /// Whether this direction accepts the pair, and in which order
    /// (`GAMV.FOR` ll. 549–582).
    #[must_use]
    pub fn accept(&self, pair: &Separation) -> Option<PairOrder> {
        let [dx, dy, dz] = pair.delta;
        let mut dxy = pair.horizontal_squared.max(0.0).sqrt();
        let dcazm = if dxy < EPSLON {
            1.0
        } else {
            dot2(dx, self.ux, dy, self.uy) / dxy
        };
        if dcazm.abs() < self.cos_atol {
            return None;
        }
        if cross2(self.ux, dy, self.uy, dx).abs() > self.bandh {
            return None;
        }
        if dcazm < 0.0 {
            dxy = -dxy;
        }
        let dcdec = if pair.is_coincident() {
            0.0
        } else {
            let c = dot2(dxy, self.uh, dz, self.uz) / pair.distance;
            if c.abs() < self.cos_dtol {
                return None;
            }
            c
        };
        if cross2(self.uh, dz, self.uz, dxy).abs() > self.bandv {
            return None;
        }
        Some(if dcazm >= 0.0 && dcdec >= 0.0 {
            PairOrder::AsGiven
        } else {
            PairOrder::Swapped
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const ORIGIN: Point3 = Point3::new(0.0, 0.0, 0.0);

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    /// Horizontal direction with wide bandwidths.
    fn horizontal(azm: f64, atol: f64) -> Direction {
        Direction::new(azm, atol, 1e9, 0.0, 90.0, 1e9).unwrap()
    }

    fn accepts(d: Direction, pi: PiConstant, head: Point3) -> Option<PairOrder> {
        d.window(pi).accept(&Separation::new(ORIGIN, head))
    }

    // --- W1: construction errors ----------------------------------------

    #[test]
    fn invalid_directions_are_errors() {
        assert!(
            Direction::new(0.0, 0.0, 1.0, 0.0, 10.0, 1.0).is_err(),
            "atol = 0"
        );
        assert!(
            Direction::new(0.0, 10.0, 1.0, 0.0, -1.0, 1.0).is_err(),
            "dtol < 0"
        );
        assert!(
            Direction::new(0.0, 10.0, -1.0, 0.0, 10.0, 1.0).is_err(),
            "bandh < 0"
        );
        assert!(
            Direction::new(0.0, 10.0, 1.0, 0.0, 10.0, -1.0).is_err(),
            "bandv < 0"
        );
        assert!(
            Direction::new(f64::NAN, 10.0, 1.0, 0.0, 10.0, 1.0).is_err(),
            "NaN"
        );
        assert!(
            Direction::new(0.0, 10.0, 0.0, 0.0, 10.0, 0.0).is_ok(),
            "zero bands"
        );
    }

    #[test]
    fn separation_runs_from_tail_to_head() {
        let s = Separation::new(p(1.0, 2.0, 3.0), p(2.0, 2.0, 3.0));
        assert_eq!(
            s.delta.map(f64::to_bits),
            [1.0, 0.0, 0.0].map(f64::to_bits),
            "delta is head − tail"
        );
        assert!((s.distance() - 1.0).abs() < 1e-15, "distance");
        assert!(Separation::new(ORIGIN, ORIGIN).is_coincident(), "h = 0");
        assert!(
            !Separation::new(ORIGIN, p(1e-15, 0.0, 0.0)).is_coincident(),
            "h > 1e-20"
        );
    }

    // --- W2: the four tests --------------------------------------------------

    #[test]
    fn azimuth_test_is_symmetric_in_sense_and_sets_order() {
        let d = horizontal(0.0, 22.5);
        let pi = PiConstant::Exact;
        assert_eq!(
            accepts(d, pi, p(0.0, 1.0, 0.0)),
            Some(PairOrder::AsGiven),
            "north"
        );
        assert_eq!(
            accepts(d, pi, p(0.0, -1.0, 0.0)),
            Some(PairOrder::Swapped),
            "south"
        );
        assert_eq!(accepts(d, pi, p(1.0, 0.0, 0.0)), None, "east");
        // 20° east of north: inside 22.5°.
        let (s20, c20) = 20.0_f64.to_radians().sin_cos();
        assert!(
            accepts(d, pi, p(s20, c20, 0.0)).is_some(),
            "inside tolerance"
        );
        let (s25, c25) = 25.0_f64.to_radians().sin_cos();
        assert!(
            accepts(d, pi, p(s25, c25, 0.0)).is_none(),
            "outside tolerance"
        );
    }

    #[test]
    fn azimuth_is_clockwise_from_north() {
        let d = horizontal(90.0, 10.0);
        assert!(
            accepts(d, PiConstant::Exact, p(1.0, 0.0, 0.0)).is_some(),
            "east"
        );
        assert!(
            accepts(d, PiConstant::Exact, p(0.0, 1.0, 0.0)).is_none(),
            "north"
        );
    }

    #[test]
    fn vertical_pairs_pass_the_azimuth_test() {
        // No horizontal extent: dcazm = 1 for any azimuth; the dip test then
        // decides.
        let d = Direction::new(37.0, 5.0, 1e9, -90.0, 10.0, 1e9).unwrap();
        assert_eq!(
            accepts(d, PiConstant::Exact, p(0.0, 0.0, 1.0)),
            Some(PairOrder::Swapped)
        );
        let narrow = Direction::new(37.0, 5.0, 1e9, 0.0, 10.0, 1e9).unwrap();
        assert!(
            accepts(narrow, PiConstant::Exact, p(0.0, 0.0, 1.0)).is_none(),
            "a horizontal direction with a narrow dip window rejects it"
        );
        // dtol = 90: cos(dtol) equals the pair's |dcdec| bit for bit (both
        // are cos(90°·π/180)), so the pair is on the edge and passes.
        assert!(
            accepts(horizontal(37.0, 5.0), PiConstant::Exact, p(0.0, 0.0, 1.0)).is_some(),
            "dtol = 90"
        );
    }

    #[test]
    fn horizontal_bandwidth_limits_perpendicular_offset() {
        let d = Direction::new(0.0, 45.0, 1.0, 0.0, 90.0, 1e9).unwrap();
        // 10 north and 1 east: offset exactly 1 (accepted), then just over.
        assert!(
            accepts(d, PiConstant::Exact, p(1.0, 10.0, 0.0)).is_some(),
            "on the band"
        );
        assert!(
            accepts(d, PiConstant::Exact, p(1.001, 10.0, 0.0)).is_none(),
            "outside"
        );
    }

    #[test]
    fn dip_is_negative_downward() {
        // Dip −45 points down-north; a pair going 1 north and 1 down matches.
        let d = Direction::new(0.0, 10.0, 1e9, -45.0, 10.0, 1e9).unwrap();
        assert_eq!(
            accepts(d, PiConstant::Exact, p(0.0, 1.0, -1.0)),
            Some(PairOrder::AsGiven)
        );
        assert!(
            accepts(d, PiConstant::Exact, p(0.0, 1.0, 1.0)).is_none(),
            "up-north"
        );
        // Reversed pair (up-south) is the same line: accepted, swapped.
        assert_eq!(
            accepts(d, PiConstant::Exact, p(0.0, -1.0, 1.0)),
            Some(PairOrder::Swapped)
        );
    }

    #[test]
    fn dip_is_measured_in_the_pairs_own_vertical_plane() {
        // Azimuth 0 ± 45, dip −30 ± 5. The pair runs 40° east of north and
        // dips 30° in its own vertical plane. A cone about the direction
        // vector would reject it (the angle between them is about 35°); GSLIB
        // accepts it.
        let d = Direction::new(0.0, 45.0, 1e9, -30.0, 5.0, 1e9).unwrap();
        let (s, c) = 40.0_f64.to_radians().sin_cos();
        let (sd, cd) = 30.0_f64.to_radians().sin_cos();
        let head = p(cd * s, cd * c, -sd);
        assert!(accepts(d, PiConstant::Exact, head).is_some());
    }

    #[test]
    fn vertical_bandwidth_applies_even_when_omnidirectional() {
        let d = Direction::new(0.0, 90.0, 1e9, 0.0, 90.0, 0.5).unwrap();
        assert!(d.is_omnidirectional(), "atol = 90");
        assert!(
            accepts(d, PiConstant::Exact, p(0.0, 10.0, 0.4)).is_some(),
            "within"
        );
        assert!(
            accepts(d, PiConstant::Exact, p(0.0, 10.0, 0.6)).is_none(),
            "outside"
        );
    }

    #[test]
    fn coincident_pairs_skip_the_dip_test() {
        let d = Direction::new(0.0, 10.0, 1e9, -60.0, 1.0, 1e9).unwrap();
        assert_eq!(
            accepts(d, PiConstant::Exact, ORIGIN),
            Some(PairOrder::AsGiven)
        );
    }

    // --- W3, F-020: boundary pairs and π -------------------------------------

    /// F-020, computed by hand from GAMV.FOR's operation order with
    /// PI = 3.14159265, azm = 0, atol = 45:
    /// - the azimuth vector is (cos a, sin a), a = 90·PI/180 ≈ π/2 − 1.8e-9, so
    ///   it leans about 1.8e-9 rad east of north;
    /// - (1, 1) is then about 1.8e-9 rad closer than 45°: |dcazm| − cos(atol)
    ///   ≈ +6.3e-10, accepted;
    /// - (−1, 1) is about 1.8e-9 rad farther: ≈ −1.9e-9, rejected.
    #[test]
    fn gslib_pi_treats_mirror_diagonals_differently() {
        let d = horizontal(0.0, 45.0);
        assert!(
            accepts(d, PiConstant::Gslib, p(1.0, 1.0, 0.0)).is_some(),
            "NE"
        );
        assert!(
            accepts(d, PiConstant::Gslib, p(-1.0, 1.0, 0.0)).is_none(),
            "NW"
        );
    }

    /// F-020 (corrected): exact π does **not** make mirror-image edge pairs
    /// symmetric in general. They lie within about 1e-16 of the edge, so the
    /// outcome is whatever IEEE arithmetic gives in GSLIB's operation order:
    /// azimuths 0 and 90 reject both diagonals, but azimuth 45 accepts east
    /// and rejects north. Recorded so a change in arithmetic is noticed.
    #[test]
    fn exact_pi_edge_pairs_are_decided_by_rounding() {
        let cases = [
            (0.0, (1.0, 1.0), false),
            (0.0, (-1.0, 1.0), false),
            (90.0, (1.0, 1.0), false),
            (90.0, (1.0, -1.0), false),
            (45.0, (1.0, 0.0), true),
            (45.0, (0.0, 1.0), false),
        ];
        for (azm, (x, y), expected) in cases {
            let got = accepts(horizontal(azm, 45.0), PiConstant::Exact, p(x, y, 0.0));
            assert_eq!(got.is_some(), expected, "azm {azm}, pair ({x}, {y})");
        }
    }

    /// F-023: `atol = 90` is omnidirectional, but cos(90°) is slightly
    /// positive in floating point, so pairs exactly perpendicular to the
    /// azimuth fail the azimuth test. With azm = 90 that is every
    /// north–south pair. Reproduced, and flagged (`docs/spec/gamv.md`).
    #[test]
    fn atol_90_rejects_exactly_perpendicular_pairs() {
        for pi in [PiConstant::Exact, PiConstant::Gslib] {
            let east = horizontal(90.0, 90.0);
            assert!(
                accepts(east, pi, p(0.0, 1.0, 0.0)).is_none(),
                "{pi:?}: N–S pair"
            );
            assert!(
                accepts(east, pi, p(1.0, 1.0, 0.0)).is_some(),
                "{pi:?}: diagonal"
            );
            // With azm = 0 the azimuth vector and the threshold are computed
            // from the same angle, so east–west pairs sit exactly on the edge
            // and pass.
            let north = horizontal(0.0, 90.0);
            assert!(
                accepts(north, pi, p(1.0, 0.0, 0.0)).is_some(),
                "{pi:?}: E–W pair"
            );
        }
    }

    /// F-024: a pair along the azimuth but against the dip is `Swapped` in
    /// both input orders, so which datum is the tail depends on the order
    /// of the data. Reproduced (W5), and flagged.
    #[test]
    fn mixed_signs_swap_in_both_orders() {
        // Direction north, dip −45 (down-north), dtol 80. The pair goes 1
        // north and 2 up: dcazm = 1 > 0, while dcdec = (cos 45 − 2 sin 45)/√5
        // ≈ −0.32 < 0, within cos 80 ≈ 0.17.
        let dir = Direction::new(0.0, 10.0, 1e9, -45.0, 80.0, 1e9).unwrap();
        let head = p(0.0, 1.0, 2.0);
        let window = dir.window(PiConstant::Exact);
        let forward = Separation::new(ORIGIN, head);
        let backward = Separation::new(head, ORIGIN);
        assert_eq!(window.accept(&forward), Some(PairOrder::Swapped), "forward");
        assert_eq!(
            window.accept(&backward),
            Some(PairOrder::Swapped),
            "backward"
        );
    }

    // --- Properties ------------------------------------------------------

    fn any_direction() -> impl Strategy<Value = Direction> {
        (
            -360.0..360.0_f64,
            0.1..89.9_f64,
            0.0..20.0_f64,
            -90.0..90.0_f64,
            0.1..89.9_f64,
            0.0..20.0_f64,
        )
            .prop_map(|(a, at, bh, d, dt, bv)| Direction::new(a, at, bh, d, dt, bv).unwrap())
    }

    fn any_delta() -> impl Strategy<Value = Point3> {
        (-10.0..10.0_f64, -10.0..10.0_f64, -10.0..10.0_f64).prop_map(|(x, y, z)| p(x, y, z))
    }

    fn any_pi() -> impl Strategy<Value = PiConstant> {
        prop_oneof![Just(PiConstant::Exact), Just(PiConstant::Gslib)]
    }

    proptest! {
        /// IEEE negation commutes with rounding, so reversing a pair negates
        /// every intermediate exactly and acceptance is identical. The order
        /// rule (tail = given tail iff `dcazm ≥ 0` and `dcdec ≥ 0`) is not
        /// antisymmetric, though: when the two signs differ, both orders are
        /// `Swapped` (F-024). So only `AsGiven` must flip.
        #[test]
        fn reversing_a_pair_keeps_acceptance(
            d in any_direction(), head in any_delta(), pi in any_pi()
        ) {
            let f = Separation::new(ORIGIN, head);
            let b = Separation::new(head, ORIGIN);
            prop_assert_eq!(f.distance().to_bits(), b.distance().to_bits());
            {
                let w = d.window(pi);
                match (w.accept(&f), w.accept(&b)) {
                    (Some(x), Some(y)) if f.is_coincident() => prop_assert_eq!(x, y),
                    (Some(PairOrder::AsGiven), Some(y)) => {
                        prop_assert_eq!(y, PairOrder::Swapped);
                    }
                    (None, None) | (Some(_), Some(_)) => {}
                    other => prop_assert!(false, "asymmetric acceptance: {:?}", other),
                }
            }
        }

        /// With both tolerances strictly above 90°, both cosine thresholds
        /// are negative, and unbounded bandwidths accept every pair.
        #[test]
        fn tolerances_above_90_accept_everything(
            azm in -360.0..360.0_f64, dip in -90.0..90.0_f64,
            atol in 90.5..180.0_f64, dtol in 90.5..180.0_f64,
            head in any_delta(), pi in any_pi()
        ) {
            let d = Direction::new(azm, atol, f64::MAX, dip, dtol, f64::MAX).unwrap();
            let s = Separation::new(ORIGIN, head);
            prop_assert!(d.window(pi).accept(&s).is_some());
        }
    }
}
