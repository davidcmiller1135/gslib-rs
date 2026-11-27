//! Points, distances, and GSLIB's anisotropy rotation (`setrot`/`sqdist`).

use crate::error::{Error, Result};

/// A location in 3-D space.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point3 {
    /// Easting.
    pub x: f64,
    /// Northing.
    pub y: f64,
    /// Elevation.
    pub z: f64,
}

impl Point3 {
    /// A point from its coordinates.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Component-wise `self - other`.
    #[must_use]
    pub fn delta(self, other: Self) -> [f64; 3] {
        [self.x - other.x, self.y - other.y, self.z - other.z]
    }

    /// Euclidean distance.
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        norm(self.delta(other))
    }

    /// True if every coordinate is finite.
    #[must_use]
    pub const fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

/// Two points closer than this are treated as the same location.
pub const COINCIDENT: f64 = 1e-10;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let [a0, a1, a2] = a;
    let [b0, b1, b2] = b;
    a0.mul_add(b0, a1.mul_add(b1, a2 * b2))
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

/// Orientation of an anisotropy ellipsoid, in degrees, using GSLIB's
/// conventions (verified against `SETROT.FOR` and `GAMV.FOR`; F-004).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Angles {
    /// Azimuth of the major axis, clockwise from north (+y).
    pub azimuth: f64,
    /// Dip of the major axis from horizontal. GSLIB convention: **negative is
    /// downward**.
    pub dip: f64,
    /// Rotation (plunge) of the minor axes about the major axis.
    pub plunge: f64,
}

/// Ranges of an anisotropy ellipsoid along its three principal axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ranges {
    /// Range along the major horizontal axis (`a_hmax`).
    pub major: f64,
    /// Range along the minor horizontal axis (`a_hmin`).
    pub minor: f64,
    /// Range along the vertical axis (`a_vert`).
    pub vertical: f64,
}

impl Ranges {
    /// The same range in every direction.
    #[must_use]
    pub const fn isotropic(range: f64) -> Self {
        Self {
            major: range,
            minor: range,
            vertical: range,
        }
    }
}

/// Maps a separation vector into a frame where the anisotropy ellipsoid
/// becomes a sphere of radius `major`.
///
/// This is GSLIB's `setrot`: `scaled_distance` returns GSLIB's
/// `sqrt(sqdist(...))`, so a point at distance `ranges.minor` along the minor
/// axis is at scaled distance `ranges.major`.
///
/// The matrix matches `vendor/gslib-2.0/src/GSLIB/SETROT.FOR` term for term
/// (F-004). Two divergences:
/// - GSLIB clamps anisotropy ratios to `EPSLON = 1e-20`; here non-positive
///   ranges are an error.
/// - GSLIB converts degrees with `π ≈ 3.141592654`; here the exact value is
///   used (≈1e-10 relative difference in angles).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Anisotropy {
    rows: [[f64; 3]; 3],
}

impl Anisotropy {
    /// No rotation and no scaling: `scaled_distance` is Euclidean distance.
    #[must_use]
    pub const fn isotropic() -> Self {
        Self {
            rows: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Build the rotation for an ellipsoid with the given orientation and
    /// ranges.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if an angle is not finite or a range is not
    /// finite and positive.
    pub fn new(angles: Angles, ranges: Ranges) -> Result<Self> {
        let Angles {
            azimuth,
            dip,
            plunge,
        } = angles;
        if !(azimuth.is_finite() && dip.is_finite() && plunge.is_finite()) {
            return Err(Error::invalid("anisotropy angles", "must be finite"));
        }
        for range in [ranges.major, ranges.minor, ranges.vertical] {
            if !(range.is_finite() && range > 0.0) {
                return Err(Error::invalid(
                    "anisotropy ranges",
                    format!("must be finite and positive, got {range}"),
                ));
            }
        }

        // GSLIB branches on azimuth < 270 to keep alpha in range; the two
        // branches differ by exactly 360 degrees, so sin/cos are identical.
        let alpha = (90.0 - azimuth).to_radians();
        let beta = (-dip).to_radians();
        let theta = plunge.to_radians();
        let (sin_az, cos_az) = alpha.sin_cos();
        let (sin_dip, cos_dip) = beta.sin_cos();
        let (sin_pl, cos_pl) = theta.sin_cos();
        let afac1 = ranges.major / ranges.minor;
        let afac2 = ranges.major / ranges.vertical;

        Ok(Self {
            rows: [
                [cos_dip * cos_az, cos_dip * sin_az, -sin_dip],
                [
                    afac1 * (sin_pl * sin_dip).mul_add(cos_az, -cos_pl * sin_az),
                    afac1 * (sin_pl * sin_dip).mul_add(sin_az, cos_pl * cos_az),
                    afac1 * sin_pl * cos_dip,
                ],
                [
                    afac2 * (cos_pl * sin_dip).mul_add(cos_az, sin_pl * sin_az),
                    afac2 * (cos_pl * sin_dip).mul_add(sin_az, -sin_pl * cos_az),
                    afac2 * cos_pl * cos_dip,
                ],
            ],
        })
    }

    /// Anisotropic distance between two points, in units of the major range.
    #[must_use]
    pub fn scaled_distance(&self, a: Point3, b: Point3) -> f64 {
        norm(self.scaled(a, b))
    }

    /// The square of [`Anisotropy::scaled_distance`], computed without a
    /// square root, as GSLIB's `sqdist` does (search limits compare it with
    /// `radius²`).
    #[must_use]
    pub fn squared_scaled_distance(&self, a: Point3, b: Point3) -> f64 {
        let v = self.scaled(a, b);
        dot(v, v)
    }

    fn scaled(&self, a: Point3, b: Point3) -> [f64; 3] {
        let d = a.delta(b);
        let [r0, r1, r2] = self.rows;
        [dot(r0, d), dot(r1, d), dot(r2, d)]
    }
}

/// A 2-D rotation of coordinates about an origin (GSLIB `rotcoord`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotation2 {
    origin: (f64, f64),
    sin: f64,
    cos: f64,
}

impl Rotation2 {
    /// Rotation by `angle_radians` (counter-clockwise) about `origin`.
    #[must_use]
    pub fn new(origin: (f64, f64), angle_radians: f64) -> Self {
        let (sin, cos) = angle_radians.sin_cos();
        Self { origin, sin, cos }
    }

    /// Original coordinates to rotated coordinates.
    #[must_use]
    pub fn to_rotated(self, (x, y): (f64, f64)) -> (f64, f64) {
        let dx = x - self.origin.0;
        let dy = y - self.origin.1;
        (
            dy.mul_add(self.sin, dx * self.cos),
            dy.mul_add(self.cos, -dx * self.sin),
        )
    }

    /// Rotated coordinates back to original coordinates.
    #[must_use]
    pub fn to_original(self, (rx, ry): (f64, f64)) -> (f64, f64) {
        (
            self.origin.0 + rx.mul_add(self.cos, -ry * self.sin),
            self.origin.1 + rx.mul_add(self.sin, ry * self.cos),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    const EPS: f64 = 1e-12;
    const ORIGIN: Point3 = Point3::new(0.0, 0.0, 0.0);

    fn ranges(major: f64, minor: f64, vertical: f64) -> Ranges {
        Ranges {
            major,
            minor,
            vertical,
        }
    }

    #[test]
    fn distance_is_euclidean() {
        assert_abs_diff_eq!(ORIGIN.distance(Point3::new(3.0, 4.0, 12.0)), 13.0);
    }

    #[test]
    fn isotropic_rotation_matches_euclidean_distance() {
        let rot = Anisotropy::new(Angles::default(), Ranges::isotropic(5.0)).unwrap();
        let p = Point3::new(1.0, -2.0, 3.0);
        assert_abs_diff_eq!(
            rot.scaled_distance(ORIGIN, p),
            ORIGIN.distance(p),
            epsilon = EPS
        );
    }

    #[test]
    fn azimuth_zero_puts_major_axis_north() {
        let rot = Anisotropy::new(Angles::default(), ranges(10.0, 2.0, 1.0)).unwrap();
        // A point at the major range due north is at scaled distance = major.
        let north = Point3::new(0.0, 10.0, 0.0);
        assert_abs_diff_eq!(rot.scaled_distance(ORIGIN, north), 10.0, epsilon = EPS);
        // A point at the minor range due east is also at scaled distance = major.
        let east = Point3::new(2.0, 0.0, 0.0);
        assert_abs_diff_eq!(rot.scaled_distance(ORIGIN, east), 10.0, epsilon = EPS);
        // And the vertical range straight down.
        let down = Point3::new(0.0, 0.0, -1.0);
        assert_abs_diff_eq!(rot.scaled_distance(ORIGIN, down), 10.0, epsilon = EPS);
    }

    #[test]
    fn azimuth_ninety_puts_major_axis_east() {
        let angles = Angles {
            azimuth: 90.0,
            ..Angles::default()
        };
        let rot = Anisotropy::new(angles, ranges(10.0, 2.0, 1.0)).unwrap();
        let east = Point3::new(10.0, 0.0, 0.0);
        assert_abs_diff_eq!(rot.scaled_distance(ORIGIN, east), 10.0, epsilon = EPS);
        let north = Point3::new(0.0, 2.0, 0.0);
        assert_abs_diff_eq!(rot.scaled_distance(ORIGIN, north), 10.0, epsilon = EPS);
    }

    #[test]
    fn negative_dip_tilts_major_axis_downward() {
        let angles = Angles {
            dip: -30.0,
            ..Angles::default()
        };
        let rot = Anisotropy::new(angles, ranges(10.0, 2.0, 1.0)).unwrap();
        // Major axis points north and 30 degrees down.
        let (s, c) = 30.0_f64.to_radians().sin_cos();
        let along_major = Point3::new(0.0, 10.0 * c, -10.0 * s);
        assert_abs_diff_eq!(
            rot.scaled_distance(ORIGIN, along_major),
            10.0,
            epsilon = 1e-9
        );
    }

    #[test]
    fn rejects_non_positive_ranges() {
        assert!(
            Anisotropy::new(Angles::default(), ranges(10.0, 0.0, 1.0)).is_err(),
            "zero minor range must be rejected"
        );
        assert!(
            Anisotropy::new(Angles::default(), ranges(f64::NAN, 1.0, 1.0)).is_err(),
            "NaN range must be rejected"
        );
    }

    #[test]
    fn rotation2_known_value() {
        let rot = Rotation2::new((0.0, 0.0), std::f64::consts::FRAC_PI_2);
        let (rx, ry) = rot.to_rotated((1.0, 0.0));
        assert_abs_diff_eq!(rx, 0.0, epsilon = EPS);
        assert_abs_diff_eq!(ry, -1.0, epsilon = EPS);
    }

    proptest! {
        #[test]
        fn scaled_distance_is_symmetric(
            az in 0.0..360.0_f64, dip in -90.0..90.0_f64, pl in -90.0..90.0_f64,
            minor in 0.1..10.0_f64, vert in 0.1..10.0_f64,
            ax in -50.0..50.0_f64, ay in -50.0..50.0_f64, az_ in -50.0..50.0_f64,
            bx in -50.0..50.0_f64, by in -50.0..50.0_f64, bz in -50.0..50.0_f64,
        ) {
            let angles = Angles { azimuth: az, dip, plunge: pl };
            let rot = Anisotropy::new(angles, ranges(10.0, minor, vert)).unwrap();
            let a = Point3::new(ax, ay, az_);
            let b = Point3::new(bx, by, bz);
            prop_assert!((rot.scaled_distance(a, b) - rot.scaled_distance(b, a)).abs() < 1e-9);
        }

        #[test]
        fn rotation_preserves_lengths_when_isotropic(
            az in 0.0..360.0_f64, dip in -90.0..90.0_f64, pl in -90.0..90.0_f64,
            x in -50.0..50.0_f64, y in -50.0..50.0_f64, z in -50.0..50.0_f64,
        ) {
            let angles = Angles { azimuth: az, dip, plunge: pl };
            let rot = Anisotropy::new(angles, Ranges::isotropic(3.0)).unwrap();
            let p = Point3::new(x, y, z);
            prop_assert!((rot.scaled_distance(ORIGIN, p) - ORIGIN.distance(p)).abs() < 1e-9);
        }

        #[test]
        fn rotation2_roundtrips(
            ox in -100.0..100.0_f64, oy in -100.0..100.0_f64, angle in -7.0..7.0_f64,
            x in -100.0..100.0_f64, y in -100.0..100.0_f64,
        ) {
            let rot = Rotation2::new((ox, oy), angle);
            let (bx, by) = rot.to_original(rot.to_rotated((x, y)));
            prop_assert!((bx - x).abs() < 1e-9 && (by - y).abs() < 1e-9);
        }
    }
}
