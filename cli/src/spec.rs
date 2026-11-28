//! Parsing variogram structure and grid axis specifications from the
//! command line.
//!
//! Format: `TYPE:CC:A_HMAX[:A_HMIN:A_VERT[:AZIMUTH:DIP:PLUNGE]]`, mirroring a
//! GSLIB parameter file's structure lines. `TYPE` is `sph`, `exp`, `gau`, or
//! `hol`. A power structure is `pow:CC:OMEGA` (unit scale, isotropic).

use gslib_rs::geometry::{Angles, Ranges};
use gslib_rs::model::{Kind, Structure};
use gslib_rs::simulation::Axis;

fn number(field: &str, what: &str) -> Result<f64, String> {
    field
        .trim()
        .parse()
        .map_err(|_| format!("{what} {field:?} is not a number"))
}

/// Parse one structure specification.
///
/// # Errors
/// A human-readable message describing the first problem.
pub(crate) fn parse_structure(spec: &str) -> Result<Structure, String> {
    let fields: Vec<&str> = spec.split(':').collect();
    let (kind_name, rest) = fields
        .split_first()
        .ok_or("empty structure specification")?;
    let numbers = rest
        .iter()
        .enumerate()
        .map(|(i, f)| number(f, &format!("field {}", i + 2)))
        .collect::<Result<Vec<f64>, String>>()?;
    let err = |e: gslib_rs::Error| e.to_string();
    if *kind_name == "pow" {
        let [cc, omega] = numbers.as_slice() else {
            return Err("power structure must be pow:CC:OMEGA".to_owned());
        };
        return Structure::isotropic(Kind::Power { exponent: *omega }, *cc, 1.0).map_err(err);
    }
    let kind = match *kind_name {
        "sph" => Kind::Spherical,
        "exp" => Kind::Exponential,
        "gau" => Kind::Gaussian,
        "hol" => Kind::HoleEffect,
        other => {
            return Err(format!(
                "unknown structure type {other:?} (expected sph, exp, gau, hol, pow)"
            ));
        }
    };
    let (cc, ranges, angles) = match *numbers.as_slice() {
        [cc, a] => (cc, Ranges::isotropic(a), Angles::default()),
        [cc, major, minor, vertical] => (
            cc,
            Ranges {
                major,
                minor,
                vertical,
            },
            Angles::default(),
        ),
        [cc, major, minor, vertical, azimuth, dip, plunge] => (
            cc,
            Ranges {
                major,
                minor,
                vertical,
            },
            Angles {
                azimuth,
                dip,
                plunge,
            },
        ),
        _ => {
            return Err(format!(
                "{spec:?}: expected TYPE:CC:A_HMAX, TYPE:CC:A_HMAX:A_HMIN:A_VERT, \
                 or TYPE:CC:A_HMAX:A_HMIN:A_VERT:AZIMUTH:DIP:PLUNGE"
            ));
        }
    };
    Structure::new(kind, cc, angles, ranges).map_err(err)
}

/// Parse a grid axis `N:MIN:SIZ`, as GSLIB's `nx xmn xsiz` parameter line:
/// `N` nodes, the first at `MIN`, `SIZ` apart.
///
/// # Errors
/// A message if a field is missing or not a number, or `N` is not a
/// positive integer.
pub(crate) fn parse_axis(spec: &str) -> Result<Axis, String> {
    let fields: Vec<&str> = spec.split(':').collect();
    let [count, first, size] = fields.as_slice() else {
        return Err(format!("expected N:MIN:SIZ, got {spec:?}"));
    };
    let count: usize = count
        .trim()
        .parse()
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| format!("node count {count:?} is not a positive integer"))?;
    Ok(Axis {
        count,
        first_center: number(first, "first node")?,
        spacing: number(size, "spacing")?,
    })
}

/// Parse search angles `AZIMUTH:DIP:PLUNGE` in degrees (GSLIB's `sang1 sang2
/// sang3`).
///
/// # Errors
/// A message if there are not three numbers.
pub(crate) fn parse_angles(spec: &str) -> Result<Angles, String> {
    let fields: Vec<&str> = spec.split(':').collect();
    let [azimuth, dip, plunge] = fields.as_slice() else {
        return Err(format!("expected AZIMUTH:DIP:PLUNGE, got {spec:?}"));
    };
    Ok(Angles {
        azimuth: number(azimuth, "azimuth")?,
        dip: number(dip, "dip")?,
        plunge: number(plunge, "plunge")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gslib_rs::geometry::Point3;
    use gslib_rs::model::VariogramModel;

    fn gamma_at(spec: &str, p: Point3) -> f64 {
        let model = VariogramModel::new(0.0, vec![parse_structure(spec).unwrap()]).unwrap();
        model.gamma(Point3::default(), p)
    }

    #[test]
    fn isotropic_spec() {
        let s = parse_structure("sph:0.9:10").unwrap();
        assert_eq!(s.kind(), Kind::Spherical);
        assert!((s.contribution() - 0.9).abs() < 1e-15, "contribution");
    }

    #[test]
    fn anisotropic_spec_orients_major_axis() {
        // Major axis east (azimuth 90): range 10 east, 2 north.
        let spec = "sph:1:10:2:1:90:0:0";
        assert!(
            (gamma_at(spec, Point3::new(10.0, 0.0, 0.0)) - 1.0).abs() < 1e-12,
            "east at range"
        );
        assert!(
            (gamma_at(spec, Point3::new(0.0, 2.0, 0.0)) - 1.0).abs() < 1e-12,
            "north at range"
        );
        assert!(
            gamma_at(spec, Point3::new(2.0, 0.0, 0.0)) < 0.5,
            "east well within range"
        );
    }

    #[test]
    fn power_spec() {
        let s = parse_structure("pow:2:1.5").unwrap();
        assert_eq!(s.kind(), Kind::Power { exponent: 1.5 });
    }

    #[test]
    fn rejects_bad_specs() {
        for bad in [
            "",
            "sph",
            "sph:1",
            "sph:1:2:3",
            "xyz:1:2",
            "sph:a:2",
            "sph:1:-2",
            "pow:1:3",
            "pow:1",
        ] {
            assert!(parse_structure(bad).is_err(), "{bad:?} should be rejected");
        }
    }

    #[test]
    fn axes_parse_like_gslib_lines() {
        let a = parse_axis("50:0.5:1").unwrap();
        assert_eq!(
            (a.count, a.first_center.to_bits(), a.spacing.to_bits()),
            (50, 0.5_f64.to_bits(), 1.0_f64.to_bits())
        );
        assert!(parse_axis("10:-5.0:2").is_ok(), "negative origin");
        assert!(
            parse_axis("0:0:1")
                .unwrap_err()
                .contains("positive integer")
        );
        assert!(parse_axis("2.5:0:1").is_err(), "fractional count");
        assert!(parse_axis("10:0").unwrap_err().contains("N:MIN:SIZ"));
        assert!(parse_axis("10:a:1").unwrap_err().contains("first node"));
    }

    #[test]
    fn angles_parse() {
        let a = parse_angles("30:-10:5").unwrap();
        assert_eq!(
            [a.azimuth, a.dip, a.plunge].map(f64::to_bits),
            [30.0_f64, -10.0, 5.0].map(f64::to_bits)
        );
        assert!(
            parse_angles("30:0")
                .unwrap_err()
                .contains("AZIMUTH:DIP:PLUNGE")
        );
        assert!(parse_angles("x:0:0").unwrap_err().contains("azimuth"));
    }
}
