//! Parsing variogram structure specifications from the command line.
//!
//! Format: `TYPE:CC:A_HMAX[:A_HMIN:A_VERT[:AZIMUTH:DIP:PLUNGE]]`, mirroring a
//! GSLIB parameter file's structure lines. `TYPE` is `sph`, `exp`, `gau`, or
//! `hol`. A power structure is `pow:CC:OMEGA` (unit scale, isotropic).

use gslib_rs::geometry::{Angles, Ranges};
use gslib_rs::model::{Kind, Structure};

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
}
