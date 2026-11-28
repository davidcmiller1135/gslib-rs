//! `kt3d` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/KT3D/KT3D.FOR` (`readparm`) and
//! `KT3D.PAR`: data file; `ixl iyl izl ivrl iextv`; `tmin tmax`;
//! `koption`; jackknife file; its columns; `idbg`; debug file; output file;
//! `nx xmn xsiz` (and y, z); `nxdis nydis nzdis`; `ndmin ndmax`; `noct`;
//! `radius radius1 radius2`; `sang1 sang2 sang3`; `ktype skmean`; nine
//! drift flags; `itrend`; external drift file; its column; `nst c0`; then
//! per structure `it cc ang1 ang2 ang3` and `a_hmax a_hmin a_vert`. Limits
//! from `KT3D.INC`: `MAXSAM = 120` (`ndmax`), `MAXNST = 4`.
//!
//! Generated so far (search spec step 3, harness step 7): grid estimation
//! (`koption = 0`) by simple or ordinary point kriging (`nxdis = nydis =
//! nzdis = 1`), no drift, with the debug file at level 2 (F-014, D-027).

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::vmodel::{MAX_STRUCTURES, Structure};

/// The estimates file the recipe asks `kt3d` to write.
pub const OUTPUT_FILE: &str = "kt3d.out";
/// The debugging file (`idbg = 2`): the data and weights of every solved
/// system.
pub const DEBUG_FILE: &str = "kt3d.dbg";
/// Most data per estimate (`MAXSAM`).
pub const MAX_DATA: u32 = 120;

/// Which kriging `kt3d` performs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kriging {
    /// Simple kriging around `sk_mean` (`ktype = 0`).
    Simple,
    /// Ordinary kriging (`ktype = 1`).
    Ordinary,
}

/// Parameters of a `kt3d` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Input data file; must be listed in the case's `inputs`.
    pub data: String,
    /// 1-based columns of x, y and z (0: the grid origin on that axis).
    pub coordinate_columns: [u32; 3],
    /// 1-based column of the variable.
    pub variable_column: u32,
    /// `[tmin, tmax]`; absent: no trimming.
    #[serde(default)]
    pub trimming_limits: Option<[f64; 2]>,
    /// Nodes along x, y and z (`nx ny nz`).
    pub grid_counts: [u32; 3],
    /// The first node's coordinates (`xmn ymn zmn`).
    pub grid_origin: [f64; 3],
    /// Node spacing (`xsiz ysiz zsiz`).
    pub grid_spacing: [f64; 3],
    /// `ndmin`.
    pub min_data: u32,
    /// `ndmax`.
    pub max_data: u32,
    /// `noct`; 0: no octant search.
    #[serde(default)]
    pub per_octant: u32,
    /// `radius radius1 radius2`.
    pub search_radii: [f64; 3],
    /// `sang1 sang2 sang3`.
    #[serde(default)]
    pub search_angles: [f64; 3],
    /// Simple or ordinary kriging.
    pub kriging: Kriging,
    /// `skmean` (used by simple kriging only).
    #[serde(default)]
    pub sk_mean: f64,
    /// Nugget effect `c0`.
    pub nugget: f64,
    /// Nested structures (at least one, as `kt3d` requires).
    pub structures: Vec<Structure>,
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema, or a value is outside
    /// what both `kt3d` and gslib-rs accept.
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let p: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        p.validate()?;
        Ok(p)
    }

    fn validate(&self) -> Result<(), String> {
        if self.variable_column == 0 {
            return Err("variable_column is 1-based".to_owned());
        }
        if self.data.trim().is_empty() || self.data.len() > 40 {
            return Err("data must be a file name of 1-40 characters (GSLIB reads a40)".to_owned());
        }
        crate::check_trimming(self.trimming_limits)?;
        if self.grid_counts.contains(&0) {
            return Err("grid_counts must be at least 1".to_owned());
        }
        let positive = |v: &[f64]| v.iter().all(|x| x.is_finite() && *x > 0.0);
        if !positive(&self.grid_spacing) || !self.grid_origin.iter().all(|x| x.is_finite()) {
            return Err("grid origin must be finite and spacing positive".to_owned());
        }
        if self.max_data == 0 || self.max_data > MAX_DATA || self.min_data > self.max_data {
            return Err(format!(
                "need 0 <= min_data <= max_data <= {MAX_DATA}, max_data >= 1"
            ));
        }
        if !positive(&self.search_radii) || !self.search_angles.iter().all(|x| x.is_finite()) {
            return Err("search radii must be positive and angles finite".to_owned());
        }
        if !(self.sk_mean.is_finite() && self.nugget.is_finite() && self.nugget >= 0.0) {
            return Err("sk_mean must be finite and the nugget finite and >= 0".to_owned());
        }
        if self.structures.is_empty() || self.structures.len() > MAX_STRUCTURES {
            return Err(format!("need 1 to {MAX_STRUCTURES} structures"));
        }
        for s in &self.structures {
            if !(1..=5).contains(&s.kind) || s.kind == 4 {
                // The power model has no covariance, which kt3d fakes (PMX);
                // gslib-rs refuses it for simple kriging (ROADMAP §2).
                return Err(format!(
                    "structure type {} is not 1, 2, 3 or 5 (power is out of scope)",
                    s.kind
                ));
            }
            let numbers = std::iter::once(s.cc).chain(s.angles).chain(s.ranges);
            if numbers.into_iter().any(|v: f64| !v.is_finite()) || s.cc < 0.0 {
                return Err("structures need finite values and cc >= 0".to_owned());
            }
            if s.ranges.iter().any(|&a| a <= 0.0) {
                return Err("ranges must be positive".to_owned());
            }
        }
        Ok(())
    }

    /// The GSLIB parameter file. Numbers use Rust's shortest round-trip
    /// form, which gfortran's list-directed READ parses to the same `f64`.
    #[must_use]
    pub fn to_par(&self) -> String {
        let [cx, cy, cz] = self.coordinate_columns;
        let [nx, ny, nz] = self.grid_counts;
        let [xmn, ymn, zmn] = self.grid_origin;
        let [xsiz, ysiz, zsiz] = self.grid_spacing;
        let [r0, r1, r2] = self.search_radii;
        let [s1, s2, s3] = self.search_angles;
        let ktype = match self.kriging {
            Kriging::Simple => 0,
            Kriging::Ordinary => 1,
        };
        let mut par = String::new();
        // Writing to a String cannot fail.
        let mut line = |text: String| {
            let _ = writeln!(par, "{text}");
        };
        line("                  Parameters for KT3D".to_owned());
        line("                  *******************".to_owned());
        line(String::new());
        line("START OF PARAMETERS:".to_owned());
        line(format!("{}   \\file with data", self.data));
        line(format!(
            "{cx} {cy} {cz} {} 0   \\columns for X, Y, Z, var, sec var",
            self.variable_column
        ));
        line(format!(
            "{}   \\trimming limits",
            crate::trimming_line(self.trimming_limits)
        ));
        line("0   \\option: 0=grid, 1=cross, 2=jackknife".to_owned());
        line("none.dat   \\file with jackknife data".to_owned());
        line("1 2 0 3 0   \\columns for X,Y,Z,vr and sec var".to_owned());
        line("2   \\debugging level: 0,1,2,3".to_owned());
        line(format!("{DEBUG_FILE}   \\file for debugging output"));
        line(format!("{OUTPUT_FILE}   \\file for kriged output"));
        line(format!("{nx} {xmn:?} {xsiz:?}   \\nx,xmn,xsiz"));
        line(format!("{ny} {ymn:?} {ysiz:?}   \\ny,ymn,ysiz"));
        line(format!("{nz} {zmn:?} {zsiz:?}   \\nz,zmn,zsiz"));
        line("1 1 1   \\x,y and z block discretization".to_owned());
        line(format!(
            "{} {}   \\min, max data for kriging",
            self.min_data, self.max_data
        ));
        line(format!(
            "{}   \\max per octant (0-> not used)",
            self.per_octant
        ));
        line(format!("{r0:?} {r1:?} {r2:?}   \\maximum search radii"));
        line(format!(
            "{s1:?} {s2:?} {s3:?}   \\angles for search ellipsoid"
        ));
        line(format!(
            "{ktype} {:?}   \\0=SK,1=OK,2=non-st SK,3=exdrift",
            self.sk_mean
        ));
        line("0 0 0 0 0 0 0 0 0   \\drift: x,y,z,xx,yy,zz,xy,xz,zy".to_owned());
        line("0   \\0, variable; 1, estimate trend".to_owned());
        line("none.dat   \\gridded file with drift/mean".to_owned());
        line("4   \\column number in gridded file".to_owned());
        line(format!(
            "{} {:?}   \\nst, nugget effect",
            self.structures.len(),
            self.nugget
        ));
        for s in &self.structures {
            let [a1, a2, a3] = s.angles;
            let [h1, h2, h3] = s.ranges;
            line(format!(
                "{} {:?} {a1:?} {a2:?} {a3:?}   \\it,cc,ang1,ang2,ang3",
                s.kind, s.cc
            ));
            line(format!("{h1:?} {h2:?} {h3:?}   \\a_hmax, a_hmin, a_vert"));
        }
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PARAMS: &str = r#"
data = "smooth.dat"
coordinate_columns = [1, 2, 0]
variable_column = 4
grid_counts = [8, 8, 1]
grid_origin = [5.0, 5.0, 0.0]
grid_spacing = [12.0, 12.0, 1.0]
min_data = 4
max_data = 12
search_radii = [40.0, 40.0, 10.0]
kriging = "ordinary"
nugget = 0.1
structures = [{ kind = 1, cc = 0.9, angles = [0.0, 0.0, 0.0], ranges = [50.0, 50.0, 50.0] }]
"#;

    #[test]
    fn par_layout() {
        let p = Params::from_table(&PARAMS.parse().unwrap()).unwrap();
        let par = p.to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert_eq!(lines[3], "START OF PARAMETERS:");
        assert!(lines[4].starts_with("smooth.dat "), "{par}");
        assert!(lines[5].starts_with("1 2 0 4 0 "), "{par}");
        assert!(lines[6].starts_with("-1.0e21 1.0e21 "), "{par}");
        assert!(lines[7].starts_with("0 "), "grid option: {par}");
        assert!(lines[10].starts_with("2 "), "debug level: {par}");
        assert!(lines[11].starts_with("kt3d.dbg "), "{par}");
        assert!(lines[12].starts_with("kt3d.out "), "{par}");
        assert!(lines[13].starts_with("8 5.0 12.0 "), "{par}");
        assert!(lines[16].starts_with("1 1 1 "), "point kriging: {par}");
        assert!(lines[17].starts_with("4 12 "), "{par}");
        assert!(lines[18].starts_with("0 "), "no octants: {par}");
        assert!(lines[19].starts_with("40.0 40.0 10.0 "), "{par}");
        assert!(lines[21].starts_with("1 0.0 "), "OK: {par}");
        assert!(
            lines[22].starts_with("0 0 0 0 0 0 0 0 0 "),
            "no drift: {par}"
        );
        assert!(lines[26].starts_with("1 0.1 "), "nst, c0: {par}");
        assert!(lines[27].starts_with("1 0.9 0.0 0.0 0.0 "), "{par}");
        assert!(lines[28].starts_with("50.0 50.0 50.0 "), "{par}");
        assert_eq!(lines.len(), 29);
    }

    #[test]
    fn rejects_out_of_range() {
        let bad = |from: &str, to: &str| {
            Params::from_table(&PARAMS.replace(from, to).parse().unwrap()).is_err()
        };
        assert!(bad("max_data = 12", "max_data = 121"), "MAXSAM");
        assert!(bad("min_data = 4", "min_data = 13"), "min above max");
        assert!(bad("kind = 1", "kind = 4"), "power model");
        assert!(
            bad("grid_counts = [8, 8, 1]", "grid_counts = [8, 0, 1]"),
            "empty axis"
        );
        assert!(
            bad(
                "search_radii = [40.0, 40.0, 10.0]",
                "search_radii = [40.0, 0.0, 10.0]"
            ),
            "zero radius"
        );
        assert!(
            bad("kriging = \"ordinary\"", "kriging = \"universal\""),
            "unknown kriging"
        );
    }
}
