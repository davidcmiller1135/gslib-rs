//! `sgsim` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/SGSIM/SGSIM.FOR` (`readparm`) and
//! `SGSIM.PAR`. Generated: simple kriging (`ktype = 0`, mean 0) with **no
//! transform** (`itrans = 0`), so data, if any, are normal scores already;
//! `sstrat = 0` (data and node searches); an isotropic search.
//!
//! - **Without data** the data file does not exist, so `sgsim` resets
//!   `ndmin`, `ndmax` and `sstrat` (lines 307–316). With a search radius
//!   below the node spacing nothing is ever found, and each node is
//!   `gmean + 1·gauinv(p)` with `gmean = 0` (lines 939–958): the RNG spec's
//!   RN7 cases.
//! - **With data**, or a larger radius: the statistical-parity cases of
//!   harness step 8 (R8, D-011).
//!
//! `itrans` is written as 0 even without data: `sgsim`'s message for a
//! missing data file says it resets `itrans`, but the code does not (F-017).

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::vmodel::{MAX_STRUCTURES, Structure};

/// The realizations file the recipe asks `sgsim` to write.
pub const OUTPUT_FILE: &str = "sgsim.out";
/// A data file name that the case never provides (unconditional).
pub const NO_DATA_FILE: &str = "nodata.dat";

/// Parameters of an `sgsim` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Nodes along x, y and z.
    pub grid_counts: [u32; 3],
    /// The first node's coordinates.
    pub grid_origin: [f64; 3],
    /// Node spacing.
    pub grid_spacing: [f64; 3],
    /// Number of realizations (`nsim`).
    pub realizations: u32,
    /// ACORN seed, `0 < seed < 2³⁰`.
    pub seed: u32,
    /// The isotropic search radius (data and nodes).
    pub search_radius: f64,
    /// Multiple-grid levels (`mults = 1, nmult`); 0: off.
    #[serde(default)]
    pub multiple_grid: u32,
    /// Conditioning data in normal-score units; absent: unconditional.
    #[serde(default)]
    pub data: Option<String>,
    /// 1-based columns of x, y and z in the data (0: not used).
    #[serde(default = "default_coordinates")]
    pub coordinate_columns: [u32; 3],
    /// 1-based column of the variable in the data.
    #[serde(default = "default_variable")]
    pub variable_column: u32,
    /// `ndmin`.
    #[serde(default)]
    pub min_data: u32,
    /// `ndmax`.
    #[serde(default = "default_max_data")]
    pub max_data: u32,
    /// `nodmax`.
    #[serde(default = "default_max_nodes")]
    pub max_nodes: u32,
    /// Nugget effect of the normal-score model.
    #[serde(default)]
    pub nugget: f64,
    /// Nested structures of the normal-score model (default: one spherical
    /// structure, sill 1, range 10).
    #[serde(default = "default_structures")]
    pub structures: Vec<Structure>,
}

const fn default_coordinates() -> [u32; 3] {
    [1, 2, 0]
}
const fn default_variable() -> u32 {
    3
}
const fn default_max_data() -> u32 {
    8
}
const fn default_max_nodes() -> u32 {
    12
}
fn default_structures() -> Vec<Structure> {
    vec![Structure {
        kind: 1,
        cc: 1.0,
        angles: [0.0; 3],
        ranges: [10.0; 3],
    }]
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema, or the case is not
    /// one this generator supports.
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let p: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        if p.grid_counts.contains(&0) || p.realizations == 0 {
            return Err("grid_counts and realizations must be at least 1".to_owned());
        }
        if !p.grid_spacing.iter().all(|s| s.is_finite() && *s > 0.0)
            || !p.grid_origin.iter().all(|o| o.is_finite())
        {
            return Err("grid origin must be finite and spacing positive".to_owned());
        }
        if p.seed == 0 || p.seed >= 1 << 30 {
            return Err("seed must satisfy 0 < seed < 2^30 (ACORN)".to_owned());
        }
        if !(p.search_radius.is_finite() && p.search_radius > 0.0) {
            return Err("search_radius must be finite and positive".to_owned());
        }
        if p.grid_counts[2] != 1 {
            return Err("the stock sgsim is 2-D only (MAXZ = 1, F-017)".to_owned());
        }
        if p.max_data == 0 || p.min_data > p.max_data || p.max_nodes == 0 {
            return Err("need 0 <= min_data <= max_data, max_data and max_nodes >= 1".to_owned());
        }
        if let Some(d) = &p.data
            && (d.trim().is_empty() || d.len() > 40 || p.variable_column == 0)
        {
            return Err(
                "data must be a 1-40 character file name with a 1-based variable column".to_owned(),
            );
        }
        if !(p.nugget.is_finite() && p.nugget >= 0.0)
            || p.structures.is_empty()
            || p.structures.len() > MAX_STRUCTURES
            || p.structures.iter().any(|s| {
                !(1..=5).contains(&s.kind)
                    || s.kind == 4
                    || s.cc < 0.0
                    || s.ranges.iter().any(|&a| !(a.is_finite() && a > 0.0))
                    || !s.angles.iter().all(|a| a.is_finite())
            })
        {
            return Err(format!(
                "need a nugget >= 0 and 1 to {MAX_STRUCTURES} structures of type 1, 2, 3 or 5"
            ));
        }
        Ok(p)
    }

    /// The GSLIB parameter file.
    #[must_use]
    pub fn to_par(&self) -> String {
        let [nx, ny, nz] = self.grid_counts;
        let [xmn, ymn, zmn] = self.grid_origin;
        let [xsiz, ysiz, zsiz] = self.grid_spacing;
        let r = self.search_radius;
        let (mults, nmult) = if self.multiple_grid > 0 {
            (1, self.multiple_grid)
        } else {
            (0, 1)
        };
        let mut par = String::new();
        // Writing to a String cannot fail.
        let mut line = |text: String| {
            let _ = writeln!(par, "{text}");
        };
        line("                  Parameters for SGSIM".to_owned());
        line("                  ********************".to_owned());
        line(String::new());
        line("START OF PARAMETERS:".to_owned());
        let [cx, cy, cz] = self.coordinate_columns;
        line(format!(
            "{}   \\file with data",
            self.data.as_deref().unwrap_or(NO_DATA_FILE)
        ));
        line(format!(
            "{cx} {cy} {cz} {} 0 0   \\columns for X,Y,Z,vr,wt,sec.var.",
            self.variable_column
        ));
        line("-1.0e21 1.0e21   \\trimming limits".to_owned());
        line("0   \\transform the data (0=no, 1=yes)".to_owned());
        line("sgsim.trn   \\file for output trans table".to_owned());
        line("0   \\consider ref. dist (0=no, 1=yes)".to_owned());
        line("none.dat   \\file with ref. dist distribution".to_owned());
        line("1 2   \\columns for vr and wt".to_owned());
        line("0.0 15.0   \\zmin,zmax(tail extrapolation)".to_owned());
        line("1 0.0   \\lower tail option, parameter".to_owned());
        line("1 15.0   \\upper tail option, parameter".to_owned());
        line("0   \\debugging level: 0,1,2,3".to_owned());
        line("sgsim.dbg   \\file for debugging output".to_owned());
        line(format!("{OUTPUT_FILE}   \\file for simulation output"));
        line(format!(
            "{}   \\number of realizations to generate",
            self.realizations
        ));
        line(format!("{nx} {xmn:?} {xsiz:?}   \\nx,xmn,xsiz"));
        line(format!("{ny} {ymn:?} {ysiz:?}   \\ny,ymn,ysiz"));
        line(format!("{nz} {zmn:?} {zsiz:?}   \\nz,zmn,zsiz"));
        line(format!("{}   \\random number seed", self.seed));
        line(format!(
            "{} {}   \\min and max original data for sim",
            self.min_data, self.max_data
        ));
        line(format!(
            "{}   \\number of simulated nodes to use",
            self.max_nodes
        ));
        // sstrat: without data sgsim forces 1 itself; with data, 0 (search
        // data and nodes).
        line(format!(
            "{}   \\assign data to nodes (0=no, 1=yes)",
            u8::from(self.data.is_none())
        ));
        line(format!(
            "{mults} {nmult}   \\multiple grid search (0=no, 1=yes),num"
        ));
        line("0   \\maximum data per octant (0=not used)".to_owned());
        line(format!("{r:?} {r:?} {r:?}   \\maximum search radii"));
        line("0.0 0.0 0.0   \\angles for search ellipsoid".to_owned());
        line("0 0.0 1.0   \\ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC".to_owned());
        line("none.dat   \\file with LVM, EXDR, or COLC variable".to_owned());
        line("4   \\column for secondary variable".to_owned());
        line(format!(
            "{} {:?}   \\nst, nugget effect",
            self.structures.len(),
            self.nugget
        ));
        for st in &self.structures {
            let [a1, a2, a3] = st.angles;
            let [h1, h2, h3] = st.ranges;
            line(format!(
                "{} {:?} {a1:?} {a2:?} {a3:?}   \\it,cc,ang1,ang2,ang3",
                st.kind, st.cc
            ));
            line(format!("{h1:?} {h2:?} {h3:?}   \\a_hmax, a_hmin, a_vert"));
        }
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PARAMS: &str = "grid_counts = [12, 10, 1]\ngrid_origin = [0.5, 0.5, 0.5]\ngrid_spacing = [1.0, 1.0, 1.0]\nrealizations = 2\nseed = 69069\nsearch_radius = 0.5";

    #[test]
    fn par_layout() {
        let par = Params::from_table(&PARAMS.parse().unwrap())
            .unwrap()
            .to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert!(lines[4].starts_with("nodata.dat "), "{par}");
        assert!(lines[7].starts_with("0 "), "no transform: {par}");
        assert!(lines[17].starts_with("sgsim.out "), "{par}");
        assert!(lines[18].starts_with("2 "), "{par}");
        assert!(lines[19].starts_with("12 0.5 1.0 "), "{par}");
        assert!(lines[22].starts_with("69069 "), "{par}");
        assert!(lines[26].starts_with("0 1 "), "no multiple grid: {par}");
        assert!(lines[28].starts_with("0.5 0.5 0.5 "), "{par}");
        assert!(lines[30].starts_with("0 "), "SK: {par}");
        assert_eq!(lines.len(), 36);
    }

    #[test]
    fn rejects_unsupported_cases() {
        let bad = |from: &str, to: &str| {
            Params::from_table(&PARAMS.replace(from, to).parse().unwrap()).is_err()
        };
        assert!(
            bad("search_radius = 0.5", "search_radius = 0.0"),
            "zero radius"
        );
        assert!(
            bad("grid_counts = [12, 10, 1]", "grid_counts = [12, 10, 2]"),
            "3-D"
        );
        assert!(bad("seed = 69069", "seed = 0"), "seed 0");
        assert!(bad("seed = 69069", "seed = 1073741824"), "seed 2^30");
        assert!(
            bad("realizations = 2", "realizations = 0"),
            "no realizations"
        );
        let multi =
            Params::from_table(&(PARAMS.to_owned() + "\nmultiple_grid = 2").parse().unwrap())
                .unwrap()
                .to_par();
        assert!(
            multi.lines().nth(26).unwrap().starts_with("1 2 "),
            "{multi}"
        );
    }
}
