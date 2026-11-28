//! `sgsim` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/SGSIM/SGSIM.FOR` (`readparm`) and
//! `SGSIM.PAR`. Generated so far (RNG spec RN7): **unconditional runs with
//! nothing to search**. The data file does not exist, so `sgsim` resets
//! `ndmin`, `ndmax` and `sstrat` (lines 307–316); the search radius is below
//! the node spacing, so no simulated node is ever found. Each node is then
//! `gmean + 1·gauinv(p)` with `gmean = 0` (lines 939–958), which exposes
//! ACORN's warm-up, the random path and every deviate.
//!
//! `itrans` is written as 0: `sgsim`'s message for a missing data file says
//! it resets `itrans` too, but the code does not (F-017).

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

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
    /// The isotropic search radius; must be below every node spacing.
    pub search_radius: f64,
    /// Multiple-grid levels (`mults = 1, nmult`); 0: off.
    #[serde(default)]
    pub multiple_grid: u32,
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
        let shortest = p.grid_spacing.iter().copied().fold(f64::INFINITY, f64::min);
        if !(p.search_radius > 0.0 && p.search_radius < shortest) {
            return Err(format!(
                "search_radius must be in (0, {shortest}): only cases with nothing to \
                 search are generated so far (RNG spec RN7)"
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
        line(format!("{NO_DATA_FILE}   \\file with data"));
        line("1 2 0 3 0 0   \\columns for X,Y,Z,vr,wt,sec.var.".to_owned());
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
        line("0 8   \\min and max original data for sim".to_owned());
        line("12   \\number of simulated nodes to use".to_owned());
        line("1   \\assign data to nodes (0=no, 1=yes)".to_owned());
        line(format!(
            "{mults} {nmult}   \\multiple grid search (0=no, 1=yes),num"
        ));
        line("0   \\maximum data per octant (0=not used)".to_owned());
        line(format!("{r:?} {r:?} {r:?}   \\maximum search radii"));
        line("0.0 0.0 0.0   \\angles for search ellipsoid".to_owned());
        line("0 0.0 1.0   \\ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC".to_owned());
        line("none.dat   \\file with LVM, EXDR, or COLC variable".to_owned());
        line("4   \\column for secondary variable".to_owned());
        line("1 0.0   \\nst, nugget effect".to_owned());
        line("1 1.0 0.0 0.0 0.0   \\it,cc,ang1,ang2,ang3".to_owned());
        line("10.0 10.0 10.0   \\a_hmax, a_hmin, a_vert".to_owned());
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
            bad("search_radius = 0.5", "search_radius = 1.0"),
            "radius reaches nodes"
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
