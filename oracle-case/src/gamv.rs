//! `gamv` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/GAMV/GAMV.FOR` (`readparm`): data file;
//! `ixl iyl izl`; `nvar ivar...`; `tmin tmax`; output file; `nlag`; `xlag`;
//! `xltol`; `ndir`; per direction `azm atol bandwh dip dtol bandwd`;
//! `isill`; `nvarg`; per variogram `ivtail ivhead ivtype`. Limits from
//! `GAMV.INC`: `MAXLAG = 200`.
//!
//! Only the omnidirectional semivariogram of one variable is generated so far
//! (plan step 6): one direction with `atol = dtol = 90` and bandwidths of
//! `1e21`, which `gamv` treats as omnidirectional (`atol >= 90`).

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// The output file the recipe asks `gamv` to write.
pub const OUTPUT_FILE: &str = "gamv.out";
/// Most lags `gamv` accepts (`MAXLAG`).
pub const MAX_LAGS: u32 = 200;

/// Parameters of a `gamv` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Input data file; must be listed in the case's `inputs`.
    pub data: String,
    /// 1-based columns of x, y and z (0 means "not used": coordinate 0).
    pub coordinate_columns: [u32; 3],
    /// 1-based column of the variable.
    pub variable_column: u32,
    /// Number of lags (`nlag`).
    pub nlag: u32,
    /// Lag spacing (`xlag`).
    pub lag_spacing: f64,
    /// Lag tolerance (`xltol`). Must be positive: GSLIB would silently
    /// replace a non-positive value with `xlag / 2`.
    pub lag_tolerance: f64,
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema, a column is invalid,
    /// `nlag` is outside `1..=200`, or the spacing or tolerance is not finite
    /// and positive.
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let p: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        if p.variable_column == 0 {
            return Err("variable_column is 1-based".to_owned());
        }
        if !(1..=MAX_LAGS).contains(&p.nlag) {
            return Err(format!("nlag must be 1 to {MAX_LAGS}"));
        }
        for (what, v) in [
            ("lag_spacing", p.lag_spacing),
            ("lag_tolerance", p.lag_tolerance),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(format!("{what} must be finite and positive"));
            }
        }
        if p.data.trim().is_empty() || p.data.len() > 40 {
            return Err("data must be a file name of 1-40 characters (GSLIB reads a40)".to_owned());
        }
        Ok(p)
    }

    /// The GSLIB parameter file: one omnidirectional semivariogram, no
    /// trimming, no sill standardization.
    #[must_use]
    pub fn to_par(&self) -> String {
        let [x, y, z] = self.coordinate_columns;
        let mut par = String::new();
        // Writing to a String cannot fail.
        let _ = write!(
            par,
            "                  Parameters for GAMV\n\
             \x20                 *******************\n\
             \n\
             START OF PARAMETERS:\n\
             {}   \\file with data\n\
             {x} {y} {z}   \\columns for X, Y, Z coordinates\n\
             1 {}   \\number of variables, column numbers\n\
             -1.0e21 1.0e21   \\trimming limits\n\
             {OUTPUT_FILE}   \\file for variogram output\n\
             {}   \\number of lags\n\
             {}   \\lag separation distance\n\
             {}   \\lag tolerance\n\
             1   \\number of directions\n\
             0.0 90.0 1.0e21 0.0 90.0 1.0e21   \\azm,atol,bandh,dip,dtol,bandv\n\
             0   \\standardize sills? (0=no, 1=yes)\n\
             1   \\number of variograms\n\
             1 1 1   \\tail var., head var., variogram type\n",
            self.data, self.variable_column, self.nlag, self.lag_spacing, self.lag_tolerance
        );
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PARAMS: &str = "data = \"s.dat\"\ncoordinate_columns = [1, 2, 0]\nvariable_column = 4\nnlag = 10\nlag_spacing = 5.0\nlag_tolerance = 2.5";

    #[test]
    fn par_layout() {
        let par = Params::from_table(&PARAMS.parse().unwrap())
            .unwrap()
            .to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert!(lines[4].starts_with("s.dat "), "{par}");
        assert!(lines[5].starts_with("1 2 0 "), "{par}");
        assert!(lines[6].starts_with("1 4 "), "{par}");
        assert!(lines[8].starts_with("gamv.out "), "{par}");
        assert!(lines[9].starts_with("10 "), "{par}");
        assert!(
            lines[10].starts_with("5 ") && lines[11].starts_with("2.5 "),
            "{par}"
        );
        assert!(
            lines[13].starts_with("0.0 90.0 1.0e21 0.0 90.0 1.0e21"),
            "{par}"
        );
        assert!(lines[16].starts_with("1 1 1 "), "{par}");
        assert_eq!(lines.len(), 17);
    }

    #[test]
    fn rejects_out_of_range() {
        let bad = |from: &str, to: &str| {
            Params::from_table(&PARAMS.replace(from, to).parse().unwrap()).is_err()
        };
        assert!(bad("nlag = 10", "nlag = 201"), "MAXLAG");
        assert!(
            bad("lag_tolerance = 2.5", "lag_tolerance = 0.0"),
            "explicit tolerance"
        );
        assert!(bad("variable_column = 4", "variable_column = 0"), "1-based");
    }
}
