//! `backtr` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/BACKTR/BACKTR.FOR`: data file; `ivr`;
//! `tmin tmax`; output file; transformation-table file; `zmin zmax`;
//! `ltail ltpar`; `utail utpar`.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// Output file with the data plus a back-transformed column.
pub const OUTPUT_FILE: &str = "backtr.out";

/// Parameters of a `backtr` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// File of normal scores to back-transform; must be in `inputs`.
    pub data: String,
    /// 1-based column of the normal scores.
    pub variable_column: u32,
    /// Transformation table (`value score` per line); must be in `inputs`.
    pub table: String,
    /// Lower bound for the lower tail.
    pub zmin: f64,
    /// Upper bound for the upper tail.
    pub zmax: f64,
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema, a column is 0, or
    /// `zmin`/`zmax` are not finite with `zmin < zmax`.
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let p: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        if p.variable_column == 0 {
            return Err("variable_column is 1-based".to_owned());
        }
        if !(p.zmin.is_finite() && p.zmax.is_finite() && p.zmin < p.zmax) {
            return Err("zmin and zmax must be finite with zmin < zmax".to_owned());
        }
        for name in [&p.data, &p.table] {
            if name.trim().is_empty() || name.len() > 40 {
                return Err("file names must be 1-40 characters (GSLIB reads a40)".to_owned());
            }
        }
        Ok(p)
    }

    /// The GSLIB parameter file. Both tails use option 1 (linear in
    /// probability), the only tail gslib-rs implements so far.
    #[must_use]
    pub fn to_par(&self) -> String {
        let mut par = String::new();
        // Writing to a String cannot fail.
        let _ = write!(
            par,
            "                  Parameters for BACKTR\n\
             \x20                 *********************\n\
             \n\
             START OF PARAMETERS:\n\
             {}   \\file with data\n\
             {}   \\column with Gaussian variable\n\
             -1.0e21 1.0e21   \\trimming limits\n\
             {OUTPUT_FILE}   \\file for output\n\
             {}   \\file with input transformation table\n\
             {} {}   \\minimum and maximum data value\n\
             1 1.0   \\lower tail option and parameter\n\
             1 1.0   \\upper tail option and parameter\n",
            self.data, self.variable_column, self.table, self.zmin, self.zmax
        );
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_layout() {
        let t: toml::Table =
            "data = \"s.dat\"\nvariable_column = 4\ntable = \"t.trn\"\nzmin = 0.0\nzmax = 50.0"
                .parse()
                .unwrap();
        let par = Params::from_table(&t).unwrap().to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert!(lines[4].starts_with("s.dat "), "{par}");
        assert!(lines[5].starts_with("4 "), "{par}");
        assert!(lines[7].starts_with("backtr.out "), "{par}");
        assert!(lines[8].starts_with("t.trn "), "{par}");
        assert!(lines[9].starts_with("0 50 "), "{par}");
        assert!(
            lines[10].starts_with("1 1.0 ") && lines[11].starts_with("1 1.0 "),
            "{par}"
        );
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn rejects_bad_bounds() {
        let t: toml::Table =
            "data = \"s.dat\"\nvariable_column = 4\ntable = \"t\"\nzmin = 5.0\nzmax = 1.0"
                .parse()
                .unwrap();
        assert!(Params::from_table(&t).is_err(), "zmin > zmax");
    }
}
