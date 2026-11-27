//! `nscore` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/NSCORE/NSCORE.FOR`: data file;
//! `ivr iwt`; `tmin tmax`; `ismooth`; reference-distribution file; `isvr
//! iswt`; output file; transformation-table file.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// Output file with the data plus a normal-score column.
pub const OUTPUT_FILE: &str = "nscore.out";
/// Output transformation table (`value score` per line).
pub const TABLE_FILE: &str = "nscore.trn";

/// Parameters of an `nscore` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Input data file; must be listed in the case's `inputs`.
    pub data: String,
    /// 1-based column of the variable.
    pub variable_column: u32,
    /// 1-based column of the declustering weight, or 0 for equal weights.
    #[serde(default)]
    pub weight_column: u32,
    /// `[tmin, tmax]`: GSLIB keeps `tmin <= v < tmax` (`NSCORE.FOR`) and
    /// writes −999 for trimmed rows. Absent: no trimming.
    #[serde(default)]
    pub trimming_limits: Option<[f64; 2]>,
}

impl Params {
    /// Read and check the parameters of a case.
    ///
    /// # Errors
    /// A message if the table does not fit the schema or a column is 0.
    pub fn from_table(table: &toml::Table) -> Result<Self, String> {
        let p: Self = toml::Value::Table(table.clone())
            .try_into()
            .map_err(|e: toml::de::Error| e.to_string())?;
        if p.variable_column == 0 {
            return Err("variable_column is 1-based".to_owned());
        }
        crate::check_trimming(p.trimming_limits)?;
        if p.data.trim().is_empty() || p.data.len() > 40 {
            return Err("data must be a file name of 1-40 characters (GSLIB reads a40)".to_owned());
        }
        Ok(p)
    }

    /// The GSLIB parameter file. No reference distribution is used.
    #[must_use]
    pub fn to_par(&self) -> String {
        let mut par = String::new();
        // Writing to a String cannot fail.
        let _ = write!(
            par,
            "                  Parameters for NSCORE\n\
             \x20                 *********************\n\
             \n\
             START OF PARAMETERS:\n\
             {}   \\file with data\n\
             {} {}   \\columns for variable and weight\n\
             {}   \\trimming limits\n\
             0   \\1=transform according to specified ref. dist.\n\
             none.dat   \\file with reference dist.\n\
             1 2   \\columns for variable and weight\n\
             {OUTPUT_FILE}   \\file for output\n\
             {TABLE_FILE}   \\file for output transformation table\n",
            self.data,
            self.variable_column,
            self.weight_column,
            crate::trimming_line(self.trimming_limits)
        );
        par
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_layout() {
        let t: toml::Table = "data = \"smooth.dat\"\nvariable_column = 4\nweight_column = 5"
            .parse()
            .unwrap();
        let par = Params::from_table(&t).unwrap().to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert_eq!(lines[3], "START OF PARAMETERS:");
        assert!(lines[4].starts_with("smooth.dat "), "{par}");
        assert!(lines[5].starts_with("4 5 "), "{par}");
        assert!(lines[10].starts_with("nscore.out "), "{par}");
        assert!(lines[11].starts_with("nscore.trn "), "{par}");
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn trimming_limits_reach_the_par_file() {
        let t: toml::Table = "data = \"a.dat\"\nvariable_column = 4".parse().unwrap();
        let line = |p: &Params| p.to_par().lines().nth(6).unwrap().to_owned();
        let mut p = Params::from_table(&t).unwrap();
        assert!(line(&p).starts_with("-1.0e21 1.0e21 "), "default");
        p.trimming_limits = Some([-998.0, 1.0e21]);
        assert!(line(&p).starts_with("-998.0 1e21 "), "{}", line(&p));
        for bad in ["[1.0, 1.0]", "[2.0, 1.0]", "[-inf, 1.0]", "[0.0, nan]"] {
            let bad_table: toml::Table =
                format!("data = \"a.dat\"\nvariable_column = 4\ntrimming_limits = {bad}")
                    .parse()
                    .unwrap();
            assert!(Params::from_table(&bad_table).is_err(), "{bad}");
        }
    }

    #[test]
    fn rejects_bad_columns() {
        let t: toml::Table = "data = \"a.dat\"\nvariable_column = 0".parse().unwrap();
        assert!(Params::from_table(&t).is_err(), "column 0");
    }
}
