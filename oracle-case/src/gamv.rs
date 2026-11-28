//! `gamv` cases: typed parameters and the GSLIB parameter file.
//!
//! Layout from `vendor/gslib-2.0/src/GAMV/GAMV.FOR` (`readparm`): data file;
//! `ixl iyl izl`; `nvar ivar...`; `tmin tmax`; output file; `nlag`; `xlag`;
//! `xltol`; `ndir`; per direction `azm atol bandwh dip dtol bandwd`;
//! `isill`; `nvarg`; per variogram `ivtail ivhead ivtype`. Limits from
//! `GAMV.INC`: `MAXLAG = 200`.
//!
//! By default a case is one omnidirectional semivariogram of one variable:
//! one direction with `atol = dtol = 90` and bandwidths of `1e21`, which
//! `gamv` treats as omnidirectional (`atol >= 90`). Directional and
//! multi-variable cases (`docs/spec/gamv.md`, plan step 6) give
//! `variable_columns`, `directions`, `variograms` and `standardize`; the
//! defaults write exactly the parameter file of the original cases. Limits
//! from `GAMV.INC`: `MAXVAR = 5` (variables plus indicator cut-offs),
//! `MAXDIR = 5`, `MXVARG = 5`.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// The output file the recipe asks `gamv` to write.
pub const OUTPUT_FILE: &str = "gamv.out";
/// Most lags `gamv` accepts (`MAXLAG`).
pub const MAX_LAGS: u32 = 200;
/// Most variables plus indicator cut-offs (`MAXVAR`).
pub const MAX_VARIABLES: usize = 5;
/// Most directions (`MAXDIR`).
pub const MAX_DIRECTIONS: usize = 5;
/// Most variograms (`MXVARG`).
pub const MAX_VARIOGRAMS: usize = 5;

/// One variogram line: `ivtail ivhead ivtype [cut]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Variogram {
    /// Tail variable: 1-based index into the case's variables.
    pub tail: u32,
    /// Head variable: 1-based index into the case's variables.
    pub head: u32,
    /// GSLIB's `ivtype`: 1 to 10, or -3.
    pub kind: i32,
    /// The cut-off of an indicator variogram (types 9 and 10 only).
    #[serde(default)]
    pub cutoff: Option<f64>,
}

/// Parameters of a `gamv` case (`[parameters]` in `case.toml`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Params {
    /// Input data file; must be listed in the case's `inputs`.
    pub data: String,
    /// 1-based columns of x, y and z (0 means "not used": coordinate 0).
    pub coordinate_columns: [u32; 3],
    /// 1-based column of the variable (one variable). Exactly one of this and
    /// `variable_columns` is given.
    #[serde(default)]
    pub variable_column: Option<u32>,
    /// 1-based columns of several variables.
    #[serde(default)]
    pub variable_columns: Option<Vec<u32>>,
    /// Directions `[azm, atol, bandh, dip, dtol, bandv]`. Absent: one
    /// omnidirectional direction.
    #[serde(default)]
    pub directions: Option<Vec<[f64; 6]>>,
    /// Variograms. Absent: the semivariogram of variable 1.
    #[serde(default)]
    pub variograms: Option<Vec<Variogram>>,
    /// `isill = 1`: standardise sills.
    #[serde(default)]
    pub standardize: bool,
    /// Number of lags (`nlag`).
    pub nlag: u32,
    /// Lag spacing (`xlag`).
    pub lag_spacing: f64,
    /// Lag tolerance (`xltol`). Must be positive: GSLIB would silently
    /// replace a non-positive value with `xlag / 2`.
    pub lag_tolerance: f64,
    /// `[tmin, tmax]`. With one variable, GSLIB keeps a record iff
    /// `tmin <= v < tmax` (`GAMV.FOR`, `readdata`); its closed pair test
    /// can then reject nothing more (F-011). Absent: no trimming.
    #[serde(default)]
    pub trimming_limits: Option<[f64; 2]>,
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
        let columns = p.columns()?;
        if columns.contains(&0) {
            return Err("variable columns are 1-based".to_owned());
        }
        p.check_directions_and_variograms(columns.len())?;
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
        crate::check_trimming(p.trimming_limits)?;
        if p.data.trim().is_empty() || p.data.len() > 40 {
            return Err("data must be a file name of 1-40 characters (GSLIB reads a40)".to_owned());
        }
        Ok(p)
    }

    /// The variable columns, from `variable_column` or `variable_columns`.
    ///
    /// # Errors
    /// A message unless exactly one of the two is given, with 1 to 5 columns.
    pub fn columns(&self) -> Result<Vec<u32>, String> {
        match (&self.variable_column, &self.variable_columns) {
            (Some(c), None) => Ok(vec![*c]),
            (None, Some(cs)) if (1..=MAX_VARIABLES).contains(&cs.len()) => Ok(cs.clone()),
            (None, Some(_)) => Err(format!("1 to {MAX_VARIABLES} variable_columns")),
            _ => Err("give exactly one of variable_column and variable_columns".to_owned()),
        }
    }

    /// The variograms, or the default semivariogram of variable 1.
    #[must_use]
    pub fn variogram_list(&self) -> Vec<Variogram> {
        self.variograms.clone().unwrap_or_else(|| {
            vec![Variogram {
                tail: 1,
                head: 1,
                kind: 1,
                cutoff: None,
            }]
        })
    }

    fn check_directions_and_variograms(&self, nvar: usize) -> Result<(), String> {
        if let Some(dirs) = &self.directions {
            if !(1..=MAX_DIRECTIONS).contains(&dirs.len()) {
                return Err(format!("1 to {MAX_DIRECTIONS} directions"));
            }
            for d in dirs {
                if d.iter().any(|v| !v.is_finite()) || d[2] < 0.0 || d[5] < 0.0 {
                    return Err(format!("direction {d:?}: finite, bandwidths >= 0"));
                }
            }
        }
        let varios = self.variogram_list();
        if !(1..=MAX_VARIOGRAMS).contains(&varios.len()) {
            return Err(format!("1 to {MAX_VARIOGRAMS} variograms"));
        }
        let mut cuts = 0;
        for v in &varios {
            let in_range = |i: u32| usize::try_from(i).is_ok_and(|i| (1..=nvar).contains(&i));
            if !(in_range(v.tail) && in_range(v.head)) {
                return Err(format!(
                    "variogram {v:?}: tail and head must be 1 to {nvar}"
                ));
            }
            let indicator = matches!(v.kind, 9 | 10);
            if !(matches!(v.kind, 1..=10 | -3)) {
                return Err(format!("variogram {v:?}: kind must be 1 to 10 or -3"));
            }
            if indicator != v.cutoff.is_some() {
                return Err(format!(
                    "variogram {v:?}: a cutoff exactly for kinds 9 and 10"
                ));
            }
            if v.cutoff.is_some_and(|c| !c.is_finite()) {
                return Err(format!("variogram {v:?}: finite cutoff"));
            }
            if indicator {
                cuts += 1;
                // GSLIB stops: 'tmin interferes with indicators!' and 'tmax ...'.
                let [tmin, tmax] = self.trimming_limits.unwrap_or([-1.0e21, 1.0e21]);
                if tmin > 0.0 || tmax <= 1.0 {
                    return Err("indicators need tmin <= 0 and tmax > 1".to_owned());
                }
            }
        }
        if nvar + cuts > MAX_VARIABLES {
            return Err(format!(
                "{nvar} variables and {cuts} cut-offs exceed MAXVAR = {MAX_VARIABLES}"
            ));
        }
        Ok(())
    }

    /// The GSLIB parameter file (layout from `readparm`).
    #[must_use]
    pub fn to_par(&self) -> String {
        let [x, y, z] = self.coordinate_columns;
        let columns = self.columns().unwrap_or_default();
        let column_list: Vec<String> = columns.iter().map(u32::to_string).collect();
        let directions = self.directions.as_ref().map_or_else(
            || "1   \\number of directions\n0.0 90.0 1.0e21 0.0 90.0 1.0e21   \\azm,atol,bandh,dip,dtol,bandv\n".to_owned(),
            |dirs| {
                let mut text = format!("{}   \\number of directions\n", dirs.len());
                for d in dirs {
                    let fields: Vec<String> = d.iter().map(|v| format!("{v:?}")).collect();
                    // Writing to a String cannot fail.
                    let _ = writeln!(text, "{}   \\azm,atol,bandh,dip,dtol,bandv", fields.join(" "));
                }
                text
            },
        );
        let varios = self.variogram_list();
        let mut vario_lines = format!("{}   \\number of variograms\n", varios.len());
        for v in &varios {
            let cut = v.cutoff.map(|c| format!(" {c:?}")).unwrap_or_default();
            // Writing to a String cannot fail.
            let _ = writeln!(
                vario_lines,
                "{} {} {}{cut}   \\tail var., head var., variogram type",
                v.tail, v.head, v.kind
            );
        }
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
             {} {}   \\number of variables, column numbers\n\
             {}   \\trimming limits\n\
             {OUTPUT_FILE}   \\file for variogram output\n\
             {}   \\number of lags\n\
             {}   \\lag separation distance\n\
             {}   \\lag tolerance\n\
             {directions}\
             {}   \\standardize sills? (0=no, 1=yes)\n\
             {vario_lines}",
            self.data,
            columns.len(),
            column_list.join(" "),
            crate::trimming_line(self.trimming_limits),
            self.nlag,
            self.lag_spacing,
            self.lag_tolerance,
            u8::from(self.standardize),
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
        assert!(
            bad(
                "variable_column = 4",
                "variable_columns = [4, 5]\nvariable_column = 4"
            ),
            "both column forms"
        );
        assert!(
            bad(
                "lag_tolerance = 2.5",
                "lag_tolerance = 2.5\ntrimming_limits = [3.0, 3.0]"
            ),
            "empty trimming interval"
        );
    }

    #[test]
    fn directional_multi_variable_par() {
        let text = PARAMS.replace("variable_column = 4", "variable_columns = [4, 5]")
            + "\nstandardize = true\ndirections = [[0.0, 22.5, 25.0, -10.0, 22.5, 25.0], [90.0, 0.0, 5.0, 0.0, 10.0, 5.0]]\n\
               variograms = [{ tail = 1, head = 2, kind = 3 }, { tail = 2, head = 2, kind = 9, cutoff = 1.5 }]";
        let par = Params::from_table(&text.parse().unwrap()).unwrap().to_par();
        let lines: Vec<&str> = par.lines().collect();
        assert!(lines[6].starts_with("2 4 5 "), "{par}");
        assert!(lines[12].starts_with("2 "), "{par}");
        assert!(
            lines[13].starts_with("0.0 22.5 25.0 -10.0 22.5 25.0 "),
            "{par}"
        );
        assert!(lines[14].starts_with("90.0 0.0 5.0 0.0 10.0 5.0 "), "{par}");
        assert!(lines[15].starts_with("1 "), "isill: {par}");
        assert!(lines[16].starts_with("2 "), "{par}");
        assert!(lines[17].starts_with("1 2 3 "), "{par}");
        assert!(lines[18].starts_with("2 2 9 1.5 "), "{par}");
    }

    /// Every committed gamv case still gets exactly its committed parameter
    /// file, so the extension changes no existing fixture.
    #[test]
    fn existing_cases_keep_their_parameter_files() {
        let root =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/gslib/gamv");
        let mut checked = 0;
        for entry in std::fs::read_dir(&root).unwrap() {
            let dir = entry.unwrap().path();
            let Ok(committed) = std::fs::read_to_string(dir.join("params.par")) else {
                continue;
            };
            let case: toml::Table = std::fs::read_to_string(dir.join("case.toml"))
                .unwrap()
                .parse()
                .unwrap();
            let params = case["parameters"].as_table().unwrap();
            assert_eq!(
                Params::from_table(params).unwrap().to_par(),
                committed,
                "{}",
                dir.display()
            );
            checked += 1;
        }
        assert!(checked >= 9, "{checked} cases");
    }

    #[test]
    fn rejects_bad_variograms() {
        let with = |extra: &str| {
            Params::from_table(&(PARAMS.to_owned() + "\n" + extra).parse().unwrap()).is_err()
        };
        assert!(
            with("variograms = [{ tail = 2, head = 1, kind = 1 }]"),
            "no variable 2"
        );
        assert!(
            with("variograms = [{ tail = 1, head = 1, kind = 11 }]"),
            "kind"
        );
        assert!(
            with("variograms = [{ tail = 1, head = 1, kind = 9 }]"),
            "cutoff"
        );
        assert!(
            with("variograms = [{ tail = 1, head = 1, kind = 1, cutoff = 2.0 }]"),
            "cutoff on a non-indicator"
        );
        assert!(
            with(
                "trimming_limits = [0.5, 9.0]\nvariograms = [{ tail = 1, head = 1, kind = 9, cutoff = 2.0 }]"
            ),
            "tmin > 0 with indicators"
        );
        assert!(with("directions = []"), "no direction");
    }

    #[test]
    fn trimming_limits_reach_the_par_file() {
        let with = PARAMS.to_owned() + "\ntrimming_limits = [0.5, 3.25]";
        let par = Params::from_table(&with.parse().unwrap()).unwrap().to_par();
        assert!(
            par.lines().nth(7).unwrap().starts_with("0.5 3.25 "),
            "{par}"
        );
        let untrimmed = Params::from_table(&PARAMS.parse().unwrap())
            .unwrap()
            .to_par();
        assert!(
            untrimmed
                .lines()
                .nth(7)
                .unwrap()
                .starts_with("-1.0e21 1.0e21 "),
            "{untrimmed}"
        );
    }
}
