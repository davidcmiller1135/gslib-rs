//! Reading `gamv` output (`vendor/gslib-2.0/src/GAMV/GAMV.FOR`, `writeout`).
//!
//! Per variogram and direction: a title line (`Semivariogram tail:... head:...
//! direction  1`), then `nlag + 2` rows `il dis gam np hm tm` (plus `hv tv`
//! for correlograms). Lag 1 is the zero-distance class; lag `il ≥ 2` is
//! centred on `(il − 2)·xlag`. Empty lags print zeros.
//!
//! [`legacy`] builds gslib-rs's request from a case the way GSLIB reads it
//! (`docs/spec/gamv.md`, plan step 6).

pub mod legacy;

/// One row of `gamv` output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    /// Lag index (1-based).
    pub il: u32,
    /// Mean separation.
    pub distance: f64,
    /// The variogram value.
    pub gamma: f64,
    /// Pair count as GSLIB reports it.
    pub pairs: u64,
    /// Mean of head values.
    pub head_mean: f64,
    /// Mean of tail values.
    pub tail_mean: f64,
    /// Correlograms only: variances of the head and tail values.
    pub variances: Option<(f64, f64)>,
}

/// One title line and its rows.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// The title line.
    pub title: String,
    /// The rows in file order.
    pub rows: Vec<Row>,
}

/// Parse a `gamv` output file.
///
/// # Errors
/// A message naming the first line that does not fit the format.
pub fn parse_output(text: &str) -> Result<Vec<Table>, String> {
    let mut tables: Vec<Table> = Vec::new();
    for (n, line) in text.lines().enumerate().map(|(i, l)| (i + 1, l)) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }
        let numeric = tokens.first().is_some_and(|t| t.parse::<u32>().is_ok());
        if !numeric {
            tables.push(Table {
                title: line.trim().to_owned(),
                rows: Vec::new(),
            });
            continue;
        }
        let (head, extra) = tokens.split_at(tokens.len().min(6));
        let Ok([il, dis, gam, np, hm, tm]) = <[&str; 6]>::try_from(head) else {
            return Err(format!(
                "line {n}: expected 6 or 8 columns, got {}",
                tokens.len()
            ));
        };
        let extra = match extra {
            [] => None,
            [hv, tv] => Some((*hv, *tv)),
            _ => {
                return Err(format!(
                    "line {n}: expected 6 or 8 columns, got {}",
                    tokens.len()
                ));
            }
        };
        let real = |t: &str| {
            t.parse::<f64>()
                .map_err(|_| format!("line {n}: {t:?} is not a number"))
        };
        let row = Row {
            il: il.parse().map_err(|_| format!("line {n}: bad lag index"))?,
            distance: real(dis)?,
            gamma: real(gam)?,
            pairs: np
                .parse()
                .map_err(|_| format!("line {n}: bad pair count"))?,
            head_mean: real(hm)?,
            tail_mean: real(tm)?,
            variances: extra
                .map(|(hv, tv)| Ok::<_, String>((real(hv)?, real(tv)?)))
                .transpose()?,
        };
        tables
            .last_mut()
            .ok_or_else(|| format!("line {n}: data before the first title"))?
            .rows
            .push(row);
    }
    Ok(tables)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_table() {
        let text = "Semivariogram  tail:V  head:V  direction  1\n   1  0.0E+000  0.0E+000  8  1.5E+000  1.5E+000\n   2  1.0E+000  2.5E-001  4  1.0E+000  2.0E+000\n";
        let t = parse_output(text).unwrap();
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].rows.len(), 2);
        assert_eq!(t[0].rows[1].pairs, 4);
        assert!((t[0].rows[1].tail_mean - 2.0).abs() < 1e-15, "{t:?}");
    }

    #[test]
    fn rejects_rows_before_a_title() {
        assert!(parse_output("   1 0 0 0 0 0\n").is_err(), "no title");
    }
}
