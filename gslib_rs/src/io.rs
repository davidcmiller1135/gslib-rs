// gslib_rs/src/io.rs

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use ndarray::Array2;

/// In‐memory representation of a GSLIB data file.
pub struct GslibData {
    pub title: String,
    pub var_names: Vec<String>,
    pub data: Array2<f64>, // shape: (n_records, nvar)
}

/// Read a simplified Geo‐EAS (GSLIB) ASCII file.
/// 
/// # Errors
/// Returns an `io::Error` if any I/O or parsing step fails.
pub fn read_gslib<P: AsRef<Path>>(path: P) -> io::Result<GslibData> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();

    // 1) Title
    let title = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "Missing title"))??
        .trim_end()
        .to_string();

    // 2) Number of variables
    let nvar: usize = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "Missing nvar"))??
        .trim()
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // 3) Variable names
    let mut var_names = Vec::with_capacity(nvar);
    for _ in 0..nvar {
        let line = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "Missing var name"))??;
        var_names.push(line.trim_end().to_string());
    }

    // 4) Data lines
    let mut records = Vec::new();
    for line in lines {
        let l = line?;
        let toks: Vec<_> = l
            .split_whitespace()
            .map(|s| s.parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if toks.len() != nvar {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Expected {} values, got {}", nvar, toks.len()),
            ));
        }
        records.push(toks);
    }

    let nrec = records.len();
    let mut arr = Array2::zeros((nrec, nvar));
    for (i, row) in records.into_iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            arr[(i, j)] = v;
        }
    }

    Ok(GslibData { title, var_names, data: arr })
}

/// Write a `GslibData` back to a simplified Geo-EAS file.
pub fn write_gslib<P: AsRef<Path>>(path: P, ds: &GslibData) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "{}", ds.title)?;
    writeln!(file, "{}", ds.var_names.len())?;
    for name in &ds.var_names {
        writeln!(file, "{}", name)?;
    }
    let (nrec, nvar) = ds.data.dim();
    for i in 0..nrec {
        for j in 0..nvar {
            write!(file, "{}{}", ds.data[(i, j)], if j + 1 < nvar { " " } else { "\n" })?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use tempfile::NamedTempFile;

    /// Helper to create a small in-memory GSLIB file.
    fn make_sample() -> (NamedTempFile, GslibData) {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "My Data Title").unwrap();
        writeln!(tmp, "3").unwrap();
        writeln!(tmp, "X Coordinate").unwrap();
        writeln!(tmp, "Y Coordinate").unwrap();
        writeln!(tmp, "Value").unwrap();
        writeln!(tmp, "0.0 1.0 10.0").unwrap();
        writeln!(tmp, "2.0 3.0 20.0").unwrap();
        writeln!(tmp, "4.0 5.0 30.0").unwrap();

        let expected = GslibData {
            title: "My Data Title".into(),
            var_names: vec!["X Coordinate".into(), "Y Coordinate".into(), "Value".into()],
            data: array![[0.0, 1.0, 10.0], [2.0, 3.0, 20.0], [4.0, 5.0, 30.0]],
        };
        (tmp, expected)
    }

    #[test]
    fn test_read_gslib_and_write_roundtrip() {
        let (tmp, expected) = make_sample();
        let path = tmp.path().to_path_buf();

        // read
        let ds = read_gslib(&path).unwrap();
        assert_eq!(ds.title, expected.title);
        assert_eq!(ds.var_names, expected.var_names);
        assert_eq!(ds.data, expected.data);

        // write out to a new file and read again
        let tmp2 = NamedTempFile::new().unwrap();
        write_gslib(tmp2.path(), &ds).unwrap();
        let ds2 = read_gslib(tmp2.path()).unwrap();
        assert_eq!(ds2.title, expected.title);
        assert_eq!(ds2.var_names, expected.var_names);
        assert_eq!(ds2.data, expected.data);
    }
}
