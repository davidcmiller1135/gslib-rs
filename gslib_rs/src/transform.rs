// gslib_rs/src/transform.rs

//! Data transformations: normal‐score, back‐transform, coordinate rotation.

use std::cmp::Ordering;
use statrs::distribution::{Normal, ContinuousCDF};


/// Map your raw data into standard‐normal scores via the “rank‐and‐invert” method
/// (Fortran `nscore`).  If `weights` is `None`, uses equal weights.
pub fn normal_score_transform(
    data: &[f64],
    weights: Option<&[f64]>,
) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }

    // 1) assign weights
    let w: Vec<f64> = match weights {
        Some(ws) => {
            assert_eq!(ws.len(), n, "weights length mismatch");
            ws.to_vec()
        }
        None => vec![1.0; n],
    };

    // 2) sort data by value, keeping original indices
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| data[i].partial_cmp(&data[j]).unwrap_or(Ordering::Equal));

    // 3) build sorted weights & accumulate positions
    let sorted_w: Vec<f64> = idx.iter().map(|&i| w[i]).collect();
    let twt: f64 = sorted_w.iter().sum();
    let mut oldcp = 0.0;
    let mut positions = vec![0.0; n];
    for (k, &wi) in sorted_w.iter().enumerate() {
        let cp = oldcp + wi / twt;
        positions[k] = (cp + oldcp) / 2.0;
        oldcp = cp;
    }

    // 4) invert via Normal(0,1).inverse_cdf
    let norm = Normal::new(0.0, 1.0).unwrap();
    let sorted_scores: Vec<f64> =
        positions.iter().map(|&p| norm.inverse_cdf(p)).collect();

    // 5) reorder back to original order
    let mut result = vec![0.0; n];
    for (k, &orig_i) in idx.iter().enumerate() {
        result[orig_i] = sorted_scores[k];
    }
    result
}

/// Simple back‐transform: given your original data‐vector `data` and
/// the normal‐scores you previously computed `scores`, interpolate
/// new normals `z` back to the original data scale.
/// 
/// * `data` and `scores` must be the same length and correspond row‐by‐row.
/// * If `z` lies outside the range of `scores`, returns the nearest endpoint.
pub fn back_transform(
    data: &[f64],
    scores: &[f64],
    z: &[f64],
) -> Vec<f64> {
    assert_eq!(data.len(), scores.len(), "data/scores length mismatch");
    let n = data.len();
    if n == 0 || z.is_empty() {
        return Vec::new();
    }

    // sort the (score, data) pairs by score
    let mut pairs: Vec<(f64, f64)> = scores.iter().cloned().zip(data.iter().cloned()).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    let scores_s: Vec<f64> = pairs.iter().map(|&(s, _)| s).collect();
    let data_s:   Vec<f64> = pairs.iter().map(|&(_, d)| d).collect();
    let last = n - 1;

    z.iter().map(|&zi| {
        if zi <= scores_s[0] {
            data_s[0]
        } else if zi >= scores_s[last] {
            data_s[last]
        } else {
            // find the bracket
            let idx = scores_s.binary_search_by(|x| x.partial_cmp(&zi).unwrap_or(Ordering::Equal));
            match idx {
                Ok(i) => data_s[i],
                Err(insert) => {
                    let j = insert - 1;
                    let t = (zi - scores_s[j]) / (scores_s[j + 1] - scores_s[j]);
                    data_s[j] + t * (data_s[j + 1] - data_s[j])
                }
            }
        }
    }).collect()
}

/// Convert a list of 2D points **from** the rotated system **to** the original system
/// (Fortran `rotcoord` with `ireverse=1`).
pub fn to_original(
    rot_pts: &[(f64, f64)],
    origin: (f64, f64),
    angle_rad: f64,
) -> Vec<(f64, f64)> {
    let (xorig, yorig) = origin;
    let (s, c) = (angle_rad.sin(), angle_rad.cos());
    rot_pts.iter().map(|&(rx, ry)| {
        let xx = xorig - ry * s + rx * c;
        let yy = yorig + ry * c + rx * s;
        (xx, yy)
    }).collect()
}

/// Convert a list of 2D points **from** the original system **to** the rotated system
/// (Fortran `rotcoord` with `ireverse≠1`).
pub fn to_rotated(
    pts: &[(f64, f64)],
    origin: (f64, f64),
    angle_rad: f64,
) -> Vec<(f64, f64)> {
    let (xorig, yorig) = origin;
    let (s, c) = (angle_rad.sin(), angle_rad.cos());
    pts.iter().map(|&(xx, yy)| {
        let dy = yy - yorig;
        let dx = xx - xorig;
        let rx = dy * s + dx * c;
        let ry = dy * c - dx * s;
        (rx, ry)
    }).collect()
}

/// Placeholder for the full cell‐declustering routine (`declus.for`).
pub fn decluster<T>(_data: &[T]) -> Vec<T> {
    unimplemented!("full declus implementation coming later")
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPS: f64 = 1e-8;

    #[test]
    fn test_normal_score_transform_basic() {
        let data = [3.0, 1.0, 2.0];
        let scores = normal_score_transform(&data, None);
        // sorted data = [1,2,3] ⇒ positions = [1/6,1/2,5/6]
        // so scores[ data==2 ] should be ~0
        let idx_mid = data.iter().position(|&v| v == 2.0).unwrap();
        assert!(scores[idx_mid].abs() < EPS);
        // first & last should be symmetric
        let idx_lo = data.iter().position(|&v| v == 1.0).unwrap();
        let idx_hi = data.iter().position(|&v| v == 3.0).unwrap();
        assert!((scores[idx_lo] + scores[idx_hi]).abs() < 1e-6);
        assert!(scores[idx_lo] < 0.0 && scores[idx_hi] > 0.0);
    }

    #[test]
    fn test_back_transform() {
        let data = [10.0, 20.0, 30.0];
        let scores = [-1.0, 0.0, 1.0];
        // back‐transform z=0 ⇒ 20
        let out = back_transform(&data, &scores, &[0.0]);
        assert_eq!(out, vec![20.0]);
        // outside range clamps
        let out2 = back_transform(&data, &scores, &[-2.0, 2.0]);
        assert_eq!(out2, vec![10.0, 30.0]);
    }

    #[test]
    fn test_rotation_roundtrip() {
        use std::f64::consts::PI;
        let pts = vec![(1.0, 0.0), (2.0, 2.0)];
        let origin = (0.5, 0.5);
        let angle = PI / 3.0;
        let rot = to_rotated(&pts, origin, angle);
        let back = to_original(&rot, origin, angle);
        for (orig, got) in pts.iter().zip(back.iter()) {
            assert!((orig.0 - got.0).abs() < EPS);
            assert!((orig.1 - got.1).abs() < EPS);
        }
    }
}
