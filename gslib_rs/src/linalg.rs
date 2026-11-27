//! Dense linear solves for kriging systems.

use nalgebra::{DMatrix, DVector};

use crate::error::{Error, Result};

/// Build the symmetric `n × n` matrix `f(items[i], items[j])`.
pub(crate) fn symmetric_matrix<T, F>(items: &[T], f: F) -> DMatrix<f64>
where
    F: Fn(&T, &T) -> f64,
{
    let n = items.len();
    // nalgebra fills column-major: column j outer, row i inner.
    let f = &f;
    DMatrix::from_iterator(
        n,
        n,
        items
            .iter()
            .flat_map(|b| items.iter().map(move |a| f(a, b))),
    )
}

/// Solve `A·x = b` for symmetric positive-definite `A` (Cholesky).
///
/// Simple-kriging covariance matrices must be SPD; failure means duplicate
/// locations or an invalid model, and is reported rather than papered over.
pub(crate) fn solve_spd(
    a: DMatrix<f64>,
    b: &DVector<f64>,
    system: &'static str,
) -> Result<DVector<f64>> {
    let x = a
        .cholesky()
        .ok_or(Error::SingularSystem { system })?
        .solve(b);
    finite(x, system)
}

/// Solve `A·x = b` for general square `A` (LU with partial pivoting). Used for
/// ordinary kriging, whose Lagrange row makes the matrix indefinite.
pub(crate) fn solve_general(
    a: DMatrix<f64>,
    b: &DVector<f64>,
    system: &'static str,
) -> Result<DVector<f64>> {
    let x = a.lu().solve(b).ok_or(Error::SingularSystem { system })?;
    finite(x, system)
}

fn finite(x: DVector<f64>, system: &'static str) -> Result<DVector<f64>> {
    if x.iter().all(|v| v.is_finite()) {
        Ok(x)
    } else {
        Err(Error::SingularSystem { system })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn symmetric_matrix_layout() {
        let m = symmetric_matrix(&[1.0, 2.0, 3.0], |a, b| a * 10.0 + b);
        assert_abs_diff_eq!(m[(0, 2)], 13.0);
        assert_abs_diff_eq!(m[(2, 0)], 31.0);
    }

    #[test]
    fn spd_solve_and_failure() {
        let a = DMatrix::from_row_slice(2, 2, &[4.0, 2.0, 2.0, 3.0]);
        let b = DVector::from_vec(vec![2.0, 1.0]);
        let x = solve_spd(a, &b, "test").unwrap();
        assert_abs_diff_eq!(x[0], 0.5, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 0.0, epsilon = 1e-12);
        let singular = DMatrix::from_element(2, 2, 1.0);
        assert!(
            solve_spd(singular, &b, "test").is_err(),
            "singular SPD must fail"
        );
    }

    #[test]
    fn general_solve_pivots() {
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 1.0, 0.0]);
        let b = DVector::from_vec(vec![1.0, 2.0]);
        let x = solve_general(a, &b, "test").unwrap();
        assert_abs_diff_eq!(x[0], 2.0);
        assert_abs_diff_eq!(x[1], 1.0);
    }
}
