// gslib_rs/src/utils.rs

use lazy_static::lazy_static;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand::distributions::{Distribution, Uniform};
use std::sync::Mutex;
use rand_distr::StandardNormal;

lazy_static! {
    static ref GLOBAL_RNG: Mutex<StdRng> = Mutex::new(StdRng::from_entropy());
}

/// Set the RNG seed (deterministic runs).
pub fn set_seed(seed: u64) {
    let mut guard = GLOBAL_RNG.lock().unwrap();
    *guard = StdRng::seed_from_u64(seed);
}

/// Equivalent to Fortran RNG: returns uniform random in [0,1)
pub fn rand_f() -> f64 {
    let mut guard = GLOBAL_RNG.lock().unwrap();
    // Uniform::new is [low, high)
    let uni = Uniform::new(0.0_f64, 1.0_f64);
    uni.sample(&mut *guard)
}

/// Compute string length as Fortran’s LEN_TRIM: trim trailing spaces.
pub fn strlen_f(s: &str) -> usize {
    s.trim_end_matches(' ').len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rand_f() {
        set_seed(42);
        let a = rand_f();
        let b = rand_f();
        assert!(a > 0.0 && a < 1.0, "a = {}", a);
        assert!(b > 0.0 && b < 1.0, "b = {}", b);

        set_seed(42);
        assert_eq!(rand_f(), a);
        assert_eq!(rand_f(), b);
    }

    #[test]
    fn test_strlen_f() {
        assert_eq!(strlen_f("HELLO   "), 5);
        assert_eq!(strlen_f("   trimmed"), 10);
        assert_eq!(strlen_f("no_spaces"), 9);
    }
}


// gslib_rs/src/utils.rs
// … above: imports, GLOBAL_RNG, set_seed, rand_f, strlen_f …

use std::cmp::Ordering;

/// In-place sort of a single slice (Fortran `sortem.for`).
pub fn sort_em(vals: &mut [f64]) {
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
}

/// In-place parallel sort of two slices by the first key (Fortran `dsortem.for`).
/// Panics if `keys.len() != vals.len()`.
pub fn dsort_em(keys: &mut [f64], vals: &mut [f64]) {
    assert_eq!(keys.len(), vals.len(), "dsort_em: lengths differ");
    // Zip, sort, unzip
    let mut combined: Vec<(f64, f64)> = keys
        .iter()
        .cloned()
        .zip(vals.iter().cloned())
        .collect();
    combined.sort_by(|(k1, _), (k2, _)| k1.partial_cmp(k2).unwrap_or(Ordering::Equal));
    for (i, (k, v)) in combined.into_iter().enumerate() {
        keys[i] = k;
        vals[i] = v;
    }
}

/// Exact-match lookup (Fortran `getindx.for`).
/// Returns `Some(i)` if `sorted[i] == target`; otherwise `None`.
pub fn get_index(sorted: &[f64], target: f64) -> Option<usize> {
    match sorted.binary_search_by(|x| x.partial_cmp(&target).unwrap_or(Ordering::Equal)) {
        Ok(i) => Some(i),
        Err(_) => None,
    }
}

/// Bracket-style locate (Fortran `locate.for`):
/// returns the largest index `i` such that `sorted[i] <= target`.
/// If `target < sorted[0]`, returns `0`; if `target >= sorted[last]`, returns `sorted.len() - 1`.
pub fn locate(sorted: &[f64], target: f64) -> usize {
    if sorted.is_empty() {
        panic!("locate: empty input slice");
    }
    match sorted.binary_search_by(|x| x.partial_cmp(&target).unwrap_or(Ordering::Equal)) {
        Ok(i) => i,
        Err(0) => 0,
        Err(n) if n >= sorted.len() => sorted.len() - 1,
        Err(n) => n - 1,
    }
}

#[cfg(test)]
mod sort_tests {
    use super::*;

    #[test]
    fn test_sort_em() {
        let mut a = [3.1, 1.2, 2.3, 2.3];
        sort_em(&mut a);
        assert_eq!(a, [1.2, 2.3, 2.3, 3.1]);
    }

    #[test]
    fn test_dsort_em() {
        let mut keys = [3.0, 1.0, 2.0];
        let mut vals = [30.0, 10.0, 20.0];
        dsort_em(&mut keys, &mut vals);
        assert_eq!(keys, [1.0, 2.0, 3.0]);
        assert_eq!(vals, [10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_get_index() {
        let data = [0.0, 1.5, 3.0, 4.5];
        assert_eq!(get_index(&data, 0.0), Some(0));
        assert_eq!(get_index(&data, 3.0), Some(2));
        assert_eq!(get_index(&data, 2.0), None);
        assert_eq!(get_index(&data, 5.0), None);
    }

    #[test]
    fn test_locate() {
        let data = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(locate(&data, 0.5), 0);
        assert_eq!(locate(&data, 1.0), 0);
        assert_eq!(locate(&data, 1.5), 0);
        assert_eq!(locate(&data, 2.0), 1);
        assert_eq!(locate(&data, 3.7), 2);
        assert_eq!(locate(&data, 4.0), 3);
        assert_eq!(locate(&data, 5.0), 3);
    }
}

// gslib_rs/src/utils.rs
// … above: imports, GLOBAL_RNG, set_seed, rand_f, strlen_f, sort_em, dsort_em, get_index, locate …

use nalgebra::{DMatrix, DVector};

/// Raises `x` to integer power `p` (Fortran `powint`).
pub fn pow_int(x: f64, p: i32) -> f64 {
    x.powi(p)
}

/// Raises `x` to integer power `p` (double‐precision variant, Fortran `dpowint`).
pub fn d_pow_int(x: f64, p: i32) -> f64 {
    x.powi(p)
}

/// Solve `A·x = b` via LU elimination **without** pivoting (Fortran `ksol`).
/// ksol: “A standard Gaussian elimination algorithm without pivot search … adequate for SK and OK applications.” :contentReference[oaicite:0]{index=0}
pub fn solve_kriging_system(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    // nalgebra’s .lu() uses partial pivoting by default, but we assume
    // the matrix ordering (OK/SK) makes pivot search unnecessary.
    a.clone().lu().solve(b)
}

/// Solve `A·x = b` via LU elimination **with** partial pivoting (Fortran `ktsol`).
/// ktsol: “A Gaussian elimination algorithm with partial pivoting is necessary for KT.” :contentReference[oaicite:1]{index=1}
pub fn solve_kriging_system_pivot(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    a.clone().lu().solve(b)
}

/// Draw a sample from N(0,1), mirroring Fortran `randn` (GAUINV after RAND).
pub fn randn() -> f64 {
    let mut guard = GLOBAL_RNG.lock().unwrap();
    // Sample() returns an f64 directly
    StandardNormal.sample(&mut *guard)
}

/// Check Fortran‐style name (chknam): 1st letter, <=8 chars, alphanumeric only.
pub fn check_name(name: &str) -> bool {
    let s = name.trim();
    if s.len() < 1 || s.len() > 8 { return false; }
    let mut chars = s.chars();
    if !chars.next().unwrap().is_ascii_alphabetic() { return false; }
    chars.all(|c| c.is_ascii_alphanumeric())
}

/// Format a title string (chktitle): trim trailing spaces, optionally enforce max length.
pub fn format_title(title: &str) -> String {
    let mut s = title.trim_end().to_string();
    if s.len() > 80 {
        s.truncate(80);
    }
    s
}

/// Convert integer to ordinal text (numtext), e.g. 1→"1st", 2→"2nd", 11→"11th", etc.
pub fn num_text(n: i32) -> String {
    let suffix = match n % 100 {
        11 | 12 | 13 => "th",
        _ => match n % 10 {
            1 => "st",
            2 => "nd",
            3 => "rd",
            _ => "th",
        }
    };
    format!("{}{}", n, suffix)
}


#[cfg(test)]
mod linear_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_pow_int() {
        assert_eq!(pow_int(2.0, 3), 8.0);
        assert_eq!(pow_int(5.0, 0), 1.0);
        assert_eq!(pow_int(2.0, -2), 0.25);
    }

    #[test]
    fn test_d_pow_int() {
        assert_eq!(d_pow_int(3.0, 4), 81.0);
        assert_eq!(d_pow_int(2.0, -1), 0.5);
    }

    #[test]
    fn test_solve_kriging_system_identity() {
        let a = DMatrix::identity(3, 3);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let x = solve_kriging_system(&a, &b).unwrap();
        assert_eq!(x, b);
    }

    #[test]
    fn test_solve_kriging_system_pivot_simple() {
        // a simple 2×2 system
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 1.0, 0.0]);
        let b = DVector::from_vec(vec![1.0, 2.0]);
        let x = solve_kriging_system_pivot(&a, &b).unwrap();
        // should solve [0 1; 1 0]·x = [1;2] ⇒ x = [2;1]
        assert_eq!(x, DVector::from_vec(vec![2.0, 1.0]));
    }
    #[test]
    fn test_randn_repeatability() {
        set_seed(123);
        let a = randn();
        let b = randn();
        assert!(a.is_finite() && b.is_finite());

        set_seed(123);
        assert_eq!(randn(), a);
        assert_eq!(randn(), b);
    }

    #[test]
    fn test_check_name() {
        assert!(check_name("VAR1"));
        assert!(check_name("Abc123"));
        assert!(!check_name("1BAD"));      // starts with digit
        assert!(!check_name("TOO_LONGNAME"));
        assert!(!check_name("INV@LID"));
    }

    #[test]
    fn test_format_title() {
        let long = "A".repeat(90);
        let formatted = format_title(&long);
        assert_eq!(formatted.len(), 80);
        assert_eq!(format_title(" Hello  "), " Hello");
    }

    #[test]
    fn test_num_text() {
        assert_eq!(num_text(1), "1st");
        assert_eq!(num_text(2), "2nd");
        assert_eq!(num_text(3), "3rd");
        assert_eq!(num_text(4), "4th");
        assert_eq!(num_text(11), "11th");
        assert_eq!(num_text(112), "112th");
        assert_eq!(num_text(123), "123rd");
    }
}




