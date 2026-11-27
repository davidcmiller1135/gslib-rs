//! The crate's error type.

use std::path::PathBuf;

/// Everything that can go wrong in `gslib_rs`.
///
/// Library functions never panic on bad input; they return one of these.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Reading or writing a file failed.
    #[error("{path}: {source}")]
    Io {
        /// The file involved.
        path: PathBuf,
        /// The underlying I/O error.
        #[source]
        source: std::io::Error,
    },

    /// A GSLIB file could not be parsed.
    #[error("line {line}: {message}")]
    Parse {
        /// 1-based line number.
        line: usize,
        /// What was wrong.
        message: String,
    },

    /// An argument was outside its valid domain (e.g. a negative range).
    #[error("invalid {what}: {reason}")]
    InvalidInput {
        /// Which argument.
        what: &'static str,
        /// Why it was rejected.
        reason: String,
    },

    /// Two inputs that must correspond element-by-element do not.
    #[error("{what}: lengths differ ({left} vs {right})")]
    LengthMismatch {
        /// Which inputs.
        what: &'static str,
        /// Length of the first input.
        left: usize,
        /// Length of the second input.
        right: usize,
    },

    /// A kriging system could not be solved, usually because two data share a
    /// location or the model is not positive definite for this configuration.
    #[error("{system} kriging system is singular or not positive definite")]
    SingularSystem {
        /// Which kind of system (e.g. "simple", "ordinary").
        system: &'static str,
    },

    /// The operation needs a covariance, but the model has no sill (power model).
    #[error("model has no finite sill, so it has no covariance; use ordinary kriging")]
    NoCovariance,

    /// Rendering a plot failed.
    #[error("plotting failed: {0}")]
    Plot(String),
}

/// `Result` specialised to [`Error`].
pub type Result<T, E = Error> = std::result::Result<T, E>;

impl Error {
    /// Shorthand for [`Error::InvalidInput`].
    pub(crate) fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        Self::InvalidInput {
            what,
            reason: reason.into(),
        }
    }
}

/// Fail with [`Error::LengthMismatch`] unless `left == right`.
pub(crate) const fn check_lengths(what: &'static str, left: usize, right: usize) -> Result<()> {
    if left == right {
        Ok(())
    } else {
        Err(Error::LengthMismatch { what, left, right })
    }
}
