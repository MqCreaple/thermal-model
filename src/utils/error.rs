use std::{error::Error, fmt::{Debug, Display}};

/// Error type for when a NaN is detected.
pub enum NanError {
    Default,
    NegativeSquareRoot,
    NegativeLogarithm,
}

impl Default for NanError {
    fn default() -> Self {
        NanError::Default
    }
}

impl Display for NanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Debug for NanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NanError::Default => write!(f, "NaN detected with unspecified reason."),
            NanError::NegativeSquareRoot => write!(f, "Nan detected when taking square root of negative number."),
            NanError::NegativeLogarithm => write!(f, "Nan detected when taking logarithm of negative number."),
        }
    }
}

impl Error for NanError {}