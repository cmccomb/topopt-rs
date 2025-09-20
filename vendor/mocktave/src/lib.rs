#![deny(warnings)]
#![warn(clippy::all, clippy::pedantic)]

use std::collections::HashMap;
use std::error::Error;
use std::fmt;

/// Error returned when attempting to access values produced by the
/// unavailable external Octave interpreter.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MocktaveError {
    kind: MocktaveErrorKind,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum MocktaveErrorKind {
    MissingValue(String),
}

impl MocktaveError {
    fn missing(name: &str) -> Self {
        Self {
            kind: MocktaveErrorKind::MissingValue(name.to_string()),
        }
    }
}

impl fmt::Display for MocktaveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            MocktaveErrorKind::MissingValue(name) => {
                write!(f, "mocktave data for '{name}' is unavailable")
            }
        }
    }
}

impl Error for MocktaveError {}

/// Placeholder results returned from the stubbed interpreter.
#[derive(Clone, Debug, Default)]
pub struct InterpreterResults {
    matrices: HashMap<String, Vec<Vec<f64>>>,
}

impl InterpreterResults {
    /// Construct an empty result set.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }

    /// Attach a matrix payload to the results. This helper exists to make
    /// it easy for tests that wish to inject deterministic data.
    #[must_use]
    pub fn with_matrix(mut self, name: impl Into<String>, values: Vec<Vec<f64>>) -> Self {
        let _ = self.matrices.insert(name.into(), values);
        self
    }

    /// Retrieve a matrix by name.
    pub fn get_matrix(&self, name: &str) -> Result<Vec<Vec<f64>>, MocktaveError> {
        self.matrices
            .get(name)
            .cloned()
            .ok_or_else(|| MocktaveError::missing(name))
    }
}

/// Evaluate an Octave/MATLAB snippet. In this stubbed implementation no
/// external interpreter is available, so the result set is always empty.
#[must_use]
pub fn eval(_input: &str) -> InterpreterResults {
    InterpreterResults::empty()
}
