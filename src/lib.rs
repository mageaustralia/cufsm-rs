//! # cufsm-rs
//!
//! An independent, dependency-free Rust port of [CUFSM](https://www.ce.jhu.edu/cufsm/), the
//! finite strip method for elastic buckling of thin-walled sections by Benjamin W. Schafer and
//! co-workers at Johns Hopkins University. Not affiliated with or endorsed by the CUFSM authors.
//!
//! The port follows CUFSM's MATLAB source function for function, so its results can be
//! checked against CUFSM's to round-off: see `tests/` and `oracle/`.
//!
//! ```
//! use cufsm::{Model, Material, Node, Element, BoundaryCondition, stripmain};
//!
//! // A 100 x 2 plate strip in uniform compression, simply supported on its long edges.
//! let nodes: Vec<Node> = (0..=10)
//!     .map(|i| {
//!         let mut n = Node::new(10.0 * i as f64, 0.0, 1.0);
//!         if i == 0 || i == 10 {
//!             n.free = [true, false, true, true]; // w fixed at the edges
//!         }
//!         n
//!     })
//!     .collect();
//! let elements = (0..10).map(|i| Element { ni: i, nj: i + 1, t: 2.0, mat: 0 }).collect();
//! let model = Model { materials: vec![Material::isotropic(200e3, 0.3)], nodes, elements, constraints: vec![], springs: vec![] };
//! let r = stripmain(&model, &[100.0], &[vec![1.0]], BoundaryCondition::SS, 1).unwrap();
//! // Plate buckling at a square half-wave: k = 4, sigma = 4 pi² E / (12 (1 - nu²)) (t / b)².
//! let k4 = 4.0 * std::f64::consts::PI.powi(2) * 200e3 / (12.0 * (1.0 - 0.09)) * (2.0f64 / 100.0).powi(2);
//! assert!((r[0].load_factors[0] / k4 - 1.0).abs() < 1e-3);
//! ```

// Index loops mirror CUFSM's MATLAB line for line, so a reader can hold the two side by side;
// iterator rewrites would hide the correspondence the parity tests depend on.
#![allow(clippy::needless_range_loop)]

pub mod analysis;
pub mod bc;
pub mod cfsm;
pub mod cutwp;
pub mod dense;
pub mod linalg;
pub mod model;
pub mod section;
pub mod strip;
pub mod template;

pub use analysis::{signature_minima, signature_ss, stripmain, LengthResult, Minimum};
pub use model::{BoundaryCondition, Constraint, Dof, Element, Material, Model, Node, Spring};
pub use section::{grosprop, stresgen, Actions, GrossProperties};
pub use template::{templatecalc, Shape, Template};

/// Why an analysis could not run.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// The model is not one the analysis can be run on; the text says what is wrong.
    InvalidModel(String),
    /// The reduced elastic stiffness is not positive definite: the section has a mechanism (a
    /// DOF with no stiffness) at this length. `dof` is the reduced DOF where it showed.
    NotPositiveDefinite { dof: usize },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidModel(s) => write!(f, "invalid model: {s}"),
            Error::NotPositiveDefinite { dof } => {
                write!(f, "the elastic stiffness is not positive definite (a mechanism at reduced DOF {dof})")
            }
        }
    }
}

impl std::error::Error for Error {}

/// The README's example, compiled and run as a doctest so it stays true.
#[doc = include_str!("../README.md")]
#[cfg(doctest)]
pub struct ReadmeDoctests;
