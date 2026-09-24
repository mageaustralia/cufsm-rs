//! The cross-section model: materials, nodes, strips (elements) and constraints.
//!
//! Mirrors CUFSM's `prop`, `node`, `elem` and `constraints` arrays, with 0-based indices in place
//! of CUFSM's 1-based node and material numbers.

/// An orthotropic material, as CUFSM's `prop` row `[matnum Ex Ey vx vy G]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Material {
    pub ex: f64,
    pub ey: f64,
    pub vx: f64,
    pub vy: f64,
    pub g: f64,
}

impl Material {
    /// An isotropic material: `G = E / (2 (1 + nu))`.
    pub fn isotropic(e: f64, nu: f64) -> Self {
        Material {
            ex: e,
            ey: e,
            vx: nu,
            vy: nu,
            g: e / (2.0 * (1.0 + nu)),
        }
    }
}

/// A nodal degree of freedom, numbered as CUFSM numbers them in `constraints`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dof {
    /// `u`, in the plane of the section along global x (CUFSM code 1).
    X = 1,
    /// `w`, in the plane of the section along global z (CUFSM code 2).
    Z = 2,
    /// `v`, along the member (CUFSM code 3).
    Y = 3,
    /// `theta`, the rotation about the member axis (CUFSM code 4).
    Theta = 4,
}

/// A node, as CUFSM's `node` row `[node# x z dofx dofz dofy dofrot stress]`.
///
/// `free` holds CUFSM's four DOF flags in that same column order - x, z, y, rotation - `true`
/// being free (CUFSM `1`) and `false` fixed (CUFSM `0`). `stress` is the reference stress at the
/// node; the buckling load factor multiplies it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Node {
    pub x: f64,
    pub z: f64,
    pub free: [bool; 4],
    pub stress: f64,
}

impl Node {
    /// A node free in all four DOFs.
    pub fn new(x: f64, z: f64, stress: f64) -> Self {
        Node {
            x,
            z,
            free: [true; 4],
            stress,
        }
    }
}

/// A strip between two nodes, as CUFSM's `elem` row `[elem# nodei nodej t matnum]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Element {
    pub ni: usize,
    pub nj: usize,
    pub t: f64,
    pub mat: usize,
}

/// A master-slave constraint, CUFSM's `constraints` row `[node#e dofe coeff node#k dofk]`:
/// `dof_e` at `node_e` is eliminated and set to `coeff` times `dof_k` at `node_k`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Constraint {
    pub node_e: usize,
    pub dof_e: Dof,
    pub coeff: f64,
    pub node_k: usize,
    pub dof_k: Dof,
}

/// A complete section model.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Model {
    pub materials: Vec<Material>,
    pub nodes: Vec<Node>,
    pub elements: Vec<Element>,
    pub constraints: Vec<Constraint>,
}

impl Model {
    /// Checks every index and every size, so the analysis never has to.
    pub fn validate(&self) -> Result<(), crate::Error> {
        use crate::Error::InvalidModel as Bad;
        if self.nodes.is_empty() {
            return Err(Bad("the model has no nodes".into()));
        }
        if self.elements.is_empty() {
            return Err(Bad("the model has no elements".into()));
        }
        let n = self.nodes.len();
        for (i, e) in self.elements.iter().enumerate() {
            if e.ni >= n || e.nj >= n {
                return Err(Bad(format!(
                    "element {i} refers to a node that does not exist"
                )));
            }
            if e.ni == e.nj {
                return Err(Bad(format!("element {i} starts and ends at the same node")));
            }
            if e.mat >= self.materials.len() {
                return Err(Bad(format!(
                    "element {i} refers to a material that does not exist"
                )));
            }
            if e.t.is_nan() || e.t <= 0.0 {
                return Err(Bad(format!(
                    "element {i} has a thickness that is not positive"
                )));
            }
            let (a, b) = (&self.nodes[e.ni], &self.nodes[e.nj]);
            if (b.x - a.x).hypot(b.z - a.z) <= 0.0 {
                return Err(Bad(format!("element {i} has zero width")));
            }
        }
        for (i, c) in self.constraints.iter().enumerate() {
            if c.node_e >= n || c.node_k >= n {
                return Err(Bad(format!(
                    "constraint {i} refers to a node that does not exist"
                )));
            }
        }
        Ok(())
    }
}

/// The boundary condition at the loaded ends of the member.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryCondition {
    /// Simply supported at both ends (`'S-S'`). The signature-curve condition.
    SS,
    /// Clamped at both ends (`'C-C'`).
    CC,
    /// Simply supported at one end, clamped at the other (`'S-C'`).
    SC,
    /// Clamped at one end, free at the other (`'C-F'`).
    CF,
    /// Clamped at one end, guided at the other (`'C-G'`).
    CG,
}

impl BoundaryCondition {
    /// Reads CUFSM's string form, `'S-S'` and the rest (either order for the mixed ones).
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "S-S" => Some(Self::SS),
            "C-C" => Some(Self::CC),
            "S-C" | "C-S" => Some(Self::SC),
            "C-F" | "F-C" => Some(Self::CF),
            "C-G" | "G-C" => Some(Self::CG),
            _ => None,
        }
    }
}
