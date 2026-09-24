//! Reading reference fixtures written by `oracle/` into the crate's own types.
#![allow(dead_code)]

use cufsm::{BoundaryCondition, Constraint, Dof, Element, Material, Model, Node};
use serde_json::Value;

pub fn load(path: &str) -> Vec<Value> {
    let full = format!("{}/tests/fixtures/{path}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&full).unwrap_or_else(|e| panic!("{full}: {e}"));
    match serde_json::from_str::<Value>(&text).expect("fixture JSON") {
        Value::Array(v) => v,
        v => vec![v],
    }
}

/// A number, a vector or `[]` as a vector (Octave's `jsonencode` writes a one-element vector as a
/// bare number, and a one-row matrix as a flat array).
pub fn vec_of(v: &Value) -> Vec<f64> {
    match v {
        Value::Number(n) => vec![n.as_f64().unwrap()],
        Value::Array(a) => a.iter().flat_map(vec_of).collect(),
        Value::Null => vec![],
        _ => panic!("not numeric: {v}"),
    }
}

/// Rows of a matrix; a flat array is one row.
pub fn rows_of(v: &Value) -> Vec<Vec<f64>> {
    match v {
        Value::Array(a) if a.iter().all(|x| x.is_array()) => a.iter().map(vec_of).collect(),
        Value::Array(a) if a.is_empty() => vec![],
        other => vec![vec_of(other)],
    }
}

/// Each entry of a list whose entries may themselves have been written as bare numbers.
pub fn list_of(v: &Value) -> Vec<Vec<f64>> {
    match v {
        Value::Array(a) => a.iter().map(vec_of).collect(),
        other => vec![vec_of(other)],
    }
}

fn dof(code: f64) -> Dof {
    match code as i64 {
        1 => Dof::X,
        2 => Dof::Z,
        3 => Dof::Y,
        4 => Dof::Theta,
        c => panic!("dof code {c}"),
    }
}

/// The model CUFSM analysed, stresses included, from an oracle record.
pub fn model_of(r: &Value) -> Model {
    let e = r["E"].as_f64().unwrap();
    let nu = r["nu"].as_f64().unwrap();
    let nodes = rows_of(&r["node"])
        .into_iter()
        .map(|n| Node {
            x: n[1],
            z: n[2],
            free: [n[3] != 0.0, n[4] != 0.0, n[5] != 0.0, n[6] != 0.0],
            stress: n[7],
        })
        .collect();
    let elements = rows_of(&r["elem"])
        .into_iter()
        .map(|el| Element {
            ni: el[1] as usize - 1,
            nj: el[2] as usize - 1,
            t: el[3],
            mat: 0,
        })
        .collect();
    let constraints = match &r["constraints"] {
        Value::Number(_) | Value::Null => vec![],
        c => rows_of(c)
            .into_iter()
            .filter(|row| row.len() >= 5)
            .map(|c| Constraint {
                node_e: c[0] as usize - 1,
                dof_e: dof(c[1]),
                coeff: c[2],
                node_k: c[3] as usize - 1,
                dof_k: dof(c[4]),
            })
            .collect(),
    };
    Model {
        materials: vec![Material::isotropic(e, nu)],
        nodes,
        elements,
        constraints,
    }
}

pub fn bc_of(r: &Value) -> BoundaryCondition {
    BoundaryCondition::parse(r["bc"].as_str().unwrap()).unwrap()
}

/// Largest entry-wise difference between two matrices, relative to the reference's largest entry.
pub fn rel_diff(got: &cufsm::dense::Mat, want: &[Vec<f64>]) -> f64 {
    assert_eq!(got.n, want.len(), "matrix size");
    let scale = want
        .iter()
        .flatten()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    let mut worst = 0.0_f64;
    for (i, row) in want.iter().enumerate() {
        for (j, w) in row.iter().enumerate() {
            worst = worst.max((got.get(i, j) - w).abs());
        }
    }
    worst / scale
}

/// Modal assurance criterion: 1 for parallel vectors, whatever their sign or scale.
pub fn mac(a: &[f64], b: &[f64]) -> f64 {
    let ab: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let aa: f64 = a.iter().map(|x| x * x).sum();
    let bb: f64 = b.iter().map(|x| x * x).sum();
    ab * ab / (aa * bb)
}

/// The reduced elastic stiffness's condition estimate at one length: (largest / smallest Cholesky
/// pivot)². A global mode at a long length rests on a near-cancellation of membrane terms, so
/// rounding moves its load factor by about eps x this in any double-precision solver.
pub fn cond_estimate(m: &Model, a: f64, bc: BoundaryCondition, m_a: &[f64]) -> f64 {
    use cufsm::analysis::{assemble, constraint_basis, msort, reduce};
    let terms = msort(m_a);
    let (k, _) = assemble(m, a, bc, &terms);
    let k = match constraint_basis(m, terms.len()) {
        Some(rb) => reduce(&k, &rb),
        None => k,
    };
    let l = cufsm::dense::cholesky(&k.symmetrised()).expect("K positive definite");
    let piv: Vec<f64> = (0..l.n).map(|j| l.get(j, j)).collect();
    (piv.iter().cloned().fold(0.0, f64::max) / piv.iter().cloned().fold(f64::INFINITY, f64::min))
        .powi(2)
}

/// Load factor `lf` of a length whose lowest is `lf1` is held to `1e-10` plus `1e-12` times the
/// condition estimate, times `lf / lf1`: a local or distortional mode to the fixed part, a global
/// mode at a long length to its rounding limit. The ratio because the solve fixes `1/λ` to an
/// absolute accuracy set by the largest `1/λ`, so a higher mode carries proportionally more
/// relative error.
pub fn load_factor_tolerance(cond: f64, lf: f64, lf1: f64) -> f64 {
    (1e-10 + 1e-12 * cond) * (lf / lf1).max(1.0)
}
