//! Parity with MATLAB CUFSM: the official compiled CUFSM v5.66 run under the MATLAB R2025b
//! Runtime (`stripmain.m`, MATLAB's own `eigs`) on a 140 x 45 x 10 x 1.2 lipped channel with
//! rounded corners, in compression and in bending, at 74 half-wavelengths, 10 modes each. The
//! fixture comes from the CufsmSharp project (MIT); see its `about` field.
//!
//! This is the check that does not go through Octave: real MATLAB, real `eigs`.
//!
//! The local modes agree to 1e-14 - 1e-12. The difference grows with length, to about 1e-4 on the
//! global mode at 300 in, because the 2.5 mm corner strips make `K` very ill-conditioned: there,
//! against the same matrices solved in 40-digit arithmetic (`oracle/high_precision.py`), MATLAB is
//! out by -3.1e-5 and this crate by +5.9e-5. Both are at the rounding limit, so the tolerance is
//! the condition-scaled one the Octave comparison uses.

mod common;
use common::*;
use cufsm::{stripmain, BoundaryCondition, Element, Material, Model, Node};

#[test]
fn load_factors_match_matlab_cufsm() {
    let d = &load("matlab_cufsm566.json")[0];
    let e = d["E"].as_f64().unwrap();
    let nu = d["nu"].as_f64().unwrap();
    let g = d["G"].as_f64().unwrap();
    let t = d["t"].as_f64().unwrap();
    let mut worst = (0.0_f64, String::new());
    let mut compared = 0;
    for c in d["cases"].as_array().unwrap() {
        let name = c["case"].as_str().unwrap();
        let nodes: Vec<Node> = rows_of(&c["nodes"])
            .into_iter()
            .map(|n| Node {
                x: n[1],
                z: n[2],
                free: [n[3] != 0.0, n[4] != 0.0, n[5] != 0.0, n[6] != 0.0],
                stress: n[7],
            })
            .collect();
        let elements = (0..nodes.len() - 1)
            .map(|i| Element {
                ni: i,
                nj: i + 1,
                t,
                mat: 0,
            })
            .collect();
        let m = Model {
            materials: vec![Material {
                ex: e,
                ey: e,
                vx: nu,
                vy: nu,
                g,
            }],
            nodes,
            elements,
            constraints: vec![],
        };
        let lengths = vec_of(&c["lengths"]);
        let m_all = vec![vec![1.0]; lengths.len()];
        let got = stripmain(&m, &lengths, &m_all, BoundaryCondition::SS, 10).unwrap();
        let want = rows_of(&c["load_factors"]);
        for (l, res) in got.iter().enumerate() {
            let tol =
                load_factor_tolerance(cond_estimate(&m, lengths[l], BoundaryCondition::SS, &[1.0]));
            for (i, w) in want[l].iter().enumerate().filter(|(_, w)| **w > 0.0) {
                let d = (res.load_factors[i] / w - 1.0).abs();
                if d > worst.0 {
                    worst = (d, format!("{name} length {} mode {}", lengths[l], i + 1));
                }
                assert!(
                    d <= tol,
                    "{name} length {} mode {}: {} vs MATLAB {w} (relative {d:e}, allowed {tol:e})",
                    lengths[l],
                    i + 1,
                    res.load_factors[i]
                );
                compared += 1;
            }
        }
    }
    eprintln!(
        "{compared} load factors compared with MATLAB CUFSM; worst relative difference {:e} ({})",
        worst.0, worst.1
    );
    assert!(compared > 1400, "only {compared} compared");
}
