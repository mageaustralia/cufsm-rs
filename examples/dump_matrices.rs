//! Writes a model's global K and Kg (simply supported, one term) at one length as JSON, for the
//! high-precision check in `oracle/high_precision.py`:
//!
//!   cargo run --example dump_matrices -- octave "<case name>" <length> > km.json
//!   cargo run --example dump_matrices -- matlab <AXIAL|BENDING> <length> > km.json
use cufsm::analysis::assemble;
use cufsm::{BoundaryCondition, Element, Material, Model, Node};
use serde_json::Value;

fn v(x: &Value) -> Vec<f64> {
    match x {
        Value::Number(n) => vec![n.as_f64().unwrap()],
        Value::Array(a) => a.iter().flat_map(v).collect(),
        _ => vec![],
    }
}

fn node(n: &[f64]) -> Node {
    Node {
        x: n[1],
        z: n[2],
        free: [true; 4],
        stress: n[7],
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (source, name, a) = (&args[0], &args[1], args[2].parse::<f64>().unwrap());
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/");
    let m = if source == "matlab" {
        let d: Value = serde_json::from_str(
            &std::fs::read_to_string(format!("{dir}matlab_cufsm566.json")).unwrap(),
        )
        .unwrap();
        let c = d["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["case"] == name.as_str())
            .unwrap();
        let nodes: Vec<Node> = c["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|n| node(&v(n)))
            .collect();
        let t = d["t"].as_f64().unwrap();
        let elements = (0..nodes.len() - 1)
            .map(|i| Element {
                ni: i,
                nj: i + 1,
                t,
                mat: 0,
            })
            .collect();
        let (e, nu, g) = (
            d["E"].as_f64().unwrap(),
            d["nu"].as_f64().unwrap(),
            d["G"].as_f64().unwrap(),
        );
        Model {
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
        }
    } else {
        let all: Vec<Value> = serde_json::from_str(
            &std::fs::read_to_string(format!("{dir}cufsm_octave.json")).unwrap(),
        )
        .unwrap();
        let r = all.iter().find(|r| r["name"] == name.as_str()).unwrap();
        let nodes = r["node"]
            .as_array()
            .unwrap()
            .iter()
            .map(|n| node(&v(n)))
            .collect();
        let elements = r["elem"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                let e = v(e);
                Element {
                    ni: e[1] as usize - 1,
                    nj: e[2] as usize - 1,
                    t: e[3],
                    mat: 0,
                }
            })
            .collect();
        Model {
            materials: vec![Material::isotropic(
                r["E"].as_f64().unwrap(),
                r["nu"].as_f64().unwrap(),
            )],
            nodes,
            elements,
            constraints: vec![],
        }
    };
    let (k, kg) = assemble(&m, a, BoundaryCondition::SS, &[1.0]);
    let rows = |x: &cufsm::dense::Mat| {
        (0..x.n)
            .map(|i| (0..x.n).map(|j| x.get(i, j)).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    };
    println!(
        "{}",
        serde_json::json!({ "K": rows(&k.symmetrised()), "Kg": rows(&kg.symmetrised()) })
    );
}
