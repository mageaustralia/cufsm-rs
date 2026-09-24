//! Writes one fixture case's global K and Kg at one length as JSON, for a high-precision check:
//! cargo run --example dump_matrices -- "<case>" <length> > out.json
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

fn main() {
    let mut args = std::env::args().skip(1);
    let name = args.next().unwrap();
    let a: f64 = args.next().unwrap().parse().unwrap();
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/cufsm_octave.json"
    ))
    .unwrap();
    let all: Vec<Value> = serde_json::from_str(&text).unwrap();
    let r = all.iter().find(|r| r["name"] == name.as_str()).unwrap();
    let nodes = r["node"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| {
            let n = v(n);
            Node {
                x: n[1],
                z: n[2],
                free: [true; 4],
                stress: n[7],
            }
        })
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
    let m = Model {
        materials: vec![Material::isotropic(
            r["E"].as_f64().unwrap(),
            r["nu"].as_f64().unwrap(),
        )],
        nodes,
        elements,
        constraints: vec![],
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
