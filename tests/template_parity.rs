//! The section template against CUFSM's own `templatecalc.m` (and `template_out_to_in.m` for the
//! cases given in outside dimensions): every node of every template case in the fixture, placed
//! where CUFSM placed it.

mod common;
use common::*;
use cufsm::template::{templatecalc, Shape, Template};
use cufsm::Material;

#[test]
fn template_nodes_match_cufsm() {
    let mut cases = 0;
    for r in load("cufsm_octave.json") {
        let Some(tp) = r.get("template") else {
            continue;
        };
        let name = r["name"].as_str().unwrap();
        let f = |k: &str| tp[k].as_f64().unwrap();
        let u = |k: &str| tp[k].as_f64().unwrap() as usize;
        let t = Template {
            shape: if f("CorZ") == 2.0 { Shape::Z } else { Shape::C },
            h: f("h"),
            b1: f("b1"),
            b2: f("b2"),
            d1: f("d1"),
            d2: f("d2"),
            r1: f("r1"),
            r2: f("r2"),
            r3: f("r3"),
            r4: f("r4"),
            q1: f("q1"),
            q2: f("q2"),
            t: f("t"),
            nh: u("nh"),
            nb1: u("nb1"),
            nb2: u("nb2"),
            nd1: u("nd1"),
            nd2: u("nd2"),
            nr1: u("nr1"),
            nr2: u("nr2"),
            nr3: u("nr3"),
            nr4: u("nr4"),
            centerline: f("center") == 1.0,
        };
        let m = templatecalc(&t, Material::isotropic(203000.0, 0.3));
        let want = rows_of(&r["node"]);
        assert_eq!(m.nodes.len(), want.len(), "{name}: node count");
        for (i, (n, w)) in m.nodes.iter().zip(&want).enumerate() {
            assert!(
                (n.x - w[1]).abs() < 1e-12 && (n.z - w[2]).abs() < 1e-12,
                "{name}: node {} at ({}, {}), CUFSM ({}, {})",
                i + 1,
                n.x,
                n.z,
                w[1],
                w[2]
            );
        }
        let el = rows_of(&r["elem"]);
        assert_eq!(m.elements.len(), el.len(), "{name}: element count");
        assert!(
            m.elements
                .iter()
                .zip(&el)
                .all(|(e, w)| e.ni + 1 == w[1] as usize
                    && e.nj + 1 == w[2] as usize
                    && (e.t - w[3]).abs() < 1e-15),
            "{name}: elements"
        );
        cases += 1;
    }
    assert!(cases >= 12, "only {cases} template cases");
}
