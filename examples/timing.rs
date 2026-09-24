//! Wall-clock for a signature curve and a cFSM classification on a meshed lipped channel.
use cufsm::cfsm::{classify, Norm, Orth};
use cufsm::template::{templatecalc, Shape, Template};
use cufsm::{grosprop, signature_ss, stresgen, Actions, BoundaryCondition, Material};
use std::time::Instant;

fn main() {
    for mesh in [8, 16, 32] {
        let mut m = templatecalc(
            &Template::outside(Shape::C, 200.0, 76.0, 15.0, 1.9, 3.0, mesh),
            Material::isotropic(203000.0, 0.3),
        );
        let p = grosprop(&m);
        stresgen(
            &mut m,
            &Actions {
                p: 1.0,
                ..Default::default()
            },
            &p,
            false,
        );
        let t = Instant::now();
        let curve = signature_ss(&m, 10).unwrap();
        let sig = t.elapsed();
        let t = Instant::now();
        let first: Vec<_> = curve.iter().step_by(10).cloned().collect();
        classify(&m, &first, BoundaryCondition::SS, Orth::Axial, Norm::Vector).unwrap();
        let cl = t.elapsed() / first.len() as u32;
        println!("{} nodes ({} DOF): signature curve, 100 lengths x 10 modes: {:.0} ms; classification: {:.1} ms per length", m.nodes.len(), 4 * m.nodes.len(), sig.as_secs_f64() * 1e3, cl.as_secs_f64() * 1e3);
    }
}
