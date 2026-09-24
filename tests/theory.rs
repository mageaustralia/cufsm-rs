//! Checks that need no oracle: closed-form buckling, and invariants any correct analysis keeps.
//! They hold whatever CUFSM itself says, so a defect CUFSM and this port shared would still show.

use cufsm::analysis::{assemble, buckling_eigen};
use cufsm::{
    grosprop, stresgen, stripmain, Actions, BoundaryCondition, Element, Material, Model, Node,
};
use std::f64::consts::PI;

const E: f64 = 200_000.0;
const NU: f64 = 0.3;

/// An open polyline of strips, each segment cut into `sub` strips, all nodes free.
fn polyline(points: &[(f64, f64)], t: f64, sub: usize) -> Model {
    let mut pts = vec![points[0]];
    for w in points.windows(2) {
        let ((x0, z0), (x1, z1)) = (w[0], w[1]);
        for k in 1..=sub {
            let f = k as f64 / sub as f64;
            pts.push((x0 + (x1 - x0) * f, z0 + (z1 - z0) * f));
        }
    }
    let nodes = pts
        .iter()
        .map(|&(x, z)| Node::new(x, z, 0.0))
        .collect::<Vec<_>>();
    let elements = (0..nodes.len() - 1)
        .map(|i| Element {
            ni: i,
            nj: i + 1,
            t,
            mat: 0,
        })
        .collect();
    Model {
        materials: vec![Material::isotropic(E, NU)],
        nodes,
        elements,
        constraints: vec![],
        springs: vec![],
    }
}

fn loaded(mut m: Model, actions: Actions) -> Model {
    let p = grosprop(&m);
    stresgen(&mut m, &actions, &p, false);
    m
}

fn lowest(m: &Model, a: f64) -> f64 {
    stripmain(m, &[a], &[vec![1.0]], BoundaryCondition::SS, 1).unwrap()[0].load_factors[0]
}

/// A plate simply supported on its long edges: k = 4 at a square half-wave, and nowhere lower.
#[test]
fn simply_supported_plate_buckles_at_k_4() {
    let (b, t) = (100.0, 2.0);
    let mut m = polyline(&[(0.0, 0.0), (b, 0.0)], t, 16);
    let last = m.nodes.len() - 1;
    m.nodes[0].free[1] = false; // w fixed at both long edges
    m.nodes[last].free[1] = false;
    for n in &mut m.nodes {
        n.stress = 1.0;
    }
    let sigma_k4 = 4.0 * PI * PI * E / (12.0 * (1.0 - NU * NU)) * (t / b).powi(2);
    let at_b = lowest(&m, b);
    assert!(
        (at_b / sigma_k4 - 1.0).abs() < 2e-4,
        "k at a = b: {}",
        4.0 * at_b / sigma_k4
    );
    // The minimum over the half-wavelength is at a = b.
    for a in [0.6 * b, 0.8 * b, 1.25 * b, 1.6 * b] {
        assert!(
            lowest(&m, a) > at_b,
            "a = {a} buckles below the square half-wave"
        );
    }
    // A plate with one long edge free (the outstand): k = 0.425 + (b/a)², approached at long a.
    let mut o = m.clone();
    o.nodes[last].free[1] = true;
    let a = 20.0 * b;
    let k_outstand = 0.425 + (b / a).powi(2);
    let got = lowest(&o, a) / sigma_k4 * 4.0;
    assert!(
        (got / k_outstand - 1.0).abs() < 5e-3,
        "outstand k {got} vs {k_outstand}"
    );
}

/// A doubly symmetric I-section: web cut into 8 strips, each half-flange into `nf`, the flanges
/// branching from the web's end nodes.
fn i_section(nf: usize) -> Model {
    let (h, bf, tw, tf) = (200.0, 100.0, 5.0, 8.0);
    let mut nodes: Vec<Node> = (0..=8)
        .map(|i| Node::new(0.0, h * i as f64 / 8.0, 0.0))
        .collect();
    let mut elements: Vec<Element> = (0..8)
        .map(|i| Element {
            ni: i,
            nj: i + 1,
            t: tw,
            mat: 0,
        })
        .collect();
    for (root, z) in [(8usize, h), (0usize, 0.0)] {
        for side in [-1.0, 1.0] {
            let mut prev = root;
            for k in 1..=nf {
                nodes.push(Node::new(side * bf / 2.0 * k as f64 / nf as f64, z, 0.0));
                let id = nodes.len() - 1;
                elements.push(Element {
                    ni: prev,
                    nj: id,
                    t: tf,
                    mat: 0,
                });
                prev = id;
            }
        }
    }
    Model {
        materials: vec![Material::isotropic(E, NU)],
        nodes,
        elements,
        constraints: vec![],
        springs: vec![],
    }
}

/// A long doubly symmetric I-section in compression buckles about its weak axis at the Euler
/// load. The flanges hang off the web's end nodes, so the model branches.
///
/// A strip's in-plane displacement across its width is linear, so a coarse flange cannot take up
/// the Poisson contraction of in-plane bending and reads slightly stiff (towards E / (1 - nu²)):
/// the excess falls as the flanges are refined, and is below 0.1% at eight strips a side.
#[test]
fn long_i_section_buckles_at_the_euler_load() {
    let mut excess = vec![];
    for nf in [2, 4, 8] {
        let m0 = i_section(nf);
        let p = grosprop(&m0);
        let m = loaded(
            m0,
            Actions {
                p: 1.0,
                ..Default::default()
            },
        );
        let iyy = p.ixx.min(p.izz);
        let a = 12000.0;
        let euler = PI * PI * E * iyy / (a * a);
        excess.push(lowest(&m, a) / euler - 1.0);
    }
    eprintln!("I-section excess over Euler at 2, 4, 8 strips: {excess:?}");
    assert!(
        excess.windows(2).all(|w| w[1] < w[0]),
        "not converging: {excess:?}"
    );
    assert!(excess[2].abs() < 1e-3, "8 strips a side: {excess:?}");
}

/// Every eigenpair returned satisfies K φ = λ Kg φ to round-off: the solve, not just its match
/// to a reference.
#[test]
fn every_eigenpair_satisfies_the_eigenproblem() {
    let base = polyline(
        &[
            (60.0, 15.0),
            (60.0, 0.0),
            (0.0, 0.0),
            (0.0, 150.0),
            (60.0, 150.0),
            (60.0, 135.0),
        ],
        1.5,
        3,
    );
    for actions in [
        Actions {
            p: 1.0,
            ..Default::default()
        },
        Actions {
            mxx: 1e6,
            ..Default::default()
        },
    ] {
        let m = loaded(base.clone(), actions);
        for a in [30.0, 300.0, 3000.0] {
            let (k, kg) = assemble(&m, a, BoundaryCondition::SS, &[1.0]);
            let (lfs, modes) = buckling_eigen(&k, &kg, 8).unwrap();
            let knorm = k.max_abs() * k.n as f64;
            for (lf, phi) in lfs.iter().zip(&modes) {
                let pn = phi.iter().map(|v| v * v).sum::<f64>().sqrt();
                let mut worst: f64 = 0.0;
                for i in 0..k.n {
                    let r: f64 = (0..k.n)
                        .map(|j| (k.get(i, j) - lf * kg.get(i, j)) * phi[j])
                        .sum();
                    worst = worst.max(r.abs());
                }
                assert!(
                    worst <= 1e-12 * knorm * pn,
                    "a = {a}, λ = {lf}: residual {worst:e}"
                );
            }
        }
    }
}

/// Load factors scale with E, inversely with the reference stress, and do not change when the
/// section is mirrored or its nodes renumbered end to end.
#[test]
fn invariants() {
    let pts = [
        (60.0, 15.0),
        (60.0, 0.0),
        (0.0, 0.0),
        (0.0, 150.0),
        (60.0, 150.0),
        (60.0, 135.0),
    ];
    let base = loaded(
        polyline(&pts, 1.5, 3),
        Actions {
            p: 1.0,
            ..Default::default()
        },
    );
    let lengths = [20.0, 120.0, 600.0, 3000.0];
    let m_all = vec![vec![1.0]; lengths.len()];
    let run = |m: &Model| -> Vec<f64> {
        stripmain(m, &lengths, &m_all, BoundaryCondition::SS, 3)
            .unwrap()
            .iter()
            .flat_map(|r| r.load_factors.clone())
            .collect()
    };
    let r0 = run(&base);
    let close = |a: &[f64], b: &[f64], f: f64, what: &str| {
        for (x, y) in a.iter().zip(b) {
            assert!((x / (y * f) - 1.0).abs() < 1e-8, "{what}: {x} vs {}", y * f);
        }
    };

    let mut stiff = base.clone();
    stiff.materials[0] = Material::isotropic(2.0 * E, NU);
    close(&run(&stiff), &r0, 2.0, "E doubled");

    let mut stressed = base.clone();
    for n in &mut stressed.nodes {
        n.stress *= 4.0;
    }
    close(&run(&stressed), &r0, 0.25, "stress x 4");

    let mut mirrored = base.clone();
    for n in &mut mirrored.nodes {
        n.x = -n.x;
    }
    close(&run(&mirrored), &r0, 1.0, "mirrored");

    let n = base.nodes.len();
    let mut reversed = base.clone();
    reversed.nodes.reverse();
    reversed.elements = base
        .elements
        .iter()
        .rev()
        .map(|e| Element {
            ni: n - 1 - e.nj,
            nj: n - 1 - e.ni,
            ..*e
        })
        .collect();
    close(&run(&reversed), &r0, 1.0, "renumbered");

    let mut moved = base.clone();
    for n in &mut moved.nodes {
        n.x += 1234.5;
        n.z -= 678.9;
    }
    close(&run(&moved), &r0, 1.0, "translated");
}

/// Local buckling converges from above as the strips are refined.
#[test]
fn refining_the_mesh_converges_from_above() {
    let pts = [
        (60.0, 15.0),
        (60.0, 0.0),
        (0.0, 0.0),
        (0.0, 150.0),
        (60.0, 150.0),
        (60.0, 135.0),
    ];
    let lfs: Vec<f64> = [1, 2, 4, 8]
        .iter()
        .map(|&s| {
            lowest(
                &loaded(
                    polyline(&pts, 1.5, s),
                    Actions {
                        p: 1.0,
                        ..Default::default()
                    },
                ),
                100.0,
            )
        })
        .collect();
    for w in lfs.windows(2) {
        assert!(
            w[1] <= w[0] * (1.0 + 1e-12),
            "refinement raised the load factor: {lfs:?}"
        );
    }
    // Richardson: the last step changes it by far less than the first.
    assert!(
        (lfs[3] - lfs[2]).abs() < 0.05 * (lfs[1] - lfs[0]).abs(),
        "{lfs:?}"
    );
}

/// A lipped channel's signature curve has a local minimum and then a distortional one, and each
/// refined minimum sits at or below the sampled points either side of it.
#[test]
fn signature_minima_of_a_lipped_channel() {
    use cufsm::template::{templatecalc, Shape, Template};
    use cufsm::{signature_minima, signature_ss};
    let m = templatecalc(
        &Template::outside(Shape::C, 200.0, 76.0, 15.0, 1.9, 3.0, 12),
        Material::isotropic(E, NU),
    );
    let m = loaded(
        m,
        Actions {
            p: 1.0,
            ..Default::default()
        },
    );
    let curve = signature_ss(&m, 1).unwrap();
    let mins = signature_minima(&curve);
    assert!(mins.len() >= 2, "{mins:?}");
    let (local, dist) = (mins[0], mins[1]);
    // Local buckling at a half-wavelength of the order of the web depth; distortional several
    // times longer and higher.
    assert!(local.length > 50.0 && local.length < 250.0, "{local:?}");
    assert!(
        dist.length > 2.0 * local.length && dist.load_factor > local.load_factor,
        "{dist:?}"
    );
    for mn in &mins {
        let at = curve
            .iter()
            .min_by(|a, b| {
                (a.length / mn.length)
                    .ln()
                    .abs()
                    .total_cmp(&(b.length / mn.length).ln().abs())
            })
            .unwrap();
        assert!(
            mn.load_factor <= at.load_factors[0] * (1.0 + 1e-12),
            "{mn:?} above the sample at {}",
            at.length
        );
        assert!(
            mn.load_factor > 0.95 * at.load_factors[0],
            "{mn:?} far below the sample at {}",
            at.length
        );
    }
}

/// Springs stiffen: every load factor with a spring added is at least the one without, and a
/// stiff enough foundation spring on a node's out-of-plane DOF approaches fixing that DOF.
#[test]
fn springs_stiffen_and_a_stiff_one_approaches_a_fixed_dof() {
    use cufsm::Spring;
    let pts = [
        (60.0, 15.0),
        (60.0, 0.0),
        (0.0, 0.0),
        (0.0, 150.0),
        (60.0, 150.0),
        (60.0, 135.0),
    ];
    let base = loaded(
        polyline(&pts, 1.5, 3),
        Actions {
            p: 1.0,
            ..Default::default()
        },
    );
    let web_mid = 9; // the node halfway up the web
    let spring = |kw: f64| Spring {
        ni: web_mid,
        nj: None,
        ku: kw,
        kv: 0.0,
        kw: 0.0,
        kq: 0.0,
        local: false,
        discrete: false,
        ys_fraction: 0.0,
    };
    let lengths = [100.0, 400.0, 1500.0];
    for a in lengths {
        let free = lowest(&base, a);
        let mut soft = base.clone();
        soft.springs.push(spring(1.0));
        let mut stiff = base.clone();
        stiff.springs.push(spring(1e8));
        let mut fixed = base.clone();
        fixed.nodes[web_mid].free[0] = false; // u, the web's out-of-plane DOF here
        let (s, k, f) = (lowest(&soft, a), lowest(&stiff, a), lowest(&fixed, a));
        assert!(
            s >= free * (1.0 - 1e-12) && k >= s * (1.0 - 1e-12),
            "a = {a}: {free} {s} {k}"
        );
        assert!(
            (k / f - 1.0).abs() < 1e-4,
            "a = {a}: stiff spring {k} vs fixed {f}"
        );
    }
}
