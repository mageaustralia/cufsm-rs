//! cFSM on sections with fewer than two corners: a flat plate (none), and an angle, a T and a
//! cruciform (one). CUFSM cannot classify these: with one corner `Kpp` is singular, and rounding
//! keeps a torsion or bending pattern that is zero at every main node, one base vector too many.
//! pyCUFSM fails the same way. So the checks here are theory, not parity.

use cufsm::cfsm::{
    base_column, base_properties, classify, constr_planar_xz, stripmain_constrained, Norm, Orth,
    Spaces,
};
use cufsm::{grosprop, stripmain, BoundaryCondition, Element, Material, Model, Node};
use std::f64::consts::PI;

const E: f64 = 200_000.0;
const NU: f64 = 0.3;
const SHORT: f64 = 30.0;
const LONG: f64 = 30_000.0;

/// Branches of strips (each a polyline cut into `sub` strips) that share nodes where their
/// points coincide, all in unit compression.
fn section(branches: &[&[(f64, f64)]], sub: usize) -> Model {
    let mut nodes: Vec<Node> = vec![];
    let mut elements = vec![];
    let at = |nodes: &mut Vec<Node>, (x, z): (f64, f64)| match nodes
        .iter()
        .position(|n| (n.x - x).abs() < 1e-9 && (n.z - z).abs() < 1e-9)
    {
        Some(i) => i,
        None => {
            nodes.push(Node::new(x, z, 1.0));
            nodes.len() - 1
        }
    };
    for pts in branches {
        for w in pts.windows(2) {
            let ((x0, z0), (x1, z1)) = (w[0], w[1]);
            let mut prev = at(&mut nodes, (x0, z0));
            for k in 1..=sub {
                let f = k as f64 / sub as f64;
                let next = at(&mut nodes, (x0 + (x1 - x0) * f, z0 + (z1 - z0) * f));
                elements.push(Element {
                    ni: prev,
                    nj: next,
                    t: 1.5,
                    mat: 0,
                });
                prev = next;
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

fn one_corner_sections() -> Vec<(&'static str, Model)> {
    vec![
        (
            "unequal angle",
            section(&[&[(0.0, 0.0), (0.0, 150.0), (60.0, 150.0)]], 4),
        ),
        (
            "equal angle",
            section(&[&[(0.0, 80.0), (0.0, 0.0), (80.0, 0.0)]], 4),
        ),
        (
            "T",
            section(
                &[
                    &[(-50.0, 0.0), (0.0, 0.0), (50.0, 0.0)],
                    &[(0.0, 0.0), (0.0, -80.0)],
                ],
                4,
            ),
        ),
        (
            "cruciform",
            section(
                &[
                    &[(-50.0, 0.0), (0.0, 0.0), (50.0, 0.0)],
                    &[(0.0, 50.0), (0.0, 0.0), (0.0, -50.0)],
                ],
                4,
            ),
        ),
    ]
}

/// The minor-axis Euler stress at length `l`.
fn euler_minor(m: &Model, l: f64) -> f64 {
    let p = grosprop(m);
    PI * PI * E * p.i11.min(p.i22) / (l * l) / p.a
}

fn g_only() -> Spaces {
    Spaces {
        global: true,
        ..Default::default()
    }
}

/// A corner translation moves the rest of a one-corner section rigidly: the edge and sub nodes
/// by the same amount, with no rotation. The rigid turn about the corner, which costs no
/// transverse energy, is left out.
#[test]
fn a_corner_translation_moves_a_one_corner_section_rigidly() {
    for (name, m) in one_corner_sections() {
        let bp = base_properties(&m);
        assert_eq!(bp.ncno, 1, "{name}");
        let rp = constr_planar_xz(&m, &bp, 1.0, 1000.0, BoundaryCondition::SS).unwrap();
        let (neno, nmno, nsno) = (bp.nmno - bp.ncno, bp.nmno, bp.nsno);
        // p block: edge x, edge z, main θ, sub x, sub z, sub θ; columns: corner x, corner z.
        let blocks = [
            (0, neno),
            (neno, neno),
            (2 * neno, nmno),
            (2 * neno + nmno, nsno),
            (2 * neno + nmno + nsno, nsno),
            (2 * neno + nmno + 2 * nsno, nsno),
        ];
        for (col, moves) in [
            (0, [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
            (1, [0.0, 1.0, 0.0, 0.0, 1.0, 0.0]),
        ] {
            for ((from, len), want) in blocks.iter().zip(moves) {
                for i in *from..from + len {
                    assert!(
                        (rp.get(i, col) - want).abs() < 1e-9,
                        "{name}: row {i}, column {col}: {}",
                        rp.get(i, col)
                    );
                }
            }
        }
    }
}

/// Short, the lowest mode is the legs buckling: local. Long, it is minor-axis flexure: global,
/// at the Euler load. Restricted to G, the load is Euler's times 1 / (1 − ν²), as for any section:
/// cFSM's global space holds the transverse membrane strain at zero.
#[test]
fn one_corner_sections_classify_by_theory() {
    let lens = [SHORT, LONG];
    let m1 = vec![vec![1.0]; 2];
    for (name, m) in one_corner_sections() {
        let (_, ngm, ndm, _) = base_column(&m, LONG, BoundaryCondition::SS, &[1.0]).unwrap();
        assert_eq!(
            (ngm, ndm),
            (3, 0),
            "{name}: axial and two bending; torsion has no warping at the main nodes"
        );
        let r = stripmain(&m, &lens, &m1, BoundaryCondition::SS, 1).unwrap();
        let c = classify(&m, &r, BoundaryCondition::SS, Orth::Axial, Norm::Vector).unwrap();
        for cls in [c[0][0], c[1][0]] {
            assert!(
                (cls.iter().sum::<f64>() - 100.0).abs() < 1e-6,
                "{name}: {cls:?}"
            );
        }
        assert!(c[0][0][2] > 99.0, "{name} short: {:?}", c[0][0]);
        assert!(c[1][0][0] > 95.0, "{name} long: {:?}", c[1][0]);
        let euler = euler_minor(&m, LONG);
        assert!(
            (r[1].load_factors[0] / euler - 1.0).abs() < 0.02,
            "{name}: {} vs Euler {euler}",
            r[1].load_factors[0]
        );
        let g = stripmain_constrained(
            &m,
            &[LONG],
            &[vec![1.0]],
            BoundaryCondition::SS,
            1,
            g_only(),
        )
        .unwrap();
        let want = euler / (1.0 - NU * NU);
        assert!(
            (g[0].load_factors[0] / want - 1.0).abs() < 0.005,
            "{name}: G only {} vs {want}",
            g[0].load_factors[0]
        );
    }
}

/// A flat plate has no corners. Its global space is axial and in-plane bending, so restricted to
/// G it buckles about its major axis. Its out-of-plane buckling, at any length, has no warping at
/// the main nodes and is local by cFSM's definitions.
#[test]
fn a_flat_plate_classifies() {
    let m = section(&[&[(0.0, 0.0), (60.0, 0.0)]], 4);
    let (_, ngm, ndm, _) = base_column(&m, LONG, BoundaryCondition::SS, &[1.0]).unwrap();
    assert_eq!((ngm, ndm), (2, 0));
    let lens = [SHORT, 300.0, LONG];
    let m1 = vec![vec![1.0]; 3];
    let r = stripmain(&m, &lens, &m1, BoundaryCondition::SS, 1).unwrap();
    let c = classify(&m, &r, BoundaryCondition::SS, Orth::Axial, Norm::Vector).unwrap();
    for (l, modes) in lens.iter().zip(&c) {
        assert!(modes[0][2] > 99.9, "{l}: {:?}", modes[0]);
    }
    let p = grosprop(&m);
    let major = PI * PI * E * p.i11.max(p.i22) / (LONG * LONG) / p.a / (1.0 - NU * NU);
    let g = stripmain_constrained(
        &m,
        &[LONG],
        &[vec![1.0]],
        BoundaryCondition::SS,
        1,
        g_only(),
    )
    .unwrap();
    assert!(
        (g[0].load_factors[0] / major - 1.0).abs() < 0.005,
        "G only {} vs {major}",
        g[0].load_factors[0]
    );
}
