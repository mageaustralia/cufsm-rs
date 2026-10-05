//! Models the analysis must refuse with an error rather than panic on.

use cufsm::cfsm::{
    base_column, classify, classify_with, stripmain_constrained, Norm, OSpace, Orth, Spaces,
};
use cufsm::{signature_ss, stripmain, BoundaryCondition, Element, Error, Material, Model, Node};

/// A 100 x 50 x 15 lipped channel, t = 1.5, in uniform compression, with an extra node at
/// `orphan` that no strip touches.
fn channel_with_orphan(orphan: usize) -> Model {
    let pts = [
        (50.0, 15.0),
        (50.0, 0.0),
        (0.0, 0.0),
        (0.0, 100.0),
        (50.0, 100.0),
        (50.0, 85.0),
    ];
    let mut nodes: Vec<Node> = pts.iter().map(|&(x, z)| Node::new(x, z, 1.0)).collect();
    nodes.insert(orphan, Node::new(25.0, 50.0, 1.0));
    let idx = |i: usize| if i >= orphan { i + 1 } else { i };
    Model {
        materials: vec![Material::isotropic(203e3, 0.3)],
        nodes,
        elements: (0..5)
            .map(|i| Element {
                ni: idx(i),
                nj: idx(i + 1),
                t: 1.5,
                mat: 0,
            })
            .collect(),
        constraints: vec![],
        springs: vec![],
    }
}

fn refused<T>(what: &str, orphan: usize, r: Result<T, Error>) {
    let want = format!("node {orphan} belongs to no element");
    match r {
        Err(Error::InvalidModel(msg)) => assert_eq!(msg, want, "{what}"),
        Err(e) => panic!("{what}: wrong error {e}"),
        Ok(_) => panic!("{what}: accepted a node that belongs to no element"),
    }
}

/// A node no strip touches, first, between the others, or last: every entry point that takes a
/// model refuses it by name. cFSM used to take it for a sub-node and index out of bounds.
#[test]
fn a_node_in_no_element_is_refused() {
    let bc = BoundaryCondition::SS;
    let all = Spaces {
        global: true,
        distortional: true,
        local: true,
        other: true,
    };
    let gd = Spaces {
        global: true,
        distortional: true,
        local: false,
        other: false,
    };
    // Modes of the connected section, to classify against the bad one: the model is checked
    // before the modes are read.
    let good = channel_with_orphan(0);
    let good = Model {
        nodes: good.nodes[1..].to_vec(),
        elements: (0..5)
            .map(|i| Element {
                ni: i,
                nj: i + 1,
                t: 1.5,
                mat: 0,
            })
            .collect(),
        ..good
    };
    let res = stripmain(&good, &[100.0], &[vec![1.0]], bc, 2).unwrap();
    for orphan in [0, 3, 6] {
        let m = channel_with_orphan(orphan);
        refused("validate", orphan, m.validate());
        refused(
            "stripmain",
            orphan,
            stripmain(&m, &[100.0], &[vec![1.0]], bc, 2),
        );
        refused("signature_ss", orphan, signature_ss(&m, 1));
        for spaces in [all, gd] {
            refused(
                "stripmain_constrained",
                orphan,
                stripmain_constrained(&m, &[100.0], &[vec![1.0, 2.0]], bc, 2, spaces),
            );
        }
        refused("base_column", orphan, base_column(&m, 100.0, bc, &[1.0]));
        refused(
            "classify",
            orphan,
            classify(&m, &res, bc, Orth::Axial, Norm::Vector),
        );
        refused(
            "classify_with",
            orphan,
            classify_with(&m, &res, bc, Orth::Natural, Norm::None, OSpace::Vector),
        );
    }
}
