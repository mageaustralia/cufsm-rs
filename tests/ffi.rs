//! The C-ABI in `src/ffi.rs`: its buffers carry exactly what the Rust API returns.
#![cfg(feature = "ffi")]

use cufsm::cfsm::{stripmain_constrained, Spaces};
use cufsm::ffi::*;
use cufsm::{stripmain, BoundaryCondition, Element, Material, Model, Node};

/// A 100 x 50 x 15 lipped channel, t = 1.5, in uniform compression: flat buffers and the same
/// model through the Rust API.
fn channel() -> (Vec<f64>, Vec<f64>, Vec<f64>, Model) {
    let pts = [
        (50.0, 15.0),
        (50.0, 0.0),
        (0.0, 0.0),
        (0.0, 100.0),
        (50.0, 100.0),
        (50.0, 85.0),
    ];
    let mat = Material::isotropic(203e3, 0.3);
    let mats = vec![mat.ex, mat.ey, mat.vx, mat.vy, mat.g];
    let nodes: Vec<f64> = pts
        .iter()
        .flat_map(|&(x, z)| [x, z, 1.0, 1.0, 1.0, 1.0, 1.0])
        .collect();
    let elems: Vec<f64> = (0..5)
        .flat_map(|i| [i as f64, (i + 1) as f64, 1.5, 0.0])
        .collect();
    let model = Model {
        materials: vec![mat],
        nodes: pts.iter().map(|&(x, z)| Node::new(x, z, 1.0)).collect(),
        elements: (0..5)
            .map(|i| Element {
                ni: i,
                nj: i + 1,
                t: 1.5,
                mat: 0,
            })
            .collect(),
        constraints: vec![],
        springs: vec![],
    };
    (mats, nodes, elems, model)
}

/// The last-error slot is one per process: tests that fail on purpose take turns.
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());
fn serial() -> std::sync::MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

const LENS: [f64; 6] = [10.0, 50.0, 100.0, 300.0, 1000.0, 3000.0];

fn signature(params: [f64; 3], bc: &str, lens: &[f64], out: &mut [f64]) -> isize {
    let (mats, nodes, elems, _) = channel();
    unsafe {
        cufsm_signature(
            params.as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            bc.as_ptr(),
            bc.len(),
            lens.as_ptr(),
            lens.len(),
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    }
}

fn modes(params: [f64; 3], bc: &str, lens: &[f64], out: &mut [f64]) -> isize {
    let (mats, nodes, elems, _) = channel();
    unsafe {
        cufsm_modes(
            params.as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            bc.as_ptr(),
            bc.len(),
            lens.as_ptr(),
            lens.len(),
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    }
}

fn last_error() -> String {
    let (p, n) = (cufsm_last_error_ptr(), cufsm_last_error_len());
    String::from_utf8(unsafe { std::slice::from_raw_parts(p as *const u8, n) }.to_vec()).unwrap()
}

/// Rows of `L, λ, G, D, L` (spaces = 7), each value what the Rust API gives, and the return is
/// the number of values written.
#[test]
fn signature_rows_match_the_rust_api() {
    let _g = serial();
    let (_, _, _, model) = channel();
    let m1 = vec![vec![1.0]; LENS.len()];
    let free = stripmain(&model, &LENS, &m1, BoundaryCondition::SS, 1).unwrap();
    let only = |g, d, l| Spaces {
        global: g,
        distortional: d,
        local: l,
        other: false,
    };
    let cons: Vec<_> = [
        only(true, false, false),
        only(false, true, false),
        only(false, false, true),
    ]
    .into_iter()
    .map(|sp| stripmain_constrained(&model, &LENS, &m1, BoundaryCondition::SS, 1, sp).unwrap())
    .collect();

    let mut out = vec![f64::MAX; 40];
    let n = signature([1.0, 7.0, 1.0], "S-S", &LENS, &mut out);
    assert_eq!(n, 30);
    for (i, row) in out[..30].chunks_exact(5).enumerate() {
        assert_eq!(row[0], LENS[i]);
        assert_eq!(row[1], free[i].load_factors[0]);
        for (c, run) in cons.iter().enumerate() {
            assert_eq!(
                row[2 + c],
                run[i].load_factors[0],
                "length {i}, column {}",
                2 + c
            );
        }
    }
    assert!(
        out[30..].iter().all(|v| *v == f64::MAX),
        "wrote past the rows"
    );

    // Only O (bit 8): one extra column, the others untouched.
    let n = signature([1.0, 8.0, 1.0], "S-S", &LENS, &mut out);
    assert_eq!(n, 18);
}

/// S-S with more than one term would be the curve's minimum over L, L/2, ..., not the signature
/// curve: refused, as are terms out of range for the other boundary conditions.
#[test]
fn terms_are_checked() {
    let _g = serial();
    let mut out = vec![0.0; 64];
    assert!(signature([6.0, 0.0, 1.0], "S-S", &LENS, &mut out) < 0);
    assert!(last_error().contains("terms must be 1"), "{}", last_error());
    assert!(signature([0.0, 0.0, 1.0], "C-C", &LENS, &mut out) < 0);
    assert!(signature([2.5, 0.0, 1.0], "C-C", &LENS, &mut out) < 0);
    // C-C at physical lengths with three terms: the Rust API's answer.
    let (_, _, _, model) = channel();
    let m = vec![vec![1.0, 2.0, 3.0]; LENS.len()];
    let r = stripmain(&model, &LENS, &m, BoundaryCondition::CC, 1).unwrap();
    assert_eq!(signature([3.0, 0.0, 1.0], "C-C", &LENS, &mut out), 12);
    for i in 0..LENS.len() {
        assert_eq!(out[2 * i + 1], r[i].load_factors[0]);
    }
}

/// Every bad input is a negative return and a message, never a read past a buffer; repeated
/// failures replace the message (the old one is freed as it was allocated).
#[test]
fn bad_inputs_fail_cleanly() {
    let _g = serial();
    let mut out = vec![0.0; 64];
    let cases: [([f64; 3], &str, &[f64], &str); 7] = [
        ([1.0, 16.0, 1.0], "S-S", &LENS, "spaces"),
        ([1.0, 1.5, 1.0], "S-S", &LENS, "spaces"),
        ([1.0, 0.0, 0.0], "S-S", &LENS, "neigs"),
        ([1.0, 0.0, 1.0], "X-X", &LENS, "unknown boundary"),
        ([1.0, 0.0, 1.0], "S-S", &[100.0, f64::NAN], "length"),
        ([1.0, 0.0, 1.0], "S-S", &[100.0, -1.0], "length"),
        ([1.0, 0.0, 99.0], "S-S", &LENS, "neigs"),
    ];
    for (params, bc, lens, why) in cases {
        assert_eq!(signature(params, bc, lens, &mut out), -1, "{why}");
        assert!(last_error().contains(why), "{why}: {}", last_error());
        assert_eq!(modes(params, bc, lens, &mut out), -1, "modes: {why}");
        assert!(last_error().contains(why), "modes {why}: {}", last_error());
    }
    // A short output buffer is refused before anything is written.
    let mut short = vec![7.0; 5];
    assert_eq!(signature([1.0, 0.0, 1.0], "S-S", &LENS, &mut short), -1);
    assert!(short.iter().all(|v| *v == 7.0));
}

/// A node's free flag and an element's node index must be exact.
#[test]
fn node_and_element_buffers_are_checked() {
    let _g = serial();
    let (mats, mut nodes, mut elems, _) = channel();
    let params = [1.0, 0.0, 1.0];
    let mut out = vec![0.0; 64];
    let call = |nodes: &[f64], elems: &[f64], out: &mut [f64]| unsafe {
        cufsm_signature(
            params.as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            b"S-S".as_ptr(),
            3,
            LENS.as_ptr(),
            LENS.len(),
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    };
    nodes[3] = 0.5;
    assert_eq!(call(&nodes, &elems, &mut out), -1);
    assert!(last_error().contains("free flag"));
    nodes[3] = 1.0;
    elems[9] = 6.0;
    assert_eq!(call(&nodes, &elems, &mut out), -1);
    assert!(last_error().contains("node_j"), "{}", last_error());
    assert_eq!(call(&nodes[..4], &elems, &mut out), -1);
    assert_eq!(
        unsafe {
            cufsm_signature(
                std::ptr::null(),
                3,
                mats.as_ptr(),
                mats.len(),
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                b"S-S".as_ptr(),
                3,
                LENS.as_ptr(),
                LENS.len(),
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
                out.as_mut_ptr(),
                out.len(),
            )
        },
        -1
    );
}

/// Springs and constraints reach the model: the same answers as the Rust API.
#[test]
fn springs_and_constraints_reach_the_model() {
    let _g = serial();
    let (_, _, _, mut model) = channel();
    model.springs.push(cufsm::model::Spring {
        ni: 0,
        nj: None,
        ku: 1e4,
        kv: 0.0,
        kw: 1e4,
        kq: 0.0,
        local: false,
        discrete: true,
        ys_fraction: 0.0,
    });
    let lens = [100.0];
    let m1 = vec![vec![1.0]];
    let want = stripmain(&model, &lens, &m1, BoundaryCondition::SS, 1).unwrap()[0].load_factors[0];
    let (mats, nodes, elems, _) = channel();
    let params = [1.0, 0.0, 1.0];
    let springs = [0.0, -1.0, 1e4, 0.0, 1e4, 0.0, 0.0, 1.0, 0.0];
    let mut out = vec![0.0; 8];
    let n = unsafe {
        cufsm_signature(
            params.as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            b"S-S".as_ptr(),
            3,
            lens.as_ptr(),
            lens.len(),
            springs.as_ptr(),
            springs.len(),
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert_eq!(n, 2, "{}", last_error());
    assert_eq!(out[1], want);
    // a bad dof code is refused with a message
    let bad = [0.0, 9.0, 1.0, 1.0, 1.0];
    assert_eq!(
        unsafe {
            cufsm_signature(
                params.as_ptr(),
                3,
                mats.as_ptr(),
                mats.len(),
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                b"S-S".as_ptr(),
                3,
                lens.as_ptr(),
                lens.len(),
                std::ptr::null(),
                0,
                bad.as_ptr(),
                bad.len(),
                out.as_mut_ptr(),
                out.len(),
            )
        },
        -1
    );
    assert!(last_error().contains("dof"), "{}", last_error());
}

/// Per length: the lowest mode's class, its terms and its shape, as the Rust API gives them.
#[test]
fn modes_match_the_rust_api() {
    let _g = serial();
    let (_, _, _, model) = channel();
    let lens = [100.0, 1000.0];
    let m1 = vec![vec![1.0]; 2];
    let r = stripmain(&model, &lens, &m1, BoundaryCondition::SS, 1).unwrap();
    let cls = cufsm::cfsm::classify(
        &model,
        &r,
        BoundaryCondition::SS,
        cufsm::cfsm::Orth::Axial,
        cufsm::cfsm::Norm::Vector,
    )
    .unwrap();
    let per = 5 + 1 + 4 * 6;
    let mut out = vec![f64::MAX; 2 * per + 3];
    assert_eq!(
        modes([1.0, 0.0, 1.0], "S-S", &lens, &mut out),
        (2 * per) as isize
    );
    for (i, row) in out[..2 * per].chunks_exact(per).enumerate() {
        assert_eq!(&row[..4], &cls[i][0]);
        assert_eq!((row[4], row[5]), (1.0, 1.0));
        assert_eq!(&row[6..], &r[i].modes[0][..]);
    }
    assert!(out[2 * per..].iter().all(|v| *v == f64::MAX));
    // Local at 100 mm, global at 1000 mm.
    assert!(
        out[2] > 50.0 && out[per] + out[per + 1] > 50.0,
        "{:?} {:?}",
        &out[..4],
        &out[per..per + 4]
    );
}

/// An angle (one corner) and a flat plate (none) classify through the C interface too.
#[test]
fn modes_classify_an_angle_and_a_plate() {
    let _g = serial();
    let params = [1.0, 0.0, 1.0];
    let mat = Material::isotropic(203e3, 0.3);
    let mats = [mat.ex, mat.ey, mat.vx, mat.vy, mat.g];
    let lens = [30.0, 30_000.0];
    for pts in [
        &[(0.0, 150.0), (0.0, 0.0), (60.0, 0.0)][..],
        &[(0.0, 0.0), (60.0, 0.0)][..],
    ] {
        let nodes: Vec<f64> = pts
            .iter()
            .flat_map(|&(x, z)| [x, z, 1.0, 1.0, 1.0, 1.0, 1.0])
            .collect();
        let elems: Vec<f64> = (0..pts.len() - 1)
            .flat_map(|i| [i as f64, (i + 1) as f64, 1.5, 0.0])
            .collect();
        let per = 5 + 1 + 4 * pts.len();
        let mut out = vec![0.0; 2 * per];
        let n = unsafe {
            cufsm_modes(
                params.as_ptr(),
                3,
                mats.as_ptr(),
                mats.len(),
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                b"S-S".as_ptr(),
                3,
                lens.as_ptr(),
                lens.len(),
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
                out.as_mut_ptr(),
                out.len(),
            )
        };
        assert_eq!(
            n,
            (2 * per) as isize,
            "{pts:?}: {}",
            if n < 0 { last_error() } else { String::new() }
        );
        for row in out.chunks_exact(per) {
            assert!(
                (row[..4].iter().sum::<f64>() - 100.0).abs() < 1e-6,
                "{:?}",
                &row[..4]
            );
        }
    }
}

#[test]
fn per_element_thickness_and_orthotropy_reach_the_model() {
    let _g = serial();
    let (_, mut nodes, mut elems, mut model) = channel();
    // web thicker, and an orthotropic second material on the flanges
    let ortho = Material {
        ex: 203e3,
        ey: 150e3,
        vx: 0.3,
        vy: 0.3 * 150e3 / 203e3,
        g: 70e3,
    };
    let iso = model.materials[0];
    let mats = vec![
        iso.ex, iso.ey, iso.vx, iso.vy, iso.g, ortho.ex, ortho.ey, ortho.vx, ortho.vy, ortho.g,
    ];
    elems[2 * 4 + 2] = 2.5; // element 2 (the web): t = 2.5
    elems[4 + 3] = 1.0; // element 1: material 1
    elems[3 * 4 + 3] = 1.0; // element 3: material 1
    model.materials.push(ortho);
    model.elements[2].t = 2.5;
    model.elements[1].mat = 1;
    model.elements[3].mat = 1;
    // pin node 2's rotation and y through the new columns
    nodes[2 * 7 + 4] = 0.0;
    nodes[2 * 7 + 5] = 0.0;
    model.nodes[2].free = [true, true, false, false];
    let lens = [50.0, 300.0];
    let mut out = vec![0.0; 4];
    let n = unsafe {
        cufsm_signature(
            [1.0, 0.0, 1.0].as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            "S-S".as_ptr(),
            3,
            lens.as_ptr(),
            2,
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert_eq!(n, 4, "{}", last_error());
    let want = stripmain(
        &model,
        &lens,
        &[vec![1.0], vec![1.0]],
        BoundaryCondition::SS,
        1,
    )
    .unwrap();
    for (row, r) in out.chunks(2).zip(&want) {
        assert_eq!(row[1], r.load_factors[0]);
    }
}

#[test]
fn bad_model_buffers_fail_cleanly() {
    let _g = serial();
    let (mats, nodes, elems, _) = channel();
    let run = |mats: &[f64], nodes: &[f64], elems: &[f64]| {
        let mut out = vec![0.0; 2];
        let n = unsafe {
            cufsm_signature(
                [1.0, 0.0, 1.0].as_ptr(),
                3,
                mats.as_ptr(),
                mats.len(),
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                "S-S".as_ptr(),
                3,
                [100.0].as_ptr(),
                1,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
                out.as_mut_ptr(),
                out.len(),
            )
        };
        (n, last_error())
    };
    let mut e = elems.clone();
    e[3] = 1.0; // element 0 refers to material 1; there is only material 0
    assert!(
        run(&mats, &nodes, &e).1.contains("material"),
        "mat out of range"
    );
    let mut m = mats.clone();
    m[2] = 1.1;
    m[3] = 1.0; // vx * vy >= 1
    assert!(run(&m, &nodes, &elems).1.contains("vx"), "vx*vy");
    let mut e = elems.clone();
    e[2] = 0.0; // t = 0
    assert!(run(&mats, &nodes, &e).1.contains("thickness"), "t");
    let mut nd = nodes.clone();
    nd[4] = 0.5; // free_y not 0 or 1
    assert!(run(&mats, &nd, &elems).1.contains("free flag"), "flag");
    assert!(
        run(&mats[..4], &nodes, &elems).1.contains("mats length"),
        "mats stride"
    );
    assert!(
        run(&mats, &nodes[..6], &elems).1.contains("nodes length"),
        "nodes stride"
    );
    let mut out = vec![0.0; 2];
    let n = unsafe {
        cufsm_signature(
            [1.0, 0.0, 99.0].as_ptr(),
            3,
            mats.as_ptr(),
            mats.len(),
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            "S-S".as_ptr(),
            3,
            [100.0].as_ptr(),
            1,
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert!(n < 0 && last_error().contains("neigs"));
}

#[test]
fn alloc_round_trips() {
    let _g = serial();
    for n in [0, 1, 7, 8, 9, 4096] {
        let p = cufsm_alloc(n);
        assert!(!p.is_null() && p as usize % 8 == 0);
        unsafe {
            std::ptr::write_bytes(p, 0xAB, n);
            cufsm_dealloc(p, n);
        }
    }
    unsafe { cufsm_dealloc(std::ptr::null_mut(), 8) };
}
