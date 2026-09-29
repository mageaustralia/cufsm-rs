//! The C-ABI in `src/ffi.rs`: its buffers carry exactly what the Rust API returns.
#![cfg(feature = "ffi")]

use cufsm::cfsm::{stripmain_constrained, Spaces};
use cufsm::ffi::*;
use cufsm::{stripmain, BoundaryCondition, Element, Material, Model, Node};

/// A 100 x 50 x 15 lipped channel, t = 1.5, in uniform compression: flat buffers and the same
/// model through the Rust API.
fn channel() -> (Vec<f64>, Vec<f64>, Model) {
    let pts = [
        (50.0, 15.0),
        (50.0, 0.0),
        (0.0, 0.0),
        (0.0, 100.0),
        (50.0, 100.0),
        (50.0, 85.0),
    ];
    let nodes: Vec<f64> = pts
        .iter()
        .flat_map(|&(x, z)| [x, z, 1.0, 1.0, 1.0])
        .collect();
    let elems: Vec<f64> = (0..5).flat_map(|i| [i as f64, (i + 1) as f64]).collect();
    let model = Model {
        materials: vec![Material::isotropic(203e3, 0.3)],
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
    (nodes, elems, model)
}

/// The last-error slot is one per process: tests that fail on purpose take turns.
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());
fn serial() -> std::sync::MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

const LENS: [f64; 6] = [10.0, 50.0, 100.0, 300.0, 1000.0, 3000.0];

fn signature(params: [f64; 5], bc: &str, lens: &[f64], out: &mut [f64]) -> isize {
    let (nodes, elems, _) = channel();
    unsafe {
        cufsm_signature(
            params.as_ptr(),
            5,
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            bc.as_ptr(),
            bc.len(),
            lens.as_ptr(),
            lens.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    }
}

fn modes(params: [f64; 5], bc: &str, lens: &[f64], out: &mut [f64]) -> isize {
    let (nodes, elems, _) = channel();
    unsafe {
        cufsm_modes(
            params.as_ptr(),
            5,
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            bc.as_ptr(),
            bc.len(),
            lens.as_ptr(),
            lens.len(),
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
    let (_, _, model) = channel();
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
    let n = signature([1.5, 203e3, 0.3, 1.0, 7.0], "S-S", &LENS, &mut out);
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
    let n = signature([1.5, 203e3, 0.3, 1.0, 8.0], "S-S", &LENS, &mut out);
    assert_eq!(n, 18);
}

/// S-S with more than one term would be the curve's minimum over L, L/2, ..., not the signature
/// curve: refused, as are terms out of range for the other boundary conditions.
#[test]
fn terms_are_checked() {
    let _g = serial();
    let mut out = vec![0.0; 64];
    assert!(signature([1.5, 203e3, 0.3, 6.0, 0.0], "S-S", &LENS, &mut out) < 0);
    assert!(last_error().contains("terms must be 1"), "{}", last_error());
    assert!(signature([1.5, 203e3, 0.3, 0.0, 0.0], "C-C", &LENS, &mut out) < 0);
    assert!(signature([1.5, 203e3, 0.3, 2.5, 0.0], "C-C", &LENS, &mut out) < 0);
    // C-C at physical lengths with three terms: the Rust API's answer.
    let (_, _, model) = channel();
    let m = vec![vec![1.0, 2.0, 3.0]; LENS.len()];
    let r = stripmain(&model, &LENS, &m, BoundaryCondition::CC, 1).unwrap();
    assert_eq!(
        signature([1.5, 203e3, 0.3, 3.0, 0.0], "C-C", &LENS, &mut out),
        12
    );
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
    let cases: [([f64; 5], &str, &[f64], &str); 7] = [
        ([1.5, 203e3, 0.3, 1.0, 16.0], "S-S", &LENS, "spaces"),
        ([1.5, 203e3, 0.3, 1.0, 1.5], "S-S", &LENS, "spaces"),
        ([1.5, 203e3, 0.7, 1.0, 0.0], "S-S", &LENS, "Poisson"),
        (
            [1.5, 203e3, 0.3, 1.0, 0.0],
            "X-X",
            &LENS,
            "unknown boundary",
        ),
        (
            [1.5, 203e3, 0.3, 1.0, 0.0],
            "S-S",
            &[100.0, f64::NAN],
            "length",
        ),
        ([1.5, 203e3, 0.3, 1.0, 0.0], "S-S", &[100.0, -1.0], "length"),
        ([0.0, 203e3, 0.3, 1.0, 0.0], "S-S", &LENS, "thickness"),
    ];
    for (params, bc, lens, why) in cases {
        assert_eq!(signature(params, bc, lens, &mut out), -1, "{why}");
        assert!(last_error().contains(why), "{why}: {}", last_error());
        assert_eq!(modes(params, bc, lens, &mut out), -1, "modes: {why}");
        assert!(last_error().contains(why), "modes {why}: {}", last_error());
    }
    // A short output buffer is refused before anything is written.
    let mut short = vec![7.0; 5];
    assert_eq!(
        signature([1.5, 203e3, 0.3, 1.0, 0.0], "S-S", &LENS, &mut short),
        -1
    );
    assert!(short.iter().all(|v| *v == 7.0));
}

/// A node's free flag and an element's node index must be exact.
#[test]
fn node_and_element_buffers_are_checked() {
    let _g = serial();
    let (mut nodes, mut elems, _) = channel();
    let params = [1.5, 203e3, 0.3, 1.0, 0.0];
    let mut out = vec![0.0; 64];
    let call = |nodes: &[f64], elems: &[f64], out: &mut [f64]| unsafe {
        cufsm_signature(
            params.as_ptr(),
            5,
            nodes.as_ptr(),
            nodes.len(),
            elems.as_ptr(),
            elems.len(),
            b"S-S".as_ptr(),
            3,
            LENS.as_ptr(),
            LENS.len(),
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
                5,
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                b"S-S".as_ptr(),
                3,
                LENS.as_ptr(),
                LENS.len(),
                out.as_mut_ptr(),
                out.len(),
            )
        },
        -1
    );
}

/// Per length: the lowest mode's class, its terms and its shape, as the Rust API gives them.
#[test]
fn modes_match_the_rust_api() {
    let _g = serial();
    let (_, _, model) = channel();
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
        modes([1.5, 203e3, 0.3, 1.0, 0.0], "S-S", &lens, &mut out),
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
