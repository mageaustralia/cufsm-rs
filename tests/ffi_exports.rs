//! The C interface's exports added in ABI 2, minor 1 (`cufsm_strip`, `cufsm_classify`,
//! `cufsm_template`, `cufsm_props_wn`, `cufsm_signature_lengths`, `cufsm_signature_minima`):
//!
//! 1. each writes exactly what the Rust API returns, value for value, and
//! 2. through the interface, each matches CUFSM's own code run under Octave: the cases of
//!    `fixtures/cufsm_octave.json`, and those of `fixtures/ffi_octave.json` (written by
//!    `oracle/octave/ffi_cases.m`) for what that fixture does not hold: `signature_ss.m`'s lengths,
//!    a different non-contiguous set of terms at each length, unions of cFSM spaces, and the
//!    classification under every orth and norm option.
//!
//! The tolerances are the crate's own parity tests'. Each oracle test prints its largest difference.
#![cfg(feature = "ffi")]

mod common;
use common::*;
use cufsm::cfsm::{classify_with, stripmain_constrained, Norm, OSpace, Orth, Spaces};
use cufsm::ffi::*;
use cufsm::template::{templatecalc, Shape, Template};
use cufsm::{
    cutwp_prop2, grosprop, signature_minima, signature_ss, signature_ss_lengths, stripmain,
    BoundaryCondition, Dof, Element, LengthResult, Material, Model, Node,
};
use serde_json::Value;

/// The last-error slot is one per process: tests that fail on purpose take turns.
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());
fn serial() -> std::sync::MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

fn last_error() -> String {
    let (p, n) = (cufsm_last_error_ptr(), cufsm_last_error_len());
    String::from_utf8(unsafe { std::slice::from_raw_parts(p as *const u8, n) }.to_vec()).unwrap()
}

/// A model as the interface's flat buffers: mats, nodes, elems, springs, constraints.
struct Bufs {
    mats: Vec<f64>,
    nodes: Vec<f64>,
    elems: Vec<f64>,
    springs: Vec<f64>,
    constraints: Vec<f64>,
}

fn bufs(m: &Model) -> Bufs {
    let f = |b: bool| if b { 1.0 } else { 0.0 };
    let code = |d: Dof| match d {
        Dof::X => 1.0,
        Dof::Z => 2.0,
        Dof::Y => 3.0,
        Dof::Theta => 4.0,
    };
    Bufs {
        mats: m
            .materials
            .iter()
            .flat_map(|t| [t.ex, t.ey, t.vx, t.vy, t.g])
            .collect(),
        nodes: m
            .nodes
            .iter()
            .flat_map(|n| {
                [
                    n.x,
                    n.z,
                    f(n.free[0]),
                    f(n.free[1]),
                    f(n.free[2]),
                    f(n.free[3]),
                    n.stress,
                ]
            })
            .collect(),
        elems: m
            .elements
            .iter()
            .flat_map(|e| [e.ni as f64, e.nj as f64, e.t, e.mat as f64])
            .collect(),
        springs: m
            .springs
            .iter()
            .flat_map(|s| {
                [
                    s.ni as f64,
                    s.nj.map_or(-1.0, |j| j as f64),
                    s.ku,
                    s.kv,
                    s.kw,
                    s.kq,
                    f(s.local),
                    f(s.discrete),
                    s.ys_fraction,
                ]
            })
            .collect(),
        constraints: m
            .constraints
            .iter()
            .flat_map(|c| {
                [
                    c.node_e as f64,
                    code(c.dof_e),
                    c.coeff,
                    c.node_k as f64,
                    code(c.dof_k),
                ]
            })
            .collect(),
    }
}

/// A 100 x 50 x 15 lipped channel, t = 1.5, in uniform compression.
fn channel() -> Model {
    let pts = [
        (50.0, 15.0),
        (50.0, 0.0),
        (0.0, 0.0),
        (0.0, 100.0),
        (50.0, 100.0),
        (50.0, 85.0),
    ];
    Model {
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
    }
}

fn bc_str(bc: BoundaryCondition) -> &'static str {
    match bc {
        BoundaryCondition::SS => "S-S",
        BoundaryCondition::CC => "C-C",
        BoundaryCondition::SC => "S-C",
        BoundaryCondition::CF => "C-F",
        BoundaryCondition::CG => "C-G",
    }
}

/// `[n_1, m.., n_2, m.., ...]`.
fn flat_m_all(m_all: &[Vec<f64>]) -> Vec<f64> {
    m_all
        .iter()
        .flat_map(|m| std::iter::once(m.len() as f64).chain(m.iter().copied()))
        .collect()
}

/// `cufsm_strip`'s raw return and output.
fn strip_raw(
    b: &Bufs,
    params: [f64; 4],
    bc: &str,
    lens: &[f64],
    m_all: &[f64],
    cap: usize,
) -> (isize, Vec<f64>) {
    let mut out = vec![f64::MAX; cap];
    let n = unsafe {
        cufsm_strip(
            params.as_ptr(),
            4,
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            bc.as_ptr(),
            bc.len(),
            lens.as_ptr(),
            lens.len(),
            if m_all.is_empty() {
                std::ptr::null()
            } else {
                m_all.as_ptr()
            },
            m_all.len(),
            b.springs.as_ptr(),
            b.springs.len(),
            b.constraints.as_ptr(),
            b.constraints.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    (n, out)
}

/// One length of `cufsm_strip`'s output.
#[derive(Debug)]
struct Row {
    found: usize,
    m: Vec<f64>,
    lf: Vec<f64>,
    cls: Vec<[f64; 4]>,
    modes: Vec<Vec<f64>>,
}

fn parse_rows(out: &[f64], nlen: usize, nn: usize, neigs: usize) -> Vec<Row> {
    let mut at = 0;
    (0..nlen)
        .map(|_| {
            let found = out[at] as usize;
            let nt = out[at + 1] as usize;
            let m = out[at + 2..at + 2 + nt].to_vec();
            let blk = 5 + 4 * nn * nt;
            let mut row = Row {
                found,
                m,
                lf: vec![],
                cls: vec![],
                modes: vec![],
            };
            for k in 0..neigs {
                let o = at + 2 + nt + k * blk;
                let b = &out[o..o + blk];
                if k < found {
                    row.lf.push(b[0]);
                    row.cls.push([b[1], b[2], b[3], b[4]]);
                    row.modes.push(b[5..].to_vec());
                } else {
                    assert!(b.iter().all(|v| v.is_nan()), "a block past found is NaN");
                }
            }
            at += 2 + nt + neigs * blk;
            row
        })
        .collect()
}

/// `cufsm_strip` parsed, sized exactly from the terms after `msort`.
fn strip(
    m: &Model,
    bc: BoundaryCondition,
    lens: &[f64],
    m_all: &[Vec<f64>],
    neigs: usize,
    spaces: u32,
    classify: bool,
) -> Vec<Row> {
    let nn = m.nodes.len();
    let cap: usize = m_all
        .iter()
        .map(|t| {
            let nt = cufsm::analysis::msort(t).len();
            2 + nt + neigs * (5 + 4 * nn * nt)
        })
        .sum();
    let (n, out) = strip_raw(
        &bufs(m),
        [
            0.0,
            spaces as f64,
            neigs as f64,
            if classify { 1.0 } else { 0.0 },
        ],
        bc_str(bc),
        lens,
        &flat_m_all(m_all),
        cap + 3,
    );
    assert_eq!(
        n,
        cap as isize,
        "{}",
        if n < 0 { last_error() } else { String::new() }
    );
    assert!(
        out[cap..].iter().all(|v| *v == f64::MAX),
        "wrote past the rows"
    );
    parse_rows(&out, lens.len(), nn, neigs)
}

fn spaces_of(bits: u32) -> Spaces {
    Spaces {
        global: bits & 1 != 0,
        distortional: bits & 2 != 0,
        local: bits & 4 != 0,
        other: bits & 8 != 0,
    }
}

/// `cufsm_classify` on given results.
fn classify(
    m: &Model,
    bc: BoundaryCondition,
    results: &[LengthResult],
    opts: [f64; 3],
) -> (isize, Vec<f64>) {
    let b = bufs(m);
    let raw: Vec<f64> = results
        .iter()
        .flat_map(|r| {
            let mut v = vec![r.length, r.m_terms.len() as f64];
            v.extend(&r.m_terms);
            v.push(r.modes.len() as f64);
            for md in &r.modes {
                v.extend(md);
            }
            v
        })
        .collect();
    let total: usize = results.iter().map(|r| r.modes.len()).sum();
    let mut out = vec![f64::MAX; 4 * total + 2];
    let bc = bc_str(bc);
    let n = unsafe {
        cufsm_classify(
            opts.as_ptr(),
            3,
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            bc.as_ptr(),
            bc.len(),
            raw.as_ptr(),
            raw.len(),
            b.springs.as_ptr(),
            b.springs.len(),
            b.constraints.as_ptr(),
            b.constraints.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    if n >= 0 {
        assert!(
            out[n as usize..].iter().all(|v| *v == f64::MAX),
            "wrote past"
        );
        out.truncate(n as usize);
    }
    (n, out)
}

fn rows_to_results(rows: &[Row], lens: &[f64]) -> Vec<LengthResult> {
    rows.iter()
        .zip(lens)
        .map(|(r, &l)| LengthResult {
            length: l,
            m_terms: r.m.clone(),
            load_factors: r.lf.clone(),
            modes: r.modes.clone(),
        })
        .collect()
}

fn assert_rows_equal(
    rows: &[Row],
    want: &[LengthResult],
    cls: Option<&[Vec<[f64; 4]>]>,
    what: &str,
) {
    assert_eq!(rows.len(), want.len());
    for (i, (row, w)) in rows.iter().zip(want).enumerate() {
        assert_eq!(row.m, w.m_terms, "{what} length {i}: terms");
        assert_eq!(row.found, w.load_factors.len(), "{what} length {i}: found");
        assert_eq!(row.lf, w.load_factors, "{what} length {i}: load factors");
        assert_eq!(row.modes, w.modes, "{what} length {i}: modes");
        match cls {
            Some(c) => assert_eq!(row.cls, c[i], "{what} length {i}: classes"),
            None => assert!(
                row.cls.iter().flatten().all(|v| v.is_nan()),
                "{what}: no classes"
            ),
        }
    }
}

// ----- 1. The interface writes what the Rust API returns -----

/// The term-count form with the default classification is `cufsm_modes`, value for value.
#[test]
fn strip_term_count_form_is_cufsm_modes() {
    let _g = serial();
    let m = channel();
    let b = bufs(&m);
    let lens = [300.0, 1200.0];
    for (bc, terms, neigs) in [("C-C", 3.0, 4.0), ("S-S", 1.0, 2.0), ("C-F", 2.0, 3.0)] {
        let per = 2 + terms as usize + neigs as usize * (5 + 4 * 6 * terms as usize);
        let mut want = vec![0.0; 2 * per];
        let n = unsafe {
            cufsm_modes(
                [terms, 0.0, neigs].as_ptr(),
                3,
                b.mats.as_ptr(),
                b.mats.len(),
                b.nodes.as_ptr(),
                b.nodes.len(),
                b.elems.as_ptr(),
                b.elems.len(),
                bc.as_ptr(),
                bc.len(),
                lens.as_ptr(),
                2,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
                want.as_mut_ptr(),
                want.len(),
            )
        };
        assert_eq!(n, (2 * per) as isize, "{}", last_error());
        let (n2, got) = strip_raw(&b, [terms, 0.0, neigs, 1.0], bc, &lens, &[], 2 * per);
        assert_eq!(n2, n);
        // Bit for bit, NaN padding included.
        assert!(
            got.iter()
                .zip(&want)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "{bc}: cufsm_strip differs from cufsm_modes"
        );
    }
}

/// A different set of terms at each length, free and restricted to unions of spaces, with and
/// without the classification: the Rust API's numbers.
#[test]
fn strip_matches_the_rust_api() {
    let _g = serial();
    let m = channel();
    let lens = [250.0, 900.0, 3000.0];
    // Unsorted, repeated and zero terms: msort's clean-up is reported in the rows.
    let m_all = vec![
        vec![3.0, 1.0, 1.0, 0.0],
        vec![2.0, 4.0],
        vec![1.0, 2.0, 3.0, 7.0],
    ];
    for bc in [
        BoundaryCondition::CC,
        BoundaryCondition::SS,
        BoundaryCondition::CF,
    ] {
        let free = stripmain(&m, &lens, &m_all, bc, 4).unwrap();
        let cls = cufsm::cfsm::classify(&m, &free, bc, Orth::Axial, Norm::Vector).unwrap();
        let rows = strip(&m, bc, &lens, &m_all, 4, 0, true);
        assert_eq!(rows[0].m, vec![1.0, 3.0]);
        assert_rows_equal(&rows, &free, Some(&cls), "free");
        assert_rows_equal(
            &strip(&m, bc, &lens, &m_all, 4, 0, false),
            &free,
            None,
            "unclassified",
        );
        for bits in [1, 2, 4, 3, 6, 12, 7, 15] {
            let want = stripmain_constrained(&m, &lens, &m_all, bc, 4, spaces_of(bits)).unwrap();
            let rows = strip(&m, bc, &lens, &m_all, 4, bits, false);
            assert_rows_equal(&rows, &want, None, &format!("{bc:?} spaces {bits}"));
        }
        // Constrained and classified.
        let want = stripmain_constrained(&m, &lens, &m_all, bc, 4, spaces_of(2)).unwrap();
        let c = cufsm::cfsm::classify(&m, &want, bc, Orth::Axial, Norm::Vector).unwrap();
        assert_rows_equal(
            &strip(&m, bc, &lens, &m_all, 4, 2, true),
            &want,
            Some(&c),
            "D classified",
        );
    }
}

/// Springs and constraints reach `cufsm_strip`, free and constrained.
#[test]
fn strip_takes_springs_and_constraints() {
    let _g = serial();
    let mut m = channel();
    m.springs.push(cufsm::model::Spring {
        ni: 0,
        nj: None,
        ku: 1e3,
        kv: 0.0,
        kw: 1e3,
        kq: 0.0,
        local: false,
        discrete: false,
        ys_fraction: 0.0,
    });
    m.constraints.push(cufsm::model::Constraint {
        node_e: 5,
        dof_e: Dof::Z,
        coeff: 1.0,
        node_k: 0,
        dof_k: Dof::Z,
    });
    let lens = [600.0, 2500.0];
    let m_all = vec![vec![1.0, 2.0], vec![1.0, 3.0]];
    let bc = BoundaryCondition::CC;
    let free = stripmain(&m, &lens, &m_all, bc, 3).unwrap();
    assert_rows_equal(
        &strip(&m, bc, &lens, &m_all, 3, 0, false),
        &free,
        None,
        "free",
    );
    let want = stripmain_constrained(&m, &lens, &m_all, bc, 3, spaces_of(3)).unwrap();
    assert_rows_equal(
        &strip(&m, bc, &lens, &m_all, 3, 3, false),
        &want,
        None,
        "GD",
    );
}

/// Every orth, norm and O space option: `classify_with`'s numbers.
#[test]
fn classify_matches_the_rust_api() {
    let _g = serial();
    let m = channel();
    let lens = [60.0, 250.0, 3000.0];
    let mut refused = 0;
    for bc in [BoundaryCondition::SS, BoundaryCondition::CC] {
        let m_all = if bc == BoundaryCondition::SS {
            vec![vec![1.0]; 3]
        } else {
            vec![vec![1.0, 2.0, 3.0]; 3]
        };
        let res = stripmain(&m, &lens, &m_all, bc, 3).unwrap();
        for (o, orth) in [(1.0, Orth::Natural), (2.0, Orth::Axial), (3.0, Orth::Load)] {
            for (nv, norm) in [
                (0.0, Norm::None),
                (1.0, Norm::Vector),
                (2.0, Norm::StrainEnergy),
                (3.0, Norm::Work),
            ] {
                for (s, os) in [
                    (1.0, OSpace::St),
                    (2.0, OSpace::K),
                    (3.0, OSpace::Kg),
                    (4.0, OSpace::Vector),
                ] {
                    let (n, got) = classify(&m, bc, &res, [o, nv, s]);
                    match classify_with(&m, &res, bc, orth, norm, os) {
                        Ok(want) => {
                            assert_eq!(
                                n,
                                4 * 9,
                                "{}",
                                if n < 0 { last_error() } else { String::new() }
                            );
                            let flat: Vec<f64> = want.iter().flatten().flatten().copied().collect();
                            assert!(
                                got.iter()
                                    .zip(&flat)
                                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                                "{bc:?} orth {o} norm {nv} ospace {s}"
                            );
                        }
                        // A combination CUFSM cannot run either: the same refusal.
                        Err(e) => {
                            assert_eq!(n, -1, "{bc:?} orth {o} norm {nv} ospace {s}");
                            assert_eq!(last_error(), e.to_string());
                            refused += 1;
                        }
                    }
                }
            }
        }
    }
    // Natural orthogonality with the null-space O basis and no energy norm reaches eig() without K
    // in CUFSM; the crate refuses those, and so does the interface.
    assert!(refused > 0 && refused < 2 * 48, "{refused} refused");
}

/// Every length may carry any number of modes, none included.
#[test]
fn classify_takes_uneven_mode_counts() {
    let _g = serial();
    let m = channel();
    let lens = [60.0, 250.0, 900.0];
    let mut res = stripmain(
        &m,
        &lens,
        &[vec![1.0], vec![1.0], vec![1.0]],
        BoundaryCondition::SS,
        3,
    )
    .unwrap();
    res[0].modes.truncate(1);
    res[1].modes.clear();
    let want = classify_with(
        &m,
        &res,
        BoundaryCondition::SS,
        Orth::Axial,
        Norm::Vector,
        OSpace::St,
    )
    .unwrap();
    let (n, got) = classify(&m, BoundaryCondition::SS, &res, [2.0, 1.0, 1.0]);
    assert_eq!(n, 4 * 4, "{}", last_error());
    let flat: Vec<f64> = want.iter().flatten().flatten().copied().collect();
    assert_eq!(got, flat);
}

fn template_params(t: &Template) -> [f64; 23] {
    [
        if t.shape == Shape::C { 1.0 } else { 2.0 },
        t.h,
        t.b1,
        t.b2,
        t.d1,
        t.d2,
        t.r1,
        t.r2,
        t.r3,
        t.r4,
        t.q1,
        t.q2,
        t.t,
        t.nh as f64,
        t.nb1 as f64,
        t.nb2 as f64,
        t.nd1 as f64,
        t.nd2 as f64,
        t.nr1 as f64,
        t.nr2 as f64,
        t.nr3 as f64,
        t.nr4 as f64,
        if t.centerline { 1.0 } else { 0.0 },
    ]
}

fn template_call(p: &[f64], cap: usize) -> (isize, Vec<f64>) {
    let mut out = vec![f64::MAX; cap];
    let n = unsafe { cufsm_template(p.as_ptr(), p.len(), out.as_mut_ptr(), out.len()) };
    (n, out)
}

/// The template's bound on the output: `2 + 7 (s + 1) + 4 s`, `s` the strip counts' sum.
fn template_bound(t: &Template) -> usize {
    let s = t.nh + t.nb1 + t.nb2 + t.nd1 + t.nd2 + t.nr1 + t.nr2 + t.nr3 + t.nr4;
    2 + 7 * (s + 1) + 4 * s
}

/// The nodes and elements `templatecalc` gives, in the interface's own layouts, ready to pass back.
#[test]
fn template_matches_the_rust_api() {
    let _g = serial();
    for t in [
        Template::outside(Shape::C, 200.0, 75.0, 20.0, 1.9, 3.0, 12),
        Template::outside(Shape::Z, 250.0, 80.0, 0.0, 2.4, 0.0, 10),
        Template::outside(Shape::C, 150.0, 60.0, 15.0, 1.5, 0.0, 8),
    ] {
        let (n, out) = template_call(&template_params(&t), template_bound(&t));
        let m = templatecalc(&t, Material::isotropic(1.0, 0.3));
        let b = bufs(&m);
        let need = 2 + b.nodes.len() + b.elems.len();
        assert_eq!(n, need as isize, "{}", last_error());
        assert_eq!(out[..2], [m.nodes.len() as f64, m.elements.len() as f64]);
        assert_eq!(out[2..2 + b.nodes.len()], b.nodes[..]);
        assert_eq!(out[2 + b.nodes.len()..need], b.elems[..]);
        // The tables go straight back into the other exports.
        let mats = [200e3, 200e3, 0.3, 0.3, 200e3 / 2.6];
        let mut props = vec![0.0; 15];
        let nodes = &out[2..2 + b.nodes.len()];
        let elems = &out[2 + b.nodes.len()..need];
        let k = unsafe {
            cufsm_props(
                mats.as_ptr(),
                5,
                nodes.as_ptr(),
                nodes.len(),
                elems.as_ptr(),
                elems.len(),
                props.as_mut_ptr(),
                15,
            )
        };
        assert_eq!(k, 15, "{}", last_error());
        assert!(props[0] > 0.0);
    }
}

#[test]
fn props_wn_match_the_rust_api() {
    let _g = serial();
    let m = channel();
    let b = bufs(&m);
    let mut out = vec![f64::MAX; 15 + 6 + 1];
    let n = unsafe {
        cufsm_props_wn(
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert_eq!(n, 21, "{}", last_error());
    let g = grosprop(&m);
    let c = cutwp_prop2(&m);
    let mut want = vec![
        g.a, g.xcg, g.zcg, g.ixx, g.izz, g.ixz, g.thetap, g.i11, g.i22, c.j, c.xs, c.zs, c.cw,
        c.b1, c.b2,
    ];
    want.extend(&c.wn);
    assert_eq!(out[..21], want[..]);
    assert_eq!(out[21], f64::MAX);
}

fn signature_lengths(m: &Model) -> (isize, Vec<f64>) {
    let b = bufs(m);
    let mut out = vec![0.0; 100];
    let n = unsafe {
        cufsm_signature_lengths(
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    (n, out)
}

fn minima(curve: &[f64], cap: usize) -> (isize, Vec<f64>) {
    let mut out = vec![f64::MAX; cap];
    let n =
        unsafe { cufsm_signature_minima(curve.as_ptr(), curve.len(), out.as_mut_ptr(), out.len()) };
    (n, out)
}

/// The lengths are `signature_ss`'s, and the minima of a curve at them `signature_minima`'s.
#[test]
fn signature_lengths_and_minima_match_the_rust_api() {
    let _g = serial();
    let m = channel();
    let (n, ls) = signature_lengths(&m);
    assert_eq!(n, 100, "{}", last_error());
    assert_eq!(ls, signature_ss_lengths(&m));
    let curve = signature_ss(&m, 1).unwrap();
    assert_eq!(
        ls,
        curve.iter().map(|r| r.length).collect::<Vec<_>>(),
        "signature_ss solves at these lengths"
    );
    let want = signature_minima(&curve);
    assert!(want.len() >= 2, "local and distortional minima");
    let flat: Vec<f64> = curve
        .iter()
        .flat_map(|r| [r.length, r.load_factors[0]])
        .collect();
    let (n, out) = minima(&flat, flat.len() - 4);
    assert_eq!(n, 2 * want.len() as isize, "{}", last_error());
    for (o, w) in out.chunks(2).zip(&want) {
        assert_eq!((o[0], o[1]), (w.length, w.load_factor));
    }
    // A NaN load factor is a length without one, as signature_minima skips an empty length.
    let mut holed = curve.clone();
    holed[40].load_factors.clear();
    let mut flat_h = flat.clone();
    flat_h[81] = f64::NAN;
    let want = signature_minima(&holed);
    let (n, out) = minima(&flat_h, flat_h.len());
    assert_eq!(n, 2 * want.len() as isize);
    for (o, w) in out.chunks(2).zip(&want) {
        assert_eq!((o[0], o[1]), (w.length, w.load_factor));
    }
}

// ----- Limits and refusals -----

/// No fixed caps: more than 50 modes and more than 100 terms are accepted (the output is sized from
/// the caller's numbers), and a problem too large to address is refused with a message.
#[test]
fn neigs_and_terms_have_no_fixed_cap() {
    let _g = serial();
    // A flat plate: 2 nodes, 8 dofs a term.
    let plate = Model {
        materials: vec![Material::isotropic(200e3, 0.3)],
        nodes: vec![Node::new(0.0, 0.0, 1.0), Node::new(100.0, 0.0, 1.0)],
        elements: vec![Element {
            ni: 0,
            nj: 1,
            t: 2.0,
            mat: 0,
        }],
        constraints: vec![],
        springs: vec![],
    };
    let m_all = vec![(1..=101).map(f64::from).collect::<Vec<_>>()];
    let want = stripmain(&plate, &[5000.0], &m_all, BoundaryCondition::CC, 60).unwrap();
    let rows = strip(
        &plate,
        BoundaryCondition::CC,
        &[5000.0],
        &m_all,
        60,
        0,
        false,
    );
    assert_eq!(rows[0].found, want[0].load_factors.len());
    assert!(rows[0].found > 50, "found {}", rows[0].found);
    assert_eq!(rows[0].lf, want[0].load_factors);
    // The term-count form, through cufsm_modes too: 101 terms and 60 modes.
    let b = bufs(&plate);
    let per = 2 + 101 + 60 * (5 + 8 * 101);
    let mut out = vec![0.0; per];
    let n = unsafe {
        cufsm_modes(
            [101.0, 0.0, 60.0].as_ptr(),
            3,
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            "C-C".as_ptr(),
            3,
            [5000.0].as_ptr(),
            1,
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert_eq!(n, per as isize, "{}", last_error());
    assert_eq!(out[0] as usize, want[0].load_factors.len());
    // A billion terms cannot be addressed: refused before anything is allocated.
    let (n, _) = strip_raw(&b, [1e9, 0.0, 1.0, 0.0], "C-C", &[5000.0], &[], 4);
    assert_eq!(n, -1);
    assert!(last_error().contains("address"), "{}", last_error());
}

#[test]
fn bad_inputs_fail_cleanly() {
    let _g = serial();
    let m = channel();
    let b = bufs(&m);
    let lens = [300.0, 900.0];
    let ok_m = [1.0, 1.0, 2.0, 1.0, 3.0];
    let cases: [([f64; 4], &str, &[f64], &str); 9] = [
        ([1.0, 0.0, 1.0, 0.0], "C-C", &ok_m, "terms must be 0"),
        (
            [0.0, 0.0, 1.0, 0.0],
            "C-C",
            &[1.0, 1.0],
            "term lists for 2 lengths",
        ),
        (
            [0.0, 0.0, 1.0, 0.0],
            "C-C",
            &[1.0, 1.0, 5.0, 1.0],
            "only 1 values follow",
        ),
        (
            [0.0, 0.0, 1.0, 0.0],
            "C-C",
            &[1.0, 0.0, 1.0, 2.0],
            "no nonzero term",
        ),
        (
            [0.0, 0.0, 1.0, 0.0],
            "C-C",
            &[1.0, f64::NAN, 1.0, 2.0],
            "not finite",
        ),
        ([0.0, 0.0, 1.0, 0.0], "C-C", &[0.5, 1.0, 1.0, 2.0], "count"),
        ([0.0, 0.0, 1.0, 2.0], "C-C", &ok_m, "classify"),
        ([0.0, 16.0, 1.0, 0.0], "C-C", &ok_m, "spaces"),
        ([0.0, 0.0, 1.0, 0.0], "X-X", &ok_m, "boundary"),
    ];
    for (params, bc, m_all, why) in cases {
        let (n, out) = strip_raw(&b, params, bc, &lens, m_all, 4096);
        assert_eq!(n, -1, "{why}");
        assert!(last_error().contains(why), "{why}: {}", last_error());
        assert!(
            out.iter().all(|v| *v == f64::MAX),
            "{why}: wrote on failure"
        );
    }
    // Lower case is accepted, as the Python package accepts it.
    let (n, _) = strip_raw(&b, [0.0, 0.0, 1.0, 0.0], "c-c", &lens, &ok_m, 4096);
    assert!(n > 0, "{}", last_error());
    // A short buffer is refused before anything is written.
    let (n, out) = strip_raw(&b, [0.0, 0.0, 1.0, 0.0], "C-C", &lens, &ok_m, 10);
    assert!(n < 0 && last_error().contains("too small"));
    assert!(out.iter().all(|v| *v == f64::MAX));

    // cufsm_classify: option codes and a malformed results buffer.
    let res = stripmain(&m, &[300.0], &[vec![1.0]], BoundaryCondition::SS, 1).unwrap();
    for (opts, why) in [
        ([0.0, 1.0, 1.0], "orth"),
        ([2.0, 4.0, 1.0], "norm"),
        ([2.0, 1.0, 5.0], "ospace"),
    ] {
        assert_eq!(classify(&m, BoundaryCondition::SS, &res, opts).0, -1);
        assert!(last_error().contains(why), "{why}: {}", last_error());
    }
    let mut short = res.clone();
    short[0].modes[0].pop();
    assert_eq!(
        classify(&m, BoundaryCondition::SS, &short, [2.0, 1.0, 1.0]).0,
        -1
    );
    assert!(last_error().contains("mode count"), "{}", last_error());

    // cufsm_template: dimensions, counts, flags, and the output bound.
    let t = Template::outside(Shape::C, 200.0, 75.0, 20.0, 1.9, 3.0, 12);
    let good = template_params(&t);
    for (k, v, why) in [
        (0, 3.0, "shape"),
        (1, 0.0, "h = 0"),
        (4, -1.0, "d1"),
        (13, 0.0, "nh"),
        (16, 1.5, "nd1"),
        (22, 2.0, "centerline"),
    ] {
        let mut p = good;
        p[k] = v;
        assert_eq!(template_call(&p, 4096).0, -1, "{why}");
        assert!(last_error().contains(why), "{why}: {}", last_error());
    }
    assert_eq!(template_call(&good[..22], 4096).0, -1);
    let (n, out) = template_call(&good, template_bound(&t) - 1);
    assert!(n < 0 && last_error().contains("too small") && out.iter().all(|v| *v == f64::MAX));
    let mut huge = good;
    huge[13] = 1e15;
    assert_eq!(
        template_call(&huge, 4096).0,
        -1,
        "an absurd count is refused, not meshed"
    );

    // cufsm_signature_minima: pairs, positive lengths.
    assert_eq!(minima(&[1.0, 2.0, 3.0], 4).0, -1);
    assert_eq!(minima(&[1.0, 2.0, -3.0, 1.0], 4).0, -1);
    assert_eq!(minima(&[1.0, 2.0, 3.0, f64::INFINITY], 4).0, -1);
}

/// A node no strip touches is refused with a message naming it, by the free and the cFSM
/// analyses, the classification and the section properties alike: -1, not a trap. The cFSM
/// paths used to index out of bounds on it, which in wasm aborts the module.
#[test]
fn a_node_in_no_element_is_refused() {
    let _g = serial();
    let good = channel();
    let mut m = good.clone();
    m.nodes.push(Node::new(25.0, 50.0, 1.0));
    let b = bufs(&m);
    let why = "node 6 belongs to no element";
    let check = |what: &str, n: isize, out: &[f64]| {
        assert_eq!(n, -1, "{what}");
        let err = last_error();
        assert!(err.contains(why) && !err.contains("panic"), "{what}: {err}");
        assert!(
            out.iter().all(|v| *v == f64::MAX),
            "{what}: wrote on failure"
        );
    };
    // cufsm_strip: free, free and classified, and restricted to G+D+L and to G alone.
    for params in [
        [1.0, 0.0, 2.0, 0.0],
        [1.0, 0.0, 2.0, 1.0],
        [2.0, 7.0, 2.0, 1.0],
        [1.0, 1.0, 1.0, 0.0],
    ] {
        let (n, out) = strip_raw(&b, params, "S-S", &[100.0, 1000.0], &[], 4096);
        check(&format!("strip {params:?}"), n, &out);
    }
    // cufsm_classify, given modes of the same model with the node in place (zeros there).
    let mut res = stripmain(&good, &[100.0], &[vec![1.0]], BoundaryCondition::SS, 2).unwrap();
    for md in &mut res[0].modes {
        md.splice(12..12, [0.0, 0.0]);
        md.extend([0.0, 0.0]);
    }
    let (n, out) = classify(&m, BoundaryCondition::SS, &res, [2.0, 1.0, 1.0]);
    check("classify", n, &out[n.max(0) as usize..]);
    // Section properties and the signature lengths read the same model.
    let mut out = vec![f64::MAX; 64];
    let n = unsafe {
        cufsm_props_wn(
            b.mats.as_ptr(),
            b.mats.len(),
            b.nodes.as_ptr(),
            b.nodes.len(),
            b.elems.as_ptr(),
            b.elems.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    check("props_wn", n, &out);
    let (n, _) = signature_lengths(&m);
    assert_eq!(n, -1);
    assert!(last_error().contains(why), "{}", last_error());
}

#[test]
fn abi_is_2_minor_1() {
    assert_eq!(cufsm_abi_version(), 2);
    assert_eq!(cufsm_abi_minor(), 1);
}

// ----- 2. Through the interface, against CUFSM under Octave -----

fn ffi_fixture() -> Value {
    load("ffi_octave.json").remove(0)
}

fn case(name: &str) -> Value {
    load("cufsm_octave.json")
        .into_iter()
        .find(|r| r["name"] == name)
        .unwrap_or_else(|| panic!("no case {name}"))
}

/// `signature_ss.m`'s 100 lengths for every model of the main fixture.
#[test]
fn signature_lengths_match_cufsm() {
    let _g = serial();
    let mut worst = (0.0_f64, String::new());
    let mut n_cases = 0;
    for rec in ffi_fixture()["signature_lengths"].as_array().unwrap() {
        let name = rec["name"].as_str().unwrap();
        let m = model_of(&case(name));
        let (n, got) = signature_lengths(&m);
        assert_eq!(n, 100, "{name}: {}", last_error());
        let want = vec_of(&rec["lengths"]);
        assert_eq!(want.len(), 100);
        for (g, w) in got.iter().zip(&want) {
            let d = (g / w - 1.0).abs();
            if d > worst.0 {
                worst = (d, name.to_string());
            }
            assert!(d <= 1e-13, "{name}: {g} vs CUFSM {w}");
        }
        n_cases += 1;
    }
    assert!(n_cases >= 37, "only {n_cases}");
    eprintln!(
        "signature lengths: {n_cases} models x 100, largest relative difference {:e} ({})",
        worst.0, worst.1
    );
}

/// Load factors with a different, non-contiguous set of terms at each length (S-S, C-C, C-F, one
/// with springs), and with the main fixture's 1..n terms, through `cufsm_strip`.
#[test]
fn strip_terms_match_cufsm() {
    let _g = serial();
    // (name, model, bc, lengths, terms, CUFSM's load factors, from the new fixture)
    type Job = (
        String,
        Model,
        BoundaryCondition,
        Vec<f64>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        bool,
    );
    let mut jobs: Vec<Job> = vec![];
    for rec in ffi_fixture()["terms"].as_array().unwrap() {
        let name = rec["name"].as_str().unwrap().to_string();
        let m = model_of(&case(&name));
        let bc = BoundaryCondition::parse(rec["bc"].as_str().unwrap()).unwrap();
        jobs.push((
            name,
            m,
            bc,
            vec_of(&rec["lengths"]),
            list_of(&rec["m_all"]),
            list_of(&rec["load_factors"]),
            true,
        ));
    }
    for r in load("cufsm_octave.json") {
        let m_all = list_of(&r["m_all"]);
        if m_all.iter().all(|t| t.len() == 1) {
            continue;
        }
        jobs.push((
            r["name"].as_str().unwrap().to_string(),
            model_of(&r),
            bc_of(&r),
            vec_of(&r["lengths"]),
            m_all,
            list_of(&r["load_factors_dense"]),
            false,
        ));
    }
    let mut worst = (0.0_f64, String::new());
    let (mut points, mut fresh) = (0, 0);
    for (name, m, bc, lens, m_all, want, new) in &jobs {
        let neigs = want.iter().map(Vec::len).max().unwrap();
        let rows = strip(m, *bc, lens, m_all, neigs, 0, false);
        for (l, (row, w)) in rows.iter().zip(want).enumerate() {
            let cond = cond_estimate(m, lens[l], *bc, &m_all[l]);
            let k = w.len().min(row.lf.len());
            assert!(k > 0, "{name} length {}", lens[l]);
            for i in 0..k {
                let tol = load_factor_tolerance(cond, w[i], w[0]);
                let d = (row.lf[i] / w[i] - 1.0).abs();
                if d > worst.0 {
                    worst = (d, format!("{name} length {} mode {}", lens[l], i + 1));
                }
                assert!(
                    d <= tol,
                    "{name} length {} mode {}: {} vs CUFSM {} ({d:e}, allowed {tol:e})",
                    lens[l],
                    i + 1,
                    row.lf[i],
                    w[i]
                );
                points += 1;
                fresh += usize::from(*new);
            }
        }
    }
    eprintln!(
        "cufsm_strip terms: {points} load factors ({fresh} with per-length non-contiguous terms), \
         largest relative difference {:e} ({})",
        worst.0, worst.1
    );
    assert!(
        fresh >= 50 && points > fresh + 50,
        "only {fresh} and {points}"
    );
}

/// Load factors restricted to cFSM spaces through `cufsm_strip`: the single spaces of the main
/// fixture and the unions of the new one.
#[test]
fn strip_spaces_match_cufsm() {
    let _g = serial();
    let mut worst = (0.0_f64, String::new());
    let mut points = 0;
    let fx = ffi_fixture();
    for r in load("cufsm_octave.json")
        .into_iter()
        .filter(|r| r["cfsm"].get("ngm").is_some())
    {
        let name = r["name"].as_str().unwrap().to_string();
        let m = model_of(&r);
        let lens = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        let union = fx["cfsm"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["name"] == r["name"]);
        let mut keys: Vec<(u32, Value)> = vec![
            (1, r["cfsm"]["lf_G"].clone()),
            (2, r["cfsm"]["lf_D"].clone()),
            (4, r["cfsm"]["lf_L"].clone()),
        ];
        if let Some(u) = union {
            for (bits, k) in [(3, "lf_GD"), (6, "lf_DL"), (12, "lf_LO"), (7, "lf_GDL")] {
                keys.push((bits, u[k].clone()));
            }
        }
        for (bits, want) in keys {
            let want = list_of(&want);
            let rows = strip(&m, bc_of(&r), &lens, &m_all, 5, bits, false);
            for (l, (row, w)) in rows.iter().zip(&want).enumerate() {
                assert_eq!(
                    row.lf.is_empty(),
                    w.is_empty(),
                    "{name} spaces {bits} length {}",
                    lens[l]
                );
                let cond = restricted_cond(&m, lens[l], bc_of(&r), spaces_of(bits));
                for (i, (g, wv)) in row.lf.iter().zip(w).enumerate() {
                    let tol = load_factor_tolerance(cond, *wv, w[0]);
                    let d = (g / wv - 1.0).abs();
                    if d > worst.0 {
                        worst = (
                            d,
                            format!("{name} spaces {bits} length {} mode {}", lens[l], i + 1),
                        );
                    }
                    assert!(
                        d <= tol,
                        "{name} spaces {bits} length {} mode {}: {g} vs CUFSM {wv} ({d:e}, allowed {tol:e})",
                        lens[l],
                        i + 1
                    );
                    points += 1;
                }
            }
        }
    }
    eprintln!(
        "cufsm_strip spaces: {points} restricted load factors, largest relative difference {:e} ({})",
        worst.0, worst.1
    );
    assert!(points > 300, "only {points}");
}

/// The classification through the interface, the way a caller gets it: modes from `cufsm_strip`,
/// then `cufsm_classify` with the options, against CUFSM's. A mode with a near-repeated load
/// factor has no unique shape, and an axial or load-orthogonal basis is not unique where its
/// subproblems have repeated eigenvalues (the doubly symmetric I-section): those are skipped, as
/// in `tests/cfsm_parity.rs`. Held to 1e-6 percentage points, and the strain-energy and work norms
/// to 1e-5: they weigh each base vector by `K` or `Kg`, which carry the rounding of a long global
/// mode (2.7e-6 for the hat's first mode at 3000 mm under the natural basis).
///
/// Not compared: axial or load orthogonality with no normalisation (orth 2 or 3, norm 0). There
/// the percentages depend on how `eig(Ksub, Kgsub)` scales its eigenvectors. CUFSM forms `Ksub`
/// and `Kgsub` as products that are not exactly symmetric, so Octave takes its general solver and
/// scales each vector to a largest entry of 1; this crate's symmetric solver returns
/// `Kg`-normalised vectors, which makes its norm 0 equal to its work norm (3), and that is
/// compared with CUFSM's. `classify_matches_the_rust_api` covers the interface for every option.
#[test]
fn classify_matches_cufsm() {
    let _g = serial();
    let fx = ffi_fixture();
    let mut worst = (0.0_f64, String::new());
    let mut compared = 0;
    let mut per_option = std::collections::BTreeMap::new();
    for r in load("cufsm_octave.json")
        .into_iter()
        .filter(|r| r["cfsm"].get("ngm").is_some())
    {
        let name = r["name"].as_str().unwrap().to_string();
        if name.contains("fixed and sprung") {
            continue;
        }
        let m = model_of(&r);
        let lens = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        let rows = strip(&m, bc_of(&r), &lens, &m_all, 3, 0, false);
        let res = rows_to_results(&rows, &lens);
        let mut keys: Vec<([f64; 3], Value)> = vec![
            ([2.0, 1.0, 1.0], r["cfsm"]["classification"].clone()),
            ([1.0, 1.0, 1.0], r["cfsm"]["classification_natural"].clone()),
            ([2.0, 1.0, 2.0], r["cfsm"]["classification_ospace2"].clone()),
            ([2.0, 1.0, 3.0], r["cfsm"]["classification_ospace3"].clone()),
            ([2.0, 1.0, 4.0], r["cfsm"]["classification_ospace4"].clone()),
        ];
        if let Some(u) = fx["cfsm"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["name"] == r["name"])
        {
            for orth in 1..=3 {
                for norm in 0..=3 {
                    let k = format!("classification_orth{orth}_norm{norm}");
                    if !u[&k].is_null() {
                        keys.push(([orth as f64, norm as f64, 1.0], u[&k].clone()));
                    }
                }
            }
        }
        for (opts, want) in keys {
            if name.contains("i-section") && opts[0] != 1.0 {
                continue;
            }
            if opts[0] != 1.0 && opts[1] == 0.0 {
                continue;
            }
            let tol = if opts[1] >= 2.0 { 1e-5 } else { 1e-6 };
            let (n, got) = classify(&m, bc_of(&r), &res, opts);
            assert!(n >= 0, "{name} {opts:?}: {}", last_error());
            let mut at = 0;
            for (l, row) in rows.iter().enumerate() {
                let wrows = rows_of(&want[l]);
                for q in 0..row.found {
                    let g = &got[at + 4 * q..at + 4 * q + 4];
                    let lfs = &row.lf;
                    let distinct = (q == 0 || (lfs[q] / lfs[q - 1] - 1.0).abs() > 1e-6)
                        && (q + 1 >= lfs.len() || (lfs[q + 1] / lfs[q] - 1.0).abs() > 1e-6);
                    if !distinct || q >= wrows.len() {
                        continue;
                    }
                    let w = &wrows[q];
                    for k in 0..4 {
                        let d = (g[k] - w[k]).abs();
                        if d > worst.0 {
                            worst = (
                                d,
                                format!("{name} {opts:?} length {} mode {}", lens[l], q + 1),
                            );
                        }
                        assert!(
                            d < tol,
                            "{name} {opts:?} length {} mode {}: {g:?} vs CUFSM {w:?}",
                            lens[l],
                            q + 1
                        );
                    }
                    compared += 1;
                    *per_option.entry(format!("{opts:?}")).or_insert(0) += 1;
                }
                at += 4 * row.found;
            }
        }
    }
    eprintln!(
        "cufsm_classify: {compared} modes over {} option sets, largest difference {:e} percentage points ({})",
        per_option.len(),
        worst.0,
        worst.1
    );
    assert_eq!(
        per_option.len(),
        5 + 12 - 2 - 2,
        "every option set compared: {per_option:?}"
    );
    assert!(per_option.values().all(|&n| n >= 20), "{per_option:?}");
}

/// Every template case of the main fixture, through `cufsm_template`.
#[test]
fn template_matches_cufsm() {
    let _g = serial();
    let mut worst = 0.0_f64;
    let mut cases = 0;
    for r in load("cufsm_octave.json") {
        let Some(tp) = r.get("template") else {
            continue;
        };
        let name = r["name"].as_str().unwrap();
        let f = |k: &str| tp[k].as_f64().unwrap();
        let p = [
            f("CorZ"),
            f("h"),
            f("b1"),
            f("b2"),
            f("d1"),
            f("d2"),
            f("r1"),
            f("r2"),
            f("r3"),
            f("r4"),
            f("q1"),
            f("q2"),
            f("t"),
            f("nh"),
            f("nb1"),
            f("nb2"),
            f("nd1"),
            f("nd2"),
            f("nr1"),
            f("nr2"),
            f("nr3"),
            f("nr4"),
            f("center"),
        ];
        let (n, out) = template_call(&p, 4096);
        assert!(n > 0, "{name}: {}", last_error());
        let want = rows_of(&r["node"]);
        let el = rows_of(&r["elem"]);
        let (nn, ne) = (out[0] as usize, out[1] as usize);
        assert_eq!((nn, ne), (want.len(), el.len()), "{name}: counts");
        for (i, w) in want.iter().enumerate() {
            let o = &out[2 + 7 * i..2 + 7 * i + 7];
            let d = (o[0] - w[1]).abs().max((o[1] - w[2]).abs());
            worst = worst.max(d);
            assert!(
                d < 1e-12,
                "{name}: node {} at ({}, {}), CUFSM ({}, {})",
                i + 1,
                o[0],
                o[1],
                w[1],
                w[2]
            );
        }
        for (i, w) in el.iter().enumerate() {
            let o = &out[2 + 7 * nn + 4 * i..2 + 7 * nn + 4 * i + 4];
            assert!(
                o[0] as usize + 1 == w[1] as usize
                    && o[1] as usize + 1 == w[2] as usize
                    && (o[2] - w[3]).abs() < 1e-15,
                "{name}: element {}",
                i + 1
            );
        }
        cases += 1;
    }
    assert!(cases >= 12, "only {cases} template cases");
    eprintln!("cufsm_template: {cases} cases, largest node coordinate difference {worst:e}");
}

/// The section properties and the warping function through `cufsm_props_wn`.
#[test]
fn props_wn_match_cufsm() {
    let _g = serial();
    let mut worst_wn = 0.0_f64;
    let mut cases = 0;
    for r in load("cufsm_octave.json")
        .into_iter()
        .filter(|r| r["cfsm"].get("cutwp").is_some())
    {
        let name = r["name"].as_str().unwrap();
        let m = model_of(&r);
        let b = bufs(&m);
        let nn = m.nodes.len();
        let mut out = vec![0.0; 15 + nn];
        let n = unsafe {
            cufsm_props_wn(
                b.mats.as_ptr(),
                b.mats.len(),
                b.nodes.as_ptr(),
                b.nodes.len(),
                b.elems.as_ptr(),
                b.elems.len(),
                out.as_mut_ptr(),
                out.len(),
            )
        };
        assert_eq!(n, (15 + nn) as isize, "{name}: {}", last_error());
        let w = &r["cfsm"]["cutwp"];
        let fw = |k: &str| w[k].as_f64().unwrap();
        for (k, i) in [
            ("J", 9),
            ("xs", 10),
            ("zs", 11),
            ("Cw", 12),
            ("B1", 13),
            ("B2", 14),
        ] {
            let scale = if k == "J" || k == "Cw" {
                fw(k).abs().max(1.0)
            } else {
                fw("A").sqrt()
            };
            assert!(
                (out[i] - fw(k)).abs() <= 1e-10 * scale,
                "{name}: {k} {} vs CUFSM {}",
                out[i],
                fw(k)
            );
        }
        let wn = vec_of(&w["wn"]);
        let wscale = wn.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1.0);
        for (i, (a, b)) in out[15..].iter().zip(&wn).enumerate() {
            let d = (a - b).abs() / wscale;
            worst_wn = worst_wn.max(d);
            assert!(d <= 1e-10, "{name}: wn at node {} {a} vs CUFSM {b}", i + 1);
        }
        let p = &r["props"];
        for (k, i) in [
            ("A", 0),
            ("xcg", 1),
            ("zcg", 2),
            ("Ixx", 3),
            ("Izz", 4),
            ("Ixz", 5),
        ] {
            let wv = p[k].as_f64().unwrap();
            assert!(
                (out[i] - wv).abs() <= 1e-12 * wv.abs().max(out[0]),
                "{name}: {k}"
            );
        }
        cases += 1;
    }
    assert!(cases >= 7, "only {cases}");
    eprintln!(
        "cufsm_props_wn: {cases} cases, largest wn difference {worst_wn:e} of the largest |wn|"
    );
}
