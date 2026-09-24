//! The AISI Direct Strength Method Design Guide (2006): the CUFSM models behind its worked
//! examples, with the load factors CUFSM computed for them then, shipped in CUFSM's
//! `examples/2006_dsm_design_guide` (MIT). Extracted by `oracle/octave/extract_dsm_guide.m`.
//!
//! Lipped and plain channels, lipped Z, hats, angles, a sigma, a rack upright, a built-up
//! section and a deck panel, in compression and in bending about either axis, some with nodes
//! pinned or constrained, compared mode for mode at every half-wavelength CUFSM analysed where the
//! elastic stiffness is conditioned below 1e8. Beyond that (members of 100 ft and more in these
//! inch models) the 2006 results have lost their digits: on the plain channel in compression they
//! drift 1-2% from this crate at 3081 in and 5415 in, and at 4487 in between them the 2006 value is
//! twice this crate's. A long member's load factor must fall as 1/L², which `tests/theory.rs`
//! holds this crate to.
//!
//! Not compared, each for a reason checked rather than assumed:
//!
//! - `zwlip_mod_webstonly.mat` and `zwlip_mod_nowebst.mat`: the saved curve is not the saved
//!   model's. Today's CUFSM, run on the file's own model under Octave, gives this crate's values
//!   (2.06341 at 1.3 in, 1.93301 at 4.1 in, where the file says 2.02255 and 0.85293).
//! - `cwlip_modified.mat`: likewise - today's CUFSM gives 1.62269 at 1.07 in where the file says
//!   83.33653 - while its sibling `cwlip_modified_P.mat` agrees with CUFSM to every digit.
//! - `panel.mat` (one length, ten modes): today's CUFSM gives this crate's 23.84080591,
//!   30.00707115, 32.77036915 where the file says 23.84076123, 29.99625615, 32.75293268; its
//!   three sibling panel runs agree to 4e-7.
//! - Runs restricted to cFSM spaces: the 2006 cFSM defined its spaces differently (three global
//!   modes where today's has four, and a distortional space that gives the sigma 8% higher).
//!   Parity with today's cFSM is `tests/cfsm_parity.rs`.

mod common;
use common::*;
use cufsm::{stripmain, BoundaryCondition, Constraint, Dof, Element, Material, Model, Node};
use serde_json::Value;

fn model(r: &Value) -> Model {
    let prop = rows_of(&r["prop"]);
    let materials: Vec<(f64, Material)> = prop
        .iter()
        .map(|p| {
            (
                p[0],
                Material {
                    ex: p[1],
                    ey: p[2],
                    vx: p[3],
                    vy: p[4],
                    g: p[5],
                },
            )
        })
        .collect();
    let nodes = rows_of(&r["node"])
        .into_iter()
        .map(|n| Node {
            x: n[1],
            z: n[2],
            free: [n[3] != 0.0, n[4] != 0.0, n[5] != 0.0, n[6] != 0.0],
            stress: n[7],
        })
        .collect();
    let elements = rows_of(&r["elem"])
        .into_iter()
        .map(|e| Element {
            ni: e[1] as usize - 1,
            nj: e[2] as usize - 1,
            t: e[3],
            mat: if e.len() > 4 {
                materials
                    .iter()
                    .position(|(num, _)| *num == e[4])
                    .unwrap_or(0)
            } else {
                0
            },
        })
        .collect();
    let dof = |c: f64| [Dof::X, Dof::Z, Dof::Y, Dof::Theta][c as usize - 1];
    let constraints = match &r["constraints"] {
        Value::Array(a) if !a.is_empty() => rows_of(&r["constraints"])
            .into_iter()
            .map(|c| Constraint {
                node_e: c[0] as usize - 1,
                dof_e: dof(c[1]),
                coeff: c[2],
                node_k: c[3] as usize - 1,
                dof_k: dof(c[4]),
            })
            .collect(),
        _ => vec![],
    };
    Model {
        materials: materials.into_iter().map(|(_, m)| m).collect(),
        nodes,
        elements,
        constraints,
        springs: vec![],
    }
}

#[test]
fn load_factors_match_the_dsm_design_guide_runs() {
    let mut compared = 0;
    let mut skipped: Vec<String> = vec![];
    let mut runs = 0;
    let mut worst = (0.0_f64, String::new());
    for r in load("dsm_guide_2006.json") {
        let name = r["name"].as_str().unwrap().to_string();
        let gbt = vec_of(&r["gbt"]);
        let sizes = if r["gbt_sizes"].is_null() {
            gbt.clone()
        } else {
            vec_of(&r["gbt_sizes"])
        };
        let restricted = gbt.iter().any(|v| *v != 0.0) && gbt != sizes;
        if STALE.contains(&name.as_str()) || restricted {
            skipped.push(name.clone());
            continue;
        }
        let m = model(&r);
        let lengths = vec_of(&r["curve_lengths"]);
        let mut want = rows_of(&r["load_factors"]);
        // One mode saved: Octave wrote the column as a flat array.
        if want.len() == 1 && lengths.len() > 1 {
            want = want[0].iter().map(|v| vec![*v]).collect();
        }
        let nmodes = want[0].len();
        let m_all = vec![vec![1.0]; lengths.len()];
        let got = stripmain(&m, &lengths, &m_all, BoundaryCondition::SS, nmodes);
        let got = got.unwrap_or_else(|e| panic!("{name}: {e}"));
        runs += 1;
        for (l, res) in got.iter().enumerate() {
            let w: Vec<f64> = want[l].iter().copied().filter(|v| *v > 0.0).collect();
            let cond = cond_estimate(&m, lengths[l], BoundaryCondition::SS, &[1.0]);
            if cond > COND_LIMIT {
                continue;
            }
            for (i, (g, wv)) in res.load_factors.iter().zip(&w).enumerate() {
                let tol =
                    load_factor_tolerance(cond, *wv, w[0]) + DSM_GUIDE_EXTRA * (wv / w[0]).max(1.0);
                let d = (g / wv - 1.0).abs();
                if d > worst.0 {
                    worst = (d, format!("{name} length {} mode {}", lengths[l], i + 1));
                }
                assert!(
                    d <= tol,
                    "{name} length {} mode {}: {g} vs CUFSM (2006) {wv} ({d:e}, allowed {tol:e})",
                    lengths[l],
                    i + 1
                );
                compared += 1;
            }
        }
    }
    eprintln!("{compared} load factors over {runs} DSM guide runs; worst {:e} ({}); not compared: {skipped:?}", worst.0, worst.1);
    assert!(
        runs > 40 && compared > 15000,
        "{runs} runs, {compared} load factors"
    );
}

/// The 2006 results carry about five digits (the lipped channel's agree to 5e-6, the sigma's to
/// 5e-7); they are held to 1e-5 on top of the condition-scaled tolerance.
const DSM_GUIDE_EXTRA: f64 = 1e-5;
const COND_LIMIT: f64 = 1e8;
const STALE: [&str; 4] = [
    "zwlip_mod_webstonly.mat",
    "zwlip_mod_nowebst.mat",
    "cwlip_modified.mat",
    "panel.mat",
];
