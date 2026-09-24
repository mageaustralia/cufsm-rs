#![allow(clippy::needless_range_loop)]
//! Parity with CUFSM itself: every case in `fixtures/cufsm_octave.json` was run through CUFSM's
//! own MATLAB code (unmodified, under GNU Octave - see `oracle/`), and each stage it computed is
//! compared here with this crate's result for the same input:
//!
//! 1. section properties and the reference stresses (`grosprop`, `stresgen`),
//! 2. strip 1's local elastic and geometric stiffness, and both rotated to global axes,
//! 3. the assembled global `K` and `Kg` (models up to 400 DOF),
//! 4. the load factors at every length, against CUFSM's reduced `K` and `Kg` solved with the full
//!    `eig()` (`load_factors_dense`), and the first mode's shape against `stripmain.m`'s own.
//!
//! Why the full `eig()` and not `stripmain.m`'s own load factors: `stripmain.m` calls
//! `eigs(K, Kg, N, 'SM')`. MATLAB's `eigs` shift-inverts on `K` and copes with an indefinite `Kg`;
//! Octave's does not, and on every bending case (part of the section in tension, so `Kg` is
//! indefinite) it returns eigenvalues that are wrong by orders of magnitude and change from run to
//! run. The oracle therefore also solves CUFSM's own reduced problem - its assembly, its constraint
//! basis - with `eig()`, the solver `stripmain.m`'s commented-out path used. `stripmain.m`'s own
//! results are still compared wherever its first load factor agrees with that solve.
//!
//! A failure names the case and the stage, so it points at the function that differs.

mod common;
use common::*;
use cufsm::analysis::{assemble, elemprop, msort, stripmain};
use cufsm::strip::{kglocal, klocal, trans};
use cufsm::{grosprop, stresgen, Actions};

const MATRIX_TOL: f64 = 1e-12;
const STRESS_TOL: f64 = 1e-12;
const MAC_TOL: f64 = 1e-8;

fn cases() -> Vec<serde_json::Value> {
    load("cufsm_octave.json")
}

#[test]
fn section_properties_and_stresses_match_cufsm() {
    for r in cases() {
        let name = r["name"].as_str().unwrap();
        let mut m = model_of(&r);
        let want = &r["props"];
        let p = grosprop(&m);
        for (key, got) in [
            ("A", p.a),
            ("xcg", p.xcg),
            ("zcg", p.zcg),
            ("Ixx", p.ixx),
            ("Izz", p.izz),
            ("Ixz", p.ixz),
            ("I11", p.i11),
            ("I22", p.i22),
            ("thetap", p.thetap),
        ] {
            let w = want[key].as_f64().unwrap();
            assert!(
                (got - w).abs() <= 1e-12 * w.abs().max(p.a),
                "{name}: {key} {got} vs CUFSM {w}"
            );
        }
        let ac = &r["actions"];
        let f = |k: &str| ac[k].as_f64().unwrap();
        let cufsm_stress: Vec<f64> = m.nodes.iter().map(|n| n.stress).collect();
        stresgen(
            &mut m,
            &Actions {
                p: f("P"),
                mxx: f("Mxx"),
                mzz: f("Mzz"),
                m11: f("M11"),
                m22: f("M22"),
            },
            &p,
            f("unsymm") != 0.0,
        );
        let scale = cufsm_stress.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for (i, (n, w)) in m.nodes.iter().zip(&cufsm_stress).enumerate() {
            assert!(
                (n.stress - w).abs() <= STRESS_TOL * scale,
                "{name}: node {} stress {} vs CUFSM {w}",
                i + 1,
                n.stress
            );
        }
    }
}

#[test]
fn strip_matrices_match_cufsm() {
    for r in cases() {
        let name = r["name"].as_str().unwrap();
        let m = model_of(&r);
        let bc = bc_of(&r);
        let a = vec_of(&r["lengths"])[0];
        let m_a = msort(&list_of(&r["m_all"])[0]);
        let e = m.elements[0];
        let (b, alpha) = elemprop(&m)[0];
        let kl = klocal(&m.materials[0], e.t, a, b, bc, &m_a);
        let kgl = kglocal(
            a,
            b,
            m.nodes[e.ni].stress * e.t,
            m.nodes[e.nj].stress * e.t,
            bc,
            &m_a,
        );
        let s = &r["strip1"];
        for (label, got, key) in [
            ("k local", kl.clone(), "k_local"),
            ("kg local", kgl.clone(), "kg_local"),
            ("k global", trans(alpha, &kl), "k_global"),
            ("kg global", trans(alpha, &kgl), "kg_global"),
        ] {
            let d = rel_diff(&got, &rows_of(&s[key]));
            assert!(
                d <= MATRIX_TOL,
                "{name}: strip 1 {label} differs from CUFSM by {d:e} of its largest entry"
            );
        }
    }
}

#[test]
fn global_matrices_match_cufsm() {
    let mut compared = 0;
    for r in cases() {
        if r["K"].as_array().map_or(true, |a| a.is_empty()) {
            continue;
        }
        let name = r["name"].as_str().unwrap();
        let m = model_of(&r);
        let a = vec_of(&r["lengths"])[0];
        let m_a = msort(&list_of(&r["m_all"])[0]);
        let (k, kg) = assemble(&m, a, bc_of(&r), &m_a);
        let dk = rel_diff(&k, &rows_of(&r["K"]));
        let dkg = rel_diff(&kg, &rows_of(&r["Kg"]));
        assert!(
            dk <= MATRIX_TOL,
            "{name}: global K differs from CUFSM by {dk:e}"
        );
        assert!(
            dkg <= MATRIX_TOL,
            "{name}: global Kg differs from CUFSM by {dkg:e}"
        );
        compared += 1;
    }
    assert!(compared > 0, "no case carried global matrices");
}

#[test]
fn load_factors_and_modes_match_cufsm() {
    let mut points = 0;
    let mut modes_compared = 0;
    let mut worst = (0.0_f64, String::new());
    for r in cases() {
        let name = r["name"].as_str().unwrap().to_string();
        let m = model_of(&r);
        let lengths = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        let neigs = r["neigs"].as_u64().unwrap() as usize;
        let got = stripmain(&m, &lengths, &m_all, bc_of(&r), neigs)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let want_lf = list_of(&r["load_factors_dense"]);
        let eigs_lf = list_of(&r["load_factors"]);
        let want_mode = list_of(&r["mode1"]);
        for (l, res) in got.iter().enumerate() {
            let w = &want_lf[l];
            let cond = cond_estimate(&m, lengths[l], bc_of(&r), &m_all[l]);

            let n = w.len().min(res.load_factors.len());
            assert!(
                n > 0,
                "{name}: no positive load factor at length {}",
                lengths[l]
            );
            for i in 0..n {
                let tol = load_factor_tolerance(cond, w[i], w[0]);
                let d = (res.load_factors[i] / w[i] - 1.0).abs();
                if d > worst.0 {
                    worst = (d, format!("{name} length {} mode {}", lengths[l], i + 1));
                }
                assert!(
                    d <= tol,
                    "{name}: length {} mode {}: load factor {} vs CUFSM {} (relative {d:e}, allowed {tol:e})",
                    lengths[l],
                    i + 1,
                    res.load_factors[i],
                    w[i]
                );
                points += 1;
            }
            // A repeated first eigenvalue has no unique first mode; only a distinct one is compared,
            // and only where stripmain.m's eigs() found the right eigenvalue.
            let distinct = w.len() < 2 || (w[1] / w[0] - 1.0).abs() > 1e-6;
            let eigs_right = eigs_lf[l]
                .first()
                .is_some_and(|e| (e / w[0] - 1.0).abs() < 1e-6);
            if distinct && eigs_right && !want_mode[l].is_empty() {
                modes_compared += 1;
                let c = mac(&res.modes[0], &want_mode[l]);
                assert!(
                    1.0 - c <= MAC_TOL,
                    "{name}: length {}: first mode MAC {c} against CUFSM",
                    lengths[l]
                );
            }
        }
    }
    eprintln!(
        "{points} load factors compared with CUFSM; worst relative difference {:e} ({}); {modes_compared} first modes compared",
        worst.0, worst.1
    );
    assert!(
        modes_compared > 100,
        "too few modes compared: {modes_compared}"
    );
}
