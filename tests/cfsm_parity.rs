//! cFSM against CUFSM's own `analysis/cFSM/` code, run under Octave (see `oracle/`), on six
//! sections in compression: sharp and rounded lipped C, lipped Z, plain channel, hat, and a
//! branched I-section - and the lipped C again with a fixed DOF, a constraint and springs, for
//! the restricted analysis's intersection with the constrained space.
//!
//! Several cFSM steps take `null()` or `eig()`, so individual base vectors are not unique. What is
//! compared is what does not depend on the basis:
//!
//! 1. `cutwp_prop2`'s properties, including the shear centre and the warping function,
//! 2. the sizes of the G, D, L and O spaces,
//! 3. the spaces themselves (every CUFSM base vector lies in this crate's space and vice versa),
//! 4. the load factors of analyses restricted to G, D or L alone,
//! 5. the G/D/L/O classification of the unconstrained modes, with CUFSM's default settings.

mod common;
use common::*;
use cufsm::cfsm::{base_column, classify, stripmain_constrained, Norm, Orth, Spaces};
use cufsm::cutwp::cutwp_prop2;
use cufsm::linalg::{solve, RMat};
use cufsm::stripmain;

fn cfsm_cases() -> Vec<serde_json::Value> {
    load("cufsm_octave.json")
        .into_iter()
        .filter(|r| r["cfsm"].get("ngm").is_some())
        .collect()
}

fn rmat(rows: Vec<Vec<f64>>) -> RMat {
    let (r, c) = (rows.len(), rows[0].len());
    RMat {
        r,
        c,
        data: rows.into_iter().flatten().collect(),
    }
}

/// Largest relative residual of projecting each column of `b` onto the span of `a`.
fn outside_span(a: &RMat, b: &RMat) -> f64 {
    if b.c == 0 {
        return 0.0;
    }
    let ata = a.t().mul(a);
    let x = solve(&ata, &a.t().mul(b)).expect("independent columns");
    let proj = a.mul(&x);
    (0..b.c)
        .map(|j| {
            let col = b.col(j);
            let n = col.iter().map(|v| v * v).sum::<f64>().sqrt();
            let r = (0..b.r)
                .map(|i| (col[i] - proj.get(i, j)).powi(2))
                .sum::<f64>()
                .sqrt();
            r / n
        })
        .fold(0.0, f64::max)
}

#[test]
fn cutwp_properties_match_cufsm() {
    let cases = cfsm_cases();
    assert!(cases.len() >= 6);
    for r in &cases {
        let name = r["name"].as_str().unwrap();
        let p = cutwp_prop2(&model_of(r));
        let w = &r["cfsm"]["cutwp"];
        let f = |k: &str| w[k].as_f64().unwrap();
        let scale_len = f("A").sqrt();
        for (k, got, scale) in [
            ("A", p.a, f("A")),
            ("xc", p.xc, scale_len),
            ("zc", p.zc, scale_len),
            ("Ix", p.ix, f("Ix").abs().max(f("Iz").abs())),
            ("Iz", p.iz, f("Ix").abs().max(f("Iz").abs())),
            ("Ixz", p.ixz, f("Ix").abs().max(f("Iz").abs())),
            ("theta", p.theta, 1.0),
            ("I1", p.i1, f("I1")),
            ("I2", p.i2, f("I1")),
            ("J", p.j, f("J")),
            ("xs", p.xs, scale_len),
            ("zs", p.zs, scale_len),
            ("Cw", p.cw, f("Cw").abs().max(1.0)),
            ("B1", p.b1, scale_len),
            ("B2", p.b2, scale_len),
        ] {
            assert!(
                (got - f(k)).abs() <= 1e-10 * scale,
                "{name}: {k} {got} vs CUFSM {}",
                f(k)
            );
        }
        let wn = vec_of(&w["wn"]);
        let wscale = wn.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1.0);
        for (i, (a, b)) in p.wn.iter().zip(&wn).enumerate() {
            assert!(
                (a - b).abs() <= 1e-10 * wscale,
                "{name}: wn at node {} {a} vs CUFSM {b}",
                i + 1
            );
        }
    }
}

#[test]
fn modal_spaces_match_cufsm() {
    for r in &cfsm_cases() {
        let name = r["name"].as_str().unwrap();
        let m = model_of(r);
        let a = vec_of(&r["lengths"])[0];
        let cf = &r["cfsm"];
        let (bv, ngm, ndm, nlm) = base_column(&m, a, bc_of(r), &[1.0]).unwrap();
        let want = (
            cf["ngm"].as_u64().unwrap() as usize,
            cf["ndm"].as_u64().unwrap() as usize,
            cf["nlm"].as_u64().unwrap() as usize,
        );
        assert_eq!((ngm, ndm, nlm), want, "{name}: space sizes");
        let theirs = rmat(rows_of(&cf["b_v_l"]));
        let ndof = bv.r;
        for (label, g0, g1) in [
            ("G", 0, ngm),
            ("D", ngm, ngm + ndm),
            ("L", ngm + ndm, ngm + ndm + nlm),
            ("O", ngm + ndm + nlm, ndof),
        ] {
            if g1 == g0 {
                continue;
            }
            // O's columns past the strips' count are zero in both: compare the nonzero ones.
            let nonzero = |b: &RMat| -> RMat {
                let keep: Vec<usize> = (0..b.c)
                    .filter(|&j| b.col(j).iter().any(|v| *v != 0.0))
                    .collect();
                let mut o = RMat::zeros(b.r, keep.len());
                for (k, &j) in keep.iter().enumerate() {
                    o.set_col(k, &b.col(j));
                }
                o
            };
            let (ours, cufsm) = (nonzero(&bv.cols(g0, g1)), nonzero(&theirs.cols(g0, g1)));
            assert_eq!(ours.c, cufsm.c, "{name}: {label} space dimension");
            let d1 = outside_span(&ours, &cufsm);
            let d2 = outside_span(&cufsm, &ours);
            assert!(
                d1 < 1e-8 && d2 < 1e-8,
                "{name}: {label} spaces differ ({d1:e}, {d2:e})"
            );
        }
    }
}

/// The condition estimate of the restricted problem's `Rᵀ K R`, as `cond_estimate` for the full.
fn restricted_cond(m: &cufsm::Model, a: f64, bc: cufsm::BoundaryCondition, spaces: Spaces) -> f64 {
    use cufsm::cfsm::mode_select;
    let (bv, ngm, ndm, nlm) = base_column(m, a, bc, &[1.0]).unwrap();
    let r = mode_select(&bv, ngm, ndm, nlm, spaces, 4 * m.nodes.len(), 1);
    if r.c == 0 {
        return 1.0;
    }
    let (k, _) = cufsm::analysis::assemble(m, a, bc, &[1.0]);
    let kff = r.t().mul(&RMat::from_square(&k)).mul(&r).to_square();
    let l = cufsm::dense::cholesky(&kff.symmetrised()).unwrap();
    let piv: Vec<f64> = (0..l.n).map(|j| l.get(j, j)).collect();
    (piv.iter().cloned().fold(0.0, f64::max) / piv.iter().cloned().fold(f64::INFINITY, f64::min))
        .powi(2)
}

#[test]
fn restricted_load_factors_match_cufsm() {
    let mut compared = 0;
    for r in &cfsm_cases() {
        let name = r["name"].as_str().unwrap();
        let m = model_of(r);
        let lengths = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        for (key, spaces) in [
            (
                "lf_G",
                Spaces {
                    global: true,
                    ..Default::default()
                },
            ),
            (
                "lf_D",
                Spaces {
                    distortional: true,
                    ..Default::default()
                },
            ),
            (
                "lf_L",
                Spaces {
                    local: true,
                    ..Default::default()
                },
            ),
        ] {
            let got = stripmain_constrained(&m, &lengths, &m_all, bc_of(r), 5, spaces).unwrap();
            let want = list_of(&r["cfsm"][key]);
            for (l, res) in got.iter().enumerate() {
                let w = &want[l];
                assert_eq!(
                    res.load_factors.is_empty(),
                    w.is_empty(),
                    "{name} {key} length {}",
                    lengths[l]
                );
                let cond = restricted_cond(&m, lengths[l], bc_of(r), spaces);
                for (i, (g, wv)) in res.load_factors.iter().zip(w).enumerate() {
                    let tol = load_factor_tolerance(cond, *wv, w[0]);
                    let d = (g / wv - 1.0).abs();
                    assert!(d <= tol, "{name} {key} length {} mode {}: {g} vs CUFSM {wv} ({d:e}, allowed {tol:e})", lengths[l], i + 1);
                    compared += 1;
                }
            }
        }
    }
    assert!(compared > 150, "only {compared}");
}

/// Classification with CUFSM's defaults (axial orthogonality, vector norm) and with the natural
/// basis. The axial basis comes out of `eig()` on each space; where those subproblems have repeated
/// eigenvalues the basis, and the vector-normalised percentages with it, are not unique - in
/// CUFSM either. The doubly symmetric I-section is such a case (at 60 mm: 4 repeated pairs in L,
/// 2 and 4 in the two halves of O), so it is held to the natural basis only, which it has
/// uniquely (no distortional space).
fn check_classification(orth: Orth, key: &str, skip_degenerate: bool) -> usize {
    let mut compared = 0;
    for r in &cfsm_cases() {
        let name = r["name"].as_str().unwrap();
        if skip_degenerate && name.contains("i-section") {
            continue;
        }
        // The oracle classifies the bare strips' modes; a case with fixities, constraints or
        // springs is covered by the restricted load factors instead.
        if name.contains("fixed and sprung") {
            continue;
        }
        let m = model_of(r);
        let lengths = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        let res = stripmain(&m, &lengths, &m_all, bc_of(r), 3).unwrap();
        let got = classify(&m, &res, bc_of(r), orth, Norm::Vector).unwrap();
        let want = &r["cfsm"][key];
        for (l, per) in got.iter().enumerate() {
            let wrows = rows_of(&want[l]);
            let lfs = &res[l].load_factors;
            for (q, (g, w)) in per.iter().zip(&wrows).enumerate() {
                // A mode with a near-repeated load factor has no unique shape.
                let distinct = (q == 0 || (lfs[q] / lfs[q - 1] - 1.0).abs() > 1e-6)
                    && (q + 1 >= lfs.len() || (lfs[q + 1] / lfs[q] - 1.0).abs() > 1e-6);
                if !distinct {
                    continue;
                }
                for k in 0..4 {
                    assert!(
                        (g[k] - w[k]).abs() < 1e-6,
                        "{name} {key} length {} mode {}: {g:?} vs CUFSM {w:?}",
                        lengths[l],
                        q + 1
                    );
                }
                compared += 1;
            }
        }
    }
    compared
}

#[test]
fn classification_matches_cufsm() {
    let axial = check_classification(Orth::Axial, "classification", true);
    let natural = check_classification(Orth::Natural, "classification_natural", false);
    assert!(
        axial > 40 && natural > 50,
        "only {axial} and {natural} compared"
    );
}
