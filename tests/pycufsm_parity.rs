//! Parity with pyCUFSM, a second, independent implementation (Python, SciPy's QZ `eig`), on every
//! case of the CUFSM fixture it can run, with exactly the inputs CUFSM ran (see
//! `oracle/run_pycufsm.py`).
//!
//! What pyCUFSM cannot be held to, and why (pycufsm 0.x as installed from PyPI):
//!
//! - Fixed DOFs and master-slave constraints: its `constr_user` keeps `r_u_matrix[:, 0:k]` with
//!   `k` the last kept column's 0-based index, dropping that column, and writes each term's block
//!   into an identity it never clears, so the identity's leftover columns free again the fixed DOFs
//!   at the high end of the numbering. A plate simply supported on both long edges comes out as an
//!   outstand (k = 0.44, not 4). CUFSM and this crate give k = 4; so does `tests/theory.rs`.
//! - More than one longitudinal term: it fails padding its mode shapes
//!   ("could not broadcast ... (5,352) into shape (5,44)").
//!
//! Those cases are listed and skipped, not compared.

mod common;

/// pyCUFSM carries noise of its own: over these cases it sits a median 3.8e-9 from CUFSM itself
/// (without this crate involved), and on a lipped channel in bending at a 10 mm half-wave, mode 8,
/// against the same matrices solved in 40-digit arithmetic, this crate is out by 1e-16, CUFSM's
/// eig() by 2e-15 and pyCUFSM by 3.6e-9. So pyCUFSM is allowed this much more.
const PYCUFSM_EXTRA: f64 = 1e-7;
use common::*;
use cufsm::stripmain;

#[test]
fn load_factors_match_pycufsm() {
    let cufsm = load("cufsm_octave.json");
    let py = load("pycufsm.json");
    let mut compared = 0;
    let mut skipped = vec![];
    let mut worst = (0.0_f64, String::new());
    for p in &py {
        let name = p["name"].as_str().unwrap();
        if !p["error"].is_null() || p["fixed_or_constrained"].as_bool().unwrap() {
            skipped.push(name.to_string());
            continue;
        }
        let r = cufsm.iter().find(|r| r["name"] == name).unwrap();
        let m = model_of(r);
        let lengths = vec_of(&r["lengths"]);
        let m_all = list_of(&r["m_all"]);
        let neigs = r["neigs"].as_u64().unwrap() as usize;
        let got = stripmain(&m, &lengths, &m_all, bc_of(r), neigs).unwrap();
        let want = list_of(&p["load_factors"]);
        for (l, res) in got.iter().enumerate() {
            let cond = cond_estimate(&m, lengths[l], bc_of(r), &m_all[l]);
            for (i, w) in want[l].iter().enumerate().take(res.load_factors.len()) {
                let tol = load_factor_tolerance(cond, *w, want[l][0])
                    + PYCUFSM_EXTRA * (w / want[l][0]).max(1.0);
                let d = (res.load_factors[i] / w - 1.0).abs();
                if d > worst.0 {
                    worst = (d, format!("{name} length {} mode {}", lengths[l], i + 1));
                }
                assert!(
                    d <= tol,
                    "{name} length {} mode {}: {} vs pyCUFSM {w} (relative {d:e}, allowed {tol:e})",
                    lengths[l],
                    i + 1,
                    res.load_factors[i]
                );
                compared += 1;
            }
        }
    }
    eprintln!(
        "{compared} load factors compared with pyCUFSM; worst {:e} ({}); not comparable: {skipped:?}",
        worst.0, worst.1
    );
    assert!(compared > 5000, "only {compared} compared");
}
