//! The finite strip analysis, CUFSM `stripmain.m` (without cFSM modal constraints and springs,
//! which are separate stages of the port).
//!
//! For each length: assemble the global elastic and geometric stiffness `K`, `Kg` over every
//! strip and longitudinal term, reduce them to the free DOFs through the constraint basis `R`,
//! solve `K φ = λ Kg φ`, keep the positive load factors smallest first, bring the mode back to
//! every DOF and scale it so its largest entry is `+1` - each step as CUFSM takes it.
//!
//! The one step done differently is the eigen-solve. CUFSM calls MATLAB `eigs(Kff, Kgff, N, 'SM')`;
//! this solves the whole problem densely instead: `K = L Lᵀ`, then every eigenpair of the symmetric
//! `L⁻¹ Kg L⁻ᵀ`, whose eigenvalues are `1 / λ`. On a positive-definite `K` that is the same set of
//! eigenpairs, and it cannot miss one.
//!
//! Accuracy. A global mode at a long half-wavelength strains the membrane almost not at all, so its
//! load factor rests on a near-cancellation between very large membrane terms. Rounding in the
//! assembled matrices alone then moves it by about `eps x cond(K)` - a few parts in 10⁷ for a
//! slender angle at 10 m - in this crate, in CUFSM and in any double-precision code. The tests
//! allow that much and no more (see `tests/cufsm_parity.rs`); local and distortional modes,
//! which are well conditioned, agree with CUFSM to 1e-12.

use crate::dense::{backward_t, cholesky, forward, sym_eigen, Mat};
use crate::model::{BoundaryCondition, Dof, Model};
use crate::strip::{kglocal, klocal, trans};
use crate::Error;

/// The result at one length.
#[derive(Clone, Debug, PartialEq)]
pub struct LengthResult {
    /// The half-wavelength (signature curve) or physical length analysed.
    pub length: f64,
    /// The longitudinal terms used, after CUFSM's clean-up (unique, non-zero, ascending).
    pub m_terms: Vec<f64>,
    /// Positive load factors, smallest first: `λ` times the reference stresses buckles the section.
    pub load_factors: Vec<f64>,
    /// One mode per load factor over every DOF (`4 x nodes x terms`, CUFSM's ordering), scaled so
    /// its largest-magnitude entry is `+1`.
    pub modes: Vec<Vec<f64>>,
}

/// CUFSM `msort.m`: unique, non-zero, ascending.
pub fn msort(m_a: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = m_a.iter().copied().filter(|&m| m != 0.0).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v.dedup();
    v
}

/// The width and angle of every strip, CUFSM `elemprop.m`.
pub fn elemprop(model: &Model) -> Vec<(f64, f64)> {
    model
        .elements
        .iter()
        .map(|e| {
            let (a, b) = (&model.nodes[e.ni], &model.nodes[e.nj]);
            let (dx, dz) = (b.x - a.x, b.z - a.z);
            (dx.hypot(dz), dz.atan2(dx))
        })
        .collect()
}

/// Global DOF index of local strip DOF `l` (0..8, `[u1 v1 u2 v2 w1 th1 w2 th2]`) for longitudinal
/// term `m`, CUFSM `assemble.m`.
#[inline]
fn gdof(l: usize, m: usize, ni: usize, nj: usize, nnodes: usize) -> usize {
    let base = 4 * nnodes * m;
    let skip = 2 * nnodes;
    match l {
        0 => base + 2 * ni,
        1 => base + 2 * ni + 1,
        2 => base + 2 * nj,
        3 => base + 2 * nj + 1,
        4 => base + skip + 2 * ni,
        5 => base + skip + 2 * ni + 1,
        6 => base + skip + 2 * nj,
        _ => base + skip + 2 * nj + 1,
    }
}

/// Assembles the global `K` and `Kg` at length `a`.
pub fn assemble(model: &Model, a: f64, bc: BoundaryCondition, m_a: &[f64]) -> (Mat, Mat) {
    let nnodes = model.nodes.len();
    let tm = m_a.len();
    let mut k = Mat::zeros(4 * nnodes * tm);
    let mut kg = Mat::zeros(4 * nnodes * tm);
    for (e, (b, alpha)) in model.elements.iter().zip(elemprop(model)) {
        let mat = &model.materials[e.mat];
        let ty1 = model.nodes[e.ni].stress * e.t;
        let ty2 = model.nodes[e.nj].stress * e.t;
        let kl = trans(alpha, &klocal(mat, e.t, a, b, bc, m_a));
        let kgl = trans(alpha, &kglocal(a, b, ty1, ty2, bc, m_a));
        for m in 0..tm {
            for p in 0..tm {
                for r in 0..8 {
                    let gr = gdof(r, m, e.ni, e.nj, nnodes);
                    for c in 0..8 {
                        let gc = gdof(c, p, e.ni, e.nj, nnodes);
                        k.add(gr, gc, kl.get(8 * m + r, 8 * p + c));
                        kg.add(gr, gc, kgl.get(8 * m + r, 8 * p + c));
                    }
                }
            }
        }
    }
    (k, kg)
}

/// Index of `dof` at `node` within one longitudinal term's `4 x nnodes` block, CUFSM
/// `constr_user.m`: `u` and `v` interleaved first, then `w` and `theta`.
fn dof_index(node: usize, dof: Dof, nnodes: usize) -> usize {
    match dof {
        Dof::X => 2 * node,
        Dof::Y => 2 * node + 1,
        Dof::Z => 2 * nnodes + 2 * node,
        Dof::Theta => 2 * nnodes + 2 * node + 1,
    }
}

/// The constraint basis `R` (every DOF = `R` x free DOF), CUFSM `constr_user.m`, as a list of
/// sparse columns `(row, coefficient)`. `None` when nothing is fixed or constrained.
///
/// CUFSM goes on to replace `R` by an orthonormal basis of the same space (`null(null(R')')`).
/// The load factors do not depend on the basis, and the modes, brought back to every DOF and
/// normalised, do not either, so the basis is used as built.
pub fn constraint_basis(model: &Model, tm: usize) -> Option<Vec<Vec<(usize, f64)>>> {
    let nnodes = model.nodes.len();
    let ndof = 4 * nnodes;
    let fixed = model.nodes.iter().any(|n| n.free.iter().any(|f| !f));
    if !fixed && model.constraints.is_empty() {
        return None;
    }
    // One term's block: identity columns, then each constraint folds its eliminated column
    // into the kept one, then the fixed and eliminated columns go.
    let mut cols: Vec<Vec<(usize, f64)>> = (0..ndof).map(|i| vec![(i, 1.0)]).collect();
    let mut keep = vec![true; ndof];
    for (i, n) in model.nodes.iter().enumerate() {
        // CUFSM's column order: x, z, y, rotation.
        for (flag, dof) in n.free.iter().zip([Dof::X, Dof::Z, Dof::Y, Dof::Theta]) {
            if !flag {
                keep[dof_index(i, dof, nnodes)] = false;
            }
        }
    }
    for c in &model.constraints {
        let e = dof_index(c.node_e, c.dof_e, nnodes);
        let k = dof_index(c.node_k, c.dof_k, nnodes);
        let add: Vec<(usize, f64)> = cols[e].iter().map(|&(r, v)| (r, v * c.coeff)).collect();
        for (r, v) in add {
            match cols[k].iter_mut().find(|(rr, _)| *rr == r) {
                Some(slot) => slot.1 += v,
                None => cols[k].push((r, v)),
            }
        }
        keep[e] = false;
    }
    let block: Vec<Vec<(usize, f64)>> = cols
        .into_iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(c, _)| c)
        .collect();
    let mut all = Vec::with_capacity(block.len() * tm);
    for m in 0..tm {
        for col in &block {
            all.push(col.iter().map(|&(r, v)| (r + m * ndof, v)).collect());
        }
    }
    Some(all)
}

/// `Rᵀ A R`.
pub fn reduce(a: &Mat, r: &[Vec<(usize, f64)>]) -> Mat {
    let nf = r.len();
    let n = a.n;
    // A R, column by column.
    let mut ar = vec![0.0; n * nf];
    for (j, col) in r.iter().enumerate() {
        for &(s, v) in col {
            for i in 0..n {
                ar[i * nf + j] += a.get(i, s) * v;
            }
        }
    }
    let mut out = Mat::zeros(nf);
    for (i, col) in r.iter().enumerate() {
        for j in 0..nf {
            let mut sum = 0.0;
            for &(s, v) in col {
                sum += v * ar[s * nf + j];
            }
            out.set(i, j, sum);
        }
    }
    out
}

/// Solves `K φ = λ Kg φ` for the `neigs` smallest positive `λ`. `K` must be positive definite.
pub fn buckling_eigen(k: &Mat, kg: &Mat, neigs: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), Error> {
    let n = k.n;
    let k = k.symmetrised();
    let kg = kg.symmetrised();
    let l = cholesky(&k).map_err(|dof| Error::NotPositiveDefinite { dof })?;
    // C = L⁻¹ Kg L⁻ᵀ: X = L⁻¹ Kg column by column, then C = L⁻¹ Xᵀ (Kg is symmetric).
    let mut x = Mat::zeros(n);
    let mut col = vec![0.0; n];
    for j in 0..n {
        for i in 0..n {
            col[i] = kg.get(i, j);
        }
        forward(&l, &mut col);
        for i in 0..n {
            x.set(i, j, col[i]);
        }
    }
    let mut c = Mat::zeros(n);
    for j in 0..n {
        for i in 0..n {
            col[i] = x.get(j, i);
        }
        forward(&l, &mut col);
        for i in 0..n {
            c.set(i, j, col[i]);
        }
    }
    let (mu, y) = sym_eigen(&c);
    // λ = 1/μ; the positive λ are the positive μ, the smallest λ the largest μ. A μ at round-off
    // relative to the largest is not a buckling mode but the numerical null of Kg.
    let mu_max = mu.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut picked: Vec<usize> = (0..n).filter(|&i| mu[i] > mu_max * 1e-14).collect();
    picked.sort_by(|&a, &b| {
        mu[b]
            .partial_cmp(&mu[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    picked.truncate(neigs);
    let mut lfs = Vec::with_capacity(picked.len());
    let mut modes = Vec::with_capacity(picked.len());
    for i in picked {
        lfs.push(1.0 / mu[i]);
        let mut phi: Vec<f64> = (0..n).map(|r| y.get(r, i)).collect();
        backward_t(&l, &mut phi);
        modes.push(phi);
    }
    Ok((lfs, modes))
}

/// Scales a mode so its largest-magnitude entry (the first, on a tie) is `+1`.
pub fn normalise_mode(mode: &mut [f64]) {
    let mut best = 0;
    for (i, v) in mode.iter().enumerate() {
        if v.abs() > mode[best].abs() {
            best = i;
        }
    }
    let s = mode[best];
    if s != 0.0 {
        for v in mode.iter_mut() {
            *v /= s;
        }
    }
}

/// A finite strip analysis at every length, CUFSM `stripmain.m`.
///
/// `m_all[i]` holds the longitudinal terms for `lengths[i]` (`[1.0]` for the signature curve).
/// `neigs` is the number of load factors kept per length (CUFSM's default is 20).
pub fn stripmain(
    model: &Model,
    lengths: &[f64],
    m_all: &[Vec<f64>],
    bc: BoundaryCondition,
    neigs: usize,
) -> Result<Vec<LengthResult>, Error> {
    model.validate()?;
    if lengths.len() != m_all.len() {
        return Err(Error::InvalidModel(format!(
            "{} lengths but {} sets of longitudinal terms",
            lengths.len(),
            m_all.len()
        )));
    }
    let mut out = Vec::with_capacity(lengths.len());
    for (&a, m_raw) in lengths.iter().zip(m_all) {
        if a.is_nan() || a <= 0.0 {
            return Err(Error::InvalidModel(format!("length {a} is not positive")));
        }
        let m_a = msort(m_raw);
        if m_a.is_empty() {
            return Err(Error::InvalidModel(format!(
                "no longitudinal terms at length {a}"
            )));
        }
        let (k, kg) = assemble(model, a, bc, &m_a);
        let r = constraint_basis(model, m_a.len());
        let (kff, kgff) = match &r {
            Some(r) => (reduce(&k, r), reduce(&kg, r)),
            None => (k, kg),
        };
        let (lfs, reduced) = buckling_eigen(&kff, &kgff, neigs)?;
        let modes = reduced
            .into_iter()
            .map(|q| {
                let mut full = match &r {
                    Some(r) => {
                        let mut f = vec![0.0; 4 * model.nodes.len() * m_a.len()];
                        for (j, col) in r.iter().enumerate() {
                            for &(row, v) in col {
                                f[row] += v * q[j];
                            }
                        }
                        f
                    }
                    None => q,
                };
                normalise_mode(&mut full);
                full
            })
            .collect();
        out.push(LengthResult {
            length: a,
            m_terms: m_a,
            load_factors: lfs,
            modes,
        });
    }
    Ok(out)
}

/// The signature curve, CUFSM `signature_ss.m`: simply supported, one longitudinal term, at 100
/// half-wavelengths spaced logarithmically from the narrowest strip's width to 1000 times the
/// widest's. Returns the lowest load factor at each.
pub fn signature_ss(model: &Model, neigs: usize) -> Result<Vec<LengthResult>, Error> {
    model.validate()?;
    let widths: Vec<f64> = elemprop(model).into_iter().map(|(w, _)| w).collect();
    let mm = widths.iter().copied().fold(f64::INFINITY, f64::min);
    let mw = widths.iter().copied().fold(0.0, f64::max);
    let lengths = logspace(mm.log10(), (1000.0 * mw).log10(), 100);
    let m_all = vec![vec![1.0]; lengths.len()];
    stripmain(model, &lengths, &m_all, BoundaryCondition::SS, neigs)
}

/// MATLAB `logspace(a, b, n)`.
pub fn logspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    if n == 1 {
        return vec![10f64.powf(b)];
    }
    (0..n)
        .map(|i| 10f64.powf(a + (b - a) * i as f64 / (n - 1) as f64))
        .collect()
}
