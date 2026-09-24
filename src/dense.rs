//! Dense symmetric linear algebra, dependency-free: a row-major matrix, Cholesky, and the
//! symmetric eigenproblem by Householder tridiagonalisation and implicit QL (the EISPACK
//! `tred2`/`tql2` pair, as in the public-domain JAMA package).
//!
//! A finite strip model is small - four DOF per node per longitudinal term, a few hundred in all -
//! so a dense solve of every eigenpair is both exact and fast, and it has none of the convergence
//! choices an iterative solver makes.

/// A dense row-major square matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Mat {
    pub n: usize,
    pub data: Vec<f64>,
}

impl Mat {
    pub fn zeros(n: usize) -> Self {
        Mat {
            n,
            data: vec![0.0; n * n],
        }
    }
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] += v;
    }
    /// `(A + Aᵀ) / 2`.
    pub fn symmetrised(&self) -> Mat {
        let n = self.n;
        let mut s = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                s.set(i, j, 0.5 * (self.get(i, j) + self.get(j, i)));
            }
        }
        s
    }
    /// The largest absolute entry.
    pub fn max_abs(&self) -> f64 {
        self.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
    }
}

/// Lower Cholesky factor `L` with `A = L Lᵀ`, or the column whose pivot was not positive.
pub fn cholesky(a: &Mat) -> Result<Mat, usize> {
    let n = a.n;
    let mut l = Mat::zeros(n);
    for j in 0..n {
        let mut d = a.get(j, j);
        for k in 0..j {
            d -= l.get(j, k) * l.get(j, k);
        }
        if d.is_nan() || d <= 0.0 || !d.is_finite() {
            return Err(j);
        }
        let d = d.sqrt();
        l.set(j, j, d);
        for i in (j + 1)..n {
            let mut s = a.get(i, j);
            for k in 0..j {
                s -= l.get(i, k) * l.get(j, k);
            }
            l.set(i, j, s / d);
        }
    }
    Ok(l)
}

/// Solves `L x = b` in place (forward substitution).
pub fn forward(l: &Mat, b: &mut [f64]) {
    for i in 0..l.n {
        let mut s = b[i];
        for k in 0..i {
            s -= l.get(i, k) * b[k];
        }
        b[i] = s / l.get(i, i);
    }
}

/// Solves `Lᵀ x = b` in place (back substitution).
pub fn backward_t(l: &Mat, b: &mut [f64]) {
    for i in (0..l.n).rev() {
        let mut s = b[i];
        for k in (i + 1)..l.n {
            s -= l.get(k, i) * b[k];
        }
        b[i] = s / l.get(i, i);
    }
}

/// Every eigenpair of a symmetric matrix, eigenvalues ascending. Column `j` of the returned
/// matrix (`v.get(i, j)` over `i`) is the unit eigenvector of eigenvalue `j`.
pub fn sym_eigen(a: &Mat) -> (Vec<f64>, Mat) {
    let n = a.n;
    if n == 0 {
        return (Vec::new(), Mat::zeros(0));
    }
    // V is worked on as a row-major n x n array, V[i][j] = v[i * n + j].
    let mut v = a.symmetrised().data;
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tred2(n, &mut v, &mut d, &mut e);
    // tql2 rotates pairs of columns of V; held transposed, each pair is two contiguous rows.
    let mut vt = transpose(n, &v);
    tql2(n, Some(&mut vt), &mut d, &mut e);
    let v = transpose(n, &vt);
    // Sort ascending, carrying the vectors.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap_or(std::cmp::Ordering::Equal));
    let vals: Vec<f64> = order.iter().map(|&i| d[i]).collect();
    let mut vecs = Mat::zeros(n);
    for (c, &src) in order.iter().enumerate() {
        for r in 0..n {
            vecs.set(r, c, v[r * n + src]);
        }
    }
    (vals, vecs)
}

/// The `k` largest eigenpairs of a symmetric matrix, largest first, and every eigenvalue
/// (ascending). The same Householder tridiagonalisation as [`sym_eigen`], then the eigenvalues
/// alone by implicit QL, then each wanted vector by inverse iteration on the tridiagonal matrix
/// (re-orthogonalised within a cluster of close eigenvalues, as LAPACK's `dstein`) and brought
/// back through the Householder transform. Accumulating all `n` vectors in QL is most of the full
/// solve's cost; this skips it. Every vector is checked against `A x = λ x`; if one misses, the
/// full solve is used instead, so the result is never worse than [`sym_eigen`]'s.
pub fn sym_eigen_top(a: &Mat, k: usize) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let n = a.n;
    let k = k.min(n);
    if n == 0 || k == 0 {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let sym = a.symmetrised();
    let mut q = sym.data.clone();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tred2(n, &mut q, &mut d, &mut e);
    // T: diagonal d, off-diagonal (i-1, i) = e[i].
    let (td, te) = (d.clone(), e.clone());
    tql2(n, None, &mut d, &mut e);
    let mut all = d.clone();
    all.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let top: Vec<f64> = all.iter().rev().take(k).copied().collect();
    let tnorm = (0..n)
        .map(|i| td[i].abs() + te[i].abs() + if i + 1 < n { te[i + 1].abs() } else { 0.0 })
        .fold(0.0, f64::max);
    let cluster = 1e-3 * tnorm;
    let mut ys: Vec<Vec<f64>> = Vec::with_capacity(k);
    for (j, &mu) in top.iter().enumerate() {
        let sigma = mu + tnorm * f64::EPSILON * 10.0 * (1 + j) as f64;
        let mut y: Vec<f64> = (0..n)
            .map(|i| 1.0 + ((i * 7919 + j * 104729) % 1000) as f64 * 1e-3)
            .collect();
        for _ in 0..4 {
            y = tridiag_solve(&td, &te, sigma, &y);
            for (pi, prev) in ys.iter().enumerate() {
                if (top[pi] - mu).abs() <= cluster {
                    let dot: f64 = y.iter().zip(prev).map(|(a, b)| a * b).sum();
                    for (yi, pv) in y.iter_mut().zip(prev) {
                        *yi -= dot * pv;
                    }
                }
            }
            let nrm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
            if nrm == 0.0 || !nrm.is_finite() {
                break;
            }
            for v in &mut y {
                *v /= nrm;
            }
        }
        ys.push(y);
    }
    // x = Q y, then check A x = λ x.
    let anorm = sym.data.iter().fold(0.0_f64, |m, v| m.max(v.abs())) * n as f64;
    let mut xs = Vec::with_capacity(k);
    for (j, y) in ys.iter().enumerate() {
        let x: Vec<f64> = (0..n)
            .map(|r| (0..n).map(|c| q[r * n + c] * y[c]).sum())
            .collect();
        let res = (0..n)
            .map(|r| ((0..n).map(|c| sym.get(r, c) * x[c]).sum::<f64>() - top[j] * x[r]).abs())
            .fold(0.0, f64::max);
        if res.is_nan() || res > 1e-11 * anorm.max(f64::MIN_POSITIVE) {
            let (vals, v) = sym_eigen(a);
            let xs = (0..k)
                .map(|t| (0..n).map(|r| v.get(r, n - 1 - t)).collect())
                .collect();
            let top = (0..k).map(|t| vals[n - 1 - t]).collect();
            return (vals, top, xs);
        }
        xs.push(x);
    }
    (all, top, xs)
}

/// Solves `(T - σ I) x = b` for the symmetric tridiagonal `T` (diagonal `d`, off-diagonal
/// `(i-1, i) = e[i]`) by Gaussian elimination with partial pivoting. A pivot that vanishes is
/// replaced by a tiny one, as inverse iteration wants.
fn tridiag_solve(d: &[f64], e: &[f64], sigma: f64, b: &[f64]) -> Vec<f64> {
    let n = d.len();
    let tiny = f64::EPSILON
        * (d.iter()
            .chain(e.iter())
            .fold(0.0_f64, |m, v| m.max(v.abs()))
            .max(1.0));
    // Rows of the band: [sub, diag, sup1, sup2] after pivoting.
    let mut diag: Vec<f64> = d.iter().map(|v| v - sigma).collect();
    let mut sup1: Vec<f64> = (0..n)
        .map(|i| if i + 1 < n { e[i + 1] } else { 0.0 })
        .collect();
    let mut sup2 = vec![0.0; n];
    let sub: Vec<f64> = (0..n)
        .map(|i| if i + 1 < n { e[i + 1] } else { 0.0 })
        .collect();
    let mut x = b.to_vec();
    for i in 0..n.saturating_sub(1) {
        // Eliminate sub[i] (row i+1, column i).
        if sub[i].abs() > diag[i].abs() {
            // Swap rows i and i+1.
            let (a0, a1, a2) = (diag[i], sup1[i], sup2[i]);
            diag[i] = sub[i];
            sup1[i] = diag[i + 1];
            sup2[i] = sup1[i + 1];
            let f = a0 / diag[i];
            diag[i + 1] = a1 - f * sup1[i];
            sup1[i + 1] = a2 - f * sup2[i];
            x.swap(i, i + 1);
            x[i + 1] -= f * x[i];
        } else {
            if diag[i] == 0.0 {
                diag[i] = tiny;
            }
            let f = sub[i] / diag[i];
            diag[i + 1] -= f * sup1[i];
            sup1[i + 1] -= f * sup2[i];
            x[i + 1] -= f * x[i];
        }
    }
    if diag[n - 1] == 0.0 {
        diag[n - 1] = tiny;
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        if i + 1 < n {
            s -= sup1[i] * x[i + 1];
        }
        if i + 2 < n {
            s -= sup2[i] * x[i + 2];
        }
        x[i] = s / if diag[i] == 0.0 { tiny } else { diag[i] };
    }
    x
}

fn transpose(n: usize, a: &[f64]) -> Vec<f64> {
    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[j * n + i] = a[i * n + j];
        }
    }
    t
}

/// Householder reduction to tridiagonal form, accumulating the transformation in `v`.
fn tred2(n: usize, v: &mut [f64], d: &mut [f64], e: &mut [f64]) {
    let at = |i: usize, j: usize| i * n + j;
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for k in 0..i {
            scale += d[k].abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
                v[at(j, i)] = 0.0;
            }
        } else {
            for k in 0..i {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for j in 0..i {
                e[j] = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[at(j, i)] = f;
                g = e[j] + v[at(j, j)] * f;
                for k in (j + 1)..i {
                    g += v[at(k, j)] * d[k];
                    e[k] += v[at(k, j)] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[at(k, j)] -= f * e[k] + g * d[k];
                }
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..(n - 1) {
        v[at(n - 1, i)] = v[at(i, i)];
        v[at(i, i)] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[at(k, i + 1)] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[at(k, i + 1)] * v[at(k, j)];
                }
                for k in 0..=i {
                    v[at(k, j)] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[at(k, i + 1)] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
        v[at(n - 1, j)] = 0.0;
    }
    v[at(n - 1, n - 1)] = 1.0;
    e[0] = 0.0;
}

/// Implicit QL on the tridiagonal form, accumulating the eigenvectors in `vt` (held transposed).
fn tql2(n: usize, mut vt: Option<&mut [f64]>, d: &mut [f64], e: &mut [f64]) {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1 = 0.0_f64;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            let mut iter = 0;
            loop {
                iter += 1;
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for di in d.iter_mut().take(n).skip(l + 2) {
                    *di -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    // `vt` is Vᵀ: row i of vt is column i of V. Values only: no vectors to rotate.
                    if let Some(vt) = vt.as_deref_mut() {
                        let (lo, hi) = vt.split_at_mut((i + 1) * n);
                        let (ri, ri1) = (&mut lo[i * n..(i + 1) * n], &mut hi[..n]);
                        for k in 0..n {
                            let hk = ri1[k];
                            ri1[k] = s * ri[k] + c * hk;
                            ri[k] = c * ri[k] - s * hk;
                        }
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 || iter > 60 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_eigenpairs_match_the_full_solve() {
        // A symmetric matrix with a repeated top eigenvalue (a block diagonal copy), so the
        // cluster re-orthogonalisation is exercised.
        let m = 5;
        let n = 2 * m;
        let mut a = Mat::zeros(n);
        for blk in 0..2 {
            for i in 0..m {
                for j in 0..m {
                    let v = 1.0 / (1.0 + i as f64 + j as f64) + if i == j { 2.0 } else { 0.0 };
                    a.set(blk * m + i, blk * m + j, v);
                }
            }
        }
        let (full, _) = sym_eigen(&a);
        let (all, top, xs) = sym_eigen_top(&a, 4);
        for (x, y) in all.iter().zip(&full) {
            assert!((x - y).abs() < 1e-12);
        }
        for (j, x) in xs.iter().enumerate() {
            assert!((top[j] - full[n - 1 - j]).abs() < 1e-12);
            for r in 0..n {
                let ax: f64 = (0..n).map(|c| a.get(r, c) * x[c]).sum();
                assert!((ax - top[j] * x[r]).abs() < 1e-11, "pair {j}");
            }
            for (k, y) in xs.iter().enumerate().take(j) {
                let dot: f64 = x.iter().zip(y).map(|(p, q)| p * q).sum();
                assert!(dot.abs() < 1e-8, "pairs {k} and {j} not orthogonal: {dot}");
            }
        }
    }

    fn mat(n: usize, v: &[f64]) -> Mat {
        Mat {
            n,
            data: v.to_vec(),
        }
    }

    #[test]
    fn cholesky_reproduces_the_matrix() {
        let a = mat(
            3,
            &[4.0, 12.0, -16.0, 12.0, 37.0, -43.0, -16.0, -43.0, 98.0],
        );
        let l = cholesky(&a).unwrap();
        assert_eq!(l.data, vec![2.0, 0.0, 0.0, 6.0, 1.0, 0.0, -8.0, 5.0, 3.0]);
        let mut b = vec![1.0, 2.0, 3.0];
        forward(&l, &mut b);
        backward_t(&l, &mut b);
        // A x = [1 2 3]
        for i in 0..3 {
            let r: f64 = (0..3).map(|j| a.get(i, j) * b[j]).sum();
            assert!((r - (i as f64 + 1.0)).abs() < 1e-10);
        }
        assert_eq!(cholesky(&mat(2, &[1.0, 2.0, 2.0, 1.0])), Err(1));
    }

    #[test]
    fn eigenpairs_of_a_known_matrix() {
        // Eigenvalues 1, 2, 4 of [[2,0,0],[0,3,1],[0,1,3]].
        let a = mat(3, &[2.0, 0.0, 0.0, 0.0, 3.0, 1.0, 0.0, 1.0, 3.0]);
        let (vals, vecs) = sym_eigen(&a);
        for (got, want) in vals.iter().zip([2.0, 2.0, 4.0]) {
            assert!((got - want).abs() < 1e-12, "{vals:?}");
        }
        // A v = lambda v for every pair.
        for c in 0..3 {
            for i in 0..3 {
                let av: f64 = (0..3).map(|j| a.get(i, j) * vecs.get(j, c)).sum();
                assert!((av - vals[c] * vecs.get(i, c)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn eigenpairs_of_an_indefinite_matrix() {
        // A 6 x 6 symmetric indefinite matrix: every residual below round-off, vectors orthonormal.
        let n = 6;
        let mut a = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                a.set(
                    i,
                    j,
                    ((i * 7 + j * 3) % 11) as f64 - 5.0 + if i == j { 0.5 * i as f64 } else { 0.0 },
                );
            }
        }
        let a = a.symmetrised();
        let (vals, vecs) = sym_eigen(&a);
        assert!(vals.windows(2).all(|w| w[0] <= w[1]));
        assert!(vals[0] < 0.0 && vals[n - 1] > 0.0);
        for c in 0..n {
            for i in 0..n {
                let av: f64 = (0..n).map(|j| a.get(i, j) * vecs.get(j, c)).sum();
                assert!((av - vals[c] * vecs.get(i, c)).abs() < 1e-11);
            }
            for c2 in 0..n {
                let dot: f64 = (0..n).map(|i| vecs.get(i, c) * vecs.get(i, c2)).sum();
                assert!((dot - if c == c2 { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}
