//! Rectangular dense linear algebra for cFSM: products, a pivoted LU solve, the singular value
//! decomposition by one-sided Jacobi (for MATLAB's `null`), and the symmetric-definite generalized
//! eigenproblem (for MATLAB's `eig(A, B)`).
//!
//! One-sided Jacobi rather than an SVD through `AᵀA`: squaring the matrix squares its condition,
//! and `null` needs the small singular values told apart from zero at `max(m, n) · eps · σmax`,
//! below where `AᵀA` can resolve them.

use crate::dense::{backward_t, cholesky, forward, sym_eigen, Mat};

/// A dense row-major `r x c` matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct RMat {
    pub r: usize,
    pub c: usize,
    pub data: Vec<f64>,
}

impl RMat {
    pub fn zeros(r: usize, c: usize) -> Self {
        RMat {
            r,
            c,
            data: vec![0.0; r * c],
        }
    }
    pub fn eye(n: usize) -> Self {
        let mut m = RMat::zeros(n, n);
        for i in 0..n {
            m.set(i, i, 1.0);
        }
        m
    }
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.c + j]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.c + j] = v;
    }
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.c + j] += v;
    }
    pub fn from_square(m: &Mat) -> Self {
        RMat {
            r: m.n,
            c: m.n,
            data: m.data.clone(),
        }
    }
    pub fn to_square(&self) -> Mat {
        assert_eq!(self.r, self.c);
        Mat {
            n: self.r,
            data: self.data.clone(),
        }
    }
    pub fn t(&self) -> RMat {
        let mut o = RMat::zeros(self.c, self.r);
        for i in 0..self.r {
            for j in 0..self.c {
                o.set(j, i, self.get(i, j));
            }
        }
        o
    }
    pub fn mul(&self, b: &RMat) -> RMat {
        assert_eq!(self.c, b.r, "inner dimensions");
        let mut o = RMat::zeros(self.r, b.c);
        for i in 0..self.r {
            for k in 0..self.c {
                let a = self.get(i, k);
                if a == 0.0 {
                    continue;
                }
                for j in 0..b.c {
                    o.data[i * b.c + j] += a * b.data[k * b.c + j];
                }
            }
        }
        o
    }
    pub fn scale(&self, f: f64) -> RMat {
        RMat {
            r: self.r,
            c: self.c,
            data: self.data.iter().map(|v| v * f).collect(),
        }
    }
    /// Columns `from..to`.
    pub fn cols(&self, from: usize, to: usize) -> RMat {
        let mut o = RMat::zeros(self.r, to - from);
        for i in 0..self.r {
            for j in from..to {
                o.set(i, j - from, self.get(i, j));
            }
        }
        o
    }
    /// Rows `from..to`.
    pub fn rows(&self, from: usize, to: usize) -> RMat {
        RMat {
            r: to - from,
            c: self.c,
            data: self.data[from * self.c..to * self.c].to_vec(),
        }
    }
    /// Writes `b` with its top-left at `(i0, j0)`.
    pub fn put(&mut self, i0: usize, j0: usize, b: &RMat) {
        for i in 0..b.r {
            for j in 0..b.c {
                self.set(i0 + i, j0 + j, b.get(i, j));
            }
        }
    }
    pub fn col(&self, j: usize) -> Vec<f64> {
        (0..self.r).map(|i| self.get(i, j)).collect()
    }
    pub fn set_col(&mut self, j: usize, v: &[f64]) {
        for (i, x) in v.iter().enumerate() {
            self.set(i, j, *x);
        }
    }
    /// Stacks `b`'s columns after this one's.
    pub fn hcat(&self, b: &RMat) -> RMat {
        assert_eq!(self.r, b.r);
        let mut o = RMat::zeros(self.r, self.c + b.c);
        o.put(0, 0, self);
        o.put(0, self.c, b);
        o
    }
}

/// `A \ B` by LU with partial pivoting (A square). `None` if `A` is singular to working precision.
pub fn solve(a: &RMat, b: &RMat) -> Option<RMat> {
    let n = a.r;
    assert_eq!(a.c, n);
    assert_eq!(b.r, n);
    let mut lu = a.clone();
    let mut x = b.clone();
    let scale = a.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for k in 0..n {
        let (mut p, mut best) = (k, lu.get(k, k).abs());
        for i in (k + 1)..n {
            if lu.get(i, k).abs() > best {
                best = lu.get(i, k).abs();
                p = i;
            }
        }
        if best <= f64::EPSILON * scale * n as f64 {
            return None;
        }
        if p != k {
            for j in 0..n {
                let t = lu.get(k, j);
                lu.set(k, j, lu.get(p, j));
                lu.set(p, j, t);
            }
            for j in 0..x.c {
                let t = x.get(k, j);
                x.set(k, j, x.get(p, j));
                x.set(p, j, t);
            }
        }
        let piv = lu.get(k, k);
        for i in (k + 1)..n {
            let f = lu.get(i, k) / piv;
            if f == 0.0 {
                continue;
            }
            lu.set(i, k, f);
            for j in (k + 1)..n {
                lu.add(i, j, -f * lu.get(k, j));
            }
            for j in 0..x.c {
                x.add(i, j, -f * x.get(k, j));
            }
        }
    }
    for j in 0..x.c {
        for i in (0..n).rev() {
            let mut s = x.get(i, j);
            for k in (i + 1)..n {
                s -= lu.get(i, k) * x.get(k, j);
            }
            x.set(i, j, s / lu.get(i, i));
        }
    }
    Some(x)
}

/// The singular values of `A` (`m x n`) and its right singular vectors, by one-sided Jacobi on
/// the columns: `A V = U Σ`. Returns `(σ, V)`, `σ` of length `n` (zeros for a wide matrix's extra
/// columns), not sorted.
pub fn svd_right(a: &RMat) -> (Vec<f64>, RMat) {
    let (m, n) = (a.r, a.c);
    // Work on columns: g[j] is column j of A V.
    let mut g: Vec<Vec<f64>> = (0..n).map(|j| a.col(j)).collect();
    let mut v = RMat::eye(n);
    let tol = f64::EPSILON;
    for _sweep in 0..80 {
        let mut rotated = false;
        for p in 0..n {
            for q in (p + 1)..n {
                let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
                for i in 0..m {
                    alpha += g[p][i] * g[p][i];
                    beta += g[q][i] * g[q][i];
                    gamma += g[p][i] * g[q][i];
                }
                if gamma == 0.0 || gamma.abs() <= tol * (alpha * beta).sqrt() {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let t = if zeta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let (x, y) = (g[p][i], g[q][i]);
                    g[p][i] = c * x - s * y;
                    g[q][i] = s * x + c * y;
                }
                for i in 0..n {
                    let (x, y) = (v.get(i, p), v.get(i, q));
                    v.set(i, p, c * x - s * y);
                    v.set(i, q, s * x + c * y);
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let sigma: Vec<f64> = g
        .iter()
        .map(|col| col.iter().map(|x| x * x).sum::<f64>().sqrt())
        .collect();
    (sigma, v)
}

/// An orthonormal basis for the null space of `A`, as MATLAB's `null(A)`: the right singular
/// vectors whose singular value is at most `max(m, n) · eps · σmax`. The basis is not unique;
/// only the space it spans is compared anywhere.
pub fn null(a: &RMat) -> RMat {
    let n = a.c;
    if a.r == 0 {
        return RMat::eye(n);
    }
    let (sigma, v) = svd_right(a);
    let smax = sigma.iter().cloned().fold(0.0, f64::max);
    let tol = (a.r.max(n) as f64) * f64::EPSILON * smax;
    let keep: Vec<usize> = (0..n).filter(|&j| sigma[j] <= tol).collect();
    let mut o = RMat::zeros(n, keep.len());
    for (k, &j) in keep.iter().enumerate() {
        o.set_col(k, &v.col(j));
    }
    o
}

/// Every eigenpair of `A x = λ B x` for symmetric `A` and symmetric positive-definite `B`, as
/// MATLAB's `eig(A, B)` on such a pair: eigenvalues ascending, vectors normalised so that
/// `Vᵀ B V = I`. `None` if `B` is not positive definite.
pub fn eig_sym_def(a: &RMat, b: &RMat) -> Option<(Vec<f64>, RMat)> {
    let n = a.r;
    let bm = b.to_square().symmetrised();
    let l = cholesky(&bm).ok()?;
    let am = a.to_square().symmetrised();
    // C = L⁻¹ A L⁻ᵀ
    let mut x = Mat::zeros(n);
    let mut col = vec![0.0; n];
    for j in 0..n {
        for i in 0..n {
            col[i] = am.get(i, j);
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
    let (vals, y) = sym_eigen(&c);
    let mut v = RMat::zeros(n, n);
    for j in 0..n {
        let mut phi: Vec<f64> = (0..n).map(|r| y.get(r, j)).collect();
        backward_t(&l, &mut phi);
        v.set_col(j, &phi);
    }
    Some((vals, v))
}

/// Upper Cholesky factor `R` with `A = Rᵀ R`, as MATLAB's `chol(A)`.
pub fn chol_upper(a: &RMat) -> Option<RMat> {
    let l = cholesky(&a.to_square().symmetrised()).ok()?;
    Some(RMat::from_square(&l).t())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_space_of_a_rank_deficient_matrix() {
        // Rows 1 and 2 span a plane in R⁴; the third is their sum. Null space dimension 2.
        let a = RMat {
            r: 3,
            c: 4,
            data: vec![1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 3.0, 1.0, 1.0],
        };
        let z = null(&a);
        assert_eq!((z.r, z.c), (4, 2));
        let az = a.mul(&z);
        assert!(az.data.iter().all(|v| v.abs() < 1e-14), "{az:?}");
        let ztz = z.t().mul(&z);
        for i in 0..2 {
            for j in 0..2 {
                assert!((ztz.get(i, j) - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn lu_solve() {
        let a = RMat {
            r: 3,
            c: 3,
            data: vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0],
        };
        let b = RMat {
            r: 3,
            c: 1,
            data: vec![3.0, 2.0, 4.0],
        };
        let x = solve(&a, &b).unwrap();
        let r = a.mul(&x);
        for i in 0..3 {
            assert!((r.get(i, 0) - b.get(i, 0)).abs() < 1e-14);
        }
    }

    #[test]
    fn generalized_eigenpairs_are_b_normalised() {
        let a = RMat {
            r: 2,
            c: 2,
            data: vec![2.0, 1.0, 1.0, 3.0],
        };
        let b = RMat {
            r: 2,
            c: 2,
            data: vec![4.0, 1.0, 1.0, 2.0],
        };
        let (vals, v) = eig_sym_def(&a, &b).unwrap();
        let vbv = v.t().mul(&b).mul(&v);
        assert!((vbv.get(0, 0) - 1.0).abs() < 1e-13 && vbv.get(0, 1).abs() < 1e-13);
        for j in 0..2 {
            let x = v.col(j);
            for i in 0..2 {
                let ax: f64 = (0..2).map(|k| a.get(i, k) * x[k]).sum();
                let bx: f64 = (0..2).map(|k| b.get(i, k) * x[k]).sum();
                assert!((ax - vals[j] * bx).abs() < 1e-13);
            }
        }
    }
}
