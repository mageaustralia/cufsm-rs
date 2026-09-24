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
    tql2(n, &mut vt, &mut d, &mut e);
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
fn tql2(n: usize, vt: &mut [f64], d: &mut [f64], e: &mut [f64]) {
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
                    // `vt` is Vᵀ: row i of vt is column i of V.
                    let (lo, hi) = vt.split_at_mut((i + 1) * n);
                    let (ri, ri1) = (&mut lo[i * n..(i + 1) * n], &mut hi[..n]);
                    for k in 0..n {
                        h = ri1[k];
                        ri1[k] = s * ri[k] + c * h;
                        ri[k] = c * ri[k] - s * h;
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
