//! The fast path for the buckling eigenproblem: a section's strips couple only neighbouring
//! nodes, so after a reverse Cuthill-McKee reordering `K` and `Kg` are narrow bands. This finds
//! the wanted load factors by Lanczos on `L⁻¹ Kg L⁻ᵀ` with `K = L Lᵀ` factored in band, and then
//! proves none was missed with a Sturm count: `K - σ Kg = L (I - σ C) Lᵀ` has as many negative
//! pivots as there are load factors in `(0, σ)` (Sylvester's law of inertia, `K` positive
//! definite). Any doubt - no narrow band, no convergence, a pivot at round-off, a count that
//! disagrees - returns `None`, and the caller solves densely.

use crate::dense::{sym_eigen, Mat};

/// A symmetric band, lower half stored: row `i` holds columns `i - b ..= i`.
struct Band {
    n: usize,
    b: usize,
    data: Vec<f64>,
}

impl Band {
    #[inline]
    fn at(&self, i: usize, j: usize) -> f64 {
        // j <= i, i - j <= b
        self.data[i * (self.b + 1) + (j + self.b - i)]
    }
    #[inline]
    fn set(&mut self, i: usize, j: usize, v: f64) {
        let b = self.b;
        self.data[i * (b + 1) + (j + b - i)] = v;
    }
    /// The permuted band of a dense symmetric matrix: new DOF `k` is old DOF `perm[k]`.
    fn from_dense(a: &Mat, perm: &[usize], b: usize) -> Band {
        let n = perm.len();
        let mut out = Band {
            n,
            b,
            data: vec![0.0; n * (b + 1)],
        };
        for i in 0..n {
            for j in i.saturating_sub(b)..=i {
                out.set(i, j, a.get(perm[i], perm[j]));
            }
        }
        out
    }
    /// Lower Cholesky factor in the same band, or `None` if a pivot is not positive.
    fn cholesky(&self) -> Option<Band> {
        let (n, b) = (self.n, self.b);
        let mut l = Band {
            n,
            b,
            data: vec![0.0; n * (b + 1)],
        };
        for j in 0..n {
            let lo = j.saturating_sub(b);
            let mut d = self.at(j, j);
            for k in lo..j {
                d -= l.at(j, k) * l.at(j, k);
            }
            if d.is_nan() || d <= 0.0 || !d.is_finite() {
                return None;
            }
            let d = d.sqrt();
            l.set(j, j, d);
            for i in (j + 1)..(j + b + 1).min(n) {
                let lo_i = i.saturating_sub(b).max(lo);
                let mut s = self.at(i, j);
                for k in lo_i..j {
                    s -= l.at(i, k) * l.at(j, k);
                }
                l.set(i, j, s / d);
            }
        }
        Some(l)
    }
    /// `L x = y` in place (this is `L`).
    fn forward(&self, x: &mut [f64]) {
        for i in 0..self.n {
            let mut s = x[i];
            for k in i.saturating_sub(self.b)..i {
                s -= self.at(i, k) * x[k];
            }
            x[i] = s / self.at(i, i);
        }
    }
    /// `Lᵀ x = y` in place.
    fn backward_t(&self, x: &mut [f64]) {
        for i in (0..self.n).rev() {
            let mut s = x[i];
            for k in (i + 1)..(i + self.b + 1).min(self.n) {
                s -= self.at(k, i) * x[k];
            }
            x[i] = s / self.at(i, i);
        }
    }
    /// `y = A x` for the symmetric band.
    fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        for i in 0..self.n {
            for j in i.saturating_sub(self.b)..=i {
                let a = self.at(i, j);
                y[i] += a * x[j];
                if i != j {
                    y[j] += a * x[i];
                }
            }
        }
        y
    }
}

/// The number of negative pivots of `K - σ Kg` by an unpivoted banded LDLᵀ, or `None` when a
/// pivot is too near zero for its sign to mean anything.
fn negative_pivots(k: &Band, kg: &Band, sigma: f64) -> Option<usize> {
    let (n, b) = (k.n, k.b);
    let mut a = Band {
        n,
        b,
        data: k
            .data
            .iter()
            .zip(&kg.data)
            .map(|(x, y)| x - sigma * y)
            .collect(),
    };
    let scale = (0..n).map(|i| a.at(i, i).abs()).fold(0.0, f64::max);
    let mut d = vec![0.0; n];
    let mut neg = 0;
    for j in 0..n {
        let lo = j.saturating_sub(b);
        // a now holds L (unit lower, below the diagonal) as it is formed.
        let mut dj = a.at(j, j);
        for k in lo..j {
            dj -= a.at(j, k) * a.at(j, k) * d[k];
        }
        if dj.abs() <= 1e-13 * scale || !dj.is_finite() {
            return None;
        }
        d[j] = dj;
        if dj < 0.0 {
            neg += 1;
        }
        for i in (j + 1)..(j + b + 1).min(n) {
            let lo_i = i.saturating_sub(b).max(lo);
            let mut s = a.at(i, j);
            for k in lo_i..j {
                s -= a.at(i, k) * a.at(j, k) * d[k];
            }
            a.set(i, j, s / dj);
        }
    }
    Some(neg)
}

/// Reverse Cuthill-McKee on the nonzero pattern of `a + b` (dense). Returns the ordering (new ->
/// old) and its half-bandwidth.
fn rcm(a: &Mat, b: &Mat) -> (Vec<usize>, usize) {
    let n = a.n;
    let adj: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && (a.get(i, j) != 0.0 || b.get(i, j) != 0.0))
                .collect()
        })
        .collect();
    let mut order = Vec::with_capacity(n);
    let mut seen = vec![false; n];
    while order.len() < n {
        // Each component from its lowest-degree unvisited node.
        let start = (0..n)
            .filter(|&i| !seen[i])
            .min_by_key(|&i| adj[i].len())
            .unwrap();
        seen[start] = true;
        let mut head = order.len();
        order.push(start);
        while head < order.len() {
            let v = order[head];
            head += 1;
            let mut next: Vec<usize> = adj[v].iter().copied().filter(|&u| !seen[u]).collect();
            next.sort_by_key(|&u| adj[u].len());
            for u in next {
                seen[u] = true;
                order.push(u);
            }
        }
    }
    order.reverse();
    let mut pos = vec![0; n];
    for (k, &o) in order.iter().enumerate() {
        pos[o] = k;
    }
    let bw = (0..n)
        .flat_map(|i| adj[i].iter().map(move |&j| (i, j)))
        .map(|(i, j)| pos[i].abs_diff(pos[j]))
        .max()
        .unwrap_or(0);
    (order, bw)
}

/// The `neigs` smallest positive load factors of `K φ = λ Kg φ` and their modes, by the band
/// path; `None` to solve densely instead.
pub fn buckling_eigen_banded(k: &Mat, kg: &Mat, neigs: usize) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = k.n;
    if n < 60 || neigs == 0 {
        return None;
    }
    let (perm, b) = rcm(k, kg);
    if 4 * b > n {
        return None;
    }
    let kb = Band::from_dense(&k.symmetrised(), &perm, b);
    let kgb = Band::from_dense(&kg.symmetrised(), &perm, b);
    let l = kb.cholesky()?;
    let op = |x: &[f64]| -> Vec<f64> {
        let mut y = x.to_vec();
        l.backward_t(&mut y);
        let mut z = kgb.matvec(&y);
        l.forward(&mut z);
        z
    };
    // Lanczos with full re-orthogonalisation, lengthened until the wanted Ritz pairs converge.
    let mut m = (2 * neigs + 30).min(n);
    loop {
        if m >= n {
            return None;
        }
        let mut q: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        let mut v: Vec<f64> = (0..n)
            .map(|i| 1.0 + ((i * 7919) % 997) as f64 * 1e-3)
            .collect();
        let nv = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.iter_mut().for_each(|x| *x /= nv);
        q.push(v);
        let (mut alpha, mut beta) = (Vec::with_capacity(m), Vec::with_capacity(m));
        let mut last_beta = 0.0;
        for j in 0..m {
            let mut w = op(&q[j]);
            let a_j: f64 = q[j].iter().zip(&w).map(|(x, y)| x * y).sum();
            for _ in 0..2 {
                for qi in &q {
                    let dot: f64 = qi.iter().zip(&w).map(|(a, c)| a * c).sum();
                    w.iter_mut().zip(qi).for_each(|(x, y)| *x -= dot * y);
                }
            }
            alpha.push(a_j);
            let b_j = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            last_beta = b_j;
            if j + 1 == m || b_j <= f64::EPSILON * a_j.abs().max(1e-300) {
                break;
            }
            beta.push(b_j);
            w.iter_mut().for_each(|x| *x /= b_j);
            q.push(w);
        }
        let mm = alpha.len();
        let mut t = Mat::zeros(mm);
        for i in 0..mm {
            t.set(i, i, alpha[i]);
            if i + 1 < mm {
                t.set(i, i + 1, beta[i]);
                t.set(i + 1, i, beta[i]);
            }
        }
        let (theta, s) = sym_eigen(&t);
        let tmax = theta.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        // Largest positive Ritz values first, above the null of Kg.
        let mut idx: Vec<usize> = (0..mm).filter(|&i| theta[i] > tmax * 1e-14).collect();
        idx.sort_by(|&a, &c| {
            theta[c]
                .partial_cmp(&theta[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        idx.truncate(neigs);
        let converged = idx.len() == neigs.min(mm)
            && idx
                .iter()
                .all(|&i| (last_beta * s.get(mm - 1, i)).abs() <= 1e-12 * tmax);
        if !converged || idx.len() < neigs {
            m *= 2;
            continue;
        }
        let lfs: Vec<f64> = idx.iter().map(|&i| 1.0 / theta[i]).collect();
        // Sturm check just above the last load factor found: exactly that many lie below.
        let sigma = lfs[lfs.len() - 1] * (1.0 + 1e-7);
        if negative_pivots(&kb, &kgb, sigma)? != lfs.len() {
            return None;
        }
        let mut modes = Vec::with_capacity(idx.len());
        for &i in &idx {
            let mut y = vec![0.0; n];
            for (r, qr) in q.iter().enumerate().take(mm) {
                let c = s.get(r, i);
                y.iter_mut().zip(qr).for_each(|(a, x)| *a += c * x);
            }
            l.backward_t(&mut y);
            let mut phi = vec![0.0; n];
            for (kk, &o) in perm.iter().enumerate() {
                phi[o] = y[kk];
            }
            modes.push(phi);
        }
        return Some((lfs, modes));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::{assemble, buckling_eigen_dense};
    use crate::template::{templatecalc, Shape, Template};
    use crate::{grosprop, stresgen, Actions, BoundaryCondition, Material};

    /// On a meshed lipped channel the band path answers (it does not decline), and its load
    /// factors and modes are the dense solve's.
    #[test]
    fn band_path_answers_and_agrees_with_the_dense_solve() {
        for actions in [
            Actions {
                p: 1.0,
                ..Default::default()
            },
            Actions {
                mxx: 1e6,
                ..Default::default()
            },
        ] {
            let mut m = templatecalc(
                &Template::outside(Shape::C, 200.0, 76.0, 15.0, 1.9, 3.0, 12),
                Material::isotropic(203000.0, 0.3),
            );
            let p = grosprop(&m);
            stresgen(&mut m, &actions, &p, false);
            for a in [40.0, 150.0, 600.0, 3000.0] {
                let (k, kg) = assemble(&m, a, BoundaryCondition::SS, &[1.0]);
                let (bl, bm) =
                    buckling_eigen_banded(&k, &kg, 10).expect("the band path should answer here");
                let (dl, dm) = buckling_eigen_dense(&k, &kg, 10).unwrap();
                assert_eq!(bl.len(), dl.len());
                // Both solves are rounding-limited on a long global mode: allow what the parity
                // tests allow, 1e-10 + 1e-12 x cond(K) (largest over smallest Cholesky pivot, squared).
                let lf = crate::dense::cholesky(&k.symmetrised()).unwrap();
                let piv: Vec<f64> = (0..lf.n).map(|j| lf.get(j, j)).collect();
                let cond = (piv.iter().cloned().fold(0.0, f64::max)
                    / piv.iter().cloned().fold(f64::INFINITY, f64::min))
                .powi(2);
                let tol0 = 1e-10 + 1e-12 * cond;
                for (i, (x, y)) in bl.iter().zip(&dl).enumerate() {
                    assert!(
                        (x / y - 1.0).abs() < tol0 * (y / dl[0]).max(1.0),
                        "a = {a}, mode {}: band {x} dense {y}",
                        i + 1
                    );
                    // Same mode up to sign and scale, where the load factor is distinct.
                    let distinct = (i == 0 || (dl[i] / dl[i - 1] - 1.0).abs() > 1e-6)
                        && (i + 1 >= dl.len() || (dl[i + 1] / dl[i] - 1.0).abs() > 1e-6);
                    if distinct {
                        let (u, v) = (&bm[i], &dm[i]);
                        let uv: f64 = u.iter().zip(v).map(|(p, q)| p * q).sum();
                        let uu: f64 = u.iter().map(|p| p * p).sum();
                        let vv: f64 = v.iter().map(|p| p * p).sum();
                        assert!(
                            1.0 - uv * uv / (uu * vv) < 1e-8,
                            "a = {a}, mode {}: MAC",
                            i + 1
                        );
                    }
                }
            }
        }
    }

    /// It declines where it should: a small model, and a count that cannot be met.
    #[test]
    fn band_path_declines_small_models() {
        let mut m = templatecalc(
            &Template::outside(Shape::C, 150.0, 50.0, 0.0, 2.0, 0.0, 4),
            Material::isotropic(203000.0, 0.3),
        );
        let p = grosprop(&m);
        stresgen(
            &mut m,
            &Actions {
                p: 1.0,
                ..Default::default()
            },
            &p,
            false,
        );
        let (k, kg) = assemble(&m, 500.0, BoundaryCondition::SS, &[1.0]);
        assert!(k.n < 60 && buckling_eigen_banded(&k, &kg, 5).is_none());
    }
}
