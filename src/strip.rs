//! One strip's elastic and geometric stiffness, CUFSM `klocal.m`, `kglocal.m` and `trans.m`.
//!
//! Each matrix is `totalm x totalm` blocks of 8 x 8, one block per pair of longitudinal terms
//! `(m, p)`, in the local DOF order `[u1 v1 u2 v2 w1 theta1 w2 theta2]`. The expressions are
//! ported one for one, including the two flexural entries CUFSM itself marks as not symmetric
//! within a block (they are symmetric across the `(m, p)` / `(p, m)` pair).

use crate::bc::{bc_i1_5, bc_i1_5_atpoint};
use crate::dense::Mat;
use crate::model::{BoundaryCondition, Material};
use std::f64::consts::PI;

/// The local elastic stiffness of a strip, CUFSM `klocal.m`.
pub fn klocal(mat: &Material, t: f64, a: f64, b: f64, bc: BoundaryCondition, m_a: &[f64]) -> Mat {
    let Material { ex, ey, vx, vy, g } = *mat;
    let e1 = ex / (1.0 - vx * vy);
    let e2 = ey / (1.0 - vx * vy);
    let dx = ex * t.powi(3) / (12.0 * (1.0 - vx * vy));
    let dy = ey * t.powi(3) / (12.0 * (1.0 - vx * vy));
    let d1 = vx * ey * t.powi(3) / (12.0 * (1.0 - vx * vy));
    let dxy = g * t.powi(3) / 12.0;
    let tm = m_a.len();
    let mut k = Mat::zeros(8 * tm);
    let (b2, b3, b4, b5, b6) = (b * b, b.powi(3), b.powi(4), b.powi(5), b.powi(6));
    for m in 0..tm {
        for p in 0..tm {
            let c1 = m_a[m] * PI / a;
            let c2 = m_a[p] * PI / a;
            let [i1, i2, i3, i4, i5] = bc_i1_5(bc, m_a[m], m_a[p], a);
            let mut km = [[0.0; 4]; 4];
            km[0][0] = e1 * i1 / b + g * b * i5 / 3.0;
            km[0][1] = e2 * vx * (-1.0 / 2.0 / c2) * i3 - g * i5 / 2.0 / c2;
            km[0][2] = -e1 * i1 / b + g * b * i5 / 6.0;
            km[0][3] = e2 * vx * (-1.0 / 2.0 / c2) * i3 + g * i5 / 2.0 / c2;

            km[1][0] = e2 * vx * (-1.0 / 2.0 / c1) * i2 - g * i5 / 2.0 / c1;
            km[1][1] = e2 * b * i4 / 3.0 / c1 / c2 + g * i5 / b / c1 / c2;
            km[1][2] = e2 * vx * (1.0 / 2.0 / c1) * i2 - g * i5 / 2.0 / c1;
            km[1][3] = e2 * b * i4 / 6.0 / c1 / c2 - g * i5 / b / c1 / c2;

            km[2][0] = -e1 * i1 / b + g * b * i5 / 6.0;
            km[2][1] = e2 * vx * (1.0 / 2.0 / c2) * i3 - g * i5 / 2.0 / c2;
            km[2][2] = e1 * i1 / b + g * b * i5 / 3.0;
            km[2][3] = e2 * vx * (1.0 / 2.0 / c2) * i3 + g * i5 / 2.0 / c2;

            km[3][0] = e2 * vx * (-1.0 / 2.0 / c1) * i2 + g * i5 / 2.0 / c1;
            km[3][1] = e2 * b * i4 / 6.0 / c1 / c2 - g * i5 / b / c1 / c2;
            km[3][2] = e2 * vx * (1.0 / 2.0 / c1) * i2 + g * i5 / 2.0 / c1;
            km[3][3] = e2 * b * i4 / 3.0 / c1 / c2 + g * i5 / b / c1 / c2;

            let mut kf = [[0.0; 4]; 4];
            let den = 420.0 * b3;
            kf[0][0] = (5040.0 * dx * i1 - 504.0 * b2 * d1 * i2 - 504.0 * b2 * d1 * i3
                + 156.0 * b4 * dy * i4
                + 2016.0 * b2 * dxy * i5)
                / den;
            kf[0][1] = (2520.0 * b * dx * i1 - 462.0 * b3 * d1 * i2 - 42.0 * b3 * d1 * i3
                + 22.0 * b5 * dy * i4
                + 168.0 * b3 * dxy * i5)
                / den;
            kf[0][2] = (-5040.0 * dx * i1
                + 504.0 * b2 * d1 * i2
                + 504.0 * b2 * d1 * i3
                + 54.0 * b4 * dy * i4
                - 2016.0 * b2 * dxy * i5)
                / den;
            kf[0][3] = (2520.0 * b * dx * i1
                - 42.0 * b3 * d1 * i2
                - 42.0 * b3 * d1 * i3
                - 13.0 * b5 * dy * i4
                + 168.0 * b3 * dxy * i5)
                / den;

            kf[1][0] = (2520.0 * b * dx * i1 - 462.0 * b3 * d1 * i3 - 42.0 * b3 * d1 * i2
                + 22.0 * b5 * dy * i4
                + 168.0 * b3 * dxy * i5)
                / den;
            kf[1][1] = (1680.0 * b2 * dx * i1 - 56.0 * b4 * d1 * i2 - 56.0 * b4 * d1 * i3
                + 4.0 * b6 * dy * i4
                + 224.0 * b4 * dxy * i5)
                / den;
            kf[1][2] = (-2520.0 * b * dx * i1
                + 42.0 * b3 * d1 * i2
                + 42.0 * b3 * d1 * i3
                + 13.0 * b5 * dy * i4
                - 168.0 * b3 * dxy * i5)
                / den;
            kf[1][3] = (840.0 * b2 * dx * i1 + 14.0 * b4 * d1 * i2 + 14.0 * b4 * d1 * i3
                - 3.0 * b6 * dy * i4
                - 56.0 * b4 * dxy * i5)
                / den;

            kf[2][0] = kf[0][2];
            kf[2][1] = kf[1][2];
            kf[2][2] = (5040.0 * dx * i1 - 504.0 * b2 * d1 * i2 - 504.0 * b2 * d1 * i3
                + 156.0 * b4 * dy * i4
                + 2016.0 * b2 * dxy * i5)
                / den;
            kf[2][3] = (-2520.0 * b * dx * i1 + 462.0 * b3 * d1 * i2 + 42.0 * b3 * d1 * i3
                - 22.0 * b5 * dy * i4
                - 168.0 * b3 * dxy * i5)
                / den;

            kf[3][0] = kf[0][3];
            kf[3][1] = kf[1][3];
            kf[3][2] = (-2520.0 * b * dx * i1 + 462.0 * b3 * d1 * i3 + 42.0 * b3 * d1 * i2
                - 22.0 * b5 * dy * i4
                - 168.0 * b3 * dxy * i5)
                / den;
            kf[3][3] = (1680.0 * b2 * dx * i1 - 56.0 * b4 * d1 * i2 - 56.0 * b4 * d1 * i3
                + 4.0 * b6 * dy * i4
                + 224.0 * b4 * dxy * i5)
                / den;

            for r in 0..4 {
                for c in 0..4 {
                    k.set(8 * m + r, 8 * p + c, km[r][c] * t);
                    k.set(8 * m + 4 + r, 8 * p + 4 + c, kf[r][c]);
                }
            }
        }
    }
    k
}

/// The local geometric stiffness of a strip under edge stresses `ty1`, `ty2` (stress times
/// thickness at its two nodes), CUFSM `kglocal.m`.
pub fn kglocal(a: f64, b: f64, ty1: f64, ty2: f64, bc: BoundaryCondition, m_a: &[f64]) -> Mat {
    let tm = m_a.len();
    let mut kg = Mat::zeros(8 * tm);
    for m in 0..tm {
        for p in 0..tm {
            let um = m_a[m] * PI;
            let up = m_a[p] * PI;
            let [_, _, _, i4, i5] = bc_i1_5(bc, m_a[m], m_a[p], a);
            let mut gm = [[0.0; 4]; 4];
            gm[0][0] = b * (3.0 * ty1 + ty2) * i5 / 12.0;
            gm[0][2] = b * (ty1 + ty2) * i5 / 12.0;
            gm[2][0] = gm[0][2];
            gm[1][1] = b * a * a * (3.0 * ty1 + ty2) * i4 / 12.0 / um / up;
            gm[1][3] = b * a * a * (ty1 + ty2) * i4 / 12.0 / um / up;
            gm[3][1] = gm[1][3];
            gm[2][2] = b * (ty1 + 3.0 * ty2) * i5 / 12.0;
            gm[3][3] = b * a * a * (ty1 + 3.0 * ty2) * i4 / 12.0 / um / up;

            let mut gf = [[0.0; 4]; 4];
            gf[0][0] = (10.0 * ty1 + 3.0 * ty2) * b * i5 / 35.0;
            gf[0][1] = (15.0 * ty1 + 7.0 * ty2) * b * b * i5 / 210.0 / 2.0;
            gf[1][0] = gf[0][1];
            gf[0][2] = 9.0 * (ty1 + ty2) * b * i5 / 140.0;
            gf[2][0] = gf[0][2];
            gf[0][3] = -(7.0 * ty1 + 6.0 * ty2) * b * b * i5 / 420.0;
            gf[3][0] = gf[0][3];
            gf[1][1] = (5.0 * ty1 + 3.0 * ty2) * b.powi(3) * i5 / 2.0 / 420.0;
            gf[1][2] = (6.0 * ty1 + 7.0 * ty2) * b * b * i5 / 420.0;
            gf[2][1] = gf[1][2];
            gf[1][3] = -(ty1 + ty2) * b.powi(3) * i5 / 140.0 / 2.0;
            gf[3][1] = gf[1][3];
            gf[2][2] = (3.0 * ty1 + 10.0 * ty2) * b * i5 / 35.0;
            gf[2][3] = -(7.0 * ty1 + 15.0 * ty2) * b * b * i5 / 420.0;
            gf[3][2] = gf[2][3];
            gf[3][3] = (3.0 * ty1 + 5.0 * ty2) * b.powi(3) * i5 / 420.0 / 2.0;

            for r in 0..4 {
                for c in 0..4 {
                    kg.set(8 * m + r, 8 * p + c, gm[r][c]);
                    kg.set(8 * m + 4 + r, 8 * p + 4 + c, gf[r][c]);
                }
            }
        }
    }
    kg
}

/// A spring's local stiffness, CUFSM `spring_klocal.m`, in the strip DOF order. A foundation
/// spring integrates the shape functions over the length; a discrete one takes them at `ys`.
#[allow(clippy::too_many_arguments)]
pub fn spring_klocal(
    ku: f64,
    kv: f64,
    kw: f64,
    kq: f64,
    a: f64,
    bc: BoundaryCondition,
    m_a: &[f64],
    discrete: bool,
    ys: f64,
) -> Mat {
    let tm = m_a.len();
    let mut k = Mat::zeros(8 * tm);
    for m in 0..tm {
        for p in 0..tm {
            let um = m_a[m] * PI;
            let up = m_a[p] * PI;
            let (i1, i5) = if discrete {
                let [i1, i5] = bc_i1_5_atpoint(bc, m_a[m], m_a[p], a, ys);
                (i1, i5)
            } else {
                let [i1, _, _, _, i5] = bc_i1_5(bc, m_a[m], m_a[p], a);
                (i1, i5)
            };
            let kvv = kv * i5 * a * a / (um * up);
            let km = [
                [ku * i1, 0.0, -ku * i1, 0.0],
                [0.0, kvv, 0.0, -kvv],
                [-ku * i1, 0.0, ku * i1, 0.0],
                [0.0, -kvv, 0.0, kvv],
            ];
            let kf = [
                [kw * i1, 0.0, -kw * i1, 0.0],
                [0.0, kq * i1, 0.0, -kq * i1],
                [-kw * i1, 0.0, kw * i1, 0.0],
                [0.0, -kq * i1, 0.0, kq * i1],
            ];
            for r in 0..4 {
                for c in 0..4 {
                    k.set(8 * m + r, 8 * p + c, km[r][c]);
                    k.set(8 * m + 4 + r, 8 * p + 4 + c, kf[r][c]);
                }
            }
        }
    }
    k
}

/// Rotates a local strip matrix into the section's global axes, CUFSM `trans.m`:
/// `gamma * k * gammaᵀ` with `gamma` block diagonal over the longitudinal terms.
pub fn trans(alpha: f64, k: &Mat) -> Mat {
    let (c, s) = (alpha.cos(), alpha.sin());
    let mut gam = [[0.0; 8]; 8];
    gam[0][0] = c;
    gam[0][4] = -s;
    gam[1][1] = 1.0;
    gam[2][2] = c;
    gam[2][6] = -s;
    gam[3][3] = 1.0;
    gam[4][0] = s;
    gam[4][4] = c;
    gam[5][5] = 1.0;
    gam[6][2] = s;
    gam[6][6] = c;
    gam[7][7] = 1.0;
    let tm = k.n / 8;
    let mut out = Mat::zeros(k.n);
    for bm in 0..tm {
        for bp in 0..tm {
            // out_block = gam * k_block * gamᵀ
            let mut tmp = [[0.0; 8]; 8];
            for r in 0..8 {
                for q in 0..8 {
                    let mut sum = 0.0;
                    for x in 0..8 {
                        sum += gam[r][x] * k.get(8 * bm + x, 8 * bp + q);
                    }
                    tmp[r][q] = sum;
                }
            }
            for r in 0..8 {
                for q in 0..8 {
                    let mut sum = 0.0;
                    for x in 0..8 {
                        sum += tmp[r][x] * gam[q][x];
                    }
                    out.set(8 * bm + r, 8 * bp + q, sum);
                }
            }
        }
    }
    out
}
