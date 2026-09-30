//! Gross section properties and reference stresses, CUFSM `grosprop.m` and `stresgen.m`.

use crate::model::Model;
use std::f64::consts::PI;

/// Gross properties of the strip model about its centroid, CUFSM `grosprop.m`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrossProperties {
    pub a: f64,
    pub xcg: f64,
    pub zcg: f64,
    pub ixx: f64,
    pub izz: f64,
    pub ixz: f64,
    /// Principal axis angle, degrees.
    pub thetap: f64,
    pub i11: f64,
    pub i22: f64,
}

/// CUFSM `grosprop.m`: each strip is a thin rectangle, with its own-axis inertia included.
pub fn grosprop(model: &Model) -> GrossProperties {
    let (mut a, mut ax, mut az, mut axx, mut azz, mut axz) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut ixx_o, mut izz_o, mut ixz_o) = (0.0, 0.0, 0.0);
    for e in &model.elements {
        let (ni, nj) = (&model.nodes[e.ni], &model.nodes[e.nj]);
        let t = e.t;
        let delx = nj.x - ni.x;
        let delz = nj.z - ni.z;
        let theta_xx = PI / 2.0 - delz.atan2(delx);
        let theta_zz = theta_xx + PI / 2.0;
        let xcg_e = 0.5 * (ni.x + nj.x);
        let zcg_e = 0.5 * (ni.z + nj.z);
        let l = delx.hypot(delz);
        let a_e = t * l;
        let ixx_e =
            1.0 / 12.0 * t * l * (t * t * theta_xx.sin().powi(2) + l * l * theta_xx.cos().powi(2));
        let izz_e =
            1.0 / 12.0 * t * l * (t * t * theta_zz.sin().powi(2) + l * l * theta_zz.cos().powi(2));
        let ixz_e = 1.0 / 12.0
            * t
            * l
            * (t * t * theta_xx.sin() * theta_xx.cos() + l * l * theta_xx.cos() * theta_xx.sin());
        a += a_e;
        ax += a_e * xcg_e;
        az += a_e * zcg_e;
        axx += a_e * xcg_e * xcg_e;
        azz += a_e * zcg_e * zcg_e;
        axz += a_e * xcg_e * zcg_e;
        ixx_o += ixx_e;
        izz_o += izz_e;
        ixz_o += ixz_e;
    }
    let xcg = ax / a;
    let zcg = az / a;
    let ixx = ixx_o + azz - a * zcg * zcg;
    let izz = izz_o + axx - a * xcg * xcg;
    let ixz = ixz_o + axz - a * xcg * zcg;
    let thetap = 180.0 / PI * 0.5 * (-2.0 * ixz).atan2(ixx - izz);
    let r = ((0.5 * (ixx - izz)).powi(2) + ixz * ixz).sqrt();
    GrossProperties {
        a,
        xcg,
        zcg,
        ixx,
        izz,
        ixz,
        thetap,
        i11: 0.5 * (ixx + izz) + r,
        i22: 0.5 * (ixx + izz) - r,
    }
}

/// The actions whose stresses load the section, as in CUFSM's `stresgen.m`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Actions {
    /// Axial force, compression positive (CUFSM's convention).
    pub p: f64,
    pub mxx: f64,
    pub mzz: f64,
    pub m11: f64,
    pub m22: f64,
}

/// Sets every node's reference stress from `actions`, CUFSM `stresgen.m`. With `unsymmetric`
/// false the product of inertia is ignored, as CUFSM's `unsymm = 0`.
///
/// CUFSM ships two `stresgen.m`, and they disagree on the sign of the `M11` term: the one in
/// `analysis/` subtracts it, the one in `helpers/` adds it (and zeroes any contribution that comes
/// out NaN). CUFSM puts `helpers/` on its path last, so that is the one its interface runs, and
/// the one ported here.
pub fn stresgen(model: &mut Model, actions: &Actions, props: &GrossProperties, unsymmetric: bool) {
    let GrossProperties {
        a,
        xcg,
        zcg,
        ixx,
        izz,
        thetap,
        i11,
        i22,
        ..
    } = *props;
    let ixz = if unsymmetric { props.ixz } else { 0.0 };
    let th = thetap * PI / 180.0;
    let (c, s) = (th.cos(), th.sin());
    let n = model.nodes.len();
    // Each action's increment over every node; an increment with any NaN counts as zero, as
    // CUFSM's `if max(isnan(stressinc)) stressinc = 0`.
    let guard = |inc: Vec<f64>| {
        if inc.iter().any(|v| v.is_nan()) {
            vec![0.0; n]
        } else {
            inc
        }
    };
    let (dx, dz): (Vec<f64>, Vec<f64>) = model
        .nodes
        .iter()
        .map(|nd| (nd.x - xcg, nd.z - zcg))
        .unzip();
    let from_p = guard(vec![actions.p / a; n]);
    let from_m = guard(
        (0..n)
            .map(|i| {
                -((actions.mzz * ixx + actions.mxx * ixz) * dx[i]
                    - (actions.mzz * ixz + actions.mxx * izz) * dz[i])
                    / (izz * ixx - ixz * ixz)
            })
            .collect(),
    );
    // inv([c -s; s c]) * [dx; dz] = [c s; -s c] * [dx; dz]
    let from_m11 = guard(
        (0..n)
            .map(|i| actions.m11 * (-s * dx[i] + c * dz[i]) / i11)
            .collect(),
    );
    let from_m22 = guard(
        (0..n)
            .map(|i| -actions.m22 * (c * dx[i] + s * dz[i]) / i22)
            .collect(),
    );
    for (i, nd) in model.nodes.iter_mut().enumerate() {
        nd.stress = 0.0 + from_p[i] + from_m[i] + from_m11[i] + from_m22[i];
    }
}

/// The actions that first yield the section, CUFSM `yieldMP.m` (the `helpers/` copy, which
/// zeroes a result that comes out NaN): the squash load `Py = fy A` and, for each moment, the
/// moment at which the most-stressed node reaches `fy`. The Direct Strength Method divides the
/// elastic buckling loads by these.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct YieldActions {
    pub py: f64,
    pub mxx: f64,
    pub mzz: f64,
    pub m11: f64,
    pub m22: f64,
}

/// Element corner points, CUFSM `elemcornergen.m`: each strip's two ends, offset by half its
/// thickness to either face, in the order i+, j+, j-, i-.
pub fn element_corners(model: &Model) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(4 * model.elements.len());
    for e in &model.elements {
        let (ni, nj) = (&model.nodes[e.ni], &model.nodes[e.nj]);
        let th = (nj.z - ni.z).atan2(nj.x - ni.x);
        let (dx, dz) = (-th.sin() * e.t / 2.0, th.cos() * e.t / 2.0);
        out.push((ni.x + dx, ni.z + dz));
        out.push((nj.x + dx, nj.z + dz));
        out.push((nj.x - dx, nj.z - dz));
        out.push((ni.x - dx, ni.z - dz));
    }
    out
}

fn yield_at(
    points: &[(f64, f64)],
    fy: f64,
    props: &GrossProperties,
    unsymmetric: bool,
) -> YieldActions {
    let GrossProperties {
        a,
        xcg,
        zcg,
        ixx,
        izz,
        thetap,
        i11,
        i22,
        ..
    } = *props;
    let ixz = if unsymmetric { props.ixz } else { 0.0 };
    let peak = |f: &dyn Fn(f64, f64) -> f64| {
        points
            .iter()
            .map(|&(x, z)| f(x - xcg, z - zcg).abs())
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let guard = |v: f64| if v.is_nan() { 0.0 } else { v };
    let bend = |mxx: f64, mzz: f64| {
        move |dx: f64, dz: f64| {
            ((mzz * ixx + mxx * ixz) * dx - (mzz * ixz + mxx * izz) * dz) / (izz * ixx - ixz * ixz)
        }
    };
    let th = thetap * PI / 180.0;
    let (c, s) = (th.cos(), th.sin());
    YieldActions {
        py: fy * a,
        mxx: guard(fy / peak(&bend(1.0, 0.0))),
        mzz: guard(fy / peak(&bend(0.0, 1.0))),
        m11: guard(fy / peak(&|dx, dz| (-s * dx + c * dz) / i11)),
        m22: guard(fy / peak(&|dx, dz| (c * dx + s * dz) / i22)),
    }
}

pub fn yield_mp(
    model: &Model,
    fy: f64,
    props: &GrossProperties,
    unsymmetric: bool,
) -> YieldActions {
    let nodes: Vec<(f64, f64)> = model.nodes.iter().map(|n| (n.x, n.z)).collect();
    yield_at(&nodes, fy, props, unsymmetric)
}

/// First yield at the element faces, CUFSM `yieldMP_extfiber.m` (June 2026): as [`yield_mp`], but
/// the most-stressed point is searched over [`element_corners`], so the face, not the midline,
/// reaches `fy`.
pub fn yield_mp_extfiber(
    model: &Model,
    fy: f64,
    props: &GrossProperties,
    unsymmetric: bool,
) -> YieldActions {
    yield_at(&element_corners(model), fy, props, unsymmetric)
}

/// The bimoment that first yields the section, CUFSM `yieldB.m`: `fy` over the peak `|w / Cw|`.
/// A section with no warping (a flat plate, `Cw = 0`) gives 0, as CUFSM's NaN trap does, and it
/// is a positive zero: `fy / -inf` would otherwise be `-0.0`.
pub fn yield_b(fy: f64, cw: f64, wn: &[f64]) -> f64 {
    let peak = wn
        .iter()
        .map(|w| (w / cw).abs())
        .fold(f64::NEG_INFINITY, f64::max);
    let by = fy / peak;
    if by.is_finite() && by > 0.0 {
        by
    } else {
        0.0
    }
}

/// Adds the warping stress of a bimoment `b` to every node, CUFSM `warp_stress.m` with
/// `Tflag = 1`: `b * w / Cw`. If any node's increment is NaN, none is added. `wn` holds one
/// warping value per node, as [`crate::cutwp_prop2`] returns them.
pub fn add_bimoment_stress(model: &mut Model, b: f64, cw: f64, wn: &[f64]) {
    debug_assert_eq!(wn.len(), model.nodes.len(), "one warping value per node");
    let inc: Vec<f64> = wn.iter().map(|w| b * w / cw).collect();
    if inc.iter().any(|v| v.is_nan()) {
        return;
    }
    for (n, s) in model.nodes.iter_mut().zip(inc) {
        n.stress += s;
    }
}

/// Member actions fitted to the nodal stresses, CUFSM `stress_to_action.m`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StressActions {
    /// Axial force.
    pub p: f64,
    /// Moment about the principal 11 axis.
    pub m11: f64,
    /// Moment about the principal 22 axis.
    pub m22: f64,
    /// Bimoment.
    pub b: f64,
    /// `|G f - s|`: how far the fitted actions' stresses are from the model's.
    pub err: f64,
}

/// CUFSM `stress_to_action.m`: the least-squares `P, M11, M22, B` whose stresses best match the
/// model's nodal stresses. A column the section cannot carry (non-finite, as `1/I22` for a flat
/// plate or `w/Cw` with no warping) is left out and its action reported as 0.
///
/// Divergence from CUFSM, deliberate: `stress_to_action.m` has no NaN trap, so one column it
/// cannot form makes its `f = G\s` return NaN for every action, not just that one. A flat plate
/// therefore comes back from CUFSM as NaN for P as well as for B, even though it carries P
/// perfectly well. Dropping the offending column keeps the rest of the fit honest: the carried
/// actions take their true values, only the uncarryable one reads 0, and `err` is the residual of
/// the reduced fit. `yieldB` and `warp_stress` do trap this NaN in CUFSM itself, so [`yield_b`]
/// and [`add_bimoment_stress`] agree with it exactly instead of diverging.
///
/// `wn` holds one warping value per node, as [`crate::cutwp_prop2`] returns them. A column that is
/// numerically dependent on the ones before it (to `n * eps` of its norm) is dropped in the same
/// way, with its action 0: see the least-squares note in the source.
pub fn stress_to_action(
    model: &Model,
    props: &GrossProperties,
    cw: f64,
    wn: &[f64],
) -> StressActions {
    debug_assert_eq!(wn.len(), model.nodes.len(), "one warping value per node");
    let th = props.thetap * PI / 180.0;
    let (c, s) = (th.cos(), th.sin());
    let n = model.nodes.len();
    // columns of G: 1/A, z1/I11, x1/-I22, w/Cw, with (x1, z1) the principal coordinates
    let mut cols: Vec<Vec<f64>> = (0..4).map(|_| Vec::with_capacity(n)).collect();
    for (i, nd) in model.nodes.iter().enumerate() {
        let (dx, dz) = (nd.x - props.xcg, nd.z - props.zcg);
        let (x1, z1) = (c * dx + s * dz, -s * dx + c * dz);
        cols[0].push(1.0 / props.a);
        cols[1].push(z1 / props.i11);
        cols[2].push(x1 / -props.i22);
        cols[3].push(wn.get(i).copied().unwrap_or(f64::NAN) / cw);
    }
    let keep: Vec<usize> = (0..4)
        .filter(|&k| cols[k].iter().all(|v| v.is_finite()) && cols[k].iter().any(|v| *v != 0.0))
        .collect();
    let s_vec: Vec<f64> = model.nodes.iter().map(|nd| nd.stress).collect();
    let f = least_squares(
        &keep.iter().map(|&k| cols[k].clone()).collect::<Vec<_>>(),
        &s_vec,
    );
    let mut out = [0.0; 4];
    for (j, &k) in keep.iter().enumerate() {
        out[k] = f[j];
    }
    let mut resid = s_vec.clone();
    for (j, &k) in keep.iter().enumerate() {
        for (r, g) in resid.iter_mut().zip(&cols[k]) {
            *r -= g * f[j];
        }
    }
    StressActions {
        p: out[0],
        m11: out[1],
        m22: out[2],
        b: out[3],
        err: resid.iter().map(|r| r * r).sum::<f64>().sqrt(),
    }
}

/// Least squares `min |G f - s|` for a tall, thin `G` given by columns, by Householder QR.
///
/// Rank-revealing without pivoting: a column whose part outside the span of the columns before it
/// is below `n * eps` of its own norm is numerically dependent on them. It gets no pivot row and
/// its coefficient is 0, and the remaining columns keep an exact triangular solve over their own
/// pivot rows. (MATLAB's `G\s` on a rank-deficient `G` likewise returns a basic solution with
/// zeros; it picks which column to drop by pivoting, this keeps the earlier column.)
fn least_squares(cols: &[Vec<f64>], s: &[f64]) -> Vec<f64> {
    let k = cols.len();
    let n = s.len();
    let mut a: Vec<Vec<f64>> = cols.to_vec(); // a[j][i]: column j, row i
    let mut b = s.to_vec();
    let mut pivot: Vec<Option<usize>> = vec![None; k]; // the row each kept column's R entry sits in
    let mut r = 0; // the next free pivot row
    for j in 0..k {
        if r >= n {
            break;
        }
        let own = a[j].iter().map(|x| x * x).sum::<f64>().sqrt();
        let norm = (r..n).map(|i| a[j][i] * a[j][i]).sum::<f64>().sqrt();
        if norm <= n as f64 * f64::EPSILON * own {
            continue; // zero, or dependent on the columns already kept
        }
        let alpha = if a[j][r] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (0..n).map(|i| if i < r { 0.0 } else { a[j][i] }).collect();
        v[r] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        for col in a.iter_mut().skip(j) {
            let d: f64 = (r..n).map(|i| v[i] * col[i]).sum::<f64>() * 2.0 / vv;
            for i in r..n {
                col[i] -= d * v[i];
            }
        }
        let d: f64 = (r..n).map(|i| v[i] * b[i]).sum::<f64>() * 2.0 / vv;
        for i in r..n {
            b[i] -= d * v[i];
        }
        pivot[j] = Some(r);
        r += 1;
    }
    let mut f = vec![0.0; k];
    for j in (0..k).rev() {
        let Some(row) = pivot[j] else { continue };
        let mut acc = b[row];
        for m in j + 1..k {
            acc -= a[m][row] * f[m];
        }
        f[j] = acc / a[j][row];
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Element, Material, Node};

    /// A flat plate along x: A = b t, Izz = t b³ / 12, Ixx = b t³ / 12.
    #[test]
    fn a_flat_plate() {
        let m = Model {
            materials: vec![Material::isotropic(200e3, 0.3)],
            nodes: vec![
                Node::new(0.0, 0.0, 0.0),
                Node::new(50.0, 0.0, 0.0),
                Node::new(100.0, 0.0, 0.0),
            ],
            elements: vec![
                Element {
                    ni: 0,
                    nj: 1,
                    t: 2.0,
                    mat: 0,
                },
                Element {
                    ni: 1,
                    nj: 2,
                    t: 2.0,
                    mat: 0,
                },
            ],
            constraints: vec![],
            springs: vec![],
        };
        let p = grosprop(&m);
        assert!((p.a - 200.0).abs() < 1e-12);
        assert!((p.xcg - 50.0).abs() < 1e-12);
        assert!((p.izz - 2.0 * 100f64.powi(3) / 12.0).abs() < 1e-6);
        assert!((p.ixx - 100.0 * 8.0 / 12.0).abs() < 1e-9);
    }

    /// The solver recovers an exact fit, and a column that is the first one plus rounding noise
    /// is dropped (coefficient 0). Before the rank check it was kept, and the fit split the
    /// constant term between the two near-copies arbitrarily (1.874 and 1.126 here), which in
    /// stress_to_action puts a spurious action into a column the stresses never asked for.
    #[test]
    fn least_squares_exact_and_near_dependent() {
        let x: Vec<f64> = (0..12).map(|i| i as f64 / 11.0).collect();
        let c0 = vec![1.0; 12];
        let c1: Vec<f64> = x.clone();
        let s: Vec<f64> = x.iter().map(|x| 3.0 - 2.0 * x).collect();
        let f = least_squares(&[c0.clone(), c1.clone()], &s);
        assert!(
            (f[0] - 3.0).abs() < 1e-12 && (f[1] + 2.0).abs() < 1e-12,
            "{f:?}"
        );

        let near: Vec<f64> = c0
            .iter()
            .enumerate()
            .map(|(i, v)| v + 5e-16 * i as f64)
            .collect();
        let f = least_squares(&[c0, near, c1], &s);
        assert_eq!(f[1], 0.0, "{f:?}");
        assert!(
            (f[0] - 3.0).abs() < 1e-12 && (f[2] + 2.0).abs() < 1e-12,
            "{f:?}"
        );
    }
}
