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
/// A section with no warping (a flat plate, `Cw = 0`) gives 0, as CUFSM's NaN trap does.
pub fn yield_b(fy: f64, cw: f64, wn: &[f64]) -> f64 {
    let peak = wn
        .iter()
        .map(|w| (w / cw).abs())
        .fold(f64::NEG_INFINITY, f64::max);
    let by = fy / peak;
    if by.is_finite() {
        by
    } else {
        0.0
    }
}

/// Adds the warping stress of a bimoment `b` to every node, CUFSM `warp_stress.m` with
/// `Tflag = 1`: `b * w / Cw`. If any node's increment is NaN, none is added.
pub fn add_bimoment_stress(model: &mut Model, b: f64, cw: f64, wn: &[f64]) {
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
pub fn stress_to_action(
    model: &Model,
    props: &GrossProperties,
    cw: f64,
    wn: &[f64],
) -> StressActions {
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
fn least_squares(cols: &[Vec<f64>], s: &[f64]) -> Vec<f64> {
    let k = cols.len();
    if k == 0 {
        return vec![];
    }
    let n = s.len();
    let mut a: Vec<Vec<f64>> = cols.to_vec(); // a[j][i]: column j, row i
    let mut b = s.to_vec();
    for j in 0..k {
        let norm = (j..n).map(|i| a[j][i] * a[j][i]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[j][j] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (0..n).map(|i| if i < j { 0.0 } else { a[j][i] }).collect();
        v[j] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        for col in a.iter_mut().skip(j) {
            let d: f64 = (j..n).map(|i| v[i] * col[i]).sum::<f64>() * 2.0 / vv;
            for i in j..n {
                col[i] -= d * v[i];
            }
        }
        let d: f64 = (j..n).map(|i| v[i] * b[i]).sum::<f64>() * 2.0 / vv;
        for i in j..n {
            b[i] -= d * v[i];
        }
    }
    let mut f = vec![0.0; k];
    for j in (0..k).rev() {
        let mut acc = b[j];
        for m in j + 1..k {
            acc -= a[m][j] * f[m];
        }
        f[j] = if a[j][j] != 0.0 { acc / a[j][j] } else { 0.0 };
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
}
