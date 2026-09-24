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

pub fn yield_mp(
    model: &Model,
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
        model
            .nodes
            .iter()
            .map(|n| f(n.x - xcg, n.z - zcg).abs())
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
