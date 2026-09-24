//! CUFSM's section template: a lipped or plain C or Z, sharp or with rounded corners, meshed into
//! strips. CUFSM `templatecalc.m` and `template_out_to_in.m`.
//!
//! The geometry is CUFSM's, including its orientation: the web lies along +z from the bottom
//! flange, flange 1 runs along +x at `z = 0` (or `z = r1` above the first corner's centre), and
//! flange 2 along `+x` for a C or `-x` for a Z at the top.

use crate::model::{Element, Material, Model, Node};
use std::f64::consts::PI;

/// C (both flanges on one side of the web) or Z (flange 2 on the other side).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    C,
    Z,
}

/// The template's dimensions and mesh, as CUFSM's `templatecalc` arguments.
///
/// With `centerline` true every length is to the mid-thickness line and the radii are centreline
/// radii, as CUFSM's `center = 1`. With it false, `h`, `b1`, `b2`, `d1`, `d2` are outside
/// dimensions and `r1`..`r4` inside radii, converted as CUFSM's `template_out_to_in.m` (after the
/// AISI Design Manual). A zero lip is a plain C or Z; zero radii are sharp corners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Template {
    pub shape: Shape,
    /// Web depth.
    pub h: f64,
    /// Flange widths: 1 at the bottom (the first node's end), 2 at the top.
    pub b1: f64,
    pub b2: f64,
    /// Lip lengths.
    pub d1: f64,
    pub d2: f64,
    /// Corner radii: 1 bottom flange-web, 2 bottom lip-flange, 3 top flange-web, 4 top lip-flange.
    pub r1: f64,
    pub r2: f64,
    pub r3: f64,
    pub r4: f64,
    /// Lip angles from the flange, degrees (90 is a lip square to the flange).
    pub q1: f64,
    pub q2: f64,
    pub t: f64,
    /// Strips in the web, each flange, each lip and each corner.
    pub nh: usize,
    pub nb1: usize,
    pub nb2: usize,
    pub nd1: usize,
    pub nd2: usize,
    pub nr1: usize,
    pub nr2: usize,
    pub nr3: usize,
    pub nr4: usize,
    pub centerline: bool,
}

impl Template {
    /// A symmetric section from outside dimensions and an inside radius, meshed as CUFSM's
    /// defaults refined: `mesh` strips across the web, half that across each flange, and two in
    /// each lip and corner. `lip = 0` makes it plain, `ri = 0` sharp.
    pub fn outside(
        shape: Shape,
        depth: f64,
        flange: f64,
        lip: f64,
        t: f64,
        ri: f64,
        mesh: usize,
    ) -> Self {
        let rounded = if ri > 0.0 { 2 } else { 0 };
        let lipped = if lip > 0.0 { 2 } else { 0 };
        Template {
            shape,
            h: depth,
            b1: flange,
            b2: flange,
            d1: lip,
            d2: lip,
            r1: ri,
            r2: ri,
            r3: ri,
            r4: ri,
            q1: 90.0,
            q2: 90.0,
            t,
            nh: mesh.max(2),
            nb1: (mesh / 2).max(2),
            nb2: (mesh / 2).max(2),
            nd1: lipped,
            nd2: lipped,
            nr1: rounded,
            nr2: rounded,
            nr3: rounded,
            nr4: rounded,
            centerline: false,
        }
    }
}

/// Centreline dimensions from outside dimensions and inside radii, CUFSM `template_out_to_in.m`.
/// Angles in radians. Returns `(h, b1, d1, b2, d2, r1, r2, r3, r4)`.
#[allow(clippy::too_many_arguments)]
pub fn template_out_to_in(
    hh: f64,
    bb1: f64,
    dd1: f64,
    q1: f64,
    bb2: f64,
    dd2: f64,
    q2: f64,
    ri1: f64,
    ri2: f64,
    ri3: f64,
    ri4: f64,
    t: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let cl = |ri: f64| if ri == 0.0 { 0.0 } else { ri + t / 2.0 };
    let (r1, r2, r3, r4) = (cl(ri1), cl(ri2), cl(ri3), cl(ri4));
    let h = hh - t / 2.0 - r1 - r3 - t / 2.0;
    let (b1, d1) = if dd1 == 0.0 {
        (bb1 - r1 - t / 2.0, 0.0)
    } else {
        (
            bb1 - r1 - t / 2.0 - (r2 + t / 2.0) * (q1 / 2.0).tan(),
            dd1 - (r2 + t / 2.0) * (q1 / 2.0).tan(),
        )
    };
    let (b2, d2) = if dd2 == 0.0 {
        (bb2 - r3 - t / 2.0, 0.0)
    } else {
        (
            bb2 - r3 - t / 2.0 - (r4 + t / 2.0) * (q2 / 2.0).tan(),
            dd2 - (r4 + t / 2.0) * (q2 / 2.0).tan(),
        )
    };
    (h, b1, d1, b2, d2, r1, r2, r3, r4)
}

/// The meshed section, CUFSM `templatecalc.m`: every node free, reference stress 1.0 (as CUFSM
/// leaves it), one element per consecutive pair of nodes, all of `material`.
pub fn templatecalc(tp: &Template, material: Material) -> Model {
    let cz = match tp.shape {
        Shape::C => 1.0,
        Shape::Z => -1.0,
    };
    let q1 = tp.q1 * PI / 180.0;
    let q2 = tp.q2 * PI / 180.0;
    let (h, b1, d1, b2, d2, r1, r2, r3, r4) = if tp.centerline {
        (tp.h, tp.b1, tp.d1, tp.b2, tp.d2, tp.r1, tp.r2, tp.r3, tp.r4)
    } else {
        template_out_to_in(
            tp.h, tp.b1, tp.d1, q1, tp.b2, tp.d2, q2, tp.r1, tp.r2, tp.r3, tp.r4, tp.t,
        )
    };
    let sharp = r1 == 0.0 && r2 == 0.0 && r3 == 0.0 && r4 == 0.0;
    let plain = d1 == 0.0 && d2 == 0.0;
    let (geom, n): (Vec<[f64; 2]>, Vec<usize>) = match (sharp, plain) {
        (true, true) => (
            vec![[b1, 0.0], [0.0, 0.0], [0.0, h], [cz * b2, h]],
            vec![tp.nb1, tp.nh, tp.nb2],
        ),
        (true, false) => (
            vec![
                [b1 + d1 * q1.cos(), d1 * q1.sin()],
                [b1, 0.0],
                [0.0, 0.0],
                [0.0, h],
                [cz * b2, h],
                [cz * (b2 + d2 * q2.cos()), h - d2 * q2.sin()],
            ],
            vec![tp.nd1, tp.nb1, tp.nh, tp.nb2, tp.nd2],
        ),
        (false, true) => (
            vec![
                [r1 + b1, 0.0],
                [r1, 0.0],
                [0.0, r1],
                [0.0, r1 + h],
                [cz * r3, r1 + h + r3],
                [cz * (r3 + b2), r1 + h + r3],
            ],
            vec![tp.nb1, tp.nr1, tp.nh, tp.nr3, tp.nb2],
        ),
        (false, false) => (
            vec![
                [
                    r1 + b1 + r2 * (PI / 2.0 - q1).cos() + d1 * q1.cos(),
                    r2 - r2 * (PI / 2.0 - q1).sin() + d1 * q1.sin(),
                ],
                [
                    r1 + b1 + r2 * (PI / 2.0 - q1).cos(),
                    r2 - r2 * (PI / 2.0 - q1).sin(),
                ],
                [r1 + b1, 0.0],
                [r1, 0.0],
                [0.0, r1],
                [0.0, r1 + h],
                [cz * r3, r1 + h + r3],
                [cz * (r3 + b2), r1 + h + r3],
                [
                    cz * (r3 + b2 + r4 * (PI / 2.0 - q2).cos()),
                    r1 + h + r3 - r4 + r4 * (PI / 2.0 - q2).sin(),
                ],
                [
                    cz * (r3 + b2 + r4 * (PI / 2.0 - q2).cos() + d2 * q2.cos()),
                    r1 + h + r3 - r4 + r4 * (PI / 2.0 - q2).sin() - d2 * q2.sin(),
                ],
            ],
            vec![
                tp.nd1, tp.nr2, tp.nb1, tp.nr1, tp.nh, tp.nr3, tp.nb2, tp.nr4, tp.nd2,
            ],
        ),
    };
    // CUFSM writes node(k, :) by index, so a later write to the same index replaces an earlier one
    // (a segment with n = 0 strips shares its start with the next segment's).
    let total: usize = n.iter().sum::<usize>() + 1;
    let mut pts: Vec<Option<[f64; 2]>> = vec![None; total];
    for i in 0..geom.len() - 1 {
        let (start, stop) = (geom[i], geom[i + 1]);
        let nstart = n[..i].iter().sum::<usize>();
        pts[nstart] = Some(start);
        let lerp = |j: usize| {
            let f = j as f64 / n[i] as f64;
            [
                start[0] + (stop[0] - start[0]) * f,
                start[1] + (stop[1] - start[1]) * f,
            ]
        };
        let corner = |j: usize| -> [f64; 2] {
            // (r, xc, zc, qstart, dq) as templatecalc.m, by segment.
            let jn = j as f64 / n[i] as f64;
            let (r, xc, zc, qs, dq) = if plain {
                match i {
                    1 => (r1, r1, r1, PI / 2.0, PI / 2.0 * jn),
                    _ => (
                        r3,
                        cz * r3,
                        r1 + h,
                        if cz == 1.0 { PI } else { 0.0 },
                        cz * PI / 2.0 * jn,
                    ),
                }
            } else {
                match i {
                    1 => (r2, r1 + b1, r2, PI / 2.0 - q1, q1 * jn),
                    3 => (r1, r1, r1, PI / 2.0, PI / 2.0 * jn),
                    5 => (
                        r3,
                        cz * r3,
                        r1 + h,
                        if cz == 1.0 { PI } else { 0.0 },
                        cz * PI / 2.0 * jn,
                    ),
                    _ => (
                        r4,
                        cz * (r3 + b2),
                        r1 + h + r3 - r4,
                        3.0 * PI / 2.0,
                        cz * q2 * jn,
                    ),
                }
            };
            [xc + r * (qs + dq).cos(), zc - r * (qs + dq).sin()]
        };
        // Straight segments are the even ones (0-based) - CUFSM's i = 1, 3, 5, ...
        let straight = sharp || i % 2 == 0;
        for j in 1..n[i] {
            pts[nstart + j] = Some(if straight { lerp(j) } else { corner(j) });
        }
    }
    pts[total - 1] = Some(geom[geom.len() - 1]);
    let nodes: Vec<Node> = pts
        .into_iter()
        .map(|p| {
            let [x, z] = p.expect("every template node is placed");
            Node::new(x, z, 1.0)
        })
        .collect();
    let elements = (0..nodes.len() - 1)
        .map(|i| Element {
            ni: i,
            nj: i + 1,
            t: tp.t,
            mat: 0,
        })
        .collect();
    Model {
        materials: vec![material],
        nodes,
        elements,
        constraints: vec![],
        springs: vec![],
    }
}
