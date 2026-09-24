//! Thin-walled section properties with torsion and warping, CUFSM `cutwp_prop2.m` (the copy in
//! `helpers/`, which CUFSM's path puts first).
//!
//! Open sections, single- or multi-branched (the latter by `helpers/`'s "brute force" walk,
//! ported as it is). A closed single cell gets its area, centroid, inertias and Bredt torsion
//! constant; its warping quantities are `NaN`, as CUFSM leaves them.

use crate::model::Model;

/// `cutwp_prop2` outputs. Coordinates are CUFSM's `(x, z)`; `theta` is the principal-axis angle
/// in radians; `wn` is the normalised unit warping at every node.
#[derive(Clone, Debug, PartialEq)]
pub struct CutwpProps {
    pub a: f64,
    pub xc: f64,
    pub zc: f64,
    pub ix: f64,
    pub iz: f64,
    pub ixz: f64,
    pub theta: f64,
    pub i1: f64,
    pub i2: f64,
    pub j: f64,
    pub xs: f64,
    pub zs: f64,
    pub cw: f64,
    pub b1: f64,
    pub b2: f64,
    pub wn: Vec<f64>,
}

pub fn cutwp_prop2(model: &Model) -> CutwpProps {
    let coord: Vec<[f64; 2]> = model.nodes.iter().map(|n| [n.x, n.z]).collect();
    let mut ends: Vec<(usize, usize, f64)> =
        model.elements.iter().map(|e| (e.ni, e.nj, e.t)).collect();
    let nele = ends.len();
    let nnode_all = coord.len();

    // Classify: count node occurrences among element ends.
    let mut count = vec![0usize; nnode_all];
    for &(s, f, _) in &ends {
        count[s] += 1;
        count[f] += 1;
    }
    let twos = count.iter().filter(|&&c| c == 2).count();
    #[derive(PartialEq)]
    enum Kind {
        Closed,
        Open,
    }
    let kind = if twos == nele {
        Kind::Closed
    } else {
        Kind::Open
    };

    if kind == Kind::Closed {
        // Reorder the elements into a loop, as cutwp_prop2 does.
        for i in 0..nele.saturating_sub(1) {
            let target = ends[i].1;
            // [m, n] = find(ends(i,2) == en(:,1:2)) with en(i,2) = 0: the first other element
            // (column-major: first column scanned first) touching ends(i,2).
            let mut found = None;
            for col in 0..2 {
                for (m, e) in ends.iter().enumerate() {
                    let v = if col == 0 { e.0 } else { e.1 };
                    let v_masked = if m == i && col == 1 { usize::MAX } else { v };
                    if v_masked == target {
                        found = Some((m, col));
                        break;
                    }
                }
                if found.is_some() {
                    break;
                }
            }
            if let Some((m, col)) = found {
                let en_m = ends[m];
                let en_next = ends[i + 1];
                if col == 0 {
                    ends[i + 1] = en_m;
                    ends[m] = en_next;
                } else {
                    ends[i + 1] = (en_m.1, en_m.0, en_m.2);
                    ends[m] = (en_next.1, en_next.0, en_next.2);
                }
            }
        }
    }

    let mut t = vec![0.0; nele];
    let mut xm = vec![0.0; nele];
    let mut ym = vec![0.0; nele];
    let mut xd = vec![0.0; nele];
    let mut yd = vec![0.0; nele];
    let mut l = vec![0.0; nele];
    for (i, &(sn, fn_, th)) in ends.iter().enumerate() {
        t[i] = th;
        xm[i] = 0.5 * (coord[sn][0] + coord[fn_][0]);
        ym[i] = 0.5 * (coord[sn][1] + coord[fn_][1]);
        xd[i] = coord[fn_][0] - coord[sn][0];
        yd[i] = coord[fn_][1] - coord[sn][1];
        l[i] = xd[i].hypot(yd[i]);
    }
    let a: f64 = (0..nele).map(|i| l[i] * t[i]).sum();
    let mut xc = (0..nele).map(|i| l[i] * t[i] * xm[i]).sum::<f64>() / a;
    let mut yc = (0..nele).map(|i| l[i] * t[i] * ym[i]).sum::<f64>() / a;
    if (xc / a.sqrt()).abs() < 1e-12 {
        xc = 0.0;
    }
    if (yc / a.sqrt()).abs() < 1e-12 {
        yc = 0.0;
    }
    let ix: f64 = (0..nele)
        .map(|i| (yd[i] * yd[i] / 12.0 + (ym[i] - yc).powi(2)) * l[i] * t[i])
        .sum();
    let iy: f64 = (0..nele)
        .map(|i| (xd[i] * xd[i] / 12.0 + (xm[i] - xc).powi(2)) * l[i] * t[i])
        .sum();
    let mut ixy: f64 = (0..nele)
        .map(|i| (xd[i] * yd[i] / 12.0 + (xm[i] - xc) * (ym[i] - yc)) * l[i] * t[i])
        .sum();
    if (ixy / (a * a)).abs() < 1e-12 {
        ixy = 0.0;
    }
    // theta = angle(Ix - Iy - 2 Ixy i) / 2. MATLAB forms the imaginary part as 0 - 2 Ixy, which
    // is +0 when Ixy is zero; -2 * 0 would be -0 and swing atan2 from +pi to -pi.
    let theta = (0.0 - 2.0 * ixy).atan2(ix - iy) / 2.0;
    let (ct, st) = (theta.cos(), theta.sin());
    let coord12: Vec<[f64; 2]> = coord
        .iter()
        .map(|p| {
            let (x, y) = (p[0] - xc, p[1] - yc);
            [ct * x + st * y, -st * x + ct * y]
        })
        .collect();
    let mut i1 = 0.0;
    let mut i2 = 0.0;
    for (i, &(sn, fn_, _)) in ends.iter().enumerate() {
        let xm12 = 0.5 * (coord12[sn][0] + coord12[fn_][0]);
        let ym12 = 0.5 * (coord12[sn][1] + coord12[fn_][1]);
        let xd12 = coord12[fn_][0] - coord12[sn][0];
        let yd12 = coord12[fn_][1] - coord12[sn][1];
        i1 += (yd12 * yd12 / 12.0 + ym12 * ym12) * l[i] * t[i];
        i2 += (xd12 * xd12 / 12.0 + xm12 * xm12) * l[i] * t[i];
    }

    if kind == Kind::Closed {
        let p: Vec<f64> = ends
            .iter()
            .enumerate()
            .map(|(i, &(sn, fn_, _))| {
                ((coord[sn][0] - xc) * (coord[fn_][1] - yc)
                    - (coord[fn_][0] - xc) * (coord[sn][1] - yc))
                    / l[i]
            })
            .collect();
        let j = 4.0 * (0..nele).map(|i| p[i] * l[i] / 2.0).sum::<f64>().powi(2)
            / (0..nele).map(|i| l[i] / t[i]).sum::<f64>();
        return CutwpProps {
            a,
            xc,
            zc: yc,
            ix,
            iz: iy,
            ixz: ixy,
            theta,
            i1,
            i2,
            j,
            xs: f64::NAN,
            zs: f64::NAN,
            cw: f64::NAN,
            b1: f64::NAN,
            b2: f64::NAN,
            wn: vec![f64::NAN; nnode_all],
        };
    }

    let j = (0..nele).map(|i| l[i] * t[i].powi(3)).sum::<f64>() / 3.0;
    // Walk the elements from the first element's start node, always taking the first element with
    // exactly one end already reached (helpers/: if none, fall back to the last element).
    let walk = |reached: &Vec<bool>| -> usize {
        let mut i = 0;
        loop {
            let (s, f, _) = ends[i];
            if reached[s] != reached[f] {
                return i;
            }
            i += 1;
            if i >= nele {
                return nele - 1;
            }
        }
    };
    let mut w = vec![0.0; nnode_all];
    let mut reached = vec![false; nnode_all];
    reached[ends[0].0] = true;
    let (mut iwx, mut iwy) = (0.0, 0.0);
    for _ in 0..nele {
        let i = walk(&reached);
        let (sn, fn_, _) = ends[i];
        let p = ((coord[sn][0] - xc) * (coord[fn_][1] - yc)
            - (coord[fn_][0] - xc) * (coord[sn][1] - yc))
            / l[i];
        if !reached[sn] {
            reached[sn] = true;
            w[sn] = w[fn_] - p * l[i];
        } else if !reached[fn_] {
            reached[fn_] = true;
            w[fn_] = w[sn] + p * l[i];
        }
        iwx += (1.0 / 3.0 * (w[sn] * (coord[sn][0] - xc) + w[fn_] * (coord[fn_][0] - xc))
            + 1.0 / 6.0 * (w[sn] * (coord[fn_][0] - xc) + w[fn_] * (coord[sn][0] - xc)))
            * t[i]
            * l[i];
        iwy += (1.0 / 3.0 * (w[sn] * (coord[sn][1] - yc) + w[fn_] * (coord[fn_][1] - yc))
            + 1.0 / 6.0 * (w[sn] * (coord[fn_][1] - yc) + w[fn_] * (coord[sn][1] - yc)))
            * t[i]
            * l[i];
    }
    let det = ix * iy - ixy * ixy;
    let (mut xs, mut ys) = if det != 0.0 {
        (
            (iy * iwy - ixy * iwx) / det + xc,
            -(ix * iwx - ixy * iwy) / det + yc,
        )
    } else {
        (xc, yc)
    };
    if (xs / a.sqrt()).abs() < 1e-12 {
        xs = 0.0;
    }
    if (ys / a.sqrt()).abs() < 1e-12 {
        ys = 0.0;
    }
    let mut wo = vec![0.0; nnode_all];
    let mut reached = vec![false; nnode_all];
    reached[ends[0].0] = true;
    let mut wno = 0.0;
    for _ in 0..nele {
        let i = walk(&reached);
        let (sn, fn_, _) = ends[i];
        let po = ((coord[sn][0] - xs) * (coord[fn_][1] - ys)
            - (coord[fn_][0] - xs) * (coord[sn][1] - ys))
            / l[i];
        if !reached[sn] {
            reached[sn] = true;
            wo[sn] = wo[fn_] - po * l[i];
        } else if !reached[fn_] {
            reached[fn_] = true;
            wo[fn_] = wo[sn] + po * l[i];
        }
        wno += 1.0 / (2.0 * a) * (wo[sn] + wo[fn_]) * t[i] * l[i];
    }
    let wn: Vec<f64> = wo.iter().map(|v| wno - v).collect();
    let cw: f64 = ends
        .iter()
        .enumerate()
        .map(|(i, &(sn, fn_, _))| {
            1.0 / 3.0 * (wn[sn] * wn[sn] + wn[sn] * wn[fn_] + wn[fn_] * wn[fn_]) * t[i] * l[i]
        })
        .sum();
    let s12 = [
        ct * (xs - xc) + st * (ys - yc),
        -st * (xs - xc) + ct * (ys - yc),
    ];
    let (mut b1, mut b2) = (0.0, 0.0);
    for (i, &(sn, fn_, _)) in ends.iter().enumerate() {
        let (x1, y1) = (coord12[sn][0], coord12[sn][1]);
        let (x2, y2) = (coord12[fn_][0], coord12[fn_][1]);
        b1 += ((y1 + y2) * (y1 * y1 + y2 * y2) / 4.0
            + (y1 * (2.0 * x1 * x1 + (x1 + x2).powi(2))
                + y2 * (2.0 * x2 * x2 + (x1 + x2).powi(2)))
                / 12.0)
            * l[i]
            * t[i];
        b2 += ((x1 + x2) * (x1 * x1 + x2 * x2) / 4.0
            + (x1 * (2.0 * y1 * y1 + (y1 + y2).powi(2))
                + x2 * (2.0 * y2 * y2 + (y1 + y2).powi(2)))
                / 12.0)
            * l[i]
            * t[i];
    }
    b1 = b1 / i1 - 2.0 * s12[1];
    b2 = b2 / i2 - 2.0 * s12[0];
    if (b1 / a.sqrt()).abs() < 1e-12 {
        b1 = 0.0;
    }
    if (b2 / a.sqrt()).abs() < 1e-12 {
        b2 = 0.0;
    }
    CutwpProps {
        a,
        xc,
        zc: yc,
        ix,
        iz: iy,
        ixz: ixy,
        theta,
        i1,
        i2,
        j,
        xs,
        zs: ys,
        cw,
        b1,
        b2,
        wn,
    }
}
