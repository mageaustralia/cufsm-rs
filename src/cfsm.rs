//! The constrained finite strip method (cFSM): the modal spaces - global (G), distortional (D),
//! local (L) and other (O) - that decompose a section's deformation, the analysis restricted to
//! any of them, and the classification of a mode into them. CUFSM `analysis/cFSM/`, by S. Ádány,
//! Z. Li and B. W. Schafer.
//!
//! Ported function for function: `meta_elems`, `node_class`, `mode_nr`, `DOF_ordering`,
//! `constr_xz_y`, `constr_planar_xz` (with `Kglobal_transv`), `constr_ys_ym`, `constr_yd_yg`,
//! `constr_yu_yd`, `yDOFs`, `base_vectors`, `base_column`, `base_update`, `mode_select`,
//! `mode_class` and `classify`, for open sections and the uncoupled basis (`couple = 1`, CUFSM's
//! default). Several steps take `null()` or `eig()`, whose bases are not unique, so what is
//! compared with CUFSM is what does not depend on the basis: the spaces themselves, the load
//! factors of an analysis restricted to them, and the classification percentages.

use crate::analysis::{
    assemble_strips, buckling_eigen, elemprop, msort, normalise_mode, LengthResult,
};
use crate::cutwp::cutwp_prop2;
use crate::linalg::{chol_upper, eig_sym_def, null, solve, RMat};
use crate::model::{BoundaryCondition, Model};
use crate::strip::{klocal_transv, trans};
use crate::Error;
use std::f64::consts::PI;

/// Which of CUFSM's node types a node is (`node_prop(:,4)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeType {
    /// A main node where two or more non-parallel strips meet.
    Corner,
    /// A main node at a free edge (one strip).
    Edge,
    /// A node between two collinear strips.
    Sub,
}

/// A main node of the meta-element model (`m_node`).
#[derive(Clone, Debug, PartialEq)]
pub struct MNode {
    pub x: f64,
    pub z: f64,
    /// The original node.
    pub orig: usize,
    /// Meta-elements meeting here.
    pub nel: usize,
    /// Those meta-elements, 1-based and signed as CUFSM keeps them: `+j` where this node is the
    /// meta-element's first node, `-j` where it is its second.
    pub elems: Vec<i64>,
}

/// A meta-element (`m_elem`): a run of collinear strips between two main nodes.
#[derive(Clone, Debug, PartialEq)]
pub struct MElem {
    /// Its end main nodes (0-based main-node indices).
    pub n1: usize,
    pub n2: usize,
    /// The sub-nodes along it (original node indices), in the order CUFSM merged them.
    pub subs: Vec<usize>,
}

/// The section's topology for cFSM (`base_properties.m`).
#[derive(Clone, Debug, PartialEq)]
pub struct BaseProperties {
    pub node_type: Vec<NodeType>,
    /// Each original node's index in the main-then-sub numbering (`node_prop(:,2)`).
    pub new_index: Vec<usize>,
    pub m_node: Vec<MNode>,
    pub m_elem: Vec<MElem>,
    pub nmno: usize,
    pub ncno: usize,
    pub nsno: usize,
    pub ndm: usize,
    pub nlm: usize,
    /// `DOFperm` as a permutation: old DOF `i` is new DOF `perm[i]`.
    pub perm: Vec<usize>,
}

/// CUFSM `meta_elems.m`, `node_class.m`, `mode_nr.m` and `DOF_ordering.m`.
pub fn base_properties(model: &Model) -> BaseProperties {
    let nnode = model.nodes.len();
    let nelem = model.elements.len();
    let node = |i: usize| (model.nodes[i].x, model.nodes[i].z);
    // node_prop(:,3) (element count, 0 for a sub node) and node_prop(:,4) (type).
    let mut count = vec![0usize; nnode];
    let mut ntype = vec![NodeType::Edge; nnode];
    for i in 0..nnode {
        let els: Vec<usize> = (0..nelem)
            .filter(|&j| model.elements[j].ni == i || model.elements[j].nj == i)
            .collect();
        let mel = els.len();
        count[i] = mel;
        if mel == 1 {
            ntype[i] = NodeType::Edge;
        }
        if mel >= 2 {
            ntype[i] = NodeType::Corner;
        }
        if mel == 2 {
            let other = |e: usize| {
                let el = &model.elements[e];
                if el.ni == i {
                    el.nj
                } else {
                    el.ni
                }
            };
            let (n2, n3) = (other(els[0]), other(els[1]));
            let a1 = (node(n2).1 - node(i).1).atan2(node(n2).0 - node(i).0);
            let a2 = (node(i).1 - node(n3).1).atan2(node(i).0 - node(n3).0);
            if (a1 - a2).abs() < 0.0000001 {
                count[i] = 0;
                ntype[i] = NodeType::Sub;
            }
        }
    }
    // m_el: [n1, n2, subs] per element, 1-based node numbers with 0 for a removed element.
    let mut m_el: Vec<(usize, usize, Vec<usize>)> = model
        .elements
        .iter()
        .map(|e| (e.ni + 1, e.nj + 1, Vec::new()))
        .collect();
    for i in 0..nnode {
        if count[i] != 0 {
            continue;
        }
        let me = i + 1;
        let els: Vec<usize> = (0..nelem)
            .filter(|&j| m_el[j].0 == me || m_el[j].1 == me)
            .collect();
        let other = |j: usize, m_el: &Vec<(usize, usize, Vec<usize>)>| {
            if m_el[j].0 == me {
                m_el[j].1
            } else {
                m_el[j].0
            }
        };
        let no1 = other(els[0], &m_el);
        let no2 = other(els[1], &m_el);
        m_el[els[0]].0 = no1;
        m_el[els[0]].1 = no2;
        m_el[els[1]].0 = 0;
        m_el[els[1]].1 = 0;
        m_el[els[0]].2.push(i);
    }
    let mut m_elem: Vec<MElem> = m_el
        .iter()
        .filter(|e| e.0 != 0 && e.1 != 0)
        .map(|e| MElem {
            n1: e.0 - 1,
            n2: e.1 - 1,
            subs: e.2.clone(),
        })
        .collect();
    let mut new_index = vec![0usize; nnode];
    let mut m_node: Vec<MNode> = vec![];
    for i in 0..nnode {
        if count[i] != 0 {
            new_index[i] = m_node.len();
            m_node.push(MNode {
                x: model.nodes[i].x,
                z: model.nodes[i].z,
                orig: i,
                nel: count[i],
                elems: vec![],
            });
        }
    }
    let nmno = m_node.len();
    for e in &mut m_elem {
        e.n1 = new_index[e.n1];
        e.n2 = new_index[e.n2];
    }
    for (i, mn) in m_node.iter_mut().enumerate() {
        for (j, e) in m_elem.iter().enumerate() {
            if e.n1 == i {
                mn.elems.push(j as i64 + 1);
            }
            if e.n2 == i {
                mn.elems.push(-(j as i64 + 1));
            }
        }
    }
    let mut nsno = 0;
    for i in 0..nnode {
        if count[i] == 0 {
            new_index[i] = nmno + nsno;
            nsno += 1;
        }
    }
    let ncno = ntype.iter().filter(|t| **t == NodeType::Corner).count();
    let neno = ntype.iter().filter(|t| **t == NodeType::Edge).count();
    // mode_nr.m
    let mut ndm = nmno as i64 - 4;
    for mn in &m_node {
        if mn.nel > 2 {
            ndm -= mn.nel as i64 - 2;
        }
    }
    let ndm = ndm.max(0) as usize;
    let nlm = nmno + 2 * nsno + neno;
    // DOF_ordering.m, as old DOF -> new DOF.
    let mut perm = vec![0usize; 4 * nnode];
    let (mut ic, mut ie, mut is) = (0, 0, 0);
    for i in 0..nnode {
        match ntype[i] {
            NodeType::Corner => {
                perm[2 * i] = nmno + ic;
                ic += 1;
            }
            NodeType::Edge => {
                perm[2 * i] = nmno + 2 * ncno + ie;
                ie += 1;
            }
            NodeType::Sub => {
                perm[2 * i] = 4 * nmno + is;
                is += 1;
            }
        }
    }
    let (mut ic, mut is) = (0, 0);
    for i in 0..nnode {
        if ntype[i] == NodeType::Sub {
            perm[2 * i + 1] = 4 * nmno + 3 * nsno + is;
            is += 1;
        } else {
            perm[2 * i + 1] = ic;
            ic += 1;
        }
    }
    let (mut ic, mut ie, mut is) = (0, 0, 0);
    for i in 0..nnode {
        match ntype[i] {
            NodeType::Corner => {
                perm[2 * nnode + 2 * i] = nmno + ncno + ic;
                ic += 1;
            }
            NodeType::Edge => {
                perm[2 * nnode + 2 * i] = nmno + 2 * ncno + neno + ie;
                ie += 1;
            }
            NodeType::Sub => {
                perm[2 * nnode + 2 * i] = 4 * nmno + nsno + is;
                is += 1;
            }
        }
    }
    let (mut ic, mut is) = (0, 0);
    for i in 0..nnode {
        if ntype[i] == NodeType::Sub {
            perm[2 * nnode + 2 * i + 1] = 4 * nmno + 2 * nsno + is;
            is += 1;
        } else {
            perm[2 * nnode + 2 * i + 1] = 3 * nmno + ic;
            ic += 1;
        }
    }
    BaseProperties {
        node_type: ntype,
        new_index,
        m_node,
        m_elem,
        nmno,
        ncno,
        nsno,
        ndm,
        nlm,
        perm,
    }
}

/// Per meta-element: width, 1/width, inclination, sine, cosine (`m_el_dat`).
fn m_el_dat(bp: &BaseProperties) -> Vec<[f64; 5]> {
    bp.m_elem
        .iter()
        .map(|e| {
            let (a, b) = (&bp.m_node[e.n1], &bp.m_node[e.n2]);
            let bi = (b.x - a.x).hypot(b.z - a.z);
            [
                bi,
                1.0 / bi,
                (b.z - a.z).atan2(b.x - a.x),
                (b.z - a.z) / bi,
                (b.x - a.x) / bi,
            ]
        })
        .collect()
}

/// The two non-parallel meta-elements CUFSM uses at a corner (`elem1`, `elem2`, signed 1-based).
fn corner_pair(mn: &MNode, dat: &[[f64; 5]]) -> (i64, i64) {
    let elem1 = mn.elems[0];
    let mut j = 1;
    while (dat[mn.elems[j].unsigned_abs() as usize - 1][2]
        - dat[elem1.unsigned_abs() as usize - 1][2])
        .sin()
        == 0.0
    {
        j += 1;
    }
    (elem1, mn.elems[j])
}

/// The main node at the far end of signed meta-element `e` from node `mnode`, CUFSM's
/// `if elem>0: m_elem(elem,3) else m_elem(-elem,2)`.
fn far_end(bp: &BaseProperties, e: i64) -> usize {
    if e > 0 {
        bp.m_elem[e as usize - 1].n2
    } else {
        bp.m_elem[(-e) as usize - 1].n1
    }
}

/// `(r, alfa, sin, cos)` of a signed meta-element, flipped as CUFSM flips it (`flip` true for
/// the orientation that subtracts pi).
fn oriented(dat: &[[f64; 5]], e: i64, flip: bool) -> (f64, f64, f64, f64) {
    let d = dat[e.unsigned_abs() as usize - 1];
    if flip {
        (d[1], d[2] - PI, -d[3], -d[4])
    } else {
        (d[1], d[2], d[3], d[4])
    }
}

/// CUFSM `constr_xz_y.m`: `Rx`, `Rz` (corners x main nodes).
pub fn constr_xz_y(bp: &BaseProperties) -> (RMat, RMat) {
    let dat = m_el_dat(bp);
    let corners: Vec<usize> = (0..bp.nmno).filter(|&i| bp.m_node[i].nel > 1).collect();
    let mut rx = RMat::zeros(corners.len(), bp.nmno);
    let mut rz = RMat::zeros(corners.len(), bp.nmno);
    for (k, &i) in corners.iter().enumerate() {
        let (elem1, elem2) = corner_pair(&bp.m_node[i], &dat);
        let (mnode1, mnode2, mnode3) = (far_end(bp, elem1), i, far_end(bp, elem2));
        let (r1, alfa1, sin1, cos1) = oriented(&dat, elem1, elem1 > 0);
        let (r2, alfa2, sin2, cos2) = oriented(&dat, elem2, elem2 < 0);
        let det = (alfa2 - alfa1).sin();
        rx.set(k, mnode1, sin2 * r1 / det);
        rx.set(k, mnode2, (-sin1 * r2 - sin2 * r1) / det);
        rx.set(k, mnode3, sin1 * r2 / det);
        rz.set(k, mnode1, -cos2 * r1 / det);
        rz.set(k, mnode2, (cos1 * r2 + cos2 * r1) / det);
        rz.set(k, mnode3, -cos1 * r2 / det);
    }
    (rx, rz)
}

/// CUFSM `Kglobal_transv.m`: the transverse-only stiffness at one term, in the original order.
fn kglobal_transv(model: &Model, m: f64, a: f64, bc: BoundaryCondition) -> RMat {
    let nnode = model.nodes.len();
    let mut k = RMat::zeros(4 * nnode, 4 * nnode);
    for (e, (b, alpha)) in model.elements.iter().zip(elemprop(model)) {
        let kl = trans(
            alpha,
            &klocal_transv(&model.materials[e.mat], e.t, a, b, m, bc),
        );
        let idx = |l: usize| -> usize {
            let skip = 2 * nnode;
            match l {
                0 => 2 * e.ni,
                1 => 2 * e.ni + 1,
                2 => 2 * e.nj,
                3 => 2 * e.nj + 1,
                4 => skip + 2 * e.ni,
                5 => skip + 2 * e.ni + 1,
                6 => skip + 2 * e.nj,
                _ => skip + 2 * e.nj + 1,
            }
        };
        for r in 0..8 {
            for c in 0..8 {
                k.add(idx(r), idx(c), kl.get(r, c));
            }
        }
    }
    k
}

/// `DOFperm' * K * DOFperm`: a square matrix from the original DOF order into cFSM's.
fn to_new_order(k: &RMat, perm: &[usize]) -> RMat {
    let n = k.r;
    let mut o = RMat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            o.set(perm[i], perm[j], k.get(i, j));
        }
    }
    o
}

/// CUFSM `constr_planar_xz.m`: `Rp = -Kpp \ Kpc`.
pub fn constr_planar_xz(
    model: &Model,
    bp: &BaseProperties,
    m: f64,
    a: f64,
    bc: BoundaryCondition,
) -> Result<RMat, Error> {
    let ndof = 4 * model.nodes.len();
    let k = to_new_order(&kglobal_transv(model, m, a, bc), &bp.perm);
    let (p0, p1) = (bp.nmno + 2 * bp.ncno, ndof - bp.nsno);
    let (c0, c1) = (bp.nmno, bp.nmno + 2 * bp.ncno);
    let mut kpp = RMat::zeros(p1 - p0, p1 - p0);
    let mut kpc = RMat::zeros(p1 - p0, c1 - c0);
    for i in p0..p1 {
        for j in p0..p1 {
            kpp.set(i - p0, j - p0, k.get(i, j));
        }
        for j in c0..c1 {
            kpc.set(i - p0, j - c0, k.get(i, j));
        }
    }
    let x = solve(&kpp, &kpc).ok_or_else(|| {
        Error::InvalidModel("cFSM: the transverse stiffness Kpp is singular".into())
    })?;
    Ok(x.scale(-1.0))
}

/// CUFSM `constr_ys_ym.m`: `Rys` (sub nodes x main nodes), sub-node warping by linear
/// interpolation along each meta-element.
pub fn constr_ys_ym(model: &Model, bp: &BaseProperties) -> RMat {
    let mut rys = RMat::zeros(bp.nsno, bp.nmno);
    for e in &bp.m_elem {
        if e.subs.is_empty() {
            continue;
        }
        let (nod1, nod3) = (bp.m_node[e.n1].orig, bp.m_node[e.n2].orig);
        let (x1, z1) = (model.nodes[nod1].x, model.nodes[nod1].z);
        let (x3, z3) = (model.nodes[nod3].x, model.nodes[nod3].z);
        let bm = (x3 - x1).hypot(z3 - z1);
        let (nnew1, nnew3) = (bp.new_index[nod1], bp.new_index[nod3]);
        for &nod2 in &e.subs {
            let bs = (model.nodes[nod2].x - x1).hypot(model.nodes[nod2].z - z1);
            let nnew2 = bp.new_index[nod2];
            rys.set(nnew2 - bp.nmno, nnew1, (bm - bs) / bm);
            rys.set(nnew2 - bp.nmno, nnew3, bs / bm);
        }
    }
    rys
}

/// CUFSM `constr_yd_yg.m`: `Ryd = Rysmᵀ A Rysm`, the area-weighted warping inner product.
pub fn constr_yd_yg(model: &Model, bp: &BaseProperties, rys: &RMat) -> RMat {
    let nnode = model.nodes.len();
    let mut a = RMat::zeros(nnode, nnode);
    for e in &model.elements {
        let (p, q) = (&model.nodes[e.ni], &model.nodes[e.nj]);
        let da = (q.x - p.x).hypot(q.z - p.z) * e.t;
        let (n1, n2) = (bp.new_index[e.ni], bp.new_index[e.nj]);
        a.add(n1, n1, 2.0 * da);
        a.add(n2, n2, 2.0 * da);
        a.add(n1, n2, da);
        a.add(n2, n1, da);
    }
    let mut rysm = RMat::zeros(nnode, bp.nmno);
    for i in 0..bp.nmno {
        rysm.set(i, i, 1.0);
    }
    rysm.put(bp.nmno, 0, rys);
    rysm.t().mul(&a).mul(&rysm)
}

/// CUFSM `constr_yu_yd.m`: `Rud`, the warping constraints at branched nodes.
pub fn constr_yu_yd(bp: &BaseProperties) -> RMat {
    let dat = m_el_dat(bp);
    let nmno = bp.nmno;
    let mut node_reg = vec![true; nmno];
    for i in 0..nmno {
        let mn = &bp.m_node[i];
        if mn.nel > 2 {
            let (_, elem2) = corner_pair(mn, &dat);
            for j in 1..mn.nel {
                let elem3 = mn.elems[j].unsigned_abs() as usize;
                if elem3 != elem2.unsigned_abs() as usize {
                    let e = &bp.m_elem[elem3 - 1];
                    if e.n1 != i {
                        node_reg[e.n1] = false;
                    } else {
                        node_reg[e.n2] = false;
                    }
                }
            }
        }
    }
    let mut rud = RMat::zeros(nmno, nmno);
    for i in 0..nmno {
        if node_reg[i] {
            rud.set(i, i, 1.0);
        }
    }
    for i in 0..nmno {
        let mn = &bp.m_node[i];
        if mn.nel <= 2 {
            continue;
        }
        let (elem1, elem2) = corner_pair(mn, &dat);
        let (mnode1, mnode2, mnode3) = (far_end(bp, elem1), i, far_end(bp, elem2));
        let (r1, alfa1, sin1, cos1) = oriented(&dat, elem1, elem1 > 0);
        let (r2, alfa2, sin2, cos2) = oriented(&dat, elem2, elem2 < 0);
        let det = (alfa2 - alfa1).sin();
        // csr = [sin2 -sin1; -cos2 cos1] * [r1 -r1 0; 0 r2 -r2] / det
        let cs = [[sin2, -sin1], [-cos2, cos1]];
        let r = [[r1, -r1, 0.0], [0.0, r2, -r2]];
        let mut csr = [[0.0; 3]; 2];
        for p in 0..2 {
            for q in 0..3 {
                csr[p][q] = (cs[p][0] * r[0][q] + cs[p][1] * r[1][q]) / det;
            }
        }
        for j in 1..mn.nel {
            let elem3 = mn.elems[j];
            if elem3.unsigned_abs() == elem2.unsigned_abs() {
                continue;
            }
            let e = &bp.m_elem[elem3.unsigned_abs() as usize - 1];
            let mnode4 = if e.n1 != i { e.n1 } else { e.n2 };
            let (r3, _alfa3, sin3, cos3) = oriented(&dat, elem3, elem3 < 0);
            let mut rudrow = [0.0; 3];
            for q in 0..3 {
                rudrow[q] = -1.0 / r3 * (cos3 * csr[0][q] + sin3 * csr[1][q]);
            }
            rudrow[1] += 1.0;
            rud.set(mnode4, mnode1, rudrow[0]);
            rud.set(mnode4, mnode2, rudrow[1]);
            rud.set(mnode4, mnode3, rudrow[2]);
        }
    }
    // Eliminate the dependent nodes' columns, as CUFSM's `while k == 1` loop.
    loop {
        let mut changed = false;
        for i in 0..nmno {
            if node_reg[i] {
                continue;
            }
            let ind: Vec<usize> = (0..nmno).filter(|&r| rud.get(r, i) != 0.0).collect();
            if ind.is_empty() {
                continue;
            }
            changed = true;
            for &r in &ind {
                let f = rud.get(r, i);
                for c in 0..nmno {
                    let v = rud.get(i, c) * f;
                    rud.add(r, c, v);
                }
                rud.set(r, i, 0.0);
            }
        }
        if !changed {
            break;
        }
    }
    rud
}

/// CUFSM `yDOFs.m`: the warping displacements of the main nodes for the global modes (axial,
/// two bending, torsion) and the distortional ones. Returns `(dy, ngm)`.
pub fn ydofs(
    model: &Model,
    bp: &BaseProperties,
    ryd: &RMat,
    rud: &RMat,
) -> Result<(RMat, usize), Error> {
    let cw = cutwp_prop2(model);
    let th = cw.theta;
    let (c, s) = (th.cos(), th.sin());
    // [x z] * [c -s; s c]
    let rot = |x: f64, z: f64| (x * c + z * s, -x * s + z * c);
    let cg = rot(cw.xc, cw.zc);
    let mut dy = RMat::zeros(bp.nmno, 4);
    for (i, mn) in bp.m_node.iter().enumerate() {
        let xz = rot(mn.x, mn.z);
        dy.set(i, 0, 1.0);
        dy.set(i, 1, xz.1 - cg.1);
        dy.set(i, 2, xz.0 - cg.0);
        dy.set(i, 3, cw.wn[mn.orig]);
    }
    let keep: Vec<usize> = (0..4)
        .filter(|&j| (0..bp.nmno).any(|i| dy.get(i, j) != 0.0))
        .collect();
    let ngm = keep.len();
    let mut out = RMat::zeros(bp.nmno, ngm + bp.ndm);
    for (k, &j) in keep.iter().enumerate() {
        out.set_col(k, &dy.col(j));
    }
    if bp.ndm > 0 {
        let ch = chol_upper(ryd)
            .ok_or_else(|| Error::InvalidModel("cFSM: Ryd is not positive definite".into()))?;
        let g = out.cols(0, ngm);
        let junk = null(&ch.mul(&g).t());
        let junk2 = solve(&ch, &junk)
            .ok_or_else(|| Error::InvalidModel("cFSM: chol(Ryd) is singular".into()))?;
        let jjunk1 = null(&junk2.t());
        let jjunk2 = null(&rud.t());
        let jjunk3 = jjunk1.hcat(&jjunk2);
        let jjunk4 = null(&jjunk3.t());
        let junk3 = jjunk4.t().mul(ryd).mul(&jjunk4);
        let (_, v) = crate::dense::sym_eigen(&junk3.to_square().symmetrised());
        let vd = jjunk4.mul(&RMat::from_square(&v));
        if vd.c != bp.ndm {
            return Err(Error::InvalidModel(format!(
                "cFSM: {} distortional warping vectors where {} were expected",
                vd.c, bp.ndm
            )));
        }
        out.put(0, ngm, &vd);
    }
    Ok((out, ngm))
}

/// CUFSM `base_vectors.m`: the natural base vectors at one longitudinal term, in the original
/// DOF order, columns G | D | L | O.
#[allow(clippy::too_many_arguments)]
pub fn base_vectors(
    model: &Model,
    bp: &BaseProperties,
    dy: &RMat,
    ngm: usize,
    m: f64,
    a: f64,
    rx: &RMat,
    rz: &RMat,
    rp: &RMat,
    rys: &RMat,
) -> RMat {
    let nno = model.nodes.len();
    let ndof = 4 * nno;
    let (nmno, ncno, nsno, ndm, nlm) = (bp.nmno, bp.ncno, bp.nsno, bp.ndm, bp.nlm);
    let neno = nmno - ncno;
    let ngdm = ngm + ndm;
    let km = m * PI / a;
    let elprop = elemprop(model);
    // New order first; G and D columns.
    let mut b = RMat::zeros(ndof, ndof);
    let y = dy.cols(0, ngdm);
    b.put(0, 0, &y);
    let bx = rx.mul(&y);
    let bz = rz.mul(&y);
    b.put(nmno, 0, &bx);
    b.put(nmno + ncno, 0, &bz);
    let xz = b.rows(nmno, nmno + 2 * ncno).cols(0, ngdm);
    b.put(nmno + 2 * ncno, 0, &rp.mul(&xz));
    b.put(ndof - nsno, 0, &rys.mul(&y));
    for i in nmno..(ndof - nsno) {
        for j in 0..ngdm {
            let v = b.get(i, j) / km;
            b.set(i, j, v);
        }
    }
    for j in 0..ngdm {
        let nrm = b.col(j).iter().map(|v| v * v).sum::<f64>().sqrt();
        for i in 0..ndof {
            let v = b.get(i, j) / nrm;
            b.set(i, j, v);
        }
    }
    // L columns.
    for k in 0..nmno {
        b.set(3 * nmno + k, ngdm + k, 1.0);
    }
    for k in 0..nsno {
        b.set(4 * nmno + 2 * nsno + k, ngdm + nmno + k, 1.0);
    }
    let adjacent = |i: usize| -> Vec<usize> {
        (0..model.elements.len())
            .filter(|&e| model.elements[e].ni == i || model.elements[e].nj == i)
            .collect()
    };
    let mut k = 0;
    for i in 0..nno {
        if bp.node_type[i] == NodeType::Edge {
            let alfa = elprop[adjacent(i)[0]].1;
            b.set(nmno + 2 * ncno + k, ngdm + nmno + nsno + k, -alfa.sin());
            b.set(
                nmno + 2 * ncno + neno + k,
                ngdm + nmno + nsno + k,
                alfa.cos(),
            );
            k += 1;
        }
    }
    let mut k = 0;
    for i in 0..nno {
        if bp.node_type[i] == NodeType::Sub {
            let alfa = elprop[adjacent(i)[0]].1;
            b.set(4 * nmno + k, ngdm + nmno + nsno + neno + k, -alfa.sin());
            b.set(
                4 * nmno + nsno + k,
                ngdm + nmno + nsno + neno + k,
                alfa.cos(),
            );
            k += 1;
        }
    }
    // Back to the original order: rows of the first ngdm + nlm columns.
    let mut out = RMat::zeros(ndof, ndof);
    for j in 0..(ngdm + nlm) {
        for old in 0..ndof {
            out.set(old, j, b.get(bp.perm[old], j));
        }
    }
    // O columns, set directly in the original order.
    let nel = model.elements.len();
    for (i, e) in model.elements.iter().enumerate() {
        let alfa = elprop[i].1;
        let (n1, n2) = (e.ni, e.nj);
        if ngdm + nlm + i < ndof {
            out.set(2 * n1 + 1, ngdm + nlm + i, 0.5);
            out.set(2 * n2 + 1, ngdm + nlm + i, -0.5);
        }
        if ngdm + nlm + nel + i < ndof {
            out.set(2 * n1, ngdm + nlm + nel + i, -0.5 * alfa.cos());
            out.set(2 * n2, ngdm + nlm + nel + i, 0.5 * alfa.cos());
            out.set(2 * nno + 2 * n1, ngdm + nlm + nel + i, 0.5 * alfa.sin());
            out.set(2 * nno + 2 * n2, ngdm + nlm + nel + i, -0.5 * alfa.sin());
        }
    }
    out
}

/// The natural basis at one length, block diagonal over the longitudinal terms, CUFSM
/// `base_column.m`. Returns `(b_v_l, ngm, ndm, nlm)`.
pub fn base_column(
    model: &Model,
    a: f64,
    bc: BoundaryCondition,
    m_a: &[f64],
) -> Result<(RMat, usize, usize, usize), Error> {
    let mut unit = model.clone();
    for n in &mut unit.nodes {
        n.stress = 1.0;
    }
    let bp = base_properties(&unit);
    let ndof_m = 4 * unit.nodes.len();
    let tm = m_a.len();
    let mut b = RMat::zeros(ndof_m * tm, ndof_m * tm);
    let (rx, rz) = constr_xz_y(&bp);
    let rys = constr_ys_ym(&unit, &bp);
    let ryd = constr_yd_yg(&unit, &bp, &rys);
    let rud = constr_yu_yd(&bp);
    let (dy, ngm) = ydofs(&unit, &bp, &ryd, &rud)?;
    for (ml, &m) in m_a.iter().enumerate() {
        let rp = constr_planar_xz(&unit, &bp, m, a, bc)?;
        let bm = base_vectors(&unit, &bp, &dy, ngm, m, a, &rx, &rz, &rp, &rys);
        b.put(ndof_m * ml, ndof_m * ml, &bm);
    }
    Ok((b, ngm, bp.ndm, bp.nlm))
}

/// How the base vectors are orthogonalised within each space (`GBTcon.orth`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orth {
    /// The natural basis, as built.
    Natural = 1,
    /// Modal basis, axial orthogonality (CUFSM's default).
    Axial = 2,
    /// Modal basis, load-dependent orthogonality (the section's own stresses).
    Load = 3,
}

/// How each base vector is normalised (`GBTcon.norm`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Norm {
    None = 0,
    /// Unit length (CUFSM's default for classification).
    Vector = 1,
    /// Unit strain energy, `vᵀ K v = 1`.
    StrainEnergy = 2,
    /// Unit work, `vᵀ Kg v = 1`.
    Work = 3,
}

/// CUFSM `base_update.m` for the ST basis (`ospace = 1`) and the uncoupled basis (`couple = 1`),
/// CUFSM's defaults.
#[allow(clippy::too_many_arguments)]
pub fn base_update(
    model: &Model,
    b_v_l: &RMat,
    a: f64,
    bc: BoundaryCondition,
    m_a: &[f64],
    ngm: usize,
    ndm: usize,
    nlm: usize,
    orth: Orth,
    norm: Norm,
) -> Result<RMat, Error> {
    let nnodes = model.nodes.len();
    let ndof_m = 4 * nnodes;
    let tm = m_a.len();
    let mut out = RMat::zeros(ndof_m * tm, ndof_m * tm);
    let needs_k =
        matches!(norm, Norm::StrainEnergy | Norm::Work) || matches!(orth, Orth::Axial | Orth::Load);
    let mut km_model = model.clone();
    if matches!(orth, Orth::Natural | Orth::Axial) {
        for n in &mut km_model.nodes {
            n.stress = 1.0;
        }
    }
    for (ml, &m) in m_a.iter().enumerate() {
        let mut b = RMat::zeros(ndof_m, ndof_m);
        for i in 0..ndof_m {
            for j in 0..ndof_m {
                b.set(i, j, b_v_l.get(ndof_m * ml + i, ndof_m * ml + j));
            }
        }
        let (k, kg) = if needs_k {
            let (k, kg) = assemble_strips(&km_model, a, bc, &[m]);
            (RMat::from_square(&k), RMat::from_square(&kg))
        } else {
            (RMat::zeros(0, 0), RMat::zeros(0, 0))
        };
        if matches!(orth, Orth::Axial | Orth::Load) {
            let groups = [
                (0, ngm),
                (ngm, ngm + ndm),
                (ngm + ndm, ngm + ndm + nlm),
                (ngm + ndm + nlm, ngm + ndm + nlm + nnodes - 1),
                (ngm + ndm + nlm + nnodes - 1, ndof_m),
            ];
            for &(g0, g1) in &groups {
                if g1 <= g0 {
                    continue;
                }
                let sub = b.cols(g0, g1);
                let ksub = sub.t().mul(&k).mul(&sub);
                let kgsub = sub.t().mul(&kg).mul(&sub);
                let (_, mut v) = eig_sym_def(&ksub, &kgsub).ok_or_else(|| {
                    Error::InvalidModel(format!(
                        "cFSM: Kg on base vectors {}..{} is not positive definite",
                        g0 + 1,
                        g1
                    ))
                })?;
                if matches!(norm, Norm::StrainEnergy | Norm::Work) {
                    let w = if norm == Norm::StrainEnergy {
                        &ksub
                    } else {
                        &kgsub
                    };
                    let s = v.t().mul(w).mul(&v);
                    for j in 0..v.c {
                        let f = s.get(j, j).sqrt();
                        for i in 0..v.r {
                            let x = v.get(i, j) / f;
                            v.set(i, j, x);
                        }
                    }
                }
                b.put(0, g0, &sub.mul(&v));
            }
        }
        match norm {
            Norm::StrainEnergy | Norm::Work => {
                let w = if norm == Norm::StrainEnergy { &k } else { &kg };
                for j in 0..ndof_m {
                    let col = RMat {
                        r: ndof_m,
                        c: 1,
                        data: b.col(j),
                    };
                    let f = col.t().mul(w).mul(&col).get(0, 0).sqrt();
                    let scaled: Vec<f64> = col.data.iter().map(|v| v / f).collect();
                    b.set_col(j, &scaled);
                }
            }
            Norm::Vector => {
                for j in 0..ndof_m {
                    let col = b.col(j);
                    let f = col.iter().map(|v| v * v).sum::<f64>().sqrt();
                    let scaled: Vec<f64> = col.iter().map(|v| v / f).collect();
                    b.set_col(j, &scaled);
                }
            }
            Norm::None => {}
        }
        out.put(ndof_m * ml, ndof_m * ml, &b);
    }
    Ok(out)
}

/// Which modal spaces to keep (`GBTcon.glob`, `.dist`, `.local`, `.other`, each whole).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Spaces {
    pub global: bool,
    pub distortional: bool,
    pub local: bool,
    pub other: bool,
}

/// CUFSM `mode_select.m`: the columns of the chosen spaces, for every longitudinal term.
pub fn mode_select(
    b_v: &RMat,
    ngm: usize,
    ndm: usize,
    nlm: usize,
    spaces: Spaces,
    ndof_m: usize,
    tm: usize,
) -> RMat {
    let nom = ndof_m - ngm - ndm - nlm;
    let mut idx = vec![];
    for ml in 0..tm {
        let base = ndof_m * ml;
        if spaces.global {
            idx.extend((0..ngm).map(|i| base + i));
        }
        if spaces.distortional {
            idx.extend((0..ndm).map(|i| base + ngm + i));
        }
        if spaces.local {
            idx.extend((0..nlm).map(|i| base + ngm + ndm + i));
        }
        if spaces.other {
            idx.extend((0..nom).map(|i| base + ngm + ndm + nlm + i));
        }
    }
    let mut o = RMat::zeros(b_v.r, idx.len());
    for (k, &j) in idx.iter().enumerate() {
        o.set_col(k, &b_v.col(j));
    }
    o
}

/// The finite strip analysis restricted to `spaces`, CUFSM `stripmain.m` with modal constraints
/// (the natural basis spans the same spaces as any orthogonalised one, and the restricted load
/// factors depend only on the space). Sections without fixed DOFs, constraints or springs.
pub fn stripmain_constrained(
    model: &Model,
    lengths: &[f64],
    m_all: &[Vec<f64>],
    bc: BoundaryCondition,
    neigs: usize,
    spaces: Spaces,
) -> Result<Vec<LengthResult>, Error> {
    model.validate()?;
    if model.nodes.iter().any(|n| n.free.iter().any(|f| !f))
        || !model.constraints.is_empty()
        || !model.springs.is_empty()
    {
        return Err(Error::InvalidModel(
            "cFSM here takes no fixed DOFs, constraints or springs".into(),
        ));
    }
    let mut out = vec![];
    for (&a, m_raw) in lengths.iter().zip(m_all) {
        let m_a = msort(m_raw);
        let (bvl, ngm, ndm, nlm) = base_column(model, a, bc, &m_a)?;
        let r = mode_select(
            &bvl,
            ngm,
            ndm,
            nlm,
            spaces,
            4 * model.nodes.len(),
            m_a.len(),
        );
        let (k, kg) = assemble_strips(model, a, bc, &m_a);
        let (k, kg) = (RMat::from_square(&k), RMat::from_square(&kg));
        let kff = r.t().mul(&k).mul(&r).to_square();
        let kgff = r.t().mul(&kg).mul(&r).to_square();
        let (lfs, reduced) = buckling_eigen(&kff, &kgff, neigs)?;
        let modes = reduced
            .into_iter()
            .map(|q| {
                let qm = RMat {
                    r: q.len(),
                    c: 1,
                    data: q,
                };
                let mut full = r.mul(&qm).data;
                normalise_mode(&mut full);
                full
            })
            .collect();
        out.push(LengthResult {
            length: a,
            m_terms: m_a,
            load_factors: lfs,
            modes,
        });
    }
    Ok(out)
}

/// A mode's participation in the four spaces, percent: `[G, D, L, O]`, CUFSM `mode_class.m`
/// (uncoupled basis).
pub fn mode_class(
    b_v: &RMat,
    displ: &[f64],
    ngm: usize,
    ndm: usize,
    nlm: usize,
    tm: usize,
    ndof_m: usize,
) -> Result<[f64; 4], Error> {
    let groups = [
        (0, ngm),
        (ngm, ngm + ndm),
        (ngm + ndm, ngm + ndm + nlm),
        (ngm + ndm + nlm, ndof_m),
    ];
    let mut sums = [0.0; 4];
    for ml in 0..tm {
        let mut bm = RMat::zeros(ndof_m, ndof_m);
        for i in 0..ndof_m {
            for j in 0..ndof_m {
                bm.set(i, j, b_v.get(ndof_m * ml + i, ndof_m * ml + j));
            }
        }
        let d = RMat {
            r: ndof_m,
            c: 1,
            data: displ[ndof_m * ml..ndof_m * (ml + 1)].to_vec(),
        };
        let clas = solve(&bm, &d)
            .ok_or_else(|| Error::InvalidModel("cFSM: the base vectors are singular".into()))?;
        for (g, &(g0, g1)) in groups.iter().enumerate() {
            for i in g0..g1 {
                sums[g] += clas.get(i, 0).powi(2);
            }
        }
    }
    let norms = sums.map(f64::sqrt);
    let total: f64 = norms.iter().sum();
    Ok(norms.map(|v| v / total * 100.0))
}

/// Classifies every mode of an analysis into G, D, L and O, CUFSM `classify.m`, with the basis
/// orthogonalised by `orth` and normalised by `norm` (CUFSM's defaults: `Orth::Axial`,
/// `Norm::Vector`).
pub fn classify(
    model: &Model,
    results: &[LengthResult],
    bc: BoundaryCondition,
    orth: Orth,
    norm: Norm,
) -> Result<Vec<Vec<[f64; 4]>>, Error> {
    let ndof_m = 4 * model.nodes.len();
    let mut out = vec![];
    for r in results {
        let (bvl, ngm, ndm, nlm) = base_column(model, r.length, bc, &r.m_terms)?;
        let bv = base_update(
            model, &bvl, r.length, bc, &r.m_terms, ngm, ndm, nlm, orth, norm,
        )?;
        let mut per = vec![];
        for mode in &r.modes {
            per.push(mode_class(
                &bv,
                mode,
                ngm,
                ndm,
                nlm,
                r.m_terms.len(),
                ndof_m,
            )?);
        }
        out.push(per);
    }
    Ok(out)
}
