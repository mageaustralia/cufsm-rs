//! The finite tube method (FTM): linear buckling of a circular tube, after Ádány and Schafer,
//! "Computationally efficient buckling analysis of wind turbine towers" (SSRC Annual Stability
//! Conference, 2023) and "Finite tube method for buckling analysis of tubular members using
//! Fourier-approximation for the displacements" (Thin-Walled Structures 206, 2025).
//!
//! One prismatic segment of radius `R` (to the mid-surface), thickness `t` and length `L`. Each of
//! the circumferential `u`, longitudinal `v` and radial `w` displacements is a product of two full
//! Fourier series: `1, sin iθ, cos iθ` (i = 1..p) around the tube and `1, sin jπy/L, cos jπy/L`
//! (j = 1..q) along it. The strains are the paper's (Sanders-type, after Silvestre 2007): the
//! surface slopes `βx = −w,y`, `βy = u/R − w,x`, `βz = −u,y`, curvatures `κxx = βy,x`,
//! `κyy = βx,y`, `κxy = βx,x + βy,y + βz/R`, membrane strains `εx = u,x + w/R`, `εy = v,y`,
//! `γ = v,x + u,y`, and second-order strains `εy = (βx² + βz²)/2`, `γ = βx βy` (hoop terms left
//! out, as in the paper). The material matrix is the plane-stress one.
//!
//! The stresses are not found by a static analysis but set from the section actions, as in the
//! finite strip method, uniform along the tube: compression `σ = N/A` and torsion `τ = T/(2πR²t)`
//! uniform around it, bending `σ = M/(πR²t) cos θ` (compression at θ = 0) and shear
//! `τ = V/(πRt) sin θ`. Hoop stress is not carried.
//!
//! Supports are constraint equations at the ends, imposed exactly: a displacement is zero all
//! round an end when, for every circumferential term, its longitudinal series sums to zero there.
//! The load factor `λ` multiplies all four stresses together.
//!
//! What this does not carry yet: several segments (a stepped or tapered tower), stresses that vary
//! along the tube, springs, pressure, and the paper's spectral and buckling-length tools.

use crate::dense::{sym_eigen, sym_eigen_top, Mat};
use crate::linalg::{null, RMat};
use crate::Error;
use std::collections::HashMap;
use std::f64::consts::PI;

/// How an end of the tube is held.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// Nothing held.
    Free,
    /// `u = w = 0` all round (hinged). The first pinned end also holds `v`, so the tube cannot
    /// slide along its axis.
    Pinned,
    /// `u = v = w = 0` and `w,y = 0` all round.
    Clamped,
}

/// The tube: mid-surface radius, thickness and length (mm), Young's modulus (MPa), Poisson's ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tube {
    pub r: f64,
    pub t: f64,
    pub l: f64,
    pub e: f64,
    pub nu: f64,
}

/// The reference stresses (MPa) the load factor multiplies: `n` uniform compression, `m` the peak
/// bending compression (at θ = 0), `t` uniform torsional shear, `v` the peak shear from a shear
/// force (at θ = π/2).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Stresses {
    pub n: f64,
    pub m: f64,
    pub t: f64,
    pub v: f64,
}

impl Stresses {
    /// From section actions (N, N·mm): compression `N`, moment `M`, torque `T`, shear `V`, by the
    /// thin-tube formulas.
    pub fn from_actions(tube: &Tube, n: f64, m: f64, t: f64, v: f64) -> Self {
        let (r, th) = (tube.r, tube.t);
        Stresses {
            n: n / (2.0 * PI * r * th),
            m: m / (PI * r * r * th),
            t: t / (2.0 * PI * r * r * th),
            v: v / (PI * r * th),
        }
    }
}

/// One buckling mode: its load factor and the Fourier coefficients, laid out `[d][c][l]` with
/// `d` = u, v, w, `c` the circumferential term (`0`, then `sin iθ`, `cos iθ` pairs) and `l` the
/// longitudinal one, scaled so the largest is 1.
#[derive(Clone, Debug, PartialEq)]
pub struct FtmMode {
    pub load_factor: f64,
    pub coef: Vec<f64>,
    /// The circumferential wave number `i` carrying most of the radial displacement.
    pub circ_waves: usize,
}

/// The analysis: the lowest modes, smallest load factor first.
#[derive(Clone, Debug, PartialEq)]
pub struct FtmResult {
    pub tube: Tube,
    /// Circumferential harmonics `p`, and the longitudinal wave numbers used (each with its sine
    /// and cosine, plus the constant).
    pub circ: usize,
    pub long: Vec<usize>,
    pub modes: Vec<FtmMode>,
    /// The number of unknowns after the supports, over every uncoupled block.
    pub dofs: usize,
}

/// Value of the `order`-th derivative of Fourier function `idx` (0: 1; 2j−1: sin jωs; 2j: cos jωs).
fn basis(idx: usize, w: f64, s: f64, order: u8) -> f64 {
    if idx == 0 {
        return if order == 0 { 1.0 } else { 0.0 };
    }
    let om = idx.div_ceil(2) as f64 * w;
    let ph = om * s + order as f64 * PI / 2.0;
    om.powi(order as i32) * if idx % 2 == 1 { ph.sin() } else { ph.cos() }
}

/// The `order`-th derivative of longitudinal function `l` (0: 1; 2k−1: sin jₖπy/L; 2k: cos jₖπy/L)
/// for wave numbers `js`.
fn lbasis(l: usize, js: &[usize], len: f64, y: f64, order: u8) -> f64 {
    if l == 0 {
        return if order == 0 { 1.0 } else { 0.0 };
    }
    let om = js[(l - 1) / 2] as f64 * PI / len;
    let ph = om * y + order as f64 * PI / 2.0;
    om.powi(order as i32) * if l % 2 == 1 { ph.sin() } else { ph.cos() }
}

/// A term of a strain quantity: coefficient × ∂ᵃ/∂xᵃ ∂ᵇ/∂yᵇ of displacement `d`.
#[derive(Clone, Copy)]
struct Term {
    d: usize,
    a: u8,
    b: u8,
    c: f64,
}

const fn tm(d: usize, a: u8, b: u8, c: f64) -> Term {
    Term { d, a, b, c }
}

/// A quadratic form ∫∫ weight(θ) · scale · A · B dA.
struct Form {
    a: Vec<Term>,
    b: Vec<Term>,
    /// 0: 1, 1: cos θ, 2: sin θ.
    weight: usize,
    scale: f64,
}

/// 8-point Gauss–Legendre on [−1, 1].
const GAUSS: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

/// The one-dimensional integrals every stiffness entry is a sum of products of.
struct Tables {
    nc: usize,
    nl: usize,
    /// ∫₀^{2πR} weight · X_c1⁽ᵃ⁾ X_c2⁽ᵇ⁾ dx, indexed [weight][a][b][c1][c2].
    ix: Vec<f64>,
    /// ∫₀ᴸ Y_l1⁽ᵃ⁾ Y_l2⁽ᵇ⁾ dy, indexed [a][b][l1][l2].
    iy: Vec<f64>,
    /// Below this an entry of `ix` is zero, per [weight][a][b] sub-table: the trapezoid rule leaves
    /// rounding where the exact integral vanishes, and derivative orders differ by `(i/R)⁴`.
    ix_zero: Vec<f64>,
}

impl Tables {
    fn new(tube: &Tube, p: usize, js: &[usize]) -> Self {
        let (nc, nl) = (2 * p + 1, 2 * js.len() + 1);
        // Around the tube: the trapezoid rule is exact for these trigonometric products.
        let nth = 4 * p + 8;
        let dx = 2.0 * PI * tube.r / nth as f64;
        let mut ix = vec![0.0; 3 * 9 * nc * nc];
        let mut vals = vec![0.0; 3 * nc];
        for k in 0..nth {
            let th = 2.0 * PI * k as f64 / nth as f64;
            let x = tube.r * th;
            for a in 0..3u8 {
                for c in 0..nc {
                    vals[a as usize * nc + c] = basis(c, 1.0 / tube.r, x, a);
                }
            }
            for (wk, wt) in [1.0, th.cos(), th.sin()].into_iter().enumerate() {
                for a in 0..3 {
                    for b in 0..3 {
                        let base = ((wk * 3 + a) * 3 + b) * nc * nc;
                        for c1 in 0..nc {
                            let v1 = wt * dx * vals[a * nc + c1];
                            if v1 == 0.0 {
                                continue;
                            }
                            for c2 in 0..nc {
                                ix[base + c1 * nc + c2] += v1 * vals[b * nc + c2];
                            }
                        }
                    }
                }
            }
        }
        // Along it: composite Gauss, one 8-point panel per longitudinal half-wave and then some.
        let panels = js.iter().max().copied().unwrap_or(1) + 2;
        let h = tube.l / panels as f64;
        let mut iy = vec![0.0; 9 * nl * nl];
        let mut yv = vec![0.0; 3 * nl];
        for pnl in 0..panels {
            for &(g, gw) in &GAUSS {
                let y = h * (pnl as f64 + 0.5 * (g + 1.0));
                let w = 0.5 * h * gw;
                for a in 0..3u8 {
                    for l in 0..nl {
                        yv[a as usize * nl + l] = lbasis(l, js, tube.l, y, a);
                    }
                }
                for a in 0..3 {
                    for b in 0..3 {
                        let base = (a * 3 + b) * nl * nl;
                        for l1 in 0..nl {
                            let v1 = w * yv[a * nl + l1];
                            for l2 in 0..nl {
                                iy[base + l1 * nl + l2] += v1 * yv[b * nl + l2];
                            }
                        }
                    }
                }
            }
        }
        let ix_zero = ix
            .chunks(nc * nc)
            .map(|sub| 1e-10 * sub.iter().fold(0.0_f64, |m, v| m.max(v.abs())))
            .collect();
        Tables {
            nc,
            nl,
            ix,
            iy,
            ix_zero,
        }
    }

    fn x(&self, w: usize, a: u8, b: u8, c1: usize, c2: usize) -> f64 {
        let sub = (w * 3 + a as usize) * 3 + b as usize;
        let v = self.ix[sub * self.nc * self.nc + c1 * self.nc + c2];
        if v.abs() > self.ix_zero[sub] {
            v
        } else {
            0.0
        }
    }

    fn y(&self, a: u8, b: u8, l1: usize, l2: usize) -> f64 {
        self.iy[(a as usize * 3 + b as usize) * self.nl * self.nl + l1 * self.nl + l2]
    }
}

/// The elastic and geometric quadratic forms, with the stresses in the geometric ones.
fn forms(tube: &Tube, s: &Stresses) -> (Vec<Form>, Vec<Form>) {
    let r = tube.r;
    let (u, v, w) = (0, 1, 2);
    let mem = [
        vec![tm(u, 1, 0, 1.0), tm(w, 0, 0, 1.0 / r)],
        vec![tm(v, 0, 1, 1.0)],
        vec![tm(v, 1, 0, 1.0), tm(u, 0, 1, 1.0)],
    ];
    let bend = [
        vec![tm(u, 1, 0, 1.0 / r), tm(w, 2, 0, -1.0)],
        vec![tm(w, 0, 2, -1.0)],
        // Sanders: κxy = βx,x + βy,y + βz/R with βz = (u,y − v,x)/2, so rigid rotations are
        // strain-free (the Donnell −2w,xy alone makes a column bending as a beam twist).
        vec![
            tm(w, 1, 1, -2.0),
            tm(u, 0, 1, 1.5 / r),
            tm(v, 1, 0, -0.5 / r),
        ],
    ];
    let c = tube.e / (1.0 - tube.nu * tube.nu);
    let d = [
        [c, c * tube.nu, 0.0],
        [c * tube.nu, c, 0.0],
        [0.0, 0.0, c * (1.0 - tube.nu) / 2.0],
    ];
    let mut ke = vec![];
    for i in 0..3 {
        for j in 0..3 {
            if d[i][j] == 0.0 {
                continue;
            }
            ke.push(Form {
                a: mem[i].clone(),
                b: mem[j].clone(),
                weight: 0,
                scale: tube.t * d[i][j],
            });
            ke.push(Form {
                a: bend[i].clone(),
                b: bend[j].clone(),
                weight: 0,
                scale: tube.t.powi(3) / 12.0 * d[i][j],
            });
        }
    }
    let bx = vec![tm(w, 0, 1, -1.0)];
    let by = vec![tm(u, 0, 0, 1.0 / r), tm(w, 1, 0, -1.0)];
    let bz = vec![tm(u, 0, 1, 0.5), tm(v, 1, 0, -0.5)];
    // Second-order energy t(σy (βx² + βz²)/2 + τ βx βy), σy tension-positive = −(n + m cos θ);
    // Kg is its negative, so Ke φ = λ Kg φ.
    let mut kg = vec![];
    for (weight, sig) in [(0, s.n), (1, s.m)] {
        if sig != 0.0 {
            kg.push(Form {
                a: bx.clone(),
                b: bx.clone(),
                weight,
                scale: tube.t * sig,
            });
            kg.push(Form {
                a: bz.clone(),
                b: bz.clone(),
                weight,
                scale: tube.t * sig,
            });
        }
    }
    for (weight, tau) in [(0, s.t), (2, s.v)] {
        if tau != 0.0 {
            kg.push(Form {
                a: bx.clone(),
                b: by.clone(),
                weight,
                scale: -tube.t * tau,
            });
            kg.push(Form {
                a: by.clone(),
                b: bx.clone(),
                weight,
                scale: -tube.t * tau,
            });
        }
    }
    (ke, kg)
}

/// The supported space of each displacement's longitudinal series, for the constant
/// circumferential term (`[d][0]`), the first harmonic (`[d][1]`) and the higher ones (`[d][2]`):
/// the null space of its end constraints. A pinned end is a hinge with a rigid end plate: `u` and
/// `w` held, the end section free only to turn as a plane.
fn constraints(tube: &Tube, js: &[usize], base: End, top: End, tabs: &Tables) -> [[RMat; 3]; 3] {
    let nl = 2 * js.len() + 1;
    let axial_at_base = base != End::Free;
    let mut rows: [[Vec<Vec<f64>>; 3]; 3] = Default::default();
    for (end, y) in [(base, 0.0), (top, tube.l)] {
        let at = |order: u8| {
            (0..nl)
                .map(|l| lbasis(l, js, tube.l, y, order))
                .collect::<Vec<_>>()
        };
        for k in 0..3 {
            match end {
                End::Free => {}
                End::Pinned => {
                    rows[0][k].push(at(0));
                    rows[2][k].push(at(0));
                    let first_held = if y == 0.0 { true } else { !axial_at_base };
                    // The end plate turns as a plane: v is free in the first harmonic (the
                    // rotation), held in the higher ones (no warping), and the constant term is
                    // held at the first pinned end so the tube cannot slide.
                    if (k == 0 && first_held) || k == 2 {
                        rows[1][k].push(at(0));
                    }
                }
                End::Clamped => {
                    rows[0][k].push(at(0));
                    rows[1][k].push(at(0));
                    rows[2][k].push(at(0));
                    rows[2][k].push(at(1));
                }
            }
        }
    }
    rows.map(|rk| {
        rk.map(|rs| {
            let mut m = RMat::zeros(rs.len(), nl);
            for (i, r) in rs.iter().enumerate() {
                for (j, v) in r.iter().enumerate() {
                    m.set(i, j, *v);
                }
            }
            orthonormal(&null(&m), tabs)
        })
    })
}

/// `Z` made L²-orthonormal along the tube, `Zᵀ G Z = I` with `G = ∫ Yᵢ Yⱼ dy`, dropping the
/// directions the sine and cosine series share: on half a period they are nearly dependent, and
/// with ten or more wave numbers the plain basis defeats a Cholesky factorisation.
fn orthonormal(z: &RMat, tabs: &Tables) -> RMat {
    let nl = tabs.nl;
    let mut g = RMat::zeros(nl, nl);
    for a in 0..nl {
        for b in 0..nl {
            g.set(a, b, tabs.y(0, 0, a, b));
        }
    }
    let m = z.t().mul(&g).mul(z);
    let (vals, vecs) = sym_eigen(&m.to_square().symmetrised());
    let top = vals.iter().fold(0.0_f64, |a, v| a.max(*v));
    let keep: Vec<usize> = (0..vals.len()).filter(|&j| vals[j] > 1e-11 * top).collect();
    let mut out = RMat::zeros(nl, keep.len());
    for (k, &j) in keep.iter().enumerate() {
        let s = 1.0 / vals[j].sqrt();
        for i in 0..nl {
            let v: f64 = (0..z.c).map(|c| z.get(i, c) * vecs.get(c, j)).sum();
            out.set(i, k, v * s);
        }
    }
    out
}

struct Dsu(Vec<usize>);
impl Dsu {
    fn find(&mut self, i: usize) -> usize {
        let mut i = i;
        while self.0[i] != i {
            self.0[i] = self.0[self.0[i]];
            i = self.0[i];
        }
        i
    }
    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.0[a] = b;
        }
    }
}

/// Longitudinal wave numbers that suit most analyses: `1..=q`, plus a band of five around the
/// classical half-wavelength `1.728 √(Rt)` (local buckling under compression and bending).
pub fn default_long_terms(tube: &Tube, q: usize) -> Vec<usize> {
    let mut js: Vec<usize> = (1..=q).collect();
    let jstar = (tube.l / (1.728 * (tube.r * tube.t).sqrt())).round() as usize;
    if jstar > q {
        js.extend(jstar.saturating_sub(2).max(q + 1)..=jstar + 2);
    }
    js.truncate(60);
    js
}

/// The finite tube analysis: the `nmodes` lowest positive load factors of the stresses `s` on
/// `tube`, with `p` circumferential harmonics and `q` longitudinal terms.
pub fn ftm_buckle(
    tube: &Tube,
    s: &Stresses,
    base: End,
    top: End,
    p: usize,
    long_terms: &[usize],
    nmodes: usize,
) -> Result<FtmResult, Error> {
    let bad = |m: &str| Error::InvalidModel(format!("FTM: {m}"));
    if !(tube.r > 0.0 && tube.t > 0.0 && tube.l > 0.0 && tube.e > 0.0)
        || !tube.r.is_finite()
        || !tube.l.is_finite()
    {
        return Err(bad("R, t, L and E must be positive"));
    }
    if !(tube.nu > -1.0 && tube.nu < 0.5) {
        return Err(bad("Poisson's ratio must be between -1 and 0.5"));
    }
    if tube.t >= tube.r {
        return Err(bad("the wall must be thinner than the radius"));
    }
    if [s.n, s.m, s.t, s.v].iter().all(|v| *v == 0.0) {
        return Err(bad("no stress: give at least one of N, M, T, V"));
    }
    if [s.n, s.m, s.t, s.v].iter().any(|v| !v.is_finite()) {
        return Err(bad("a stress is not finite"));
    }
    // A tube held by one pin, or not at all, can swing as a rigid body. The Fourier series cannot
    // draw that rigid rotation exactly (a straight line along the tube), so the factorisation
    // would not always see it: refuse it outright.
    let held = |e: End| e != End::Free;
    if !(base == End::Clamped || top == End::Clamped || (held(base) && held(top))) {
        return Err(bad(
            "the tube is a mechanism under these end conditions: clamp one end, or pin both",
        ));
    }
    let mut js: Vec<usize> = long_terms.iter().copied().filter(|&j| j > 0).collect();
    js.sort_unstable();
    js.dedup();
    if p < 1 || js.is_empty() || nmodes < 1 {
        return Err(bad("at least one term each way, and one mode"));
    }
    if p > 80 || js.len() > 60 {
        return Err(bad("at most 80 circumferential and 60 longitudinal terms"));
    }
    let tabs = Tables::new(tube, p, &js);
    let (nc, nl) = (tabs.nc, tabs.nl);
    let (ke_forms, kg_forms) = forms(tube, s);
    let z = constraints(tube, &js, base, top, &tabs);
    let group = |d: usize, c: usize| d * nc + c;
    let ng = 3 * nc;

    // Which (displacement, circumferential term) groups couple: the uncoupled blocks solve apart.
    let mut dsu = Dsu((0..ng).collect());
    for f in ke_forms.iter().chain(&kg_forms) {
        for ta in &f.a {
            for tb in &f.b {
                for c1 in 0..nc {
                    for c2 in 0..nc {
                        if tabs.x(f.weight, ta.a, tb.a, c1, c2) != 0.0 {
                            dsu.union(group(ta.d, c1), group(tb.d, c2));
                        }
                    }
                }
            }
        }
    }
    let mut blocks: HashMap<usize, Vec<usize>> = HashMap::new();
    for g in 0..ng {
        let root = dsu.find(g);
        blocks.entry(root).or_default().push(g);
    }
    let mut blocks: Vec<Vec<usize>> = blocks.into_values().collect();
    // Harmonic-major order inside a block: Ke couples only terms of one harmonic, so it is block
    // diagonal, and Kg (bending, shear) only neighbouring harmonics, so it is banded.
    for b in &mut blocks {
        b.sort_by_key(|&g| ((g % nc).div_ceil(2), g % nc, g / nc));
    }
    blocks.sort();

    let mut found: Vec<(f64, Vec<f64>)> = vec![];
    let mut dofs = 0;
    for groups in &blocks {
        let pos: HashMap<usize, usize> = groups.iter().enumerate().map(|(i, &g)| (g, i)).collect();
        let n_full = groups.len() * nl;
        let assemble = |fs: &[Form]| {
            let mut k = vec![0.0; n_full * n_full];
            for f in fs {
                for ta in &f.a {
                    for tb in &f.b {
                        for c1 in 0..nc {
                            let Some(&i1) = pos.get(&group(ta.d, c1)) else {
                                continue;
                            };
                            for c2 in 0..nc {
                                let Some(&i2) = pos.get(&group(tb.d, c2)) else {
                                    continue;
                                };
                                let xv = tabs.x(f.weight, ta.a, tb.a, c1, c2);
                                if xv == 0.0 {
                                    continue;
                                }
                                let sc = f.scale * ta.c * tb.c * xv;
                                for l1 in 0..nl {
                                    let row = (i1 * nl + l1) * n_full + i2 * nl;
                                    for l2 in 0..nl {
                                        k[row + l2] += sc * tabs.y(ta.b, tb.b, l1, l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            k
        };
        let (ke, kg) = (assemble(&ke_forms), assemble(&kg_forms));
        // Into the supported space: block-diagonal Z, one block per group.
        let zs: Vec<&RMat> = groups
            .iter()
            .map(|g| &z[g / nc][(g % nc).div_ceil(2).min(2)])
            .collect();
        let offs: Vec<usize> = zs
            .iter()
            .scan(0, |acc, zb| {
                let o = *acc;
                *acc += zb.c;
                Some(o)
            })
            .collect();
        let m = offs.last().map_or(0, |o| o + zs[zs.len() - 1].c);
        if m == 0 {
            continue;
        }
        dofs += m;
        let reduce = |k: &[f64]| {
            let mut out = Mat::zeros(m);
            for (i1, z1) in zs.iter().enumerate() {
                for (i2, z2) in zs.iter().enumerate() {
                    // Z1ᵀ K12 Z2, skipping the (many) groups that do not couple.
                    let blank = (0..nl).all(|l1| {
                        let row = (i1 * nl + l1) * n_full + i2 * nl;
                        k[row..row + nl].iter().all(|v| *v == 0.0)
                    });
                    if blank {
                        continue;
                    }
                    let mut tmp = vec![0.0; nl * z2.c];
                    for l1 in 0..nl {
                        let row = (i1 * nl + l1) * n_full + i2 * nl;
                        for c in 0..z2.c {
                            let mut acc = 0.0;
                            for l2 in 0..nl {
                                acc += k[row + l2] * z2.get(l2, c);
                            }
                            tmp[l1 * z2.c + c] = acc;
                        }
                    }
                    for a in 0..z1.c {
                        for c in 0..z2.c {
                            let mut acc = 0.0;
                            for l1 in 0..nl {
                                acc += z1.get(l1, a) * tmp[l1 * z2.c + c];
                            }
                            out.add(offs[i1] + a, offs[i2] + c, acc);
                        }
                    }
                }
            }
            out.symmetrised()
        };
        let (ker, kgr) = (reduce(&ke), reduce(&kg));
        if kgr.max_abs() == 0.0 {
            continue;
        }
        let (l, start) = chol(&ker).ok_or_else(|| {
            bad("the tube is a mechanism under these end conditions (hold at least one end, or both ends pinned)")
        })?;
        let tp = top_positive(&l, &start, &kgr, nmodes);
        for (mu, y) in tp {
            // x = L⁻ᵀ y in the supported space, then back through Z.
            let mut x = y;
            solve_lt(&l, &start, &mut x);
            let mut coef = vec![0.0; 3 * nc * nl];
            for (i, &g) in groups.iter().enumerate() {
                let zb = zs[i];
                for li in 0..nl {
                    let mut acc = 0.0;
                    for c in 0..zb.c {
                        acc += zb.get(li, c) * x[offs[i] + c];
                    }
                    coef[g * nl + li] = acc;
                }
            }
            found.push((1.0 / mu, coef));
        }
    }
    if found.is_empty() {
        return Err(bad(
            "no buckling under these stresses (every load factor is negative or infinite)",
        ));
    }
    found.sort_by(|a, b| a.0.total_cmp(&b.0));
    found.truncate(nmodes);
    let modes = found
        .into_iter()
        .map(|(lf, mut coef)| {
            let big = coef.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            coef.iter_mut().for_each(|v| *v /= big);
            let mut share = vec![0.0; p + 1];
            for c in 0..nc {
                for l in 0..nl {
                    share[c.div_ceil(2)] += coef[(2 * nc + c) * nl + l].powi(2);
                }
            }
            let circ_waves = (0..=p)
                .max_by(|a, b| share[*a].total_cmp(&share[*b]))
                .unwrap_or(0);
            FtmMode {
                load_factor: lf,
                coef,
                circ_waves,
            }
        })
        .collect();
    Ok(FtmResult {
        tube: *tube,
        circ: p,
        long: js,
        modes,
        dofs,
    })
}

/// The first column of each row holding a non-zero at or left of the diagonal: the profile a
/// Cholesky factor keeps.
fn profile(a: &Mat) -> Vec<usize> {
    let n = a.n;
    (0..n)
        .map(|i| (0..=i).find(|&j| a.data[i * n + j] != 0.0).unwrap_or(i))
        .collect()
}

/// Lower Cholesky factor (row-major) of a symmetric positive-definite matrix, with its profile,
/// or `None`. Profile (skyline) elimination: fill stays inside each row's first non-zero.
fn chol(a: &Mat) -> Option<(Vec<f64>, Vec<usize>)> {
    let n = a.n;
    let start = profile(a);
    let mut l = a.data.clone();
    for i in 0..n {
        for j in start[i]..=i {
            let (ri, rj) = (i * n, j * n);
            let k0 = start[i].max(start[j]);
            let dot: f64 = if k0 < j {
                l[ri + k0..ri + j]
                    .iter()
                    .zip(&l[rj + k0..rj + j])
                    .map(|(x, y)| x * y)
                    .sum()
            } else {
                0.0
            };
            let v = l[ri + j] - dot;
            if i == j {
                if v.is_nan() || v <= 0.0 {
                    return None;
                }
                l[ri + i] = v.sqrt();
            } else {
                l[ri + j] = v / l[rj + j];
            }
        }
        l[i * n + i + 1..(i + 1) * n]
            .iter_mut()
            .for_each(|v| *v = 0.0);
    }
    Some((l, start))
}

/// `L x = b` in place.
fn solve_l(l: &[f64], start: &[usize], b: &mut [f64]) {
    let n = b.len();
    for i in 0..n {
        let s0 = start[i];
        let dot: f64 = l[i * n + s0..i * n + i]
            .iter()
            .zip(&b[s0..i])
            .map(|(x, y)| x * y)
            .sum();
        b[i] = (b[i] - dot) / l[i * n + i];
    }
}

/// `Lᵀ x = b` in place.
fn solve_lt(l: &[f64], start: &[usize], b: &mut [f64]) {
    let n = b.len();
    for i in (0..n).rev() {
        b[i] /= l[i * n + i];
        let bi = b[i];
        let s0 = start[i];
        for (bk, lk) in b[s0..i].iter_mut().zip(&l[i * n + s0..i * n + i]) {
            *bk -= lk * bi;
        }
    }
}

/// The `k` largest positive eigenvalues `μ` of `L⁻¹ Kg L⁻ᵀ`, with unit eigenvectors: dense below a
/// few hundred unknowns, Lanczos with full reorthogonalisation above.
fn top_positive(l: &[f64], start: &[usize], kg: &Mat, k: usize) -> Vec<(f64, Vec<f64>)> {
    let n = kg.n;
    // Kg's non-zero span in each row, for a banded product.
    let span: Vec<(usize, usize)> = kg
        .data
        .chunks_exact(n)
        .map(|row| {
            let lo = row.iter().position(|v| *v != 0.0).unwrap_or(0);
            let hi = row.iter().rposition(|v| *v != 0.0).map_or(0, |h| h + 1);
            (lo, hi.max(lo))
        })
        .collect();
    let op = |x: &[f64]| {
        let mut y = x.to_vec();
        solve_lt(l, start, &mut y);
        let mut z: Vec<f64> = kg
            .data
            .chunks_exact(n)
            .zip(&span)
            .map(|(row, &(lo, hi))| row[lo..hi].iter().zip(&y[lo..hi]).map(|(a, b)| a * b).sum())
            .collect();
        solve_l(l, start, &mut z);
        z
    };
    if n <= 400 {
        let mut a = Mat::zeros(n);
        for j in 0..n {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let col = op(&e);
            for (i, v) in col.iter().enumerate() {
                a.set(i, j, *v);
            }
        }
        let (_, vals, vecs) = sym_eigen_top(&a.symmetrised(), k.min(n));
        return vals
            .into_iter()
            .zip(vecs)
            .filter(|(mu, _)| *mu > 0.0)
            .collect();
    }
    let steps = n.min((4 * k + 80).max(120));
    let mut qv: Vec<Vec<f64>> = Vec::with_capacity(steps);
    let (mut alpha, mut beta) = (vec![], vec![]);
    let mut v: Vec<f64> = (0..n)
        .map(|i| 1.0 + ((i * 7919) % 101) as f64 / 101.0)
        .collect();
    let nv = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    v.iter_mut().for_each(|x| *x /= nv);
    for _ in 0..steps {
        let mut w = op(&v);
        let a: f64 = w.iter().zip(&v).map(|(x, y)| x * y).sum();
        qv.push(v.clone());
        alpha.push(a);
        // Full reorthogonalisation, twice.
        for _ in 0..2 {
            for qq in &qv {
                let c: f64 = w.iter().zip(qq).map(|(x, y)| x * y).sum();
                w.iter_mut().zip(qq).for_each(|(x, y)| *x -= c * y);
            }
        }
        let b = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if b < 1e-12 * a.abs().max(1e-300) {
            break;
        }
        beta.push(b);
        v = w.into_iter().map(|x| x / b).collect();
    }
    let m = alpha.len();
    let mut t = Mat::zeros(m);
    for i in 0..m {
        t.set(i, i, alpha[i]);
        if i + 1 < m {
            t.set(i, i + 1, beta[i]);
            t.set(i + 1, i, beta[i]);
        }
    }
    let (vals, vecs) = sym_eigen(&t);
    (0..m)
        .rev()
        .filter(|&j| vals[j] > 0.0)
        .take(k)
        .map(|j| {
            let mut y = vec![0.0; n];
            for (i, qq) in qv.iter().enumerate() {
                let s = vecs.get(i, j);
                y.iter_mut().zip(qq).for_each(|(a, b)| *a += s * b);
            }
            let ny = y.iter().map(|x| x * x).sum::<f64>().sqrt();
            y.iter_mut().for_each(|x| *x /= ny);
            (vals[j], y)
        })
        .collect()
}

/// The displacements `(u, v, w)` of a mode on an `nth × ny` grid, θ = 2πi/nth, y = L j/(ny − 1),
/// row-major over θ then y.
pub fn ftm_field(res: &FtmResult, mode: &FtmMode, nth: usize, ny: usize) -> Vec<[f64; 3]> {
    let (nc, nl) = (2 * res.circ + 1, 2 * res.long.len() + 1);
    let (r, l) = (res.tube.r, res.tube.l);
    let mut out = Vec::with_capacity(nth * ny);
    for i in 0..nth {
        let x = r * 2.0 * PI * i as f64 / nth as f64;
        let xs: Vec<f64> = (0..nc).map(|c| basis(c, 1.0 / r, x, 0)).collect();
        for j in 0..ny {
            let y = if ny > 1 {
                l * j as f64 / (ny - 1) as f64
            } else {
                0.0
            };
            let ys: Vec<f64> = (0..nl).map(|k| lbasis(k, &res.long, l, y, 0)).collect();
            let mut d = [0.0; 3];
            for (di, dv) in d.iter_mut().enumerate() {
                for (c, xv) in xs.iter().enumerate() {
                    let base = (di * nc + c) * nl;
                    let s: f64 = ys
                        .iter()
                        .enumerate()
                        .map(|(k, yv)| mode.coef[base + k] * yv)
                        .sum();
                    *dv += xv * s;
                }
            }
            out.push(d);
        }
    }
    out
}
