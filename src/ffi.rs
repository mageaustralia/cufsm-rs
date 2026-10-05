//! A dependency-free C-ABI so the same crate compiles to WebAssembly for the GUI. Behind the
//! `ffi` feature: off by default, so a crate that links cufsm-rs does not export these symbols.
//!
//! Build the module with `cargo rustc --release --lib --features ffi --crate-type cdylib
//! --target wasm32-unknown-unknown`.
//!
//! Everything crosses the boundary as flat `f64` buffers in linear memory: the caller
//! allocates with [`cufsm_alloc`], fills, calls, reads, then frees with [`cufsm_dealloc`]. No
//! `wasm-bindgen`, so the crate keeps its one promise (no dependencies), and the glue in the
//! page is a dozen lines of `WebAssembly.instantiate`.
//!
//! Layouts (all little-endian `f64` unless said otherwise):
//!
//! * `params`: `[terms, spaces, neigs]`, with `neigs` 1 or more. `neigs` is checked for both
//!   exports but only `cufsm_modes` uses it: `cufsm_signature` reports the lowest mode, so it
//!   solves for one.
//! * `mats`: 5 per material: `ex, ey, vx, vy, g`; an element refers to a material by 0-based row
//! * `nodes`: 7 per node, in CUFSM's column order: `x, z, free_x, free_z, free_y, free_q, stress`
//!   (the four free flags 0 or 1)
//! * `elems`: 4 per element: `node_i, node_j, t, mat`
//! * `springs`: 9 per spring (0 for none): `node_i, node_j` (-1 for ground), `k_u, k_v, k_w, k_q`,
//!   `local` (0 or 1), `discrete` (0 or 1), `ys_fraction` - CUFSM's `springs` row, v4.3 form
//! * `constraints`: 5 per constraint (0 for none): `node_e, dof_e, coeff, node_k, dof_k`, with
//!   dof codes 1 = x, 2 = z, 3 = y along the member, 4 = theta, and the row reading
//!   `u_e = coeff * u_k`
//! * `bc`: UTF-8 bytes of CUFSM's boundary condition string, `S-S`, `C-C`, `S-C`, `C-F`, `C-G`
//! * `lengths`: with the longitudinal terms `1..=terms` (1 or more) at each. For `S-S` with
//!   `terms` = 1 they are half-wavelengths: the signature curve. Otherwise they are physical member
//!   lengths, CUFSM's general boundary condition solution; with `S-S` the terms are uncoupled, so a
//!   length reports its lowest mode over 1 to `terms` half-waves (the signature curve's minimum over
//!   L, L/2, ...).
//! * `spaces`: the constrained curves wanted, bits 1 = G, 2 = D, 4 = L, 8 = O (0 to 15).
//! * `cufsm_signature` output: one row of `2 + popcount(spaces)` values per length, in input
//!   order: `L, λ`, then one column per requested space in G, D, L, O order. A length with no
//!   positive load factor writes `NaN`.
//! * `cufsm_ftm`: the finite tube method (`crate::ftm`) on one tube. `params`:
//!   `[R, t, L, E, nu, N, M, T, V, base, top, p, nmodes, nth, ny]`, the actions in N and N·mm
//!   (compression and the moment's compression side at θ = 0 positive), ends 0 free edge,
//!   1 pinned, 2 clamped, 3 free with a stiff ring, `p` circumferential harmonics; `terms` the longitudinal wave numbers. Output:
//!   `[modes found, unknowns, σN, σM, τT, τV]` (the reference stresses), then per mode
//!   `[λ, circumferential waves]` and its `(u, v, w)` on an `nth × ny` grid (θ = 2πi/nth,
//!   y = L j/(ny − 1), θ outer), scaled so the largest coefficient is 1.
//! * `cufsm_modes` output: per length, a row of `2 + nterms + neigs * (5 + 4 * nodes * nterms)`
//!   values: `found, nterms`, the `nterms` longitudinal terms `m_1..m_nterms`, then `neigs`
//!   blocks of `(λ, G, D, L, O, 4 * nodes * nterms dofs)` (the dofs in CUFSM's order, per
//!   term: `u`/`v` interleaved, then `w`/`θ` interleaved). Blocks past `found` are `NaN`, so a
//!   length short of `neigs` positive eigenvalues reports how many it found in `found`.
//! * `cufsm_props`: gross and warping section properties, no `args`. Output 15 values:
//!   `A, xcg, zcg, Ixx, Izz, Ixz, thetap_deg, I11, I22, J, xs, zs, Cw, B1, B2`.
//! * `cufsm_stresgen`: `args` is `[P, Mxx, Mzz, M11, M22, B, restrained]`, with `restrained` 0
//!   or 1 (`unsymmetric = !restrained`). Writes one reference stress per node; the nodes' own
//!   stress column in the input is ignored.
//! * `cufsm_yield`: `args` is `[fy, restrained, extreme_fibre]`, both flags 0 or 1. Writes 6
//!   values: `Py, Mxx, Mzz, M11, M22, By`.
//! * `cufsm_stress_to_action`: no `args`; reads the nodes' stress column. Writes 5 values:
//!   `P, M11, M22, B, err`.
//!
//! Added in ABI 2, minor 1 ([`cufsm_abi_minor`]); every layout above is unchanged:
//!
//! * `cufsm_strip`: `params` is `[terms, spaces, neigs, classify]`, plus an `m_all` buffer after
//!   `lengths`. With `m_all` empty (length 0) every length takes the terms `1..=terms`, as
//!   `cufsm_modes`; otherwise `terms` must be 0 and `m_all` is one list per length, each its count
//!   then its terms: `[n_1, m.., n_2, m.., ...]`, any finite values (CUFSM's `msort` sorts them,
//!   drops zeros and repeats, and the row reports the terms it used). `spaces` 0 is the free
//!   analysis; otherwise its bits (1 = G, 2 = D, 4 = L, 8 = O) together are the one cFSM space the
//!   analysis is restricted to (`stripmain_constrained`; unlike `cufsm_signature`, where each bit
//!   is its own column). `classify` 1 writes CUFSM's default classification of each mode (as
//!   `cufsm_modes` does), 0 writes `NaN` there and skips the work. Output: one `cufsm_modes` row
//!   per length, back to back; a row's size, `2 + nterms + neigs * (5 + 4 * nodes * nterms)`,
//!   follows its own `nterms` (after `msort`). With `[terms, 0, neigs, 1]` and no `m_all` the
//!   output is `cufsm_modes`', value for value. The bc string may be lower case here.
//! * `cufsm_classify`: `params` is `[orth, norm, ospace]` in CUFSM's codes (orth 1 natural,
//!   2 axial, 3 load; norm 0 none, 1 vector, 2 strain energy, 3 work; ospace 1 ST, 2 `K⁻¹`,
//!   3 `Kg⁻¹`, 4 the null space; CUFSM's defaults are 2, 1, 1), then the model (orth 3 reads the
//!   nodes' stress column), the bc, a `results` buffer, springs and constraints. `results` holds per
//!   length `L, nterms, m_1..m_nterms, k`, then `k` modes of `4 * nodes * nterms` dofs each (the
//!   terms used as given, as `cufsm_strip` reports them). Writes `G, D, L, O` percent per mode, in
//!   order: `4 * (total modes)` values.
//! * `cufsm_template`: `params` is the 23 values of CUFSM's `templatecalc`: `shape` (1 C, 2 Z),
//!   `h, b1, b2, d1, d2, r1, r2, r3, r4, q1, q2, t, nh, nb1, nb2, nd1, nd2, nr1, nr2, nr3, nr4,
//!   centerline` (q in degrees; the strip counts whole, `nh`, `nb1`, `nb2` at least 1; centerline
//!   0 for outside dimensions and inside radii). Writes `nnodes, nelems`, then the nodes and
//!   elements in this interface's own `nodes` and `elems` layouts (0-based, every node free,
//!   stress 1, material 0), ready to pass back in. `out` must hold the bound
//!   `2 + 7 * (s + 1) + 4 * s` with `s` the sum of the nine strip counts, checked before meshing.
//! * `cufsm_props_wn`: `cufsm_props`' 15 values, then `cutwp_prop2`'s warping function `wn` at each
//!   node: `15 + nodes` values.
//! * `cufsm_signature_lengths`: the 100 half-wavelengths of CUFSM's `signature_ss.m` for the model.
//! * `cufsm_signature_minima`: `curve` is `L, λ` pairs (a `NaN` λ marks a length with no positive
//!   load factor); writes `L, λ` per interior local minimum, shortest first, refined as
//!   [`crate::signature_minima`]: at most `curve_len - 4` values.
//!
//! No export caps `terms` or `neigs` (the old [`MAX_TERMS`] and [`MAX_NEIGS`]): every output is
//! sized from the caller's numbers and checked against `out_cap` before any work. The limit left
//! is memory: the analysis holds dense matrices of `(4 * nodes * nterms)²` values, and an
//! allocation the platform cannot make aborts (a trap on wasm32), as it would in any caller of the
//! crate; a size this platform cannot even address is refused with a message.
//!
//! Every export returns the number of `f64` written, or a negative number on failure, when
//! [`cufsm_last_error_ptr`]/[`cufsm_last_error_len`] hold the message in UTF-8 (valid until the
//! next failing call; the slot is one per process, so native callers on several threads must
//! serialise). A panic inside the analysis is caught and reported the same way on native
//! targets; on `wasm32-unknown-unknown`, which aborts on panic, it traps the instance.

#![allow(unsafe_code)]

use std::alloc::{alloc, dealloc, Layout};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::cfsm::{self, stripmain_constrained, Spaces};
use crate::model::{Constraint, Dof, Spring};
use crate::{stripmain, BoundaryCondition, Element, Material, Model, Node};

const N_PARAMS: usize = 3;
const N_MAT: usize = 5;
const N_NODE: usize = 7;
const N_ELEM: usize = 4;
const N_SPRING: usize = 9;
const N_CONSTRAINT: usize = 5;
/// The most longitudinal terms the interface once accepted. No longer enforced: the output is
/// sized from the caller's own numbers, so the only limit left is memory (see the module docs).
pub const MAX_TERMS: usize = 100;
/// The most eigenvalues per length the interface once accepted. No longer enforced, as
/// [`MAX_TERMS`].
pub const MAX_NEIGS: usize = 50;

static ERR_PTR: AtomicUsize = AtomicUsize::new(0);
static ERR_LEN: AtomicUsize = AtomicUsize::new(0);

fn set_error(msg: &str) {
    let old = ERR_PTR.swap(0, Ordering::Relaxed);
    let old_len = ERR_LEN.swap(0, Ordering::Relaxed);
    if old != 0 {
        // Freed as what it was allocated as: a boxed byte slice.
        drop(unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(old as *mut u8, old_len)) });
    }
    let boxed: Box<[u8]> = msg.as_bytes().into();
    let len = boxed.len();
    let ptr = Box::into_raw(boxed) as *mut u8;
    ERR_PTR.store(ptr as usize, Ordering::Relaxed);
    ERR_LEN.store(len, Ordering::Relaxed);
}

/// Runs an export's body, turning an error or a panic into a negative return and a message.
fn guarded(body: impl FnOnce() -> Result<usize, String>) -> isize {
    let msg = match catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(n)) => return n as isize,
        Ok(Err(msg)) => msg,
        Err(p) => {
            let what = p
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| p.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown".into());
            format!("internal error (panic): {what}")
        }
    };
    set_error(&msg);
    -1
}

/// The layout [`cufsm_alloc`] and [`cufsm_dealloc`] agree on: 8-aligned, at least 8 bytes.
fn block(n: usize) -> Option<Layout> {
    Layout::from_size_align(n.max(1).div_ceil(8) * 8, 8).ok()
}

/// `n` bytes of 8-aligned memory for the caller to fill, freed with [`cufsm_dealloc`] and the
/// same `n`. Null when the memory cannot be had.
#[no_mangle]
pub extern "C" fn cufsm_alloc(n: usize) -> *mut u8 {
    match block(n) {
        Some(l) => unsafe { alloc(l) },
        None => std::ptr::null_mut(),
    }
}

/// Frees a block from [`cufsm_alloc`]. `n` is the byte count passed to it.
///
/// # Safety
/// `p` must come from `cufsm_alloc(n)` and not be freed twice.
#[no_mangle]
pub unsafe extern "C" fn cufsm_dealloc(p: *mut u8, n: usize) {
    if let (false, Some(l)) = (p.is_null(), block(n)) {
        // SAFETY: the caller passes a block from cufsm_alloc(n), which allocated this layout.
        unsafe { dealloc(p, l) };
    }
}

/// Pointer to the last error's UTF-8 bytes (0 when there has been no error).
#[no_mangle]
pub extern "C" fn cufsm_last_error_ptr() -> usize {
    ERR_PTR.load(Ordering::Relaxed)
}

/// Length of the last error's UTF-8 bytes.
#[no_mangle]
pub extern "C" fn cufsm_last_error_len() -> usize {
    ERR_LEN.load(Ordering::Relaxed)
}

/// The buffer layouts' version. 2 = materials buffer, 7-value nodes, 4-value elements,
/// `[terms, spaces, neigs]` params, multi-mode `cufsm_modes`, and the section exports. It changes
/// only when an existing export's layout changes; exports added alongside raise
/// [`cufsm_abi_minor`] instead, so a page built for version 2 keeps working.
#[no_mangle]
pub extern "C" fn cufsm_abi_version() -> u32 {
    2
}

/// The exports added within [`cufsm_abi_version`] 2, which leave every earlier layout as it was.
/// 1 = `cufsm_strip`, `cufsm_classify`, `cufsm_template`, `cufsm_props_wn`,
/// `cufsm_signature_lengths` and `cufsm_signature_minima`, and no fixed caps on `terms` or `neigs`.
/// A page that needs those checks for a minor of at least 1 (a module without this export is 0).
#[no_mangle]
pub extern "C" fn cufsm_abi_minor() -> u32 {
    1
}

/// `len` values at `p`, refusing a null pointer (reading through one is undefined behaviour).
unsafe fn input<'a>(p: *const f64, len: usize, what: &str) -> Result<&'a [f64], String> {
    if p.is_null() {
        return Err(format!("no {what} given"));
    }
    // SAFETY: not null (checked above); the caller guarantees `len` readable values at `p`.
    Ok(unsafe { slice::from_raw_parts(p, len) })
}

/// A whole number in `lo..=hi`, or why not.
fn whole(v: f64, lo: usize, hi: usize, what: &str) -> Result<usize, String> {
    if v.is_finite() && v.fract() == 0.0 && v >= lo as f64 && v <= hi as f64 {
        Ok(v as usize)
    } else {
        Err(format!(
            "{what} must be a whole number from {lo} to {hi}, got {v}"
        ))
    }
}

/// `a * b`, refusing an overflow: every size here comes from the caller's numbers.
fn mul(a: usize, b: usize) -> Result<usize, String> {
    a.checked_mul(b)
        .ok_or_else(|| "the requested output is larger than this platform can address".into())
}

/// `a + b`, refusing an overflow, as [`mul`].
fn add(a: usize, b: usize) -> Result<usize, String> {
    a.checked_add(b)
        .ok_or_else(|| "the requested output is larger than this platform can address".into())
}

/// The values one length's row of modes takes: `found, nterms`, the terms, then `neigs` blocks of
/// `(λ, G, D, L, O, 4 * nodes * nterms dofs)`.
fn modes_row(nodes: usize, nterms: usize, neigs: usize) -> Result<usize, String> {
    let blk = add(5, mul(4, mul(nodes, nterms)?)?)?;
    add(add(2, nterms)?, mul(neigs, blk)?)
}

/// Refuses a problem whose dense `n x n` matrix (`n = 4 * nodes * nterms`) cannot even be
/// addressed on this platform (wasm32 above about 23,000 unknowns), rather than letting the
/// allocation abort. Below that the limit is the memory actually available, as in any caller.
fn dense_fits(nodes: usize, nterms: usize) -> Result<(), String> {
    let n = mul(4, mul(nodes, nterms)?)?;
    match n.checked_mul(n).and_then(|v| v.checked_mul(8)) {
        Some(_) => Ok(()),
        None => Err(format!(
            "{nterms} longitudinal terms on {nodes} nodes is {n} unknowns: a dense {n} x {n} \
             matrix is larger than this platform can address"
        )),
    }
}

/// The lengths, every one a positive finite number.
unsafe fn read_lengths(lengths: *const f64, lengths_len: usize) -> Result<Vec<f64>, String> {
    if lengths_len == 0 {
        return Err("no lengths given".into());
    }
    // SAFETY: the caller guarantees `lengths_len` readable values at `lengths`.
    let lens = unsafe { input(lengths, lengths_len, "lengths") }?;
    if lens.iter().any(|l| !(l.is_finite() && *l > 0.0)) {
        return Err("a length is not a positive finite number".into());
    }
    Ok(lens.to_vec())
}

/// CUFSM's boundary condition string from UTF-8 bytes. `upper` accepts lower case too, as the
/// Python package does (the original exports keep their exact-case parse).
unsafe fn read_bc(bc: *const u8, bc_len: usize, upper: bool) -> Result<BoundaryCondition, String> {
    if bc.is_null() || bc_len == 0 {
        return Err("no boundary condition given".into());
    }
    // SAFETY: not null and non-empty (checked above); the caller guarantees `bc_len` bytes at `bc`.
    let s = std::str::from_utf8(unsafe { slice::from_raw_parts(bc, bc_len) })
        .map_err(|_| "boundary condition is not UTF-8".to_string())?;
    let parsed = if upper {
        BoundaryCondition::parse(&s.to_ascii_uppercase())
    } else {
        BoundaryCondition::parse(s)
    };
    parsed.ok_or_else(|| format!("unknown boundary condition {s:?}"))
}

/// Everything an export needs, read and checked from the flat buffers. Every length is checked
/// first: a short buffer is a caller bug, and reading past it would be undefined behaviour.
struct Input {
    model: Model,
    bc: BoundaryCondition,
    lengths: Vec<f64>,
    m_all: Vec<Vec<f64>>,
    spaces: u32,
    neigs: usize,
}

/// The model from its flat buffers: materials, 7-value nodes, 4-value elements, springs and
/// constraints. Every stride is checked before anything is read.
#[allow(clippy::too_many_arguments)]
unsafe fn read_model(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
) -> Result<Model, String> {
    if mats_len == 0 || mats_len % N_MAT != 0 {
        return Err(format!(
            "mats length {mats_len} is not a positive multiple of {N_MAT}"
        ));
    }
    if nodes_len == 0 || nodes_len % N_NODE != 0 {
        return Err(format!(
            "nodes length {nodes_len} is not a positive multiple of {N_NODE}"
        ));
    }
    if elems_len == 0 || elems_len % N_ELEM != 0 {
        return Err(format!(
            "elems length {elems_len} is not a positive multiple of {N_ELEM}"
        ));
    }
    // SAFETY (these three): the caller guarantees each pointer is valid for its length.
    let raw_mats = unsafe { input(mats, mats_len, "mats") }?;
    let raw_nodes = unsafe { input(nodes, nodes_len, "nodes") }?;
    let raw_elems = unsafe { input(elems, elems_len, "elems") }?;

    let mut materials = Vec::with_capacity(mats_len / N_MAT);
    for (k, c) in raw_mats.chunks_exact(N_MAT).enumerate() {
        if !c.iter().all(|v| v.is_finite()) || c[0] <= 0.0 || c[1] <= 0.0 || c[4] <= 0.0 {
            return Err(format!(
                "material {k}: Ex, Ey and G must be positive finite numbers"
            ));
        }
        if c[2] * c[3] >= 1.0 {
            return Err(format!(
                "material {k}: vx * vy must be below 1, got {}",
                c[2] * c[3]
            ));
        }
        materials.push(Material {
            ex: c[0],
            ey: c[1],
            vx: c[2],
            vy: c[3],
            g: c[4],
        });
    }
    let flag = |v: f64| {
        if v == 0.0 || v == 1.0 {
            Ok(v == 1.0)
        } else {
            Err(format!("a node's free flag must be 0 or 1, got {v}"))
        }
    };
    let mut out_nodes = Vec::with_capacity(nodes_len / N_NODE);
    for c in raw_nodes.chunks_exact(N_NODE) {
        if !(c[0].is_finite() && c[1].is_finite() && c[6].is_finite()) {
            return Err("a node coordinate or stress is not finite".into());
        }
        let mut n = Node::new(c[0], c[1], c[6]);
        n.free = [flag(c[2])?, flag(c[3])?, flag(c[4])?, flag(c[5])?];
        out_nodes.push(n);
    }
    let n_nodes = out_nodes.len();
    let mut out_elems = Vec::with_capacity(elems_len / N_ELEM);
    for (k, c) in raw_elems.chunks_exact(N_ELEM).enumerate() {
        if !(c[2].is_finite() && c[2] > 0.0) {
            return Err(format!(
                "element {k}: thickness must be a positive finite number"
            ));
        }
        out_elems.push(Element {
            ni: whole(c[0], 0, n_nodes - 1, "an element's node_i")?,
            nj: whole(c[1], 0, n_nodes - 1, "an element's node_j")?,
            t: c[2],
            mat: whole(c[3], 0, materials.len() - 1, "an element's material")?,
        });
    }
    let springs_raw = if springs_len == 0 {
        &[][..]
    } else {
        unsafe { input(springs, springs_len, "springs") }?
    };
    let constraints_raw = if constraints_len == 0 {
        &[][..]
    } else {
        unsafe { input(constraints, constraints_len, "constraints") }?
    };
    let model = Model {
        materials,
        nodes: out_nodes,
        elements: out_elems,
        constraints: read_constraints(constraints_raw, n_nodes)?,
        springs: read_springs(springs_raw, n_nodes)?,
    };
    model.validate().map_err(|err| err.to_string())?;
    Ok(model)
}

#[allow(clippy::too_many_arguments)]
unsafe fn read_input(
    params: *const f64,
    params_len: usize,
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    bc: *const u8,
    bc_len: usize,
    lengths: *const f64,
    lengths_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
) -> Result<Input, String> {
    if params_len < N_PARAMS {
        return Err(format!("params needs {N_PARAMS} values, got {params_len}"));
    }
    // SAFETY (this and the two readers below): the caller guarantees each pointer is valid for
    // its length.
    let p = unsafe { input(params, N_PARAMS, "params") }?;
    let lens = unsafe { read_lengths(lengths, lengths_len) }?;
    let terms = p[0];
    let spaces = whole(p[1], 0, 15, "spaces")? as u32;
    let neigs = whole(p[2], 1, usize::MAX, "neigs")?;
    let bc = unsafe { read_bc(bc, bc_len, false) }?;
    // Any end condition takes one or more longitudinal terms. With S-S and one term the lengths
    // are half-wavelengths (the signature curve); with S-S and m = 1..n they are member lengths and
    // the terms are uncoupled, so each length reports its lowest mode over 1 to n half-waves, as
    // CUFSM's general boundary condition solution does.
    let terms = whole(terms, 1, usize::MAX, "terms")?;

    // SAFETY: the export's contract: every pointer valid for its length.
    let model = unsafe {
        read_model(
            mats,
            mats_len,
            nodes,
            nodes_len,
            elems,
            elems_len,
            springs,
            springs_len,
            constraints,
            constraints_len,
        )
    }?;
    dense_fits(model.nodes.len(), terms)?;
    let m: Vec<f64> = (1..=terms).map(|k| k as f64).collect();
    Ok(Input {
        model,
        bc,
        m_all: vec![m; lens.len()],
        lengths: lens,
        spaces,
        neigs,
    })
}

/// A dof code for a constraint row: 1 = x, 2 = z, 3 = y along the member, 4 = theta.
fn dof_of(v: f64, what: &str) -> Result<Dof, String> {
    match v {
        1.0 => Ok(Dof::X),
        2.0 => Ok(Dof::Z),
        3.0 => Ok(Dof::Y),
        4.0 => Ok(Dof::Theta),
        _ => Err(format!("{what} must be a dof code from 1 to 4, got {v}")),
    }
}

/// Springs from their 9-value rows: `node_i, node_j` (-1 for ground), stiffnesses, flags.
fn read_springs(raw: &[f64], n_nodes: usize) -> Result<Vec<Spring>, String> {
    if raw.len() % N_SPRING != 0 {
        return Err(format!(
            "springs length {} is not a multiple of {N_SPRING}",
            raw.len()
        ));
    }
    let flag = |v: f64, what: &str| match v {
        0.0 => Ok(false),
        1.0 => Ok(true),
        _ => Err(format!("{what} must be 0 or 1, got {v}")),
    };
    let mut out = Vec::with_capacity(raw.len() / N_SPRING);
    for c in raw.chunks_exact(N_SPRING) {
        if !c[2..=5].iter().all(|v| v.is_finite()) {
            return Err("a spring stiffness is not finite".into());
        }
        if !c[8].is_finite() {
            return Err("a spring's ys_fraction is not finite".into());
        }
        out.push(Spring {
            ni: whole(c[0], 0, n_nodes - 1, "a spring's node_i")?,
            nj: if c[1] == -1.0 {
                None
            } else {
                Some(whole(c[1], 0, n_nodes - 1, "a spring's node_j")?)
            },
            ku: c[2],
            kv: c[3],
            kw: c[4],
            kq: c[5],
            local: flag(c[6], "a spring's local flag")?,
            discrete: flag(c[7], "a spring's discrete flag")?,
            ys_fraction: c[8],
        });
    }
    Ok(out)
}

/// Equation constraints from their 5-value rows: `node_e, dof_e, coeff, node_k, dof_k`,
/// reading `u_e = coeff * u_k`.
fn read_constraints(raw: &[f64], n_nodes: usize) -> Result<Vec<Constraint>, String> {
    if raw.len() % N_CONSTRAINT != 0 {
        return Err(format!(
            "constraints length {} is not a multiple of {N_CONSTRAINT}",
            raw.len()
        ));
    }
    let mut out = Vec::with_capacity(raw.len() / N_CONSTRAINT);
    for c in raw.chunks_exact(N_CONSTRAINT) {
        if !c[2].is_finite() {
            return Err("a constraint's coeff is not finite".into());
        }
        out.push(Constraint {
            node_e: whole(c[0], 0, n_nodes - 1, "a constraint's node_e")?,
            dof_e: dof_of(c[1], "a constraint's dof_e")?,
            coeff: c[2],
            node_k: whole(c[3], 0, n_nodes - 1, "a constraint's node_k")?,
            dof_k: dof_of(c[4], "a constraint's dof_k")?,
        });
    }
    Ok(out)
}

/// `out` as a mutable slice of `need` values, refusing a null or short buffer.
unsafe fn output<'a>(out: *mut f64, out_cap: usize, need: usize) -> Result<&'a mut [f64], String> {
    if out.is_null() {
        return Err("no output buffer".into());
    }
    if out_cap < need {
        return Err(format!(
            "output buffer is too small: {need} values needed, {out_cap} given"
        ));
    }
    // SAFETY: not null and `need <= out_cap` (checked above); the caller guarantees `out_cap`
    // writable values at `out`, not aliased while the export runs.
    Ok(unsafe { slice::from_raw_parts_mut(out, need) })
}

fn first_lf(r: &crate::LengthResult) -> f64 {
    r.load_factors.first().copied().unwrap_or(f64::NAN)
}

/// The curve at the caller's lengths: always the free analysis, plus one column per space
/// asked for in `params[1]` from `stripmain_constrained`. See the module docs for the layouts.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_signature(
    params: *const f64,
    params_len: usize,
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    bc: *const u8,
    bc_len: usize,
    lengths: *const f64,
    lengths_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let i = unsafe {
            read_input(
                params,
                params_len,
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                bc,
                bc_len,
                lengths,
                lengths_len,
                springs,
                springs_len,
                constraints,
                constraints_len,
            )
        }?;
        let stride = 2 + i.spaces.count_ones() as usize;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, i.lengths.len() * stride) }?;
        let free =
            // One mode per length: the curve reports only the lowest, so neigs (for cufsm_modes)
            // would only multiply the eigen work here.
            stripmain(&i.model, &i.lengths, &i.m_all, i.bc, 1).map_err(|e| e.to_string())?;
        for (row, r) in buf.chunks_exact_mut(stride).zip(&free) {
            row[0] = r.length;
            row[1] = first_lf(r);
        }
        let mut col = 2;
        for bit in [1u32, 2, 4, 8] {
            if i.spaces & bit == 0 {
                continue;
            }
            let sp = Spaces {
                global: bit == 1,
                distortional: bit == 2,
                local: bit == 4,
                other: bit == 8,
            };
            let run = stripmain_constrained(&i.model, &i.lengths, &i.m_all, i.bc, 1, sp)
                .map_err(|e| e.to_string())?;
            for (row, r) in buf.chunks_exact_mut(stride).zip(&run) {
                row[col] = first_lf(r);
            }
            col += 1;
        }
        Ok(buf.len())
    })
}

/// `neigs` modes at each requested length, each with its load factor, its G/D/L/O classification
/// percentages and its shape.
///
/// # Safety
/// Every pointer must be valid for its length; `out` needs, per length,
/// `2 + nterms + neigs * (5 + 4 * nodes * nterms)` values.
#[no_mangle]
pub unsafe extern "C" fn cufsm_modes(
    params: *const f64,
    params_len: usize,
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    bc: *const u8,
    bc_len: usize,
    lengths: *const f64,
    lengths_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let i = unsafe {
            read_input(
                params,
                params_len,
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                bc,
                bc_len,
                lengths,
                lengths_len,
                springs,
                springs_len,
                constraints,
                constraints_len,
            )
        }?;
        let nn = i.model.nodes.len();
        let nt = i.m_all[0].len();
        let per = modes_row(nn, nt, i.neigs)?;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, mul(i.lengths.len(), per)?) }?;
        let results =
            stripmain(&i.model, &i.lengths, &i.m_all, i.bc, i.neigs).map_err(|e| e.to_string())?;
        let classes = default_classes(&i.model, &results, i.bc)?;
        for ((row, r), cls) in buf.chunks_exact_mut(per).zip(&results).zip(&classes) {
            write_modes_row(row, r, Some(cls), i.neigs);
        }
        Ok(buf.len())
    })
}

/// CUFSM's default classification (axial orthogonality, vector norm, the ST O space) of every
/// mode, as `cufsm_modes` reports it.
fn default_classes(
    model: &Model,
    results: &[crate::LengthResult],
    bc: BoundaryCondition,
) -> Result<Vec<Vec<[f64; 4]>>, String> {
    cfsm::classify(model, results, bc, cfsm::Orth::Axial, cfsm::Norm::Vector)
        .map_err(|e| e.to_string())
}

/// One length's row of modes, sized by [`modes_row`] for `r.m_terms.len()` terms: `found, nterms`,
/// the terms, then `neigs` blocks of `(λ, G, D, L, O, dofs)`, `NaN` past `found` (and in place of
/// G, D, L, O without a classification).
fn write_modes_row(
    row: &mut [f64],
    r: &crate::LengthResult,
    cls: Option<&Vec<[f64; 4]>>,
    neigs: usize,
) {
    let nt = r.m_terms.len();
    let blk = (row.len() - 2 - nt) / neigs;
    let found = r.load_factors.len().min(neigs);
    row[0] = found as f64;
    row[1] = nt as f64;
    row[2..2 + nt].copy_from_slice(&r.m_terms);
    for (k, b) in row[2 + nt..].chunks_exact_mut(blk).enumerate() {
        if k < found {
            b[0] = r.load_factors[k];
            let c = cls.and_then(|c| c.get(k).copied());
            b[1..5].copy_from_slice(&c.unwrap_or([f64::NAN; 4]));
            b[5..].copy_from_slice(&r.modes[k]);
        } else {
            b.fill(f64::NAN);
        }
    }
}

/// Gross and warping section properties; see the module docs for the layout.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_props(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, 15) }?;
        let g = crate::grosprop(&m);
        let c = crate::cutwp_prop2(&m);
        buf.copy_from_slice(&[
            g.a, g.xcg, g.zcg, g.ixx, g.izz, g.ixz, g.thetap, g.i11, g.i22, c.j, c.xs, c.zs, c.cw,
            c.b1, c.b2,
        ]);
        Ok(15)
    })
}

/// Reference stresses from member actions, CUFSM's loading panel; see the module docs.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_stresgen(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    actions: *const f64,
    actions_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let mut m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        if actions_len != 7 {
            return Err(format!(
                "actions needs 7 values (P, Mxx, Mzz, M11, M22, B, restrained), got {actions_len}"
            ));
        }
        // SAFETY: the export's contract: `actions` valid for 7 values (checked just above).
        let a = unsafe { input(actions, 7, "actions") }?;
        if !a[..6].iter().all(|v| v.is_finite()) {
            return Err("an action is not finite".into());
        }
        let restrained = match a[6] {
            0.0 => false,
            1.0 => true,
            v => return Err(format!("restrained must be 0 or 1, got {v}")),
        };
        let n = m.nodes.len();
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, n) }?;
        let g = crate::grosprop(&m);
        crate::stresgen(
            &mut m,
            &crate::Actions {
                p: a[0],
                mxx: a[1],
                mzz: a[2],
                m11: a[3],
                m22: a[4],
            },
            &g,
            !restrained,
        );
        if a[5] != 0.0 {
            let c = crate::cutwp_prop2(&m);
            crate::add_bimoment_stress(&mut m, a[5], c.cw, &c.wn);
        }
        for (o, nd) in buf.iter_mut().zip(&m.nodes) {
            *o = nd.stress;
        }
        Ok(n)
    })
}

/// First-yield actions; see the module docs.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_yield(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    args: *const f64,
    args_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        if args_len != 3 {
            return Err(format!(
                "yield args needs 3 values (fy, restrained, extreme_fibre), got {args_len}"
            ));
        }
        // SAFETY: the export's contract: `args` valid for 3 values (checked just above).
        let a = unsafe { input(args, 3, "yield args") }?;
        if !(a[0].is_finite() && a[0] > 0.0) {
            return Err("fy must be a positive finite number".into());
        }
        let bit = |v: f64, what: &str| match v {
            0.0 => Ok(false),
            1.0 => Ok(true),
            _ => Err(format!("{what} must be 0 or 1, got {v}")),
        };
        let restrained = bit(a[1], "restrained")?;
        let ext = bit(a[2], "extreme_fibre")?;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, 6) }?;
        let g = crate::grosprop(&m);
        let y = if ext {
            crate::yield_mp_extfiber(&m, a[0], &g, !restrained)
        } else {
            crate::yield_mp(&m, a[0], &g, !restrained)
        };
        let c = crate::cutwp_prop2(&m);
        buf.copy_from_slice(&[
            y.py,
            y.mxx,
            y.mzz,
            y.m11,
            y.m22,
            crate::yield_b(a[0], c.cw, &c.wn),
        ]);
        Ok(6)
    })
}

/// Generate from Stress: actions fitted to the nodes' stress column; see the module docs.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_stress_to_action(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, 5) }?;
        let g = crate::grosprop(&m);
        let c = crate::cutwp_prop2(&m);
        let s = crate::stress_to_action(&m, &g, c.cw, &c.wn);
        buf.copy_from_slice(&[s.p, s.m11, s.m22, s.b, s.err]);
        Ok(5)
    })
}

/// The finite tube method on one tube; see the module docs for the layouts.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_ftm(
    params: *const f64,
    params_len: usize,
    terms: *const f64,
    terms_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    use crate::ftm::{ftm_buckle, ftm_field, End, Stresses, Tube};
    guarded(|| {
        if params_len < 15 {
            return Err(format!("params needs 15 values, got {params_len}"));
        }
        // SAFETY: the export's contract: every pointer valid for its length.
        let p = unsafe { input(params, 15, "params") }?;
        if terms_len == 0 {
            return Err("no longitudinal terms given".into());
        }
        // SAFETY: as above.
        let js = unsafe { input(terms, terms_len, "terms") }?
            .iter()
            .map(|&j| whole(j, 1, 10_000, "a longitudinal wave number"))
            .collect::<Result<Vec<_>, _>>()?;
        if p[..9].iter().any(|v| !v.is_finite()) {
            return Err("a tube dimension, material value or action is not finite".into());
        }
        let end = |v: f64, what: &str| -> Result<End, String> {
            Ok(match whole(v, 0, 3, what)? {
                0 => End::Free,
                1 => End::Pinned,
                2 => End::Clamped,
                _ => End::Ring,
            })
        };
        let tube = Tube {
            r: p[0],
            t: p[1],
            l: p[2],
            e: p[3],
            nu: p[4],
        };
        let (base, top) = (end(p[9], "the base end")?, end(p[10], "the top end")?);
        let harmonics = whole(p[11], 1, 80, "circumferential harmonics")?;
        let nmodes = whole(p[12], 1, 20, "modes")?;
        let (nth, ny) = (
            whole(p[13], 4, 720, "θ points")?,
            whole(p[14], 2, 720, "y points")?,
        );
        let s = Stresses::from_actions(&tube, p[5], p[6], p[7], p[8]);
        let res =
            ftm_buckle(&tube, &s, base, top, harmonics, &js, nmodes).map_err(|e| e.to_string())?;
        let per = 2 + 3 * nth * ny;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, 6 + nmodes * per) }?;
        buf[..6].copy_from_slice(&[res.modes.len() as f64, res.dofs as f64, s.n, s.m, s.t, s.v]);
        let mut k = 6;
        for mode in &res.modes {
            buf[k] = mode.load_factor;
            buf[k + 1] = mode.circ_waves as f64;
            k += 2;
            for d in ftm_field(&res, mode, nth, ny) {
                buf[k..k + 3].copy_from_slice(&d);
                k += 3;
            }
        }
        Ok(k)
    })
}

/// Per-length term lists from `[n_1, m_1..m_n1, n_2, ...]`: exactly one list per length, each
/// term finite, and at least one nonzero term in each list (CUFSM's `msort` drops zeros).
fn read_m_all(raw: &[f64], nlen: usize) -> Result<Vec<Vec<f64>>, String> {
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(nlen);
    let mut k = 0;
    while k < raw.len() {
        let left = raw.len() - k - 1;
        let n = whole(raw[k], 1, usize::MAX, "a term list's count")?;
        if n > left {
            return Err(format!(
                "term list {} has count {n} but only {left} values follow",
                out.len()
            ));
        }
        let terms = &raw[k + 1..k + 1 + n];
        if !terms.iter().all(|m| m.is_finite()) {
            return Err(format!("term list {}: a term is not finite", out.len()));
        }
        if terms.iter().all(|m| *m == 0.0) {
            return Err(format!("term list {}: no nonzero term", out.len()));
        }
        out.push(terms.to_vec());
        k += 1 + n;
    }
    if out.len() != nlen {
        return Err(format!(
            "m_all holds {} term lists for {nlen} lengths",
            out.len()
        ));
    }
    Ok(out)
}

/// `stripmain` (or `stripmain_constrained`) at each length with its own longitudinal terms: the
/// modes, as `cufsm_modes` writes them, optionally classified and optionally restricted to cFSM
/// spaces. What the Python package's `strip` returns. See the module docs for the layouts.
///
/// # Safety
/// Every pointer must be valid for its length; `m_all` may be null when `m_all_len` is 0.
#[no_mangle]
pub unsafe extern "C" fn cufsm_strip(
    params: *const f64,
    params_len: usize,
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    bc: *const u8,
    bc_len: usize,
    lengths: *const f64,
    lengths_len: usize,
    m_all: *const f64,
    m_all_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        if params_len != 4 {
            return Err(format!(
                "params needs 4 values (terms, spaces, neigs, classify), got {params_len}"
            ));
        }
        // SAFETY (all the reads below): the export's contract, every pointer valid for its length.
        let p = unsafe { input(params, 4, "params") }?;
        let spaces = whole(p[1], 0, 15, "spaces")? as u32;
        let neigs = whole(p[2], 1, usize::MAX, "neigs")?;
        let classify = match p[3] {
            0.0 => false,
            1.0 => true,
            v => return Err(format!("classify must be 0 or 1, got {v}")),
        };
        let lens = unsafe { read_lengths(lengths, lengths_len) }?;
        let bc = unsafe { read_bc(bc, bc_len, true) }?;
        let model = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                springs,
                springs_len,
                constraints,
                constraints_len,
            )
        }?;
        let nn = model.nodes.len();
        let m_all = if m_all_len == 0 {
            let terms = whole(p[0], 1, usize::MAX, "terms")?;
            dense_fits(nn, terms)?;
            vec![(1..=terms).map(|k| k as f64).collect::<Vec<f64>>(); lens.len()]
        } else {
            if p[0] != 0.0 {
                return Err(format!(
                    "terms must be 0 when m_all gives the terms, got {}",
                    p[0]
                ));
            }
            read_m_all(unsafe { input(m_all, m_all_len, "m_all") }?, lens.len())?
        };
        // Each row's size from its terms as the analysis will use them (CUFSM's msort).
        let mut rows = Vec::with_capacity(lens.len());
        let mut need = 0usize;
        for m in &m_all {
            let nt = crate::analysis::msort(m).len();
            dense_fits(nn, nt)?;
            let per = modes_row(nn, nt, neigs)?;
            rows.push(per);
            need = add(need, per)?;
        }
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, need) }?;
        let results = if spaces == 0 {
            stripmain(&model, &lens, &m_all, bc, neigs)
        } else {
            let sp = Spaces {
                global: spaces & 1 != 0,
                distortional: spaces & 2 != 0,
                local: spaces & 4 != 0,
                other: spaces & 8 != 0,
            };
            stripmain_constrained(&model, &lens, &m_all, bc, neigs, sp)
        }
        .map_err(|e| e.to_string())?;
        let classes = if classify {
            Some(default_classes(&model, &results, bc)?)
        } else {
            None
        };
        let mut at = 0;
        for (i, (r, per)) in results.iter().zip(&rows).enumerate() {
            if *per != modes_row(nn, r.m_terms.len(), neigs)? {
                return Err("internal error: a row's terms changed during the analysis".into());
            }
            write_modes_row(
                &mut buf[at..at + per],
                r,
                classes.as_ref().map(|c| &c[i]),
                neigs,
            );
            at += per;
        }
        Ok(at)
    })
}

/// cFSM classification of given modes, CUFSM `classify.m` (uncoupled basis) with its orth, norm
/// and O space options: what the Python package's `classify` returns. See the module docs.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_classify(
    params: *const f64,
    params_len: usize,
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    bc: *const u8,
    bc_len: usize,
    results: *const f64,
    results_len: usize,
    springs: *const f64,
    springs_len: usize,
    constraints: *const f64,
    constraints_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    use cfsm::{Norm, OSpace, Orth};
    guarded(|| {
        if params_len != 3 {
            return Err(format!(
                "params needs 3 values (orth, norm, ospace), got {params_len}"
            ));
        }
        // SAFETY (all the reads below): the export's contract, every pointer valid for its length.
        let p = unsafe { input(params, 3, "params") }?;
        let orth = match whole(p[0], 1, 3, "orth")? {
            1 => Orth::Natural,
            2 => Orth::Axial,
            _ => Orth::Load,
        };
        let norm = match whole(p[1], 0, 3, "norm")? {
            0 => Norm::None,
            1 => Norm::Vector,
            2 => Norm::StrainEnergy,
            _ => Norm::Work,
        };
        let ospace = match whole(p[2], 1, 4, "ospace")? {
            1 => OSpace::St,
            2 => OSpace::K,
            3 => OSpace::Kg,
            _ => OSpace::Vector,
        };
        let bc = unsafe { read_bc(bc, bc_len, true) }?;
        let model = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                springs,
                springs_len,
                constraints,
                constraints_len,
            )
        }?;
        if results_len == 0 {
            return Err("no results given".into());
        }
        let raw = unsafe { input(results, results_len, "results") }?;
        let ndof = 4 * model.nodes.len();
        let mut lrs = vec![];
        let mut modes_total = 0usize;
        let mut k = 0;
        while k < raw.len() {
            let at = lrs.len();
            let left = raw.len() - k;
            if left < 2 {
                return Err(format!("results: length {at} is cut short"));
            }
            let length = raw[k];
            if !(length.is_finite() && length > 0.0) {
                return Err(format!(
                    "results: length {at} is not a positive finite number"
                ));
            }
            let nt = whole(raw[k + 1], 1, (left - 2).max(1), "a result's term count")?;
            if left < 3 + nt {
                return Err(format!("results: length {at} is cut short"));
            }
            let m_terms = raw[k + 2..k + 2 + nt].to_vec();
            if !m_terms.iter().all(|m| m.is_finite()) {
                return Err(format!(
                    "results: length {at} has a term that is not finite"
                ));
            }
            dense_fits(model.nodes.len(), nt)?;
            let per_mode = ndof * nt;
            let rest = left - 3 - nt;
            let nm = whole(raw[k + 2 + nt], 0, rest / per_mode, "a result's mode count")?;
            let start = k + 3 + nt;
            let modes = raw[start..start + nm * per_mode]
                .chunks_exact(per_mode)
                .map(<[f64]>::to_vec)
                .collect();
            modes_total = add(modes_total, nm)?;
            lrs.push(crate::LengthResult {
                length,
                m_terms,
                load_factors: vec![],
                modes,
            });
            k = start + nm * per_mode;
        }
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, mul(4, modes_total)?) }?;
        let classes =
            cfsm::classify_with(&model, &lrs, bc, orth, norm, ospace).map_err(|e| e.to_string())?;
        for (o, c) in buf.chunks_exact_mut(4).zip(classes.iter().flatten()) {
            o.copy_from_slice(c);
        }
        Ok(buf.len())
    })
}

/// CUFSM's C or Z section template, `templatecalc.m`: the meshed nodes and elements in this
/// interface's own input layouts. See the module docs.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_template(
    params: *const f64,
    params_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    use crate::template::{templatecalc, Shape, Template};
    guarded(|| {
        if params_len != 23 {
            return Err(format!("params needs 23 values, got {params_len}"));
        }
        // SAFETY: the export's contract: `params` valid for 23 values (checked just above).
        let p = unsafe { input(params, 23, "params") }?;
        let shape = match whole(p[0], 1, 2, "shape (1 = C, 2 = Z)")? {
            1 => Shape::C,
            _ => Shape::Z,
        };
        let names = [
            "h", "b1", "b2", "d1", "d2", "r1", "r2", "r3", "r4", "q1", "q2", "t",
        ];
        for (k, name) in names.iter().enumerate() {
            let v = p[1 + k];
            let ok = match *name {
                "h" | "b1" | "b2" | "t" => v.is_finite() && v > 0.0,
                "q1" | "q2" => v.is_finite(),
                _ => v.is_finite() && v >= 0.0,
            };
            if !ok {
                return Err(format!(
                    "{name} = {v} must be {}",
                    match *name {
                        "h" | "b1" | "b2" | "t" => "positive",
                        "q1" | "q2" => "finite",
                        _ => "zero or positive",
                    }
                ));
            }
        }
        let counts = ["nh", "nb1", "nb2", "nd1", "nd2", "nr1", "nr2", "nr3", "nr4"];
        let mut n = [0usize; 9];
        for (k, name) in counts.iter().enumerate() {
            let lo = usize::from(k < 3);
            n[k] = whole(p[13 + k], lo, usize::MAX, name)?;
        }
        let centerline = match p[22] {
            0.0 => false,
            1.0 => true,
            v => return Err(format!("centerline must be 0 or 1, got {v}")),
        };
        // At most sum(n) strips and one node more. Checked before meshing, so an absurd count is
        // refused here rather than allocated.
        let strips = n.iter().try_fold(0usize, |a, &b| add(a, b))?;
        let bound = add(add(2, mul(N_NODE, add(strips, 1)?)?)?, mul(N_ELEM, strips)?)?;
        if out.is_null() {
            return Err("no output buffer".into());
        }
        if out_cap < bound {
            return Err(format!(
                "output buffer is too small: {bound} values needed, {out_cap} given"
            ));
        }
        let tp = Template {
            shape,
            h: p[1],
            b1: p[2],
            b2: p[3],
            d1: p[4],
            d2: p[5],
            r1: p[6],
            r2: p[7],
            r3: p[8],
            r4: p[9],
            q1: p[10],
            q2: p[11],
            t: p[12],
            nh: n[0],
            nb1: n[1],
            nb2: n[2],
            nd1: n[3],
            nd2: n[4],
            nr1: n[5],
            nr2: n[6],
            nr3: n[7],
            nr4: n[8],
            centerline,
        };
        // The material is a placeholder, as in the Python package: the caller supplies its own.
        let m = templatecalc(&tp, Material::isotropic(1.0, 0.3));
        let need = 2 + N_NODE * m.nodes.len() + N_ELEM * m.elements.len();
        // SAFETY: the export's contract: `out` valid for `out_cap` values (and need <= bound).
        let buf = unsafe { output(out, out_cap, need) }?;
        buf[0] = m.nodes.len() as f64;
        buf[1] = m.elements.len() as f64;
        let f = |b: bool| if b { 1.0 } else { 0.0 };
        let mut k = 2;
        for nd in &m.nodes {
            buf[k..k + N_NODE].copy_from_slice(&[
                nd.x,
                nd.z,
                f(nd.free[0]),
                f(nd.free[1]),
                f(nd.free[2]),
                f(nd.free[3]),
                nd.stress,
            ]);
            k += N_NODE;
        }
        for e in &m.elements {
            buf[k..k + N_ELEM].copy_from_slice(&[e.ni as f64, e.nj as f64, e.t, e.mat as f64]);
            k += N_ELEM;
        }
        Ok(k)
    })
}

/// `cufsm_props`' 15 values followed by the warping function `wn` at every node (`cutwp_prop2`),
/// the section properties the Python package returns.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_props_wn(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        let need = 15 + m.nodes.len();
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, need) }?;
        let g = crate::grosprop(&m);
        let c = crate::cutwp_prop2(&m);
        buf[..15].copy_from_slice(&[
            g.a, g.xcg, g.zcg, g.ixx, g.izz, g.ixz, g.thetap, g.i11, g.i22, c.j, c.xs, c.zs, c.cw,
            c.b1, c.b2,
        ]);
        buf[15..].copy_from_slice(&c.wn);
        Ok(need)
    })
}

/// The 100 half-wavelengths of CUFSM's `signature_ss.m` for this model
/// ([`crate::signature_ss_lengths`]).
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn cufsm_signature_lengths(
    mats: *const f64,
    mats_len: usize,
    nodes: *const f64,
    nodes_len: usize,
    elems: *const f64,
    elems_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        // SAFETY: the export's own contract: every pointer valid for its length.
        let m = unsafe {
            read_model(
                mats,
                mats_len,
                nodes,
                nodes_len,
                elems,
                elems_len,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
            )
        }?;
        let ls = crate::signature_ss_lengths(&m);
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, ls.len()) }?;
        buf.copy_from_slice(&ls);
        Ok(ls.len())
    })
}

/// The interior local minima of a signature curve ([`crate::signature_minima`]): `curve` holds
/// `L, λ` pairs (a `NaN` λ, as `cufsm_signature` writes for a length with no positive load factor,
/// is a length without one); writes `L, λ` per minimum, shortest first.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_signature_minima(
    curve: *const f64,
    curve_len: usize,
    out: *mut f64,
    out_cap: usize,
) -> isize {
    guarded(|| {
        if curve_len == 0 || curve_len % 2 != 0 {
            return Err(format!(
                "curve length {curve_len} is not a positive multiple of 2"
            ));
        }
        // SAFETY: the export's contract: `curve` valid for `curve_len` values.
        let raw = unsafe { input(curve, curve_len, "curve") }?;
        let mut pts = Vec::with_capacity(curve_len / 2);
        for c in raw.chunks_exact(2) {
            if !(c[0].is_finite() && c[0] > 0.0) {
                return Err("a curve length is not a positive finite number".into());
            }
            if c[1].is_infinite() {
                return Err("a curve load factor is infinite".into());
            }
            pts.push(crate::LengthResult {
                length: c[0],
                m_terms: vec![],
                load_factors: if c[1].is_nan() { vec![] } else { vec![c[1]] },
                modes: vec![],
            });
        }
        let mins = crate::signature_minima(&pts);
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, 2 * mins.len()) }?;
        for (o, m) in buf.chunks_exact_mut(2).zip(&mins) {
            o[0] = m.length;
            o[1] = m.load_factor;
        }
        Ok(buf.len())
    })
}
