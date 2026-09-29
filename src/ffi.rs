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
//! * `params`: `[t, E, nu, terms, spaces]`
//! * `nodes`: 5 per node: `x, z, stress, free_x, free_z` (the last two 0 or 1)
//! * `elems`: 2 per element: `node_i, node_j`; every element takes `t` from `params`
//! * `springs`: 9 per spring (0 for none): `node_i, node_j` (-1 for ground), `k_u, k_v, k_w, k_q`,
//!   `local` (0 or 1), `discrete` (0 or 1), `ys_fraction` - CUFSM's `springs` row, v4.3 form
//! * `constraints`: 5 per constraint (0 for none): `node_e, dof_e, coeff, node_k, dof_k`, with
//!   dof codes 1 = x, 2 = z, 3 = y along the member, 4 = theta, and the row reading
//!   `u_e = coeff * u_k`
//! * `bc`: UTF-8 bytes of CUFSM's boundary condition string, `S-S`, `C-C`, `S-C`, `C-F`, `C-G`
//! * `lengths`: for `S-S`, half-wavelengths, with `terms` = 1 (the signature curve takes the
//!   single term m = 1; more terms at a half-wavelength would give the curve's minimum over
//!   L, L/2, … instead, so it is refused). For the other boundary conditions, physical member
//!   lengths, with the longitudinal terms `1..=terms` (1 to [`MAX_TERMS`]).
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
//! * `cufsm_modes` output: per length: `G, D, L, O` (the lowest mode's class percentages),
//!   `nterms`, the `nterms` longitudinal terms, then the lowest mode's `4 * nodes * nterms`
//!   entries in CUFSM's order (per term: `u`/`v` interleaved, then `w`/`θ` interleaved).
//!
//! Both return the number of `f64` written, or a negative number on failure, when
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

const N_PARAMS: usize = 5;
const N_NODE: usize = 5;
const N_ELEM: usize = 2;
const N_SPRING: usize = 9;
const N_CONSTRAINT: usize = 5;
const NEIGS: usize = 1;
/// The most longitudinal terms accepted for a boundary condition other than S-S.
pub const MAX_TERMS: usize = 100;

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

/// Everything an export needs, read and checked from the flat buffers. Every length is checked
/// first: a short buffer is a caller bug, and reading past it would be undefined behaviour.
struct Input {
    model: Model,
    bc: BoundaryCondition,
    lengths: Vec<f64>,
    m_all: Vec<Vec<f64>>,
    spaces: u32,
}

#[allow(clippy::too_many_arguments)]
unsafe fn read_input(
    params: *const f64,
    params_len: usize,
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
    if lengths_len == 0 {
        return Err("no lengths given".into());
    }
    // SAFETY (this and the three below): the caller guarantees each pointer is valid for its length.
    let p = unsafe { input(params, N_PARAMS, "params") }?;
    let raw_nodes = unsafe { input(nodes, nodes_len, "nodes") }?;
    let raw_elems = unsafe { input(elems, elems_len, "elems") }?;
    let lens = unsafe { input(lengths, lengths_len, "lengths") }?;

    let (t, e, nu) = (p[0], p[1], p[2]);
    if !(t.is_finite() && t > 0.0) {
        return Err("thickness must be a positive finite number".into());
    }
    if !(e.is_finite() && e > 0.0) {
        return Err("Young's modulus must be a positive finite number".into());
    }
    if !(nu.is_finite() && nu > -1.0 && nu < 0.5) {
        return Err(format!(
            "Poisson's ratio must be between -1 and 0.5, got {nu}"
        ));
    }
    let spaces = whole(p[4], 0, 15, "spaces")? as u32;
    if lens.iter().any(|l| !(l.is_finite() && *l > 0.0)) {
        return Err("a length is not a positive finite number".into());
    }

    if bc.is_null() || bc_len == 0 {
        return Err("no boundary condition given".into());
    }
    // SAFETY: not null and non-empty (checked above); the caller guarantees `bc_len` bytes at `bc`.
    let s = std::str::from_utf8(unsafe { slice::from_raw_parts(bc, bc_len) })
        .map_err(|_| "boundary condition is not UTF-8".to_string())?;
    let bc =
        BoundaryCondition::parse(s).ok_or_else(|| format!("unknown boundary condition {s:?}"))?;
    let terms = if bc == BoundaryCondition::SS {
        if p[3] != 1.0 {
            return Err(format!(
                "S-S is the signature curve, which takes the single term m = 1 at each half-wavelength; terms must be 1, got {}",
                p[3]
            ));
        }
        1
    } else {
        whole(p[3], 1, MAX_TERMS, "terms")?
    };

    let mut out_nodes = Vec::with_capacity(nodes_len / N_NODE);
    for c in raw_nodes.chunks_exact(N_NODE) {
        if !c[..3].iter().all(|v| v.is_finite()) {
            return Err("a node coordinate or stress is not finite".into());
        }
        let flag = |v: f64| {
            if v == 0.0 || v == 1.0 {
                Ok(v == 1.0)
            } else {
                Err(format!("a node's free flag must be 0 or 1, got {v}"))
            }
        };
        let mut n = Node::new(c[0], c[1], c[2]);
        n.free = [flag(c[3])?, flag(c[4])?, true, true];
        out_nodes.push(n);
    }
    let n_nodes = out_nodes.len();
    let mut out_elems = Vec::with_capacity(elems_len / N_ELEM);
    for c in raw_elems.chunks_exact(N_ELEM) {
        let ni = whole(c[0], 0, n_nodes - 1, "an element's node_i")?;
        let nj = whole(c[1], 0, n_nodes - 1, "an element's node_j")?;
        out_elems.push(Element { ni, nj, t, mat: 0 });
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
        materials: vec![Material::isotropic(e, nu)],
        nodes: out_nodes,
        elements: out_elems,
        constraints: read_constraints(constraints_raw, n_nodes)?,
        springs: read_springs(springs_raw, n_nodes)?,
    };
    model.validate().map_err(|err| err.to_string())?;
    let m: Vec<f64> = (1..=terms).map(|k| k as f64).collect();
    Ok(Input {
        model,
        bc,
        m_all: vec![m; lens.len()],
        lengths: lens.to_vec(),
        spaces,
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
/// asked for in `params[4]` from `stripmain_constrained`. See the module docs for the layouts.
///
/// # Safety
/// Every pointer must be valid for its length.
#[no_mangle]
pub unsafe extern "C" fn cufsm_signature(
    params: *const f64,
    params_len: usize,
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
            stripmain(&i.model, &i.lengths, &i.m_all, i.bc, NEIGS).map_err(|e| e.to_string())?;
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
            let run = stripmain_constrained(&i.model, &i.lengths, &i.m_all, i.bc, NEIGS, sp)
                .map_err(|e| e.to_string())?;
            for (row, r) in buf.chunks_exact_mut(stride).zip(&run) {
                row[col] = first_lf(r);
            }
            col += 1;
        }
        Ok(buf.len())
    })
}

/// The lowest mode at each requested length, with its G/D/L/O classification percentages.
///
/// # Safety
/// Every pointer must be valid for its length; `out` needs, per length,
/// `5 + nterms + 4 * nodes * nterms` values.
#[no_mangle]
pub unsafe extern "C" fn cufsm_modes(
    params: *const f64,
    params_len: usize,
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
        let per = 5 + nt + 4 * nn * nt;
        // SAFETY: the export's contract: `out` valid for `out_cap` values.
        let buf = unsafe { output(out, out_cap, i.lengths.len() * per) }?;
        let results =
            stripmain(&i.model, &i.lengths, &i.m_all, i.bc, NEIGS).map_err(|e| e.to_string())?;
        let classes = cfsm::classify(
            &i.model,
            &results,
            i.bc,
            cfsm::Orth::Axial,
            cfsm::Norm::Vector,
        )
        .map_err(|e| e.to_string())?;
        for ((row, r), modes) in buf.chunks_exact_mut(per).zip(&results).zip(&classes) {
            let cls = modes.first().copied().unwrap_or([f64::NAN; 4]);
            row[..4].copy_from_slice(&cls);
            row[4] = nt as f64;
            row[5..5 + nt].copy_from_slice(&r.m_terms);
            // Lowest mode only: the GUI draws and classifies mode 1.
            match r.modes.first() {
                Some(mode) => row[5 + nt..].copy_from_slice(mode),
                None => row[5 + nt..].fill(f64::NAN),
            }
        }
        Ok(buf.len())
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
