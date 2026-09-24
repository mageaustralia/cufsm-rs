# cufsm-rs

An independent, dependency-free Rust port of [CUFSM](https://www.ce.jhu.edu/cufsm/), the
finite strip method for the elastic buckling of thin-walled sections, by Benjamin W. Schafer
and co-workers at Johns Hopkins University.

**Not affiliated with or endorsed by the CUFSM authors.** CUFSM is MIT-licensed; its copyright
notice is kept in [`LICENSE`](LICENSE). If you use this in published work, cite CUFSM itself, as its
[`Citation.cff`](https://github.com/thinwalled/cufsm-git/blob/main/Citation.cff) asks: CUFSM
v5.66, Schafer, B.W., Ádány, S., Li, Z., Jin, S. (DOI 10.5281/zenodo.17771486), and for the
general end conditions, Schafer, B.W., Li, Z., "Buckling analysis of cold-formed steel members
with general boundary conditions using CUFSM: conventional and constrained finite strip methods",
20th International Specialty Conference on Cold-Formed Steel Structures, 2010, pp. 17-32.

## What it does

- Strip elastic and geometric stiffness, for all five end conditions CUFSM offers
  (S-S, C-C, S-C, C-F, C-G) and any set of longitudinal terms.
- Fixed nodal DOFs, master-slave constraints, and springs (foundation or discrete, to ground or
  between nodes, CUFSM's v4.3 form).
- Gross section properties and reference stresses from P, Mxx, Mzz, M11, M22.
- The signature curve, and its local minima (the inputs to the Direct Strength Method).
- CUFSM's C and Z template, lipped or plain, sharp or with rounded corners, from centreline
  dimensions or from outside dimensions and inside radii.

Not yet ported: cFSM (the constrained finite strip modal decomposition that labels
modes as local, distortional or global).

```rust
use cufsm::{grosprop, stresgen, signature_ss, signature_minima, Actions, Material};
use cufsm::template::{templatecalc, Shape, Template};

// A 200 x 76 x 15 x 1.9 lipped channel, 3 mm inside radii, in compression.
let mut m = templatecalc(&Template::outside(Shape::C, 200.0, 76.0, 15.0, 1.9, 3.0, 12),
                         Material::isotropic(203_000.0, 0.3));
let props = grosprop(&m);
stresgen(&mut m, &Actions { p: 1.0, ..Default::default() }, &props, false);
let curve = signature_ss(&m, 1).unwrap();
for min in signature_minima(&curve) {
    println!("half-wavelength {:.0} mm: Pcr = {:.1} kN", min.length, min.load_factor / 1e3);
}
```

## How it is checked

The point of a port is that it gives CUFSM's answers. Every claim below is a test in `tests/`.

| Reference | What is compared | Result |
|---|---|---|
| **CUFSM itself**, its MATLAB source run unmodified under GNU Octave (`oracle/`) | 28 cases: lipped, rounded, unequal and plain channels, Z in four loadings, hat, angle, plate, outside-dimension templates; all five end conditions with up to 8 terms; fixities, constraints and springs. Stage by stage: section properties, stresses, each strip's local and global matrices, the assembled `K` and `Kg`, then the load factors and first modes. | Matrices, properties and stresses to 1e-12. 7,000+ load factors and 450+ mode shapes (MAC to 1 - 1e-8). |
| **MATLAB CUFSM** v5.66, compiled, under the MATLAB R2025b Runtime (from the CufsmSharp project) | 1,480 load factors, a lipped channel with 2.5 mm corner radii, compression and bending. | Local modes to 1e-14 - 1e-12. |
| **pyCUFSM**, an independent Python port | 5,900+ load factors on every fixture case it can run. | Within its own noise (below). |
| **Theory**, no oracle (`tests/theory.rs`) | A simply supported plate at k = 4 with its minimum at a square half-wave; an outstand at k = 0.425 + (b/a)²; a long I-section at the Euler load, converging with the mesh; every eigenpair satisfying K φ = λ Kg φ to round-off; invariance to E, stress scale, mirroring, renumbering and translation; convergence from above under mesh refinement; springs that only stiffen, and a stiff one that approaches a fixed DOF. | All hold. |
| **CUFSM's template** | Every node of 14 template cases. | To 1e-12. |

### Accuracy, honestly

CUFSM solves `K φ = λ Kg φ` with MATLAB's `eigs`. This crate solves it densely instead
(Cholesky of `K`, then every eigenpair of `L⁻¹ Kg L⁻ᵀ` by Householder and implicit QL), which cannot
miss a mode and has no convergence tolerance.

For local and distortional modes the four implementations agree to 12+ digits. Checked against
the same matrices solved in 40-digit arithmetic, on a well-conditioned mode this crate was out by
1e-16, CUFSM's dense `eig` by 2e-15 and pyCUFSM by 3.6e-9.

A **global** mode at a long half-wavelength is a different matter. It barely strains the
membrane, so its load factor rests on a near-cancellation between very large membrane terms,
and any double-precision code loses about `eps × cond(K)` of it. Section models with very narrow
strips (small corner radii) are the worst case. On the MATLAB reference's longest global mode,
against 40-digit arithmetic, MATLAB CUFSM was out by -3.1e-5 and this crate by +5.9e-5. That is
the rounding limit of the method, not a porting error, and it is irrelevant at engineering
precision. The parity tests allow `1e-10 + 1e-12 × cond(K)`, scaled by `λᵢ / λ₁` for higher modes,
and no more.

### What the comparison turned up

- CUFSM ships two `stresgen.m`, in `analysis/` and `helpers/`, with opposite signs on the M11
  term. CUFSM's interface puts `helpers/` last on the path, so that version is the one ported.
- Octave's `eigs` cannot handle an indefinite `Kg` the way MATLAB's does. On bending cases it
  returns eigenvalues wrong by orders of magnitude, and different on every run. The oracle
  therefore also solves CUFSM's own reduced matrices with `eig()`.
- pyCUFSM (as released on PyPI) mishandles fixed DOFs and constraints. Its `constr_user` drops a
  column and leaves stale identity columns, so fixed DOFs at the high end of the numbering come
  back free: a plate simply supported on both long edges buckles as an outstand. It also cannot
  run more than one longitudinal term. Those cases are skipped, not compared.

## Regenerating the references

`oracle/` has everything needed; nothing in it ships with the crate.

```sh
python3 oracle/cases.py > oracle/cases.json
CUFSM_ROOT=<cufsm-git checkout> oracle/run_octave.sh oracle/cases.json tests/fixtures/cufsm_octave.json
python3 oracle/run_pycufsm.py tests/fixtures/cufsm_octave.json tests/fixtures/pycufsm.json   # NumPy < 2
cargo run --example dump_matrices -- matlab AXIAL 300 > km.json && python3 oracle/high_precision.py km.json
```

Octave's own `eigs.m` (GPL) is copied into a temporary folder at run time for the `eigs` shim to
call; it is never kept in this repository.

## Licence

MIT. See [`LICENSE`](LICENSE), which carries CUFSM's notice as well as this port's. The MATLAB
reference values in `tests/fixtures/matlab_cufsm566.json` come from the
[CufsmSharp](https://github.com/BizimGri/CufsmSharp) project (MIT). Its `about` field
records their provenance.
