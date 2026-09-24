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
- Gross section properties, reference stresses from P, Mxx, Mzz, M11, M22, and the first-yield
  actions (Py, My) the Direct Strength Method divides by.
- The signature curve, and its local minima (the inputs to the Direct Strength Method).
- CUFSM's C and Z template, lipped or plain, sharp or with rounded corners, from centreline
  dimensions or from outside dimensions and inside radii.
- cFSM, the constrained finite strip method: the global, distortional, local and other modal
  spaces, analysis restricted to any of them (pure distortional buckling, for one), and the
  G/D/L/O classification of any mode, for open sections (single- or multi-branched) with the
  uncoupled basis, CUFSM's default. A restricted analysis takes fixed DOFs, constraints and
  springs as CUFSM does, intersecting the modal space with the constrained one. `cutwp_prop2` comes with it: shear centre, torsion and
  warping constants, and the warping function.

cFSM's uncoupled basis (`couple = 1`, CUFSM's default) takes every O-space choice (`ospace` 1
to 4); the coupled basis (`couple = 2`, for several longitudinal terms) is ported as CUFSM's
`base_update.m` has it, including that branch's own numbering of the O-space choices.

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
| **AISI Direct Strength Method Design Guide (2006)**: the CUFSM models behind its worked examples, with the results CUFSM saved then (`examples/2006_dsm_design_guide` in CUFSM) | 43 runs: lipped and plain channels, Z, hats, angles, a sigma, a rack upright, a built-up section, deck panels; compression and bending about either axis; some pinned or constrained. | 20,000+ load factors to 1e-5, the 2006 results' own precision, where cond(K) < 1e8. |
| **pyCUFSM**, an independent Python port | 5,900+ load factors on every fixture case it can run. | Within its own noise (below). |
| **Theory**, no oracle (`tests/theory.rs`) | A simply supported plate at k = 4 with its minimum at a square half-wave; an outstand at k = 0.425 + (b/a)²; a long I-section at the Euler load, converging with the mesh; every eigenpair satisfying K φ = λ Kg φ to round-off; invariance to E, stress scale, mirroring, renumbering and translation; convergence from above under mesh refinement; springs that only stiffen, and a stiff one that approaches a fixed DOF. | All hold. |
| **CUFSM's template** | Every node of 14 template cases. | To 1e-12. |
| **CUFSM's cFSM** (`analysis/cFSM/`, under Octave) | Six sections: sharp and rounded lipped C, lipped Z, plain channel, hat, branched I-section, and a lipped C with a fixity, a constraint and springs. `cutwp_prop2`'s properties and warping function; the sizes of the four spaces and the spaces themselves; load factors restricted to G, D or L alone; the classification of every distinct mode, with CUFSM's defaults, with the natural basis, and with each of the other O spaces. | Properties to 1e-10, spaces to 1e-8, restricted load factors to rounding, classifications to 1e-6 percentage points. |

### Accuracy, honestly

CUFSM solves `K φ = λ Kg φ` with MATLAB's `eigs`. This crate solves it one of two ways, and
neither can miss a mode:

- **The band path** (models of 60 DOF and more). The strips couple only neighbouring nodes,
  so after a reverse Cuthill-McKee reordering `K` and `Kg` are narrow bands. `K` is factored in
  band, and Lanczos with full re-orthogonalisation on `L⁻¹ Kg L⁻ᵀ` finds the wanted modes. The
  answer is then proved with a Sturm count: `K - σ Kg` has as many negative pivots as there are
  load factors below `σ` (Sylvester's law of inertia). If the count disagrees, Lanczos does not
  converge, or the band is not narrow, the dense path takes over.
- **The dense path.** Cholesky of `K`, every eigenvalue of `L⁻¹ Kg L⁻ᵀ` by Householder and implicit
  QL, and inverse iteration for the wanted modes.

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
- CUFSM also ships two `cutwp_prop2.m`: the `helpers/` one (first on the path) runs a
  multi-branched section through the open-section walk where the `analysis/` one does not.
  `helpers/`'s is ported. Its principal angle comes from `angle(Ix - Iy - 2 Ixy i)`, whose
  imaginary part MATLAB forms as `+0` when `Ixy` is zero, which decides between +π/2 and -π/2 for a
  section like a hat. The port matches.
- CUFSM's `mode_select` fails on an empty space (selecting the distortional modes of a plain
  channel, which has none). Here that is simply no load factors.
- On a doubly symmetric section, cFSM's axial orthogonalisation meets repeated eigenvalues, so
  its modal basis, and the vector-normalised classification with it, is not unique in CUFSM
  either. The natural basis is, and that is what such a section is compared on. The coupled G
  space of a clamped channel meets exact repeats too. How a mode splits among the spaces is
  still unique, so there the D : L : O proportions are compared, to 1e-8.
- The coupled branch of `base_update.m` numbers its O-space choices one higher than the
  uncoupled branch (3, 4, 5 for `K⁻¹`, `Kg⁻¹`, null space). With the natural basis and the ST
  O space it fills in no vectors, so that combination is refused here.
- Four of the DSM Design Guide files carry a saved curve that is not their saved model's:
  today's CUFSM, run on each file's own model, gives this crate's values, not the file's
  (`cwlip_modified.mat`: 1.62269 at 1.07 in, where the file says 83.33653). They are named in
  `tests/dsm_guide_parity.rs` and left out, as are the runs restricted by the 2006 cFSM, whose
  spaces were defined differently.
- A few of those models have strips far narrower than their neighbours (0.02 in beside 0.68 in),
  which leaves K too ill-conditioned for a plain Cholesky factorisation. This crate then retries
  with K scaled to a unit diagonal. It keeps the unscaled factorisation otherwise, because that
  one keeps more digits of long global modes.
- pyCUFSM (as released on PyPI) mishandles fixed DOFs and constraints. Its `constr_user` drops a
  column and leaves stale identity columns, so fixed DOFs at the high end of the numbering come
  back free: a plate simply supported on both long edges buckles as an outstand. It also cannot
  run more than one longitudinal term. Those cases are skipped, not compared.

## Cost

Each length is solved on its own thread (single-threaded on wasm). On an Apple M-series laptop
(`cargo run --release --example timing`), a 200 x 76 x 15 x 1.9 lipped channel's signature
curve, 100 half-wavelengths x 10 modes:

| Mesh | DOF | Signature curve | cFSM classification, per length |
|---|---|---|---|
| 29 nodes | 116 | 0.10 s | 1.2 ms |
| 45 nodes | 180 | 0.22 s | 2.8 ms |
| 77 nodes | 308 | 0.58 s | 14 ms |

For comparison, on one 148-DOF section with the same 100 lengths: this crate takes 0.14 s
multi-threaded (0.27 s on one thread). CUFSM's own assembly and a LAPACK dense solve take 3.5 s
under GNU Octave, and pyCUFSM takes 7.7 s. MATLAB itself was not available to time.

## Regenerating the references

`oracle/` has everything needed; nothing in it ships with the crate.

```sh
python3 oracle/cases.py > oracle/cases.json
CUFSM_ROOT=<cufsm-git checkout> oracle/run_octave.sh oracle/cases.json tests/fixtures/cufsm_octave.json
octave-cli oracle/octave/extract_dsm_guide.m <cufsm-git>/examples/2006_dsm_design_guide/files_and_scripts tests/fixtures/dsm_guide_2006.json
python3 oracle/run_pycufsm.py tests/fixtures/cufsm_octave.json tests/fixtures/pycufsm.json   # NumPy < 2
cargo run --example dump_matrices -- matlab AXIAL 300 > km.json && python3 oracle/high_precision.py km.json
```

Octave's own `eigs.m` (GPL) is copied into a temporary folder at run time for the `eigs` shim to
call; it is never kept in this repository.

## Licence

MIT. See [`LICENSE`](LICENSE), which carries CUFSM's notice as well as this port's. The DSM Design
Guide models and results in `tests/fixtures/dsm_guide_2006.json` are extracted from CUFSM's own
repository (MIT). The MATLAB reference values in `tests/fixtures/matlab_cufsm566.json` come from the
[CufsmSharp](https://github.com/BizimGri/CufsmSharp) project (MIT). Its `about` field
records their provenance.
