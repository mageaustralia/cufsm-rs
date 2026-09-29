# cufsm-rs

A Rust port of [CUFSM](https://www.ce.jhu.edu/cufsm/), the finite strip method for the elastic
buckling of thin-walled sections, by Benjamin W. Schafer and co-workers at Johns Hopkins
University. It has no dependencies.

This is an independent port, not affiliated with or endorsed by the CUFSM authors. CUFSM is
MIT-licensed, and its copyright notice is kept in [`LICENSE`](LICENSE).

If you use this in published work, cite CUFSM as its
[`Citation.cff`](https://github.com/thinwalled/cufsm-git/blob/main/Citation.cff) asks:

- Schafer, B.W., Ádány, S., Li, Z., Jin, S. CUFSM v5.66. DOI 10.5281/zenodo.17771486.
- For general end conditions: Schafer, B.W., Li, Z. "Buckling analysis of cold-formed steel
  members with general boundary conditions using CUFSM: conventional and constrained finite
  strip methods." 20th International Specialty Conference on Cold-Formed Steel Structures,
  2010, pp. 17-32.

## Features

- Elastic and geometric strip stiffness for all five of CUFSM's end conditions (S-S, C-C, S-C,
  C-F, C-G), with any number of longitudinal terms.
- Fixed DOFs, master-slave constraints, and springs (foundation or discrete, to ground or
  between nodes), in CUFSM's v4.3 form.
- Gross section properties; reference stresses from P, Mxx, Mzz, M11 and M22; and the
  first-yield actions Py and My that the Direct Strength Method uses.
- The signature curve and its local minima, which are the inputs to the Direct Strength Method.
- CUFSM's C and Z templates: lipped or plain, sharp or rounded corners, from centreline
  dimensions or from outside dimensions and inside radii.
- cFSM, the constrained finite strip method:
  - the global, distortional, local and other (G, D, L, O) modal spaces;
  - analysis restricted to any of those spaces, for example pure distortional buckling;
  - classification of any mode into G, D, L and O.

  It covers open sections, single- or multi-branched. Restricted analyses honour fixed DOFs,
  constraints and springs the way CUFSM does. The uncoupled basis (CUFSM's default) supports
  all four O-space options. The coupled basis, used with several longitudinal terms, follows
  CUFSM's `base_update.m`, including that branch's different numbering of the O-space options.
- `cutwp_prop2`: shear centre, torsion and warping constants, and the warping function.
- An optional interface for calling the crate from JavaScript as a WebAssembly module (see
  [Calling it from a web page](#calling-it-from-a-web-page)).

## Example

```rust
use cufsm::{grosprop, stresgen, signature_ss, signature_minima, Actions, Material};
use cufsm::template::{templatecalc, Shape, Template};

// A 200 x 76 x 15 x 1.9 lipped channel with 3 mm inside radii, in compression.
let mut m = templatecalc(
    &Template::outside(Shape::C, 200.0, 76.0, 15.0, 1.9, 3.0, 12),
    Material::isotropic(203_000.0, 0.3),
);
let props = grosprop(&m);
stresgen(&mut m, &Actions { p: 1.0, ..Default::default() }, &props, false);

let curve = signature_ss(&m, 1).unwrap();
for min in signature_minima(&curve) {
    println!("half-wavelength {:.0} mm: Pcr = {:.1} kN", min.length, min.load_factor / 1e3);
}
```

The crate is named `cufsm-rs` on crates.io and imported as `cufsm`.

## Calling it from a web page

The `ffi` feature, which is off by default, adds a small C interface. With it, the crate builds
as a WebAssembly module that JavaScript can call directly, without a binding library:

```sh
cargo rustc --release --lib --features ffi --crate-type cdylib --target wasm32-unknown-unknown
```

The module has two functions:

- `cufsm_signature` returns the buckling load factor at each length. It can also return the
  load factors restricted to the G, D, L and O spaces.
- `cufsm_modes` returns the lowest mode at each length, with its G, D, L and O percentages.

Inputs and outputs are arrays of numbers in the module's memory. `src/ffi.rs` describes their
layout. For S-S, the lengths are half-wavelengths and there must be one longitudinal term,
because that is what defines the signature curve. For the other end conditions, the lengths are
member lengths and any number of terms can be used.

A crate that depends on cufsm-rs does not get these functions unless it turns the feature on.

## Validation

Each comparison below is a test in the [GitHub repository](https://github.com/mageaustralia/cufsm-rs).
The tests and the reference data are not included in the crates.io package.

| Reference | Cases | Agreement |
|---|---|---|
| CUFSM's MATLAB source, run unmodified in GNU Octave | 28 cases covering channels, Z sections, a hat, an angle, a plate and outside-dimension templates. All five end conditions, up to 8 terms, with fixities, constraints and springs. Every stage is compared: section properties, stresses, strip matrices, assembled `K` and `Kg`, load factors and mode shapes. | Matrices, properties and stresses to 1e-12. Over 7,000 load factors, and over 450 mode shapes with MAC above 1 − 1e-8. |
| Compiled MATLAB CUFSM v5.66, run under the MATLAB R2025b Runtime (values from the CufsmSharp project) | 1,480 load factors for a lipped channel with 2.5 mm corner radii, in compression and bending. | Local modes to 1e-12 or better. |
| The models behind the AISI *Direct Strength Method Design Guide* (2006), with the results CUFSM saved at the time | 43 runs: channels, Z sections, hats, angles, a sigma, a rack upright, a built-up section and deck panels, in compression and in bending about either axis. | Over 20,000 load factors to 1e-5, the precision of the saved results, where cond(K) < 1e8. |
| pyCUFSM, an independent Python port | Over 5,900 load factors, on every case it can run. | Within pyCUFSM's own precision (see below). |
| Theory (`tests/theory.rs`) | A plate at k = 4 with its minimum at a square half-wave. An outstand at k = 0.425 + (b/a)². A long I-section converging to the Euler load. Every eigenpair satisfying K φ = λ Kg φ. Invariance to E, stress scale, mirroring, renumbering and translation. Convergence from above as the mesh is refined. Springs that only stiffen, and a stiff spring approaching a fixed DOF. | All hold. |
| CUFSM's template generator | Every node of 14 template cases. | To 1e-12. |
| CUFSM's cFSM code, run in Octave | Six sections: sharp and rounded lipped C, lipped Z, plain channel, hat, branched I-section, and a lipped C with a fixity, a constraint and springs. Compared: `cutwp_prop2` properties and warping function, the four modal spaces, load factors restricted to G, D or L, and the classification of every distinct mode under each basis and O-space option. | Properties to 1e-10, spaces to 1e-8, restricted load factors to rounding, classifications to 1e-6 percentage points. |

### How eigenvalues are solved

CUFSM solves `K φ = λ Kg φ` with MATLAB's `eigs`. This crate uses one of two methods, neither of
which can skip a mode:

- **Banded (60 DOF and above).** Strips only couple neighbouring nodes, so after reverse
  Cuthill-McKee reordering `K` and `Kg` are narrow bands. `K` is factored in band form, and
  Lanczos with full reorthogonalisation finds the wanted modes of `L⁻¹ Kg L⁻ᵀ`. A Sturm count
  then confirms the result: `K − σ Kg` has one negative pivot per load factor below `σ`. If the
  count disagrees, Lanczos fails to converge, or the band is too wide, the dense method is used
  instead.
- **Dense.** Cholesky factorisation of `K`, all eigenvalues of `L⁻¹ Kg L⁻ᵀ` by Householder
  reduction and implicit QL, then inverse iteration for the wanted modes.

### Precision

For local and distortional modes, this crate, CUFSM (MATLAB and Octave) and pyCUFSM agree to
at least 12 significant digits. Against the same matrices solved in 40-digit arithmetic, on a
well-conditioned mode, the errors were 1e-16 for this crate, 2e-15 for CUFSM's dense `eig`, and
3.6e-9 for pyCUFSM.

Global modes at long half-wavelengths are less precise in every implementation. Their load
factor is a small difference between large membrane terms, so double precision loses about
`eps × cond(K)` of it. Sections with very narrow strips, such as small corner radii, are the
worst case. On the MATLAB reference's longest global mode, compared with 40-digit arithmetic,
MATLAB CUFSM was off by −3.1e-5 and this crate by +5.9e-5. That is the limit of the method in
double precision, and it is well below engineering precision. The parity tests allow
`1e-10 + 1e-12 × cond(K)`, scaled by `λᵢ / λ₁` for higher modes.

### Differences found in the references

- CUFSM ships two `stresgen.m` files, in `analysis/` and `helpers/`, with opposite signs on the
  M11 term. CUFSM's interface puts `helpers/` last on the path, so that is the version ported.
- Octave's `eigs` mishandles an indefinite `Kg`. On bending cases it returns eigenvalues that
  are wrong by orders of magnitude and change from run to run. The Octave reference therefore
  also solves CUFSM's reduced matrices with `eig()`.
- CUFSM also ships two `cutwp_prop2.m` files. The `helpers/` one, first on the path, handles
  multi-branched sections through the open-section walk and the `analysis/` one does not, so
  the `helpers/` one is ported. Its principal angle comes from `angle(Ix − Iy − 2 Ixy i)`.
  MATLAB forms the imaginary part as +0 when `Ixy` is zero, which decides between +π/2 and −π/2
  for a section such as a hat. The port does the same.
- CUFSM's `mode_select` fails on an empty space, for example the distortional modes of a plain
  channel, which has none. This crate returns no load factors instead.
- On a doubly symmetric section, cFSM meets repeated eigenvalues, so its modal basis, and the
  vector-normalised classification that depends on it, is not unique in CUFSM either. Those
  sections are compared with the natural basis, which is unique. The coupled G space of a
  clamped channel also has exact repeats. There, the D : L : O proportions of each mode are
  unique and are compared to 1e-8.
- The coupled branch of `base_update.m` numbers its O-space options one higher than the
  uncoupled branch (3, 4 and 5 for `K⁻¹`, `Kg⁻¹` and the null space). With the natural basis
  and the ST O space it produces no vectors, so this crate refuses that combination.
- Four DSM Design Guide files store a curve that does not match the model saved with them.
  Current CUFSM, run on each file's model, gives this crate's values, not the stored ones. For
  example, `cwlip_modified.mat` stores 83.33653 at 1.07 in, where both give 1.62269. These files
  are listed in `tests/dsm_guide_parity.rs` and excluded, as are runs that used the 2006 version
  of cFSM, which defined the spaces differently.
- Some of those models have strips much narrower than their neighbours (0.02 in next to
  0.68 in), which makes `K` too ill-conditioned for a plain Cholesky factorisation. This crate
  then retries with `K` scaled to a unit diagonal. It does not scale otherwise, because the
  unscaled factorisation keeps more digits in long global modes.
- pyCUFSM, as released on PyPI, mishandles fixed DOFs and constraints. Its `constr_user` drops a
  column and leaves stale identity columns, so fixed DOFs at the high end of the numbering come
  back free: a plate simply supported on both long edges buckles as an outstand. It also cannot
  run more than one longitudinal term. Those cases are skipped rather than compared.

## Performance

Each half-wavelength is solved on its own thread (single-threaded on wasm). Timings on an Apple
M-series laptop (`cargo run --release --example timing`), for the 200 x 76 x 15 x 1.9 lipped
channel above, 100 half-wavelengths and 10 modes each:

| Mesh | DOF | Signature curve | cFSM classification, per length |
|---|---|---|---|
| 29 nodes | 116 | 0.10 s | 1.2 ms |
| 45 nodes | 180 | 0.22 s | 2.8 ms |
| 77 nodes | 308 | 0.58 s | 14 ms |

For a 148-DOF section with the same 100 lengths, this crate took 0.14 s (0.27 s on one thread).
CUFSM's assembly with a LAPACK dense solve took 3.5 s under GNU Octave, and pyCUFSM took 7.7 s.
MATLAB was not available for timing.

## Regenerating the reference data

The `oracle/` directory in the GitHub repository contains the scripts. None of it is included in
the crate.

```sh
python3 oracle/cases.py > oracle/cases.json
CUFSM_ROOT=<cufsm-git checkout> oracle/run_octave.sh oracle/cases.json tests/fixtures/cufsm_octave.json
octave-cli oracle/octave/extract_dsm_guide.m <cufsm-git>/examples/2006_dsm_design_guide/files_and_scripts tests/fixtures/dsm_guide_2006.json
python3 oracle/run_pycufsm.py tests/fixtures/cufsm_octave.json tests/fixtures/pycufsm.json   # needs NumPy < 2
cargo run --example dump_matrices -- matlab AXIAL 300 > km.json && python3 oracle/high_precision.py km.json
```

Octave's `eigs.m` is GPL-licensed. It is copied into a temporary folder when the scripts run and
is never stored in this repository.

## Licence

MIT. See [`LICENSE`](LICENSE), which includes CUFSM's copyright notice as well as this port's.

The DSM Design Guide models and results in `tests/fixtures/dsm_guide_2006.json` are extracted
from CUFSM's own repository (MIT). The MATLAB reference values in
`tests/fixtures/matlab_cufsm566.json` come from the
[CufsmSharp](https://github.com/BizimGri/CufsmSharp) project (MIT); the file's `about` field
records where they came from.
