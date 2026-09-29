# FTM ↔ CalculiX oracle

An independent check on `cufsm-rs`'s finite tube method (`src/ftm.rs`): each case
is run twice: once through the Fourier-discretised tube, once through a
continuum shell model in [CalculiX](http://www.calculix.de/), and the lowest
buckling factors are compared. The two discretisations share nothing but the
geometry and the load.

## Layout

| File | What it is |
|---|---|
| `tube_inp.py` | The case list and the CalculiX deck generator (S4R cylinder, `*BUCKLE`). |
| `run_oracle.py` | Orchestrator: generate -> ccx -> parse -> `cargo run --example ftm_oracle_case` -> table. |
| `../examples/ftm_oracle_case.rs` | The same cases through the FTM, printing `case lambda` lines. |

Run (ccx from the local micromamba env, CalculiX 2.23):

```sh
export MAMBA=<tools>/bin/micromamba MAMBA_ROOT_PREFIX=<tools>/mamba
CCX_CMD="$MAMBA run -p <tools>/ccx-env ccx" python3 oracle/ftm_calculix/run_oracle.py
MESH_SCALE=2 ... run_oracle.py thin_N        # convergence knob, and per-case filter
```

## First results (2026-09-29)

| Case | ccx λ | FTM λ | FTM/ccx | |
|---|---|---|---|---|
| `long_N` - Euler column | 153.7 | 156.25 | **1.016** | ok |
| `thin_N` - R/t = 500 | 238.7 | 245.71 | **1.029** | ok, converged from below (218 -> 239 with the mesh) |
| `med_M` - stub, bending | 2582.7 | 2515.05 | **0.974** | ok |
| `med_N` - stub, compression | 2000.0 | 2401.64 | **1.201** | open: converged on both sides (ccx 3×->6× meshes agree; FTM flat over p = 8...32, q = 6...16) - the difference lives in the method's kinematics, not truncation |

Reading: the long and thin cases validate the pipeline (the `long_N` spectrum shows the
Euler doublet, and the shell value climbs to classical as the wave resolves). The stub
compression gap (`L/R = 4`, about 20%) is the one substantive difference - a likely
interactive/diamond mode that the paper's Sanders-type strain set (hoop terms left out)
cannot express. Bending of the same stub agrees to 3%.

## Deck notes learned the hard way

* `*NODE FILE` must not name an element set.
* The anti-slide of `End::Pinned` is `sum(v) = 0` around the end ring (an `*EQUATION`);
  pinning `v` at every end node clamps the bending slope (gave exactly 2× Euler), and a
  single node dims up a local mode.
* Loads are unit *stress* (force per node = its share of `2πRt`), so λ compares as MPa
  with the FTM convention.
* S4, not S4R: the reduced-integration shell is too stiff for these thin walls.

## Conventions

* **Loads are nodal `*CLOAD` forces.** An internal plate suite of ours
  (an internal plate suite of ours) found `*DLOAD ... P` on
  expanded S8R shells doubles the load; this suite never uses it. That suite's
  runner and FRD parser are the same graceful-skip / compare-the-in-house-solver
  pattern; reuse them there if displacement-field comparison is ever wanted
  (`parse_frd.py` reads the `*NODE FILE` output).
* **Ends mirror `ftm::End::Pinned`**: `u = w = 0` all round through cylindrical
  `*TRANSFORM`s (dofs 1 = radial, 2 = tangential), with the bottom end also
  holding axial `v`: Euler k = 1.
* **Cases** (`tube_inp.CASES`, kept in step with the example): a long tube that
  must land on Euler, a medium tube in the local/global transition, the same tube
  in bending, and an R/t = 500 thin tube (the convergence case).
* **Nothing is gated.** No reference numbers are committed until the oracle has
  been run on a machine with ccx; the tolerance in `run_oracle.py` is ±5% and can
  be tightened once the spread is known.

## Getting ccx

No homebrew formula exists. The working route is conda-forge's native Apple Silicon
build (CalculiX 2.23), which is what the recipe above uses; a container works just the
same via `CCX_CMD`:

```sh
CCX_CMD="docker run --rm -v $PWD/oracle/ftm_calculix/work:/work -w /work <image> ccx" \
  python3 oracle/ftm_calculix/run_oracle.py
```

Without ccx the suite prints the FTM factors and skips - it is still useful as a
regression snapshot.
