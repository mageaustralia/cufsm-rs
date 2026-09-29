# FTM ↔ CalculiX oracle

An independent check on `cufsm-rs`'s finite tube method (`src/ftm.rs`): each case
is run twice — once through the Fourier-discretised tube, once through a
continuum shell model in [CalculiX](http://www.calculix.de/) — and the lowest
buckling factors are compared. The two discretisations share nothing but the
geometry and the load.

## Layout

| File | What it is |
|---|---|
| `tube_inp.py` | The case list and the CalculiX deck generator (S4R cylinder, `*BUCKLE`). |
| `run_oracle.py` | Orchestrator: generate → ccx → parse → `cargo run --example ftm_oracle_case` → table. |
| `../examples/ftm_oracle_case.rs` | The same cases through the FTM, printing `case lambda` lines. |

Run:

```sh
python3 oracle/ftm_calculix/run_oracle.py
```

## Conventions

* **Loads are nodal `*CLOAD` forces.** An internal plate suite of ours
  (an internal plate suite of ours) found `*DLOAD ... P` on
  expanded S8R shells doubles the load; this suite never uses it. That suite's
  runner and FRD parser are the same graceful-skip / compare-the-in-house-solver
  pattern; reuse them there if displacement-field comparison is ever wanted
  (`parse_frd.py` reads the `*NODE FILE` output).
* **Ends mirror `ftm::End::Pinned`**: `u = w = 0` all round through cylindrical
  `*TRANSFORM`s (dofs 1 = radial, 2 = tangential), with the bottom end also
  holding axial `v` — Euler k = 1.
* **Cases** (`tube_inp.CASES`, kept in step with the example): a long tube that
  must land on Euler, a medium tube in the local/global transition, the same tube
  in bending, and an R/t = 500 thin tube (the convergence case).
* **Nothing is gated.** No reference numbers are committed until the oracle has
  been run on a machine with ccx; the tolerance in `run_oracle.py` is ±5% and can
  be tightened once the spread is known.

## Getting ccx

No homebrew formula exists (an old `brew install calculix-ccx` instruction
never worked). Either `conda install -c conda-forge calculix`, or point `CCX_CMD`
at a container:

```sh
CCX_CMD="docker run --rm -v $PWD/oracle/ftm_calculix/work:/work -w /work <image> ccx" \
  python3 oracle/ftm_calculix/run_oracle.py
```

Without ccx the suite prints the FTM factors and skips — it is still useful as a
regression snapshot.
