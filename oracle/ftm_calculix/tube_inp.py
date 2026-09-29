"""Generate a CalculiX input deck for the buckling of one circular tube.

This is the independent oracle for `cufsm-rs`'s finite tube method: a continuum
shell model (S4R), so nothing is shared with the Fourier-discretised tube except
the geometry and the load. The two conventions here are deliberate:

* Loads are applied as nodal ``*CLOAD`` forces. An internal plate suite of ours
  (an internal plate suite of ours) found that ``*DLOAD ...
  P`` on expanded S8R shells doubles the applied pressure; buckling with nodal
  forces avoids that path entirely.
* End conditions mirror `ftm::End::Pinned`: `u = w = 0` all round (through a
  cylindrical `*TRANSFORM`, so radial/tangential are exact per node), and the
  bottom ring's axial displacements sum to zero (an `*EQUATION`): the
  anti-slide of `End::Pinned` without restraining the end warping - pinning `v`
  all round would clamp the bending slope, a single node would dimple.

The tube axis is global Z. Nodes sit on the mid-surface radius.
"""

import math
import os
from typing import Tuple

# The cases, kept in step with examples/ftm_oracle_case.rs.
# Meshes: near-square elements of about 2 sqrt(R t), so a shell wave is resolved
# and the ends are not dominated by load-introduction dimples.
CASES = {
    "long_N": dict(r=100.0, t=2.0, l=8000.0, load="N", ntheta=32, nz=160),
    "med_N": dict(r=100.0, t=2.0, l=400.0, load="N", ntheta=32, nz=16),
    "med_M": dict(r=100.0, t=2.0, l=400.0, load="M", ntheta=32, nz=16),
    "thin_N": dict(r=500.0, t=1.0, l=2000.0, load="N", ntheta=160, nz=100),
}
E = 203_000.0
NU = 0.3


def mesh(r: float, l: float, ntheta: int, nz: int) -> Tuple[list, list]:
    """Nodes (id, x, y, z) and S4R quads (id, n1..n4) on the cylinder."""
    nodes = []
    for m in range(nz + 1):
        for k in range(ntheta):
            th = 2.0 * math.pi * k / ntheta
            nodes.append((m * ntheta + k + 1, r * math.cos(th), r * math.sin(th), l * m / nz))
    elems = []
    e = 0
    for m in range(nz):
        for k in range(ntheta):
            kn = (k + 1) % ntheta
            a = m * ntheta + k + 1
            b = m * ntheta + kn + 1
            elems.append((e + 1, a, b, b + ntheta, a + ntheta))
            e += 1
    return nodes, elems


def input_deck(name: str, r: float, t: float, l: float, load: str,
               ntheta: int = 64, nz: int = 40) -> str:
    scale = float(os.environ.get("MESH_SCALE", "1"))   # convergence knob
    ntheta, nz = int(ntheta * scale), int(nz * scale)
    nodes, elems = mesh(r, l, ntheta, nz)
    # Reference load: unit total compression, or unit peak bending stress.
    if load == "N":
        area = 2.0 * math.pi * r * t / ntheta
        top_force = lambda th: -area                        # unit stress: λ compares as MPa
    else:                                                       # peak stress = 1 MPa at θ = 0
        s0 = 1.0
        area = 2.0 * math.pi * r * t / ntheta
        top_force = lambda th: -s0 * math.cos(th) * area
    out = [f"*HEADING\n{name}: FTM oracle, {load} on a r={r} t={t} l={l} tube"]
    out.append("*NODE")
    for i, x, y, z in nodes:
        out.append(f"{i}, {x:.6f}, {y:.6f}, {z:.6f}")
    out.append("*ELEMENT, TYPE=S4, ELSET=Shell")
    for i, a, b, c, d in elems:
        out.append(f"{i}, {a}, {b}, {c}, {d}")
    out.append(f"*SHELL SECTION, ELSET=Shell, MATERIAL=STEEL\n{t}")
    out.append(f"*MATERIAL, NAME=STEEL\n*ELASTIC\n{E}, {NU}")
    # End rings: cylindrical transforms so dof 1 = radial, 2 = tangential, 3 = axial.
    out.append("*NSET, NSET=Bot, GENERATE\n1, %d, %d" % (ntheta, 1))
    out.append("*NSET, NSET=Top, GENERATE\n%d, %d, %d"
               % (nz * ntheta + 1, (nz + 1) * ntheta, 1))
    out.append("*TRANSFORM, NSET=Bot, TYPE=C\n0, 0, 0, 0, 0, 1")
    out.append("*TRANSFORM, NSET=Top, TYPE=C\n0, 0, 0, 0, 0, 1")
    out.append("*BOUNDARY")
    out.append("Bot, 1, 2")       # u = w = 0 all round
    out.append("Top, 1, 2")
    # Anti-slide as the tube's Pinned end means it: the bottom ring's axial
    # displacements sum to zero - no warping restraint, no single-point dimple.
    out.append("*EQUATION")
    out.append(str(ntheta))
    ring = [k + 1 for k in range(ntheta)]
    for i in range(0, ntheta, 4):
        terms = ring[i:i + 4]
        out.append(", ".join(f"{n}, 3, 1.0" for n in terms))
    out.append("*STEP")
    out.append("*BUCKLE\n10,")
    out.append("*CLOAD")
    for m in (nz,):
        for k in range(ntheta):
            th = 2.0 * math.pi * k / ntheta
            out.append(f"{m * ntheta + k + 1}, 3, {top_force(th):.6e}")
    out.append("*NODE FILE\nU")
    out.append("*END STEP")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "med_N"
    c = CASES[name]
    sys.stdout.write(input_deck(name, **c))
