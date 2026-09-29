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
  first (bottom) end also holds the axial `v`.

The tube axis is global Z. Nodes sit on the mid-surface radius.
"""

import math
from typing import Tuple

# The cases, kept in step with examples/ftm_oracle_case.rs.
CASES = {
    "long_N": dict(r=100.0, t=2.0, l=8000.0, load="N"),
    "med_N": dict(r=100.0, t=2.0, l=400.0, load="N"),
    "med_M": dict(r=100.0, t=2.0, l=400.0, load="M"),
    "thin_N": dict(r=500.0, t=1.0, l=2000.0, load="N"),
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
    nodes, elems = mesh(r, l, ntheta, nz)
    # Reference load: unit total compression, or unit peak bending stress.
    if load == "N":
        per = -1.0 / ntheta                                    # total N = 1 (compression)
        top_force = lambda th: per
    else:                                                       # peak stress = 1 MPa at θ = 0
        s0 = 1.0
        area = 2.0 * math.pi * r * t / ntheta
        top_force = lambda th: -s0 * math.cos(th) * area
    out = [f"*HEADING\n{name}: FTM oracle, {load} on a r={r} t={t} l={l} tube"]
    out.append("*NODE")
    for i, x, y, z in nodes:
        out.append(f"{i}, {x:.6f}, {y:.6f}, {z:.6f}")
    out.append("*ELEMENT, TYPE=S4R, ELSET=Shell")
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
    out.append("Bot, 3, 3")       # first end holds v (pinned-pinned, Euler k = 1)
    out.append("Top, 1, 2")
    out.append("*STEP")
    out.append("*BUCKLE\n10,")
    out.append("*CLOAD")
    for m in (nz,):
        for k in range(ntheta):
            th = 2.0 * math.pi * k / ntheta
            out.append(f"{m * ntheta + k + 1}, 3, {top_force(th):.6e}")
    out.append("*NODE FILE, NSET=Shell\nU")
    out.append("*END STEP")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "med_N"
    c = CASES[name]
    sys.stdout.write(input_deck(name, **c))
