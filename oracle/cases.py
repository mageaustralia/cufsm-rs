#!/usr/bin/env python3
"""Writes the reference case list the oracles run: python3 cases.py > cases.json

Standard library only. Units N and mm throughout; E = 203 000 MPa, nu = 0.3.

Sections come either from CUFSM's own template (`templatecalc.m`, run inside the Octave oracle,
so the geometry is CUFSM's) or as explicit node and element lists for shapes the template does
not make (a hat, an angle, a flat plate).
"""
import json
import math

E, NU = 203000.0, 0.3


def logspace(a, b, n):
    return [10 ** (math.log10(a) + (math.log10(b) - math.log10(a)) * i / (n - 1)) for i in range(n)]


def template(corz, h, b1, b2, d1, d2, r, t, nh=4, nb=2, nd=1, nr=0):
    return {
        "CorZ": corz, "h": h, "b1": b1, "b2": b2, "d1": d1, "d2": d2,
        "r1": r, "r2": r, "r3": r, "r4": r, "q1": 90, "q2": 90, "t": t,
        "nh": nh, "nb1": nb, "nb2": nb, "nd1": nd, "nd2": nd,
        "nr1": nr, "nr2": nr, "nr3": nr, "nr4": nr, "kipin": 1, "center": 1,
    }


def actions(P=0.0, Mxx=0.0, Mzz=0.0, M11=0.0, M22=0.0, unsymm=0):
    return {"P": P, "Mxx": Mxx, "Mzz": Mzz, "M11": M11, "M22": M22, "unsymm": unsymm}


def polyline(points, t, subdiv):
    """Nodes and elements for an open polyline, each segment cut into `subdiv` strips."""
    pts = [points[0]]
    for (x0, z0), (x1, z1) in zip(points, points[1:]):
        for k in range(1, subdiv + 1):
            pts.append((x0 + (x1 - x0) * k / subdiv, z0 + (z1 - z0) * k / subdiv))
    node = [[i + 1, x, z, 1, 1, 1, 1, 0] for i, (x, z) in enumerate(pts)]
    elem = [[i + 1, i + 1, i + 2, t, 100] for i in range(len(pts) - 1)]
    return node, elem


def signature(name, geom, act, n_lengths=40, neigs=10, lmin=10.0, lmax=10000.0, **extra):
    lengths = logspace(lmin, lmax, n_lengths)
    c = {"name": name, "E": E, "nu": NU, "actions": act, "bc": "S-S",
         "lengths": lengths, "m_all": [[1]] * len(lengths), "neigs": neigs}
    c.update(geom)
    c.update(extra)
    return c


def general(name, geom, act, bc, lengths, nterms, neigs=5, **extra):
    c = {"name": name, "E": E, "nu": NU, "actions": act, "bc": bc,
         "lengths": lengths, "m_all": [list(range(1, nterms + 1))] * len(lengths), "neigs": neigs}
    c.update(geom)
    c.update(extra)
    return c


cases = []

# Lipped channel 200 x 75 x 20 x 1.5: square corners, then rounded (r = 3, two strips per corner).
lc = {"template": template(1, 200, 75, 75, 20, 20, 0, 1.5)}
lc_r = {"template": template(1, 200, 75, 75, 20, 20, 3.0, 1.5, nh=6, nb=3, nd=2, nr=2)}
cases.append(signature("lipped-c compression", lc, actions(P=1000.0)))
cases.append(signature("lipped-c major bending", lc, actions(Mxx=1e6)))
cases.append(signature("lipped-c minor bending", lc, actions(Mzz=1e6)))
cases.append(signature("lipped-c rounded compression", lc_r, actions(P=1000.0)))
cases.append(signature("lipped-c rounded major bending", lc_r, actions(Mxx=1e6)))

# Unequal flanges and lips, and a plain (unlipped) channel.
cases.append(signature("lipped-c unequal compression", {"template": template(1, 250, 80, 60, 25, 15, 0, 2.0)}, actions(P=1000.0)))
cases.append(signature("plain channel compression", {"template": template(1, 150, 50, 50, 0, 0, 0, 2.0)}, actions(P=1000.0)))

# Lipped Z: restrained bending about the geometric axis (unsymm 0), unrestrained (unsymm 1),
# and about the principal axis.
lz = {"template": template(2, 200, 70, 70, 20, 20, 0, 1.8)}
cases.append(signature("lipped-z compression", lz, actions(P=1000.0)))
cases.append(signature("lipped-z restrained bending", lz, actions(Mxx=1e6, unsymm=0)))
cases.append(signature("lipped-z unrestrained bending", lz, actions(Mxx=1e6, unsymm=1)))
cases.append(signature("lipped-z principal bending", lz, actions(M11=1e6, unsymm=1)))

# Shapes the template does not make.
hat_node, hat_elem = polyline([(0, 0), (20, 0), (20, 60), (80, 60), (80, 0), (100, 0)], 1.2, 3)
cases.append(signature("hat compression", {"node": hat_node, "elem": hat_elem}, actions(P=1000.0)))
cases.append(signature("hat bending", {"node": hat_node, "elem": hat_elem}, actions(Mxx=1e6)))
ang_node, ang_elem = polyline([(60, 0), (0, 0), (0, 60)], 3.0, 6)
cases.append(signature("equal angle compression", {"node": ang_node, "elem": ang_elem}, actions(P=1000.0)))

# A flat plate with its long edges' out-of-plane displacement fixed: the textbook k = 4 plate.
pl_node, pl_elem = polyline([(0, 0), (100, 0)], 2.0, 10)
cases.append(signature("ss plate compression", {"node": pl_node, "elem": pl_elem}, actions(P=1000.0),
                       lmin=20, lmax=1000, fix=[[1, 5], [len(pl_node), 5]]))

# General end conditions with several longitudinal terms, on a small lipped C so the global
# matrices are carried too (4 x 11 nodes x terms <= 400 DOF).
small = {"template": template(1, 150, 60, 60, 15, 15, 0, 1.5, nh=4, nb=2, nd=1)}
for bc, terms in [("C-C", 8), ("S-C", 8), ("C-F", 8), ("C-G", 8), ("S-S", 6)]:
    cases.append(general(f"lipped-c {bc} compression", small, actions(P=1000.0), bc, [300.0, 1200.0, 4000.0], terms))
cases.append(general("lipped-c C-C major bending", small, actions(Mxx=1e6), "C-C", [600.0, 2500.0], 8))

# Master-slave constraints: the two lip tips' z displacement tied together, and the web
# mid-node's rotation fixed.
cases.append(signature("lipped-c constrained compression", small, actions(P=1000.0),
                       n_lengths=20, constraints=[[1, 2, 1.0, 11, 2]], fix=[[6, 7]]))

# Outside dimensions and inside radii, converted by CUFSM's template_out_to_in.m: a 200 x 76 x 15
# x 1.9 lipped C and a Z, with 3 mm inside radii, and a sharp-cornered plain channel.
out_c = template(1, 200, 76, 76, 15, 15, 3.0, 1.9, nh=8, nb=4, nd=2, nr=2)
out_c["center"] = 0
out_z = template(2, 200, 76, 70, 15, 15, 3.0, 1.9, nh=8, nb=4, nd=2, nr=2)
out_z["center"] = 0
out_p = template(1, 150, 50, 50, 0, 0, 0, 2.0)
out_p["center"] = 0
cases.append(signature("outside-dims lipped-c compression", {"template": out_c}, actions(P=1000.0), n_lengths=20))
cases.append(signature("outside-dims lipped-z compression", {"template": out_z}, actions(P=1000.0), n_lengths=20))
cases.append(signature("outside-dims plain channel compression", {"template": out_p}, actions(P=1000.0), n_lengths=20))

# Springs (CUFSM v4.3 form: [# nodei nodej ku kv kw kq local discrete ys]; nodej 0 = ground).
# A foundation spring holding a lip tip, a discrete spring at mid-length restraining the other lip,
# a global spring and a local one between the two lip tips, on the small lipped C, under S-S and
# under C-C with several terms.
springs_a = [[1, 1, 0, 0, 0, 5.0, 0, 0, 0, 0],
             [2, 11, 0, 0, 0, 50.0, 0, 0, 1, 0.5]]
springs_b = [[1, 1, 11, 2.0, 0, 3.0, 0, 1, 0, 0],
             [2, 1, 11, 0, 0, 1.0, 500.0, 0, 0, 0]]
cases.append(signature("lipped-c grounded springs", small, actions(P=1000.0), n_lengths=20, springs=springs_a))
cases.append(signature("lipped-c node springs", small, actions(P=1000.0), n_lengths=20, springs=springs_b))
cases.append(general("lipped-c C-C springs", small, actions(P=1000.0), "C-C", [600.0, 2500.0], 6, springs=springs_a + springs_b[:1]))

# cFSM: compression on sections the modal decomposition covers - a sharp lipped C, a rounded one
# (many small corners), a lipped Z, a plain channel, a hat, and a branched I-section (a node with
# three strips, which exercises the warping constraints at a branch).
def i_section(h, bf, tw, tf, nw, nf):
    node = [[1 + k, 0.0, h * k / nw, 1, 1, 1, 1, 0] for k in range(nw + 1)]
    elem = [[k + 1, k + 1, k + 2, tw, 100] for k in range(nw)]
    for root, z in [(nw + 1, h), (1, 0.0)]:
        for side in (-1.0, 1.0):
            prev = root
            for k in range(1, nf + 1):
                node.append([len(node) + 1, side * bf / 2 * k / nf, z, 1, 1, 1, 1, 0])
                elem.append([len(elem) + 1, prev, len(node), tf, 100])
                prev = len(node)
    return {"node": node, "elem": elem}

cfsm_lengths = [60.0, 250.0, 900.0, 3000.0]
def cf(name, geom):
    c = signature(name, geom, actions(P=1000.0), n_lengths=2, neigs=5)
    c["lengths"] = cfsm_lengths
    c["m_all"] = [[1]] * len(cfsm_lengths)
    c["cfsm"] = True
    return c
cases.append(cf("cfsm lipped-c", lc))
cases.append(cf("cfsm lipped-c rounded", lc_r))
cases.append(cf("cfsm lipped-z", lz))
cases.append(cf("cfsm plain channel", {"template": template(1, 150, 50, 50, 0, 0, 0, 2.0)}))
cases.append(cf("cfsm hat", {"node": hat_node, "elem": hat_elem}))
cases.append(cf("cfsm i-section", i_section(200.0, 100.0, 5.0, 8.0, 4, 2)))
# cFSM with a fixed DOF, a master-slave constraint and springs: CUFSM intersects the modal space
# with the constrained one, R = null([null(Rmode') null(Ruser')]').
c_fix = cf("cfsm lipped-c fixed and sprung", small)
c_fix["fix"] = [[6, 7]]
c_fix["constraints"] = [[1, 2, 1.0, 11, 2]]
c_fix["springs"] = springs_a
cases.append(c_fix)
# cFSM with several longitudinal terms (clamped ends), including the coupled basis (couple 2).
c_cc = cf("cfsm lipped-c C-C coupled", small)
c_cc["bc"] = "C-C"
c_cc["lengths"] = [600.0, 2500.0]
c_cc["m_all"] = [[1, 2, 3]] * 2
c_cc["coupled"] = True
cases.append(c_cc)

print(json.dumps(cases, indent=1))
