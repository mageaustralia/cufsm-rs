#!/usr/bin/env python3
"""Runs every case of the CUFSM (Octave) fixture through pyCUFSM, as a second, independent
implementation: python3 run_pycufsm.py ../tests/fixtures/cufsm_octave.json ../tests/fixtures/pycufsm.json

Needs pycufsm (https://github.com/ClearCalcs/pyCUFSM) and NumPy < 2, which its compiled core
requires. The inputs - nodes with CUFSM's stresses, elements, lengths, terms - are exactly the
ones CUFSM ran, so the three implementations are compared on identical models.
"""
import json
import sys

import numpy as np
from pycufsm.fsm import strip


def vec(v):
    return v if isinstance(v, list) else [v]


def rows(v):
    return v if v and isinstance(v[0], list) else [v]


src, dst = sys.argv[1], sys.argv[2]
out = []
for r in json.load(open(src)):
    E, nu = r["E"], r["nu"]
    props = np.array([[100, E, E, nu, nu, E / (2 * (1 + nu))]])
    node = np.array(rows(r["node"]), dtype=float)
    node[:, 0] -= 1  # pyCUFSM numbers nodes from 0
    elem = np.array([[e[0] - 1, e[1] - 1, e[2] - 1, e[3], 100] for e in rows(r["elem"])], dtype=float)
    cons = r["constraints"]
    if isinstance(cons, list) and cons:
        # pyCUFSM numbers nodes and DOF codes from 0 (CUFSM: 1=x 2=z 3=y 4=rotation)
        cons = np.array([[c[0] - 1, c[1] - 1, c[2], c[3] - 1, c[4] - 1] for c in rows(cons)], dtype=float)
    else:
        cons = np.array([])
    lengths = np.array(vec(r["lengths"]), dtype=float)
    m_all = np.array([vec(m) for m in r["m_all"]], dtype=float)  # every length here has the same terms
    gbt = {"glob": [0], "dist": [0], "local": [0], "other": [0], "o_space": 1, "couple": 1, "orth": 1, "norm": 0}
    # One length per call: pyCUFSM fails to pad its curve when two lengths return different
    # numbers of positive modes ("could not broadcast ... shape (9,) into shape (10,)").
    lfs = []
    error = None
    try:
        for i in range(len(lengths)):
            _, curve, _ = strip(
                props=props, nodes=node, elements=elem, lengths=lengths[i:i + 1], springs=np.array([]),
                constraints=cons, GBT_con=gbt, B_C=r["bc"], m_all=m_all[i:i + 1], n_eigs=int(r["neigs"]),
                sect_props={},
            )
            lfs.append([float(x) for x in np.asarray(curve[0]).ravel() if x > 0])
    except Exception as e:  # recorded, not fatal: the test says which cases pyCUFSM could not run
        error = f"{type(e).__name__}: {e}"
        lfs = []
    fixed = bool((node[:, 3:7] == 0).any()) or len(cons) > 0
    out.append({"name": r["name"], "load_factors": lfs, "fixed_or_constrained": fixed, "error": error})
    print(f"{r['name']}: {error or f'{len(lfs)} lengths, lowest {min(min(l) for l in lfs if l):.6g}'}")
json.dump(out, open(dst, "w"))
