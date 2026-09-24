#!/usr/bin/env python3
"""The lowest load factors of K phi = lambda Kg phi in high-precision arithmetic (mpmath).

    cargo run --example dump_matrices -- matlab AXIAL 300 > km.json
    python3 oracle/high_precision.py km.json [digits]

The matrices are the double-precision ones the analysis assembles, so this is the exact answer
for them: it separates what the eigen-solver loses from what the matrices already carry.
"""
import json
import sys

import mpmath as mp

mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
d = json.load(open(sys.argv[1]))
K = mp.matrix(d["K"])
Kg = mp.matrix(d["Kg"])
L = mp.cholesky(K)
Li = mp.inverse(L)
C = Li * Kg * Li.T
C = (C + C.T) / 2
E = mp.eigsy(C, eigvals_only=True)
mus = sorted([e for e in E if e > 0], reverse=True)[:3]
for mu in mus:
    print(mp.nstr(1 / mu, 20))
