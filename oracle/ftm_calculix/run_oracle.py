"""Compare the finite tube method (cufsm-rs) with a CalculiX shell model.

Independent oracle for `cufsm-rs`'s FTM: for each case in ``tube_inp.CASES`` this
generates a continuum-shell CalculiX input (S4R cylinder, nodal *CLOAD), extracts
the first buckling factor from ccx's ``.dat``, runs the same case through
``cargo run --example ftm_oracle_case``, and prints a comparison table.

The CalculiX side is skipped gracefully when ccx is not installed. Install with
conda (`conda install -c conda-forge calculix`) or point CCX_CMD at anything that
behaves like ``ccx -i <job>`` run in the deck's directory — for example

    CCX_CMD="docker run --rm -v $PWD:/work -w /work <image> ccx" python3 run_oracle.py

Conventions reused from an internal plate suite of ours
(`an internal plate suite of ours`): graceful skip, and
comparison against the in-house solver as the thing being validated. Its
``*DLOAD P`` shell caveat does not apply here — this suite loads with *CLOAD.
"""

import math
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tube_inp import CASES, E, NU, input_deck  # noqa: E402

ORACLE_DIR = os.path.dirname(os.path.abspath(__file__))
CRATE = os.path.dirname(os.path.dirname(ORACLE_DIR))


def parse_first_factor(dat_path: str):
    """The lowest buckling factor from a ccx .dat eigenvalue table.

    Tolerant of CalculiX versions: look for an 'E I G E N V A L U E' / 'EIGENVALUE'
    marker, then the first numeric row of `mode factor`."""
    if not os.path.exists(dat_path):
        return None
    text = open(dat_path, errors="replace").read()
    rows = []
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+([+\-0-9.Ee]+)\s*$", line)
        if m:
            rows.append(float(m.group(2)))
    if rows:
        return min(rows)                      # the lowest factor is the one wanted
    # eigenvalue tables may use packed formats; scan for the marker's numbers
    for m in re.finditer(r"E\s*I\s*G\s*E\s*N\s*V\s*A\s*L\s*U\s*E.{0,200}?([0-9.]+(?:[Ee][+\-]\d+)?)",
                         text, re.S):
        try:
            return float(m.group(1))
        except ValueError:
            continue
    return None


def run_ccx(job_dir: str, job: str):
    """Run ccx in job_dir; returns the .dat path or None when ccx is unavailable."""
    cmd = os.environ.get("CCX_CMD", "ccx")
    exe = cmd.split()[0]
    if exe == "ccx" and shutil.which("ccx") is None:
        print("  [SKIP] CalculiX not found; install with `conda install -c conda-forge calculix`")
        print("         or set CCX_CMD to a ccx equivalent (see the module docstring)")
        return None
    full = cmd.split() + ["-i", job]
    proc = subprocess.run(full, cwd=job_dir, capture_output=True, text=True)
    dat = os.path.join(job_dir, job + ".dat")
    if not os.path.exists(dat):
        print(f"  [ERROR] ccx produced no .dat for {job}: {proc.stdout[-300:]}{proc.stderr[-300:]}")
        return None
    return dat


def ftm_factors():
    """`{case: lambda}` from the crate's own engine."""
    proc = subprocess.run(
        ["cargo", "run", "--quiet", "--example", "ftm_oracle_case"],
        cwd=CRATE, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        print(proc.stderr[-600:])
        sys.exit(1)
    out = {}
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[0]] = float(parts[1])
    return out


def closed_form(name, r, t, l):
    """Long-tube Euler and classical shell stress, for orientation (not gating)."""
    if name.endswith("_N") and l > 100 * r:
        i = math.pi * r**3 * t
        a = 2 * math.pi * r * t
        return math.pi**2 * E * i / (a * l * l)
    return None


def main():
    os.makedirs(os.path.join(ORACLE_DIR, "work"), exist_ok=True)
    ftm = ftm_factors()
    print(f"{'case':8} {'ccx λ':>12} {'FTM λ':>12} {'ratio':>8}  note")
    ok = True
    for name, c in CASES.items():
        deck = input_deck(name, **c)
        job_dir = os.path.join(ORACLE_DIR, "work")
        path = os.path.join(job_dir, name + ".inp")
        open(path, "w").write(deck)
        dat = run_ccx(job_dir, name)
        lam_ftm = ftm.get(name)
        if dat is None:
            print(f"{name:8} {'-':>12} {lam_ftm:12.4f} {'-':>8}  (ccx skipped)")
            continue
        lam_ccx = parse_first_factor(dat)
        cf = closed_form(name, c["r"], c["t"], c["l"])
        note = f"closed form {cf:.1f}" if cf else ""
        if lam_ccx is None:
            print(f"{name:8} {'?':>12} {lam_ftm:12.4f} {'-':>8}  no eigenvalue in .dat {note}")
            ok = False
            continue
        ratio = lam_ftm / lam_ccx
        verdict = "ok" if 0.95 <= ratio <= 1.05 else "CHECK"
        if verdict != "ok":
            ok = False
        print(f"{name:8} {lam_ccx:12.4f} {lam_ftm:12.4f} {ratio:8.3f}  {verdict} {note}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
