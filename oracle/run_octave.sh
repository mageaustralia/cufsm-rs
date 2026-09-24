#!/usr/bin/env bash
# Runs the CUFSM reference cases through CUFSM's own MATLAB code under GNU Octave.
#
#   CUFSM_ROOT=<cufsm-git checkout> OCTAVE_ENV=<conda env with octave> ./run_octave.sh cases.json out.json
#
# Octave from conda-forge must run inside its activated environment (it segfaults without it),
# so this goes through micromamba when MAMBA is set, and plain octave-cli otherwise.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
: "${CUFSM_ROOT:?set CUFSM_ROOT to a checkout of https://github.com/thinwalled/cufsm-git}"
if [ -n "${MAMBA:-}" ]; then
  exec "$MAMBA" run -p "$OCTAVE_ENV" octave-cli --no-gui --norc "$here/octave/run_cases.m" "$CUFSM_ROOT" "$1" "$2"
else
  exec octave-cli --no-gui --norc "$here/octave/run_cases.m" "$CUFSM_ROOT" "$1" "$2"
fi
