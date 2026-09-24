#!/bin/bash
# Thin wrapper. Harbor's Oracle copies the whole solution/ directory to
# /solution, so the reference script is resolved as a sibling of this file
# rather than by an absolute path.
set -uo pipefail
python3 "$(dirname "$0")/reference_solution.py"
