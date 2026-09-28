#!/bin/bash
# Reproduce every number in both papers, from this repository alone.
#
#   bash checks/reproduce.sh
#
# It runs the script named for each claim in checks/CLAIMS.tsv, then checks that every
# number the manuscripts quote beside that claim appears in that script's output at the
# precision the paper states. A claim whose anchor has moved out from under its script is
# reported too.
#
# This is the gate that matters to a reader. It is not the whole authoring pass, which
# also measures prose and needs a writing toolkit not redistributed here.
set -eu
cd "$(dirname "$0")/.."

command -v Rscript >/dev/null || { echo "Rscript not found" >&2; exit 2; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 2; }

# Say what is missing up front. Without this a reader with no numpy watched the pass
# report seventy-two numbers "not in the output", which is true and tells them nothing.
miss=$(python3 - <<'PYX'
import importlib.util
print(" ".join(m for m in ("numpy", "scipy", "sympy", "mpmath")
                if importlib.util.find_spec(m) is None))
PYX
)
if [ -n "$miss" ]; then
  echo "The Python scripts need these and this interpreter has none of them: $miss" >&2
  echo "  python3 -m pip install -r checks/requirements.txt" >&2
  echo "The R scripts need nothing and would run; the pass checks both, so it stops here." >&2
  exit 2
fi

export CLAIMS_TSV=checks/CLAIMS.tsv
export CLAIMS_PAPERS_DIR=paper
export CLAIMS_SCRIPT_MAP='checks/calc/=checks/calc/,checks/=checks/,tangents/=tangents/'
export PYTHONPATH=checks${PYTHONPATH:+:$PYTHONPATH}

python3 checks/claims_check.py paper/PAPER2_v4_draft.md
