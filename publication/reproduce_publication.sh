#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PRISM_PYTHON:-$ROOT/.venv/bin/python}"

if [[ "$PYTHON" != */* ]]; then
  PYTHON="$(command -v "$PYTHON" || true)"
fi
if [[ ! -x "$PYTHON" ]]; then
  echo "Python environment not found at $PYTHON. Set PRISM_PYTHON to a Python with pandas, numpy, scipy, matplotlib, and statsmodels." >&2
  exit 2
fi

cd "$ROOT"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"
MPLCONFIGDIR="$(mktemp -d "${TMPDIR:-/tmp}/prism-matplotlib.XXXXXX")"
export MPLCONFIGDIR
trap 'rm -rf -- "$MPLCONFIGDIR"' EXIT
"$PYTHON" publication/manuscript_results_audit.py
"$PYTHON" publication/scripts/generate_feature_definitions.py
"$PYTHON" publication/scripts/generate_supplementary.py
"$PYTHON" publication/scripts/generate_figures.py
"$PYTHON" publication/scripts/feature_correlation.py
"$PYTHON" publication/scripts/generate_candidate_evidence.py
"$PYTHON" publication/covariate_sensitivity/run_covariate_sensitivity.py
"$PYTHON" publication/scripts/generate_provenance.py
"$PYTHON" publication/publication_lint.py
