#!/usr/bin/env bash
# Reproduce every number and figure in the paper from scratch.
#
#   ./reproduce.sh            full run (~45 min, the Monte Carlo dominates)
#   ./reproduce.sh --quick    tests and the worked examples only (~30 s)
set -euo pipefail
cd "$(dirname "$0")"

echo "== invariants =="
python3 tests/test_qcell.py

echo
echo "== Supplementary Note 4: each equation worked through numerically =="
(cd paper && python3 worked_examples.py)

if [[ "${1:-}" == "--quick" ]]; then
    echo; echo "quick run complete (scenarios and figures skipped)"
    exit 0
fi

echo
echo "== scenarios (this is the long one) =="
python3 run_cases.py

echo
echo "== figures =="
python3 make_figs.py

echo
echo "== every number in the paper, checked against those results =="
python3 paper/verify_paper.py

echo
echo "done:"
echo "  results/results.json"
echo "  figures/fig1_framework.png .. fig5_design.png"
