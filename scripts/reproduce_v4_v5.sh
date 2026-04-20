#!/usr/bin/env bash
set -euo pipefail

python -m pip install -e ".[structural,dev]"

pytest -q tests/test_selectivity_v4_data.py tests/test_selectivity_v4_features.py
python -m kira.experiments.selectivity_v4_data
python -m kira.experiments.run_selectivity_v4

pytest -q tests/test_selectivity_v5_expand_data.py tests/test_selectivity_v5_exact_core.py
python -m kira.experiments.selectivity_v5_exact_core \
  --candidate-csv data/processed/selectivity_v5_candidate_rows.csv \
  --output-csv data/processed/selectivity_v5_exact_core_rows.csv \
  --output-summary-json data/processed/selectivity_v5_exact_core_summary.json
