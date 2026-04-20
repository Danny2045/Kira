# Kira Engine

**Computational selectivity analysis for parasite-vs-human target pairs under sparse medicinal-chemistry data.**

Kira is a research repository for retrospective translational analysis, scaffold-aware selectivity benchmarking, assay-aware data curation, and residue-level mechanistic hypothesis generation. The current scientific focus is narrow by design: parasite-vs-human selectivity under sparse public assay data, not a universal drug-discovery engine.

## Current scientific scope

Kira currently supports:

- curated parasite-vs-human comparator target pairs
- retrospective selectivity analysis from public assay evidence
- compound-conditioned v4 benchmarking
- scaffold-aware leakage control using `target_pair_id::murcko_scaffold`
- assay-aware v5 expansion from ChEMBL
- exact-core aggregation with replicate handling and conflict detection

Kira does **not** currently provide:

- a production selectivity predictor for arbitrary new compounds
- wet-lab validation
- a full protein-ligand binding physics engine
- a universal multiscale biology simulator

## Key results so far

### v4: compound-conditioned scaffold-aware benchmark

v4 repaired the core v3 limitation by making the benchmark genuinely compound-conditioned.

| Quantity | Value |
|---|---:|
| Total mapped rows | 310 |
| Primary labeled candidate rows | 237 |
| Primary trainable rows | 235 |
| Full v4 feature width | 379 |
| Cross-validation | 5-fold StratifiedGroupKFold |
| Grouping | `target_pair_id::murcko_scaffold` |

v4 ablation summary:

| Ablation | Pooled AUROC | Macro pair AUROC | Macro pair Spearman | Interpretation |
|---|---:|---:|---:|---|
| A0_pair_only | 0.794 | 0.366 | -0.206 | Pair-only fails on the main within-pair metric. |
| A1_compound_only | 0.835 | 0.705 | 0.343 | Compound chemistry carries most signal. |
| A2_compound_plus_pair | 0.900 | 0.708 | 0.357 | Best classifier by macro pair AUROC. |
| A2b_pair_plus_side | 0.898 | 0.697 | 0.348 | Side summaries do not improve classification. |
| A2c_pair_plus_compat | 0.896 | 0.686 | 0.345 | Compatibility block is benchmarkable but not yet beneficial. |
| A3_full_v4 | 0.896 | 0.695 | 0.361 | Best ranking signal only, not best classifier. |

**Honest v4 claim:** compound chemotype explains most of the current predictive signal. The coarse pocket-compatibility features are scientifically motivated and benchmarkable, but they have not yet shown clear within-pair classification gain over compound + pair.

### v5: assay-aware expansion and exact core

v5 expands the data substrate before adding a more complex model.

| Quantity | Value |
|---|---:|
| Curated ChEMBL activity rows | 11,703 |
| Candidate evidence rows | 12,091 |
| Exact matched-ratio candidate rows | 1,227 |
| Exact-core rows after replicate collapse | 114 |
| Trainable exact-core rows | 110 |
| Conflicting exact-core rows | 4 |

Trainable v5 exact-core class balance:

| Pair | Trainable rows | Positive | Negative | Unique scaffolds |
|---|---:|---:|---:|---:|
| LmDHFR | 11 | 9 | 2 | 10 |
| LmPTR1 | 1 | 0 | 1 | 1 |
| SmDHODH | 3 | 0 | 3 | 3 |
| SmHDAC8 | 73 | 1 | 72 | 37 |
| TbCathB | 6 | 1 | 5 | 4 |
| TbPDEB1 | 16 | 0 | 16 | 15 |

**Honest v5 claim:** v5 successfully expands and cleans the assay-aware evidence substrate, but the strict exact core is too class-degenerate for a strong six-pair classifier benchmark. v5 does **not** yet prove a new compound-dominance shift. That requires the upcoming B0/B1/B2 benchmark on an exact or tiered core.

## Installation

### Option A: conda/mamba

```bash
mamba env create -f environment.yml
mamba activate kira
python -m pip install -e ".[structural,dev]"
```

### Option B: Docker

```bash
docker build -t kira-engine .
docker run --rm -it -v "$PWD":/app kira-engine
```

### Option C: existing Python environment

```bash
python -m pip install -e ".[structural,dev]"
```

RDKit is required for v4/v5 structure resolution and feature generation.

## Reproduce v4

```bash
pytest -q tests/test_selectivity_v4_data.py tests/test_selectivity_v4_features.py

python -m kira.experiments.selectivity_v4_data
python -m kira.experiments.run_selectivity_v4
```

Expected committed v4 state:

- `data/processed/selectivity_v4_summary.json`
- `results/selectivity_v4/summary.json`
- 235 trainable rows
- 379 full features
- A2 compound + pair is the best classifier by macro pair AUROC

## Reproduce v5 expansion

```bash
pytest -q tests/test_selectivity_v5_expand_data.py tests/test_selectivity_v5_exact_core.py

python -m kira.experiments.selectivity_v5_expand_data \
  --pair-config data/reference/selectivity_v5_target_pairs.csv \
  --output-csv data/processed/selectivity_v5_candidate_rows.csv \
  --output-summary-json data/processed/selectivity_v5_expansion_summary.json \
  --set-chunk-size 10 \
  --timeout-seconds 25

python -m kira.experiments.selectivity_v5_exact_core \
  --candidate-csv data/processed/selectivity_v5_candidate_rows.csv \
  --output-csv data/processed/selectivity_v5_exact_core_rows.csv \
  --output-summary-json data/processed/selectivity_v5_exact_core_summary.json
```

## Notebooks

The notebooks in `notebooks/` are meant as readable walkthroughs:

- `00_v0_to_v5_evolution.ipynb` — conceptual and result evolution
- `01_reproduce_v4_benchmark.ipynb` — v4 data and ablation reproduction
- `02_v5_expansion_exact_core.ipynb` — v5 expansion and exact-core aggregation

## Repository layout

```text
src/kira/experiments/
  selectivity_v4_data.py
  selectivity_v4_features.py
  run_selectivity_v4.py
  selectivity_v5_expand_data.py
  selectivity_v5_exact_core.py

data/reference/
  selectivity_v5_target_pairs.csv

data/processed/
  selectivity_v4_summary.json
  selectivity_v5_expansion_summary.json
  selectivity_v5_exact_core_summary.json
  selectivity_v5_exact_core_rows.csv
  selectivity_v5_candidate_rows.csv

results/selectivity_v4/
  summary.json
  per_pair_metrics.csv

notebooks/
  00_v0_to_v5_evolution.ipynb
  01_reproduce_v4_benchmark.ipynb
  02_v5_expansion_exact_core.ipynb
```

## Scientific claim discipline

Use:

> Kira v4 is a scaffold-aware, compound-conditioned selectivity benchmark. Its first result shows that compound chemotype currently explains most predictive signal. Kira v5 adds assay-aware data expansion and exact-core aggregation; the strict exact core is cleaner but class-degenerate, so the next benchmark must test a tiered evidence core before making stronger modeling claims.

Do not use:

> Kira is a universal causal drug-discovery engine.

Do not use:

> v5 proves compound dominance.

The v4 benchmark supports compound dominance on the current v4 substrate. The v5 exact core supports a data-quality conclusion, not yet a modeling conclusion.
