# Results Provenance

All results generated on the machine and date below.
To reproduce, install kira-engine and run the commands listed.

## Selectivity v3 Experiment
- Command: `python -m kira.experiments.run_selectivity_v3`
- Output: `selectivity_v3_results.txt`
- Key result: Mean LODO 0.519 vs 0.429 ESM-2 baseline (+21%)

## Case Study: CHEMBL155771
- Command: `python -m kira.experiments.case_study_chembl155771`
- Output: `case_study_chembl155771.txt`
- Key result: 4 selectivity-driving residue positions identified

## Data Sources
- Leishmaniasis: `data/leishmania/leish_selectivity.csv` (116 compounds)
- Trypanosomiasis: `data/trypanosoma/tryp_selectivity_expanded.csv` (104 compounds)
- Schistosomiasis: `data/publication/discovery_candidates.csv` (17 compounds)
- Pocket sequences: Mori et al. 2021 (PDB 6UY4), PMC9486128 (PDB 6HSF)
Generated: Tue Apr  7 10:27:41 CDT 2026
Python: Python 3.11.14
Platform: arm64
