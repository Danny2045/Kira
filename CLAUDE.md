# Kira - Selectivity-First Drug Repurposing Platform

## What this is
3 diseases (schistosomiasis, trypanosomiasis, leishmaniasis), 2,699 compounds,
50 targets, 311 selectivity labels, molecular docking, ESM-2 ML prediction.
21 scripts + clean rerun. Public repo: github.com/Danny2045/Kira

## Environment
- conda activate bio-builder
- Python 3.11, RDKit, pandas, scikit-learn, PyTorch, ESM-2, AutoDock Vina, meeko
- All data in data/ directory

## Key data files
- data/processed/kira_selectivity_analysis.csv (schistosomiasis, 91 records)
- data/trypanosoma/tryp_selectivity_expanded.csv (trypanosomiasis, 104 records)
- data/leishmania/leish_selectivity.csv (leishmaniasis, 116 records)
- data/processed/smtgr_docking_results.csv (43 docking results)
- data/models/esm2_embeddings_v3.npz (11 protein embeddings)
- data/models/clean_model_results.json (ML evaluation results)

## Key findings
- 61.7% of antiparasitic chemical matter is non-selective across 311 comparisons
- PTR1 is the best target (68.2x median selectivity, 34 compounds >10x)
- SmHDAC8 is a selectivity trap (90.4% non-selective)
- ML: AUROC 0.895 on 5-fold CV but LODO fails (0.35-0.55)
- ESM-2 cosine 0.99 for SmDHODH-HsDHODH yet 30.8x selectivity exists
- Protein pair features: 8.6% importance, do not transfer across diseases

## Coding rules
- Run code after writing to verify it works
- Use existing data files, don't re-query ChEMBL unless needed
- Commit with descriptive messages
- Document limitations honestly
- Always test on held-out data, not just random CV
- Scientific claims must be supported by the numbers
