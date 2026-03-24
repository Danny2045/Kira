# Kira — Computational Drug Repurposing for Schistosomiasis

A multi-pathway drug repurposing pipeline that integrates target-based bioactivity, structural similarity, whole-organism evidence, ADMET filtering, cross-species selectivity analysis, and supply chain assessment to identify candidates for schistosomiasis treatment.

## Core Finding

**Systematic selectivity analysis of 91 anti-schistosomal compounds reveals that 90.4% of SmHDAC8 chemical matter is non-selective against the human orthologue, while SmDHODH harbors compounds with >10-fold parasite selectivity.**

Four compounds exceed 10x selectivity for the parasite enzyme over the human version. The top candidate — CHEMBL155771 — achieves 23 nM potency against *S. mansoni* DHODH with 30.8x selectivity and near-perfect drug-likeness (QED = 0.89). Atovaquone, an approved antimalarial on the WHO Essential Medicines List and available in sub-Saharan Africa, shows 6-fold SmDHODH selectivity, supporting its prioritization for anti-schistosomal evaluation.

## Why This Matters

Schistosomiasis infects ~200 million people, overwhelmingly in sub-Saharan Africa. Global treatment depends almost entirely on a single drug (praziquantel). Drug repurposing — finding new uses for existing approved drugs — offers the fastest path to alternatives. But most computational repurposing studies rank compounds by target potency alone, without checking whether the compound also hits the equivalent human protein. This pipeline adds that check systematically and shows it changes the landscape dramatically.

## Pipeline Architecture

```
PrimeKG Knowledge Graph          → Disease context (human biology, known drugs)
ChEMBL Target-Based Activity     → Parasite protein IC50 data (226 measurements, 10 targets)
ChEMBL Whole-Organism Activity   → Phenotypic worm-killing data (494 records)
RDKit Structural Similarity      → Morgan fingerprints, Tanimoto to known actives
Composite Ranking (7 signals)    → Potency, target essentiality, confidence, drug stage,
                                    multi-target, similarity, whole-organism
ADMET Filtering                  → Lipinski Rule of Five, TPSA, QED drug-likeness
Human Orthologue Selectivity     → Cross-species IC50 ratios (91 compounds, 2 targets)
Literature Novelty (PubMed)      → Publication count per compound
WHO EML + Supply Chain           → African availability, cost tier, route of administration
```

## Key Results

| Target | Compounds Tested | With Selectivity Data | Non-Selective (%) | Best Selectivity |
|--------|:---:|:---:|:---:|:---:|
| SmHDAC8 | 100 | 73 | 90.4% | 11.0x (1 compound) |
| SmDHODH | 21 | 18 | 38.9% | 30.8x (3 compounds >10x) |
| SmTGR | 43 | 0 | Unknown | Experimental gap |
| Sirtuin | 5 | 0 | Unknown | — |
| SmVKR2 | 7 | 0 | Unknown | — |

### Top Selective Compounds (≥10x parasite selectivity)

| Compound | Target | Parasite IC50 (nM) | Human IC50 (nM) | Selectivity | QED |
|----------|--------|:---:|:---:|:---:|:---:|
| CHEMBL155771 | SmDHODH | 23 | 709 | 30.8x | 0.89 |
| CHEMBL4474026 | SmDHODH | 227 | 4,600 | 20.3x | 0.88 |
| CHEMBL4855490 | SmHDAC8 | 100 | 1,100 | 11.0x | 0.44 |
| CHEMBL4452960 | SmDHODH | 78 | 808 | 10.4x | 0.89 |

### Top Translational Candidate

**Atovaquone** — Approved antimalarial. IC50 = 430 nM vs SmDHODH. 6.0x selective over human DHODH. WHO Essential Medicines List. Available in sub-Saharan Africa as Malarone (oral). Known safety profile in target population.

## How It Was Built

13 scripts, each building on the last. The git history shows the complete progression from empty folder to selectivity finding.

| Script | What It Does |
|--------|-------------|
| `01_explore_primekg.py` | Download and query PrimeKG knowledge graph for schistosomiasis |
| `02_query_chembl.py` | Find *S. mansoni* protein targets and bioactivity data in ChEMBL |
| `03_build_eval_set.py` | Build ground-truth evaluation set (228 compounds, 3 classes) |
| `04_rank_and_evaluate.py` | First composite ranking algorithm (5 signals) + AUROC evaluation |
| `05_structural_similarity.py` | Add structural similarity pathway (RDKit Morgan fingerprints) |
| `06_whole_organism.py` | Add whole-organism activity pathway from ChEMBL phenotypic data |
| `07_admet_and_report.py` | ADMET filtering (Lipinski, TPSA, QED) + first candidate report |
| `08_harden_benchmark.py` | Expand negatives (16→50), train/test split, bootstrap CIs |
| `09_novelty_filter.py` | PubMed literature search for each top compound |
| `10_selectivity_analysis.py` | **Cross-species selectivity** against human orthologues |
| `11_selectivity_rerank.py` | Selectivity-adjusted re-ranking + definitive shortlist |
| `12_publication_analysis.py` | Publication tables, per-target statistics, hard negatives |
| `13_supply_chain_and_final.py` | WHO EML, African availability, deployment-adjusted scoring |

## Project Structure

```
kira/
├── data/
│   ├── raw/                          # PrimeKG (~250 MB, not in git)
│   ├── processed/                    # Pipeline outputs (CSVs)
│   ├── eval/                         # Evaluation sets (v1, v2, dev/test split)
│   ├── publication/                  # Publication-ready tables
│   └── reports/                      # Human-readable reports
├── src/kira/                         # Package structure (future)
├── tests/                            # Tests (future)
├── docs/                             # Documentation
│   └── kira-ground-to-god.md         # Complete scientific breakdown
├── mission.md                        # Why Kira exists
├── plan.md                           # Phase 1 milestones
├── agents.md                         # Coding rules
├── architecture.md                   # Design decisions
└── 01-13 scripts                     # The pipeline
```

## Requirements

- Python 3.11 (conda environment recommended)
- RDKit (cheminformatics)
- pandas, NumPy, NetworkX
- chembl_webresource_client (ChEMBL API)
- scikit-learn (evaluation metrics)
- requests (PubMed API)
- Internet connection (for ChEMBL and PubMed queries; cached after first run)

### Setup

```bash
conda create -n bio-builder python=3.11
conda activate bio-builder
conda install -c conda-forge rdkit numpy pandas matplotlib jupyterlab biopython networkx
pip install chembl_webresource_client scikit-learn pyarrow tqdm requests pytest
```

### Running the Pipeline

```bash
conda activate bio-builder
cd ~/kira

# Run scripts in order (each builds on the previous)
python 01_explore_primekg.py      # ~5 min (downloads PrimeKG)
python 02_query_chembl.py         # ~10 min (queries ChEMBL)
python 03_build_eval_set.py       # <1 min
python 04_rank_and_evaluate.py    # <1 min
python 05_structural_similarity.py # <1 min
python 06_whole_organism.py       # ~5 min (queries ChEMBL)
python 07_admet_and_report.py     # <1 min
python 08_harden_benchmark.py     # <1 min
python 09_novelty_filter.py       # ~1 min (queries PubMed)
python 10_selectivity_analysis.py # ~10 min (queries ChEMBL)
python 11_selectivity_rerank.py   # ~1 min
python 12_publication_analysis.py # <1 min
python 13_supply_chain_and_final.py # ~3 min (queries ChEMBL)
```

Total pipeline runtime: ~35 minutes (mostly network queries, cached after first run).

## Limitations

This is a prioritization tool, not a clinical recommendation engine.

1. Selectivity assessed for SmHDAC8 and SmDHODH only; SmTGR (the most biologically compelling target) lacks human orthologue data
2. Selectivity ratios computed from ChEMBL data across different labs and assay conditions
3. Only direct orthologue selectivity assessed; broader off-target pharmacology not evaluated
4. Whole-organism dose-response data available for 6 compounds only
5. No molecular docking or structural modeling performed
6. AUROC 1.0 reflects easy discrimination task (known actives vs unrelated drugs)
7. PubMed novelty screen uses exact string matching
8. Supply chain data is manually curated (WHO EML 2023 edition)

## Next Steps

1. **Experimental validation:** Test CHEMBL155771 and CHEMBL4452960 in *S. mansoni* whole-worm killing assay
2. **Atovaquone validation:** Test against adult *S. mansoni* at clinically relevant concentrations
3. **SmTGR selectivity:** Test top SmTGR inhibitors against human thioredoxin reductase 1
4. **Structural modeling:** AlphaFold 3 docking against data-desert targets
5. **Benchmark hardening:** Add structurally similar experimentally confirmed inactive compounds

## Context

Built by Daniel Ngabonziza in Kigali, Rwanda. The pipeline addresses diseases that pharmaceutical markets systematically underserve. The supply chain layer encodes deployment constraints specific to East African clinical settings — knowledge that cannot be replicated from public databases alone.

## License

MIT

## Citation

If you use this pipeline or its findings, please cite:

```
Ngabonziza, D. (2026). Systematic cross-species selectivity analysis reveals SmDHODH
as the most tractable target for schistosomiasis drug repurposing. bioRxiv [preprint].
```
