# Kira Engine — Computational Selectivity Analysis for Neglected Tropical Diseases

A computational repository for retrospective selectivity analysis, approximate structure checks, and residue-level mechanistic hypothesis generation across neglected tropical diseases.

[![CI](https://github.com/Danny2045/Kira/actions/workflows/ci.yml/badge.svg)](https://github.com/Danny2045/Kira/actions/workflows/ci.yml)

## Current Scientific Scope

- Benchmarks target-pair pocket descriptors against historical selectivity labels across neglected tropical disease target pairs.
- Generates residue-level mechanistic hypotheses from curated pockets and approximate residue-energy heuristics.
- Performs approximate structure checks using currently implemented steric, Lennard-Jones, and backbone-dihedral routines.
- Preserves archived translational and docking analyses as auditable scientific evidence.
- Does not currently provide compound-conditioned selectivity prediction, explicit protein-ligand binding energetics, a full active docking engine, or an eight-check active validation workflow.

## Core Finding

**Binding-site divergence features show a modest cross-disease transfer signal where global protein similarity is often uninformative.**

ESM-2 global embeddings can place parasite and human orthologs extremely close in representation space while remaining weakly informative for this small leave-one-disease-out benchmark. In the current proof-of-concept setting, curated pocket divergence features outperform the ESM-2 baseline on mean LODO AUROC (0.519 vs 0.429). These features are target-pair-level descriptors rather than compound-conditioned predictors, so the result should be interpreted as benchmark association rather than production selectivity prediction.

| Metric | ESM-2 (Script 19) | Pocket Features | Delta |
|--------|:-:|:-:|:-:|
| 5-fold CV AUROC | 0.895 | 0.824 | -0.071 |
| **LODO Schistosomiasis** | 0.382 | **0.606** | **+0.224** |
| **LODO Leishmaniasis** | 0.547 | **0.637** | **+0.090** |
| LODO Trypanosomiasis | 0.359 | 0.314 | -0.045 |
| **Mean LODO** | **0.429** | **0.519** | **+0.090** |

## What This Does

Kira Engine is three things in one package:

**1. Retrospective selectivity analysis** — Systematic analysis of historical parasite-versus-human selectivity evidence across three neglected tropical diseases (schistosomiasis, trypanosomiasis, leishmaniasis). 237 compounds, 6 target pairs, 311 selectivity comparisons. Summarizes observed parasite-versus-human selectivity patterns.

**2. Approximate structure-checking toolkit** — Runs the currently implemented structure checks on protein models and experimental structures. The active validator reports steric clashes, Lennard-Jones summaries, and backbone-dihedral counts. These outputs are heuristic structure-quality summaries, not full physical validation.

**3. Mechanistic hypothesis layer** — Proposes residue-level rationales for selectivity. Extracts centroid-defined or curated pockets, computes per-position physicochemical divergence (hydrophobicity, charge, volume), and highlights residue differences associated with the heuristic score.

## Quick Start

```bash
# Install
conda activate bio-builder  # or any Python 3.11+ environment
pip install -e ".[dev]"

# Run approximate structure checks on a protein structure
kira validate structure.pdb

# Run the target-pair pocket-feature benchmark
python -m kira.experiments.run_selectivity_v3

# Run tests
pytest tests/ -v
```

## Key Selectivity Results (Three-Disease Platform)

### Selectivity landscape: 311 comparisons across 3 diseases

| Target | Disease | Compounds | Non-Selective | Median Ratio | Best Ratio |
|--------|---------|:-:|:-:|:-:|:-:|
| PTR1 | Leishmaniasis | 45 | 24.4% | 68.2x | >100x |
| SmDHODH | Schistosomiasis | 18 | 38.9% | 4.7x | 30.8x |
| DHFR-TS | Leishmaniasis | 69 | 58.0% | 3.6x | 16.1x |
| PDEB1 | Trypanosomiasis | 68 | 89.7% | 1.6x | 7.0x |
| SmHDAC8 | Schistosomiasis | 73 | 90.4% | 1.0x | 11.0x |
| Cathepsin B | Trypanosomiasis | 36 | 94.4% | 0.8x | 4.1x |
| SmTGR | Schistosomiasis | 43 | 88.0% | ~1x | ~1.3x |

### Top selective compounds

| Compound | Target | Parasite IC50 (nM) | Selectivity | QED |
|----------|--------|:-:|:-:|:-:|
| CHEMBL155771 | SmDHODH | 23 | 30.8x | 0.89 |
| Atovaquone | SmDHODH | ~140 | 6.0x | — |
| CHEMBL4474026 | SmDHODH | 227 | 20.3x | 0.88 |
| CHEMBL4855490 | SmHDAC8 | 100 | 11.0x | 0.44 |

### The ESM-2 blind spot

SmDHODH and HsDHODH have ESM-2 cosine similarity of **0.9897** — nearly identical by global representation. Yet compounds achieve **30.8x selectivity** between them. Global protein embeddings can miss pocket-level divergence. The few residues that differ at the ubiquinone site (Ser53→Leu59 hydrophobicity flip, Val358→Pro364 flexibility change) may contribute to a selectivity window and illustrate the kind of local divergence the benchmark is designed to capture.

## Feature Importance: Signals Associated With Selectivity

| Feature | Importance | Interpretation |
|---------|:-:|---|
| pocket_divergence | 28% | How different the binding pocket is |
| std_physicochemical_distance | 20% | Whether divergence is concentrated at specific positions |
| n_volume_changes | 12% | Positions with large sidechain size changes |
| mean_physicochemical_distance | 9% | Average property difference across pocket |

## CLI Commands

```bash
kira validate structure.pdb
kira info structure.pdb
kira evaluate --predictions predictions.csv --ground-truth ground_truth.csv
kira selectivity --parasite parasite.pdb --human human.pdb --ligand-x 0.0 --ligand-y 0.0 --ligand-z 0.0
```

Note: the active selectivity CLI uses parasite and human PDB files plus a ligand centroid, not target names or an explicit docked ligand.

## Project Structure

```
src/kira/
├── physics/              # Approximate structure-checking utilities
│   ├── core/             # PDB parser, topology, geometry, LJ energy kernels
│   ├── checks/           # Active steric clash logic and validation scaffolding
│   └── config.py         # All thresholds (YAML-overridable)
├── causality/            # Mechanistic hypothesis modules
│   ├── binding_site.py   # Pocket extraction + comparison
│   ├── divergence.py     # Local vs global divergence profiling
│   ├── energy_decomp.py  # Residue-energy heuristic decomposition
│   └── selectivity_map.py # Residue-level heuristic attribution
├── experiments/          # Scientific experiments
│   ├── run_selectivity_v3.py     # Target-pair selectivity benchmark
│   ├── selectivity_features.py   # Pocket divergence feature extraction
│   └── validate_physics.py       # Plausibility benchmark on selected PDB structures
├── scoring.py            # Compound ranking functions
├── targets.py            # Target essentiality + ortholog map
├── chemistry.py          # ADMET calculations
├── drugs.py              # Drug identifiers
└── cli.py                # Unified CLI (5 commands)

tests/                    # 194 tests
data/                     # Selectivity CSVs, docking structures, models
```

## Limitations

1. **Selectivity ratios are from cross-lab IC50 comparisons.** Different labs, assay conditions, and years. Systematic error is not quantified.
2. **Physics engine uses element-level LJ parameters**, not full atom-type assignment. This is a ~30% approximation for some atom pairs.
3. **Pocket sequences are from literature analysis**, not automated structural alignment. Some target pairs have approximate pocket definitions.
4. **237 compounds across 6 target pairs.** The dataset is small by ML standards. Results should be interpreted as proof-of-concept, not production models.
5. **Active benchmark features are target-pair-level, not compound-conditioned.** Current benchmark results do not establish a true compound-specific selectivity predictor.
6. **Residue-level attribution is not explicit protein-ligand binding energetics.** The active attribution path uses curated pockets plus protein-only residue-energy heuristics.
7. **No experimental validation.** All results are computational. The selectivity claims need biochemical confirmation.

## License

MIT — see [LICENSE](LICENSE).

## Citation

```bibtex
@software{ngabonziza2026kira,
  author = {Ngabonziza, Daniel},
  title = {Kira Engine: Computational Selectivity Analysis for Neglected Tropical Diseases},
  year = {2026},
  url = {https://github.com/Danny2045/Kira}
}
```
