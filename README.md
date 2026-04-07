# Kira Engine — Open Causal Discovery for Neglected Tropical Diseases

A unified computational engine for drug selectivity analysis, physics validation, and mechanistic explanation across neglected tropical diseases.

[![CI](https://github.com/Danny2045/Kira/actions/workflows/ci.yml/badge.svg)](https://github.com/Danny2045/Kira/actions/workflows/ci.yml)

## Core Finding

**Binding-site divergence predicts cross-disease drug selectivity where protein language models fail.**

ESM-2 global embeddings achieve 0.99 cosine similarity between parasite and human orthologs but fail at cross-disease selectivity prediction (LODO AUROC 0.38–0.55). Per-residue pocket divergence features — physicochemical property differences at binding-site positions — improve cross-disease transfer by 21% (mean LODO 0.519 vs 0.429).

| Metric | ESM-2 (Script 19) | Pocket Features | Delta |
|--------|:-:|:-:|:-:|
| 5-fold CV AUROC | 0.895 | 0.824 | -0.071 |
| **LODO Schistosomiasis** | 0.382 | **0.606** | **+0.224** |
| **LODO Leishmaniasis** | 0.547 | **0.637** | **+0.090** |
| LODO Trypanosomiasis | 0.359 | 0.314 | -0.045 |
| **Mean LODO** | **0.429** | **0.519** | **+0.090** |

## What This Does

Kira Engine is three things in one package:

**1. Drug selectivity platform** — Systematic selectivity analysis across three neglected tropical diseases (schistosomiasis, trypanosomiasis, leishmaniasis). 237 compounds, 6 target pairs, 311 selectivity comparisons. Identifies which compounds prefer the parasite target over the human ortholog.

**2. Physics validation engine** — Validates AI-generated protein structures (from Boltz-2, AlphaFold 3, Chai-1, etc.) with 8 physics checks: steric clashes, bond geometry, Ramachandran, peptide planarity, chirality, rotamers, Lennard-Jones energy, disulfide geometry. Produces a Trust Report with accept/relax/discard recommendation.

**3. Mechanistic explanation layer** — Explains WHY selectivity exists at the residue level. Extracts binding pockets, computes per-position physicochemical divergence (hydrophobicity, charge, volume), and attributes selectivity to specific pocket differences.

## Quick Start

```bash
# Install
conda activate bio-builder  # or any Python 3.11+ environment
pip install -e ".[dev]"

# Validate a protein structure
kira validate structure.pdb

# Run the selectivity experiment
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

SmDHODH and HsDHODH have ESM-2 cosine similarity of **0.9897** — nearly identical by global representation. Yet compounds achieve **30.8x selectivity** between them. Global protein embeddings miss binding-site-level divergence. The few residues that differ at the ubiquinone binding site (Ser53→Leu59 hydrophobicity flip, Val358→Pro364 flexibility change) drive the selectivity window that ESM-2 cannot see.

## Feature Importance: What Drives Selectivity

| Feature | Importance | Interpretation |
|---------|:-:|---|
| pocket_divergence | 28% | How different the binding pocket is |
| std_physicochemical_distance | 20% | Whether divergence is concentrated at specific positions |
| n_volume_changes | 12% | Positions with large sidechain size changes |
| mean_physicochemical_distance | 9% | Average property difference across pocket |

## CLI Commands

```bash
kira validate structure.pdb      # Physics validation with Trust Report
kira info structure.pdb          # Structure summary
kira query --target SmTGR        # Target essentiality lookup
kira evaluate --predictions r.csv --ground-truth g.csv  # Benchmark scoring
kira selectivity --parasite SmDHODH --human HsDHODH     # Selectivity analysis
```

## Project Structure

```
src/kira/
├── physics/              # Physics validation engine (JAX-accelerated)
│   ├── core/             # PDB parser, topology, geometry, LJ energy kernels
│   ├── checks/           # 8 physics checks + composite scorer
│   └── config.py         # All thresholds (YAML-overridable)
├── causality/            # Mechanistic explanation layer
│   ├── binding_site.py   # Pocket extraction + comparison
│   ├── divergence.py     # Local vs global divergence profiling
│   ├── energy_decomp.py  # Per-residue energy decomposition
│   └── selectivity_map.py # Selectivity attribution
├── experiments/          # Scientific experiments
│   ├── run_selectivity_v3.py     # Main selectivity experiment
│   ├── selectivity_features.py   # Pocket divergence feature extraction
│   └── validate_physics.py       # Real PDB validation
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
5. **No experimental validation.** All results are computational. The selectivity claims need biochemical confirmation.

## License

MIT — see [LICENSE](LICENSE).

## Citation

```bibtex
@software{ngabonziza2026kira,
  author = {Ngabonziza, Daniel},
  title = {Kira Engine: Open Causal Discovery for Neglected Tropical Diseases},
  year = {2026},
  url = {https://github.com/Danny2045/Kira}
}
```
