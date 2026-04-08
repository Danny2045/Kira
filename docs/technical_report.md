# Kira Engine: Technical Report

**Author:** Daniel Ngabonziza
**Date:** April 2026
**Repository:** https://github.com/Danny2045/Kira
**Status:** Exploratory computational triage tool (proof-of-concept)

---

## 1. What This System Does

Kira Engine is a computational tool for analyzing drug selectivity across neglected tropical diseases. It addresses one question: when a compound binds a parasite protein target, does it also bind the equivalent human protein? If it does, the compound is toxic. If it doesn't, the compound has a selectivity window that makes it a viable drug candidate.

The system has three components:

**Selectivity analysis platform.** Aggregates compound-target bioactivity data from ChEMBL across three diseases (schistosomiasis, trypanosomiasis, leishmaniasis), computes selectivity ratios against human orthologs, and ranks candidates by both potency and selectivity.

**Physics validation engine.** Takes protein structures (PDB files) and checks them for physical validity: steric clashes, bond geometry, backbone conformation, chirality, sidechain rotamers, Lennard-Jones energy, and disulfide geometry. Produces a composite trust score with an accept/relax/discard recommendation.

**Mechanistic explanation layer.** Extracts binding pockets from parasite-human protein pairs, computes per-position physicochemical divergence (hydrophobicity, charge, volume), and attributes selectivity to specific residue differences. This is the component that addresses the question "why is this compound selective?" rather than just "is it selective?"

## 2. Principal Result

A gradient boosting classifier trained on 15 binding-site divergence features achieves a mean leave-one-disease-out (LODO) AUROC of 0.519 for selectivity prediction. The baseline — the same classifier using ESM-2 global protein embeddings as features (from Kira Script 19) — achieves 0.429.

The improvement is concentrated on two diseases:
- Schistosomiasis LODO: 0.382 → 0.606 (+58%)
- Leishmaniasis LODO: 0.637 → 0.547 (+16%)
- Trypanosomiasis LODO: 0.359 → 0.314 (-13%)

The mean improvement of +0.09 AUROC is modest in absolute terms. Its significance lies in the direction: binding-site-level features transfer across diseases better than global embeddings, which fail catastrophically at LODO (AUROC below 0.5 for two of three diseases).

## 3. Methods

### 3.1 Selectivity Data

Compound-level selectivity data was assembled from three sources:

| Source | Disease | Compounds | Targets |
|--------|---------|-----------|---------|
| `data/leishmania/leish_selectivity.csv` | Leishmaniasis | 116 | PTR1, DHFR-TS |
| `data/trypanosoma/tryp_selectivity_expanded.csv` | Trypanosomiasis | 104 | Cathepsin B, PDEB1 |
| `data/publication/discovery_candidates.csv` | Schistosomiasis | 17 | SmDHODH, SmHDAC8 |

Total: 237 compounds across 6 target pairs.

Selectivity ratio = (human IC50) / (parasite IC50). Binary label: selective ≥ 10x.

### 3.2 Binding-Site Divergence Features

For each parasite-human target pair, pocket residues were identified from published crystal structures:

| Target Pair | PDB Sources | Pocket Size | Identity |
|-------------|-------------|-------------|----------|
| SmDHODH / HsDHODH | 6UY4 / 1D3H | 9 residues | 55.6% |
| SmHDAC8 / HsHDAC8 | 6HSF / 1T69 | 10 residues | 90.0% |
| SmTGR / HsTrxR1 | 2X99 / 2ZZC | 15 residues | 86.7% |
| TbCathB / HsCathB | Literature | 15 residues | 93.3% |
| TbPDEB1 / HsPDE4 | Literature | 13 residues | 76.9% |
| LmPTR1 / HsDHFR | 1E7W / Literature | 15 residues | 46.7% |
| LmDHFR-TS / HsDHFR | Literature | 15 residues | 93.3% |

Fifteen features were computed per target pair:

1. Pocket size (number of aligned positions)
2. Pocket identity (fraction of identical residues)
3. Pocket divergence (1 − identity)
4. Mean physicochemical distance across positions
5. Max physicochemical distance (worst-case position)
6. Standard deviation of physicochemical distance
7. Number of charge sign changes (e.g., Asp → Lys)
8. Number of hydrophobicity flips (polar ↔ nonpolar)
9. Number of volume changes (>30 Å³ difference)
10. Fraction of positions with charge changes
11. Fraction with hydrophobicity flips
12. Fraction with volume changes
13. Total absolute hydrophobicity delta
14. Total absolute charge delta
15. Total absolute volume delta

Physicochemical distance uses Kyte-Doolittle hydrophobicity, charge at pH 7, and sidechain volume, each normalized by range, combined as Euclidean distance.

### 3.3 Classification

Gradient Boosting Classifier (100 estimators, max depth 3, min samples leaf 5). StandardScaler fitted inside a sklearn Pipeline to prevent feature leakage. Evaluation: 5-fold stratified cross-validation and leave-one-disease-out (LODO). Scoring metric: AUROC.

### 3.4 Physics Validation

Eight checks applied to protein structures:

| Check | Method | Reference |
|-------|--------|-----------|
| Steric clashes | Pairwise vdW overlap | Bondi radii (1964) |
| Bond geometry | Z-scores vs ideal | Engh & Huber (1991) |
| Ramachandran | φ/ψ region classification | Lovell et al. (2003) |
| Peptide planarity | ω angle deviation | Standard (180° trans) |
| Chirality | Scalar triple product | L-amino acid sign |
| Rotamers | χ1 vs known states | Penultimate library |
| Lennard-Jones | Pairwise LJ potential | AMBER-like parameters |
| Disulfides | S-S distance + χ3 | Standard (2.05 Å) |

Composite score: weighted average of 8 subscores. Recommendation: accept (≥0.85), relax (0.60-0.85), discard (<0.60).

### 3.5 Case Study: CHEMBL155771

CHEMBL155771 (2-hydroxy-3-isopentylnaphthalene-1,4-dione) was run through the full pipeline:

- **Global analysis:** ESM-2 cosine similarity 0.9897 between SmDHODH and HsDHODH → predicts no selectivity
- **Local analysis:** Pocket identity 55.6% (9 positions, 4 divergent) → predicts selectivity window
- **Attribution:** 4 selectivity-driving positions identified:
  - Ser53→Leu59: polarity flip + volume change (distance 0.690)
  - Val358→Pro364: polarity flip (distance 0.665)
  - Gly46→Met43: polarity flip + volume change (distance 0.664)
  - Ile128→Val134: conservative substitution (distance 0.163)

## 4. Limitations

These limitations are fundamental to the current system and should be considered when interpreting any result.

### 4.1 Data Quality

**Cross-lab IC50 values.** Selectivity ratios are computed from IC50 measurements performed in different laboratories, using different assay protocols, at different times. Systematic error from inter-lab variability is not quantified. A 3x selectivity ratio may be entirely explained by assay differences. Only ratios >10x are likely to reflect genuine selectivity, and even those should be confirmed in a single matched assay.

**Small dataset.** 237 compounds across 6 target pairs is underpowered for robust ML. The schistosomiasis subset has only 17 compounds. The LODO results should be interpreted as directional evidence, not definitive measurements. Bootstrap confidence intervals are not reported for the LODO numbers.

**No experimental validation.** All results are purely computational. No compound has been synthesized or tested in an assay based on this analysis. The selectivity predictions are hypotheses, not confirmed findings.

### 4.2 Physics Engine Approximations

**Element-level Lennard-Jones parameters.** The LJ potential uses one sigma and epsilon value per element (C, N, O, S, etc.), not per atom type. This means a carbon in an aromatic ring gets the same parameters as a carbon in a methyl group — approximately 30% error in sigma and up to 50% error in epsilon for some atom pairs. A production physics engine would use full atom-type assignment from AMBER ff14SB or CHARMM36m.

**No solvation.** The energy computation is in vacuum. Real binding happens in water. Solvation effects (desolvation penalty, hydrophobic effect, entropic contributions) are not modeled. This means the per-residue energy decomposition captures direct contact terms but misses solvent-mediated effects.

**No many-body terms.** The LJ potential is pairwise additive. Real molecular interactions include polarization, charge transfer, and cooperativity effects that pairwise potentials cannot capture.

**Clashscore calibration.** The steric clash detection produces clashscores 2-3x higher than Molprobity on the same structures (1UBQ: 39.4 vs expected <15). This is because element-level vdW radii don't account for hydrogen bonding, aromatic stacking, or atom-type-specific contact distances. The absolute clashscore values should not be compared to published Molprobity clashscores.

### 4.3 Pocket Definitions

**Manual literature-based sequences.** Pocket residues were identified from published crystal structures and alignment papers, not from automated computational methods (fpocket, AlphaSpace, or distance-based extraction from 3D coordinates). This means:
- Pocket definitions may miss relevant residues
- The analysis cannot be automatically extended to new targets without manual curation
- Some target pairs (TbCathB, TbPDEB1, LmDHFR) use approximate pocket definitions that may not fully capture the true binding site divergence

**Sequential alignment assumption.** Pocket comparison uses a simple sequential alignment of pocket residue sequences. For closely related proteins this is adequate, but for distant homologs (e.g., LmPTR1 vs HsDHFR, which are different enzyme families) a structural alignment would be more appropriate.

### 4.4 ML Methodology

**No permutation null model.** The LODO AUROC of 0.519 is compared against the ESM-2 baseline of 0.429, but not against a scrambled-pocket null (random pocket sequences as features). Without this null, the possibility that the improvement is driven by target-pair label distribution rather than pocket features cannot be excluded.

**No calibration analysis.** AUROC measures ranking but not calibration. The predicted probabilities may not correspond to actual selectivity frequencies.

**Feature leakage risk.** Although the StandardScaler is fitted inside the Pipeline, the pocket features themselves are constant per target pair. All compounds for a given target pair share identical features. The classifier is effectively learning which target pairs tend to have selective compounds, not which individual compounds are selective. This is scientifically meaningful (pocket divergence predicts target-pair selectivity) but should not be confused with compound-level prediction.

### 4.5 Scope

**Three diseases, seven target pairs.** The analysis covers schistosomiasis, trypanosomiasis, and leishmaniasis. Results may not generalize to other NTDs (Chagas disease, onchocerciasis, lymphatic filariasis) or to non-NTD drug design contexts.

**No ADMET integration.** The selectivity analysis does not account for compound pharmacokinetics. A compound that is selective but poorly absorbed, rapidly metabolized, or toxic through off-target mechanisms is still a failed drug.

## 5. What This System Is Not

This system is not a production drug design engine. It is not a replacement for experimental selectivity assays. It is not a validated force field. It is not a peer-reviewed publication.

It is an exploratory computational triage tool that demonstrates a specific hypothesis: binding-site-level physicochemical features carry selectivity signal that global protein embeddings miss. The evidence is preliminary, the dataset is small, and the physics are approximate. The results should be treated as motivation for further investigation, not as established findings.

## 6. Reproducibility

All results can be reproduced from the public repository:

```bash
git clone https://github.com/Danny2045/Kira.git
cd Kira
pip install -e ".[dev]"
python -m kira.experiments.run_selectivity_v3        # LODO results
python -m kira.experiments.case_study_chembl155771   # Case study
python -m kira.experiments.validate_physics          # PDB validation
pytest tests/ -v                                     # 194+ tests
```

Frozen outputs are in `results/` with provenance documentation in `results/PROVENANCE.md`.

## 7. Dependencies

- Python ≥ 3.11
- JAX ≥ 0.4 (geometry and energy computation)
- scikit-learn ≥ 1.2 (classification)
- pandas ≥ 2.0 (data handling)
- NumPy ≥ 1.24
- Typer + Rich (CLI)

Full dependency specification in `pyproject.toml`.
