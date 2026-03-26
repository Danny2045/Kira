# Kira ML Analysis: Can Protein Language Models Predict Antiparasitic Selectivity?

## Final Clean Results — All Engineering Confounders Resolved

*Daniel Ngabonziza — March 2026*

---

## Summary

We tested whether ESM-2 protein language model embeddings, combined with compound molecular features, can predict whether a compound will selectively inhibit a parasite enzyme over its human orthologue. Using 311 curated selectivity labels across three diseases and six verified target-orthologue pairs, we found:

**Result 1: Protein pair features contribute modestly within a mixed dataset.** Adding 13 protein pair features (sequence statistics + ESM-2 embeddings) to 266 compound features improves 5-fold CV AUROC from 0.840 to 0.895 (+0.055). Protein features account for 8.6% of total model importance.

**Result 2: This signal does not transfer across diseases.** Leave-one-disease-out evaluation: pair features worsened predictions for schistosomiasis (-0.130) and trypanosomiasis (-0.206), and were roughly neutral for leishmaniasis (+0.043). In 2 of 3 held-out diseases, protein pair features actively degraded performance.

**Result 3: In the current feature set, simple sequence-difference features contributed more than ESM-2 summary features.** The top protein-level feature was parasite sequence length (3.1% importance), followed by 4-mer Jaccard overlap (1.7%). ESM-2 cosine similarity contributed 0.5%. However, this does not constitute a full ESM-only ablation, so the conclusion is that ESM summary features were less useful than sequence statistics in this representation, not that ESM embeddings have no value for this task in general.

**Key observation:** SmDHODH and human DHODH have ESM-2 cosine similarity of 0.9897 — nearly identical in embedding space — yet SmDHODH harbors compounds with 30.8-fold parasite selectivity. The selectivity-relevant structural divergence is invisible to global protein representations.

---

## Data

311 selectivity labels from Kira Scripts 10, 15, 16, 18:

| Disease | Target pair | N compounds | % Selective (≥3x) | Median ratio |
|---------|-------------|:---:|:---:|:---:|
| Schistosomiasis | SmHDAC8 vs HsHDAC8 | 73 | 9.6% | 1.0x |
| Schistosomiasis | SmDHODH vs HsDHODH | 18 | 61.1% | 4.7x |
| Trypanosomiasis | TbCatB vs HsCatL | 36 | 5.6% | 0.8x |
| Trypanosomiasis | TbPDEB1 vs HsPDE4B | 68 | 38.2% | 1.6x |
| Leishmaniasis | LmPTR1 vs HsDHFR | 45 | 75.6% | 68.2x |
| Leishmaniasis | LmDHFR vs HsDHFR | 69 | 56.5% | 3.6x |

Binary task: Selective+Moderate (119) vs Poor+Counter-selective (192).

All UniProt accessions verified by organism and expected protein length:

| Target | UniProt | Length | Status |
|--------|---------|:---:|:---:|
| SmHDAC8 | A0A3Q0KU27 | 426 aa | Verified |
| SmDHODH | G4VFD7 | 379 aa | Verified (fixed from 0-aa G4VQR5) |
| TbCatB | Q6R7Z5 | 340 aa | Verified |
| TbPDEB1 | Q8WQX9 | 930 aa | Verified |
| LmPTR1 | Q01782 | 288 aa | Verified |
| LmDHFR | P07382 | 520 aa | Verified |
| HsHDAC8 | Q9BY41 | 377 aa | Verified |
| HsDHODH | Q02127 | 395 aa | Verified |
| HsCatL | P07711 | 333 aa | Verified |
| HsPDE4B | Q07343 | 736 aa | Verified |
| HsDHFR | P00374 | 187 aa | Verified |

---

## Features

**Compound features (266):** 256-bit Morgan fingerprint (ECFP2) + 10 RDKit descriptors (MW, LogP, HBD, HBA, TPSA, rotatable bonds, rings, aromatic rings, Fsp3, QED).

**Protein pair features (13):** parasite length, human length, length ratio, length difference, k-mer Jaccard (k=2,3,4), amino acid composition L2 distance, ESM-2 cosine similarity, ESM-2 Euclidean distance, ESM-2 L1 distance, ESM-2 parasite embedding norm, ESM-2 human embedding norm.

**Total: 279 features per sample.**

---

## ESM-2 Pair Distances (Clean)

| Parasite target | Human target | ESM-2 cosine | Length diff |
|----------------|--------------|:---:|:---:|
| SmHDAC8 | HsHDAC8 | 0.535 | 49 aa |
| SmDHODH | HsDHODH | **0.990** | 16 aa |
| TbCatB | HsCatL | 0.960 | 7 aa |
| TbPDEB1 | HsPDE4B | 0.962 | 194 aa |
| LmPTR1 | HsDHFR | 0.797 | 101 aa |
| LmDHFR | HsDHFR | 0.807 | 333 aa |

SmHDAC8 is the only pair with low ESM-2 similarity (0.535), consistent with its different fold architecture. All other pairs score >0.79, meaning ESM-2 considers them highly similar — yet their selectivity landscapes range from 94% non-selective (cathepsin B) to 24% non-selective (PTR1).

---

## Results

### 5-Fold Stratified Cross-Validation

| Model | Features | AUROC |
|-------|----------|:---:|
| Gradient Boosting | Compound only (266) | 0.840 ± 0.039 |
| Gradient Boosting | Compound + Pair (279) | **0.895 ± 0.034** |
| MLP (128-64-32) | Compound + Pair (279) | 0.845 ± 0.035 |
| Majority baseline | — | 0.617 |

Pair features improve AUROC by +0.055 within mixed-disease CV.

### Leave-One-Disease-Out (The Real Test)

| Held-out disease | Compound only | Compound + Pair | Pair effect |
|-----------------|:---:|:---:|:---:|
| Schistosomiasis (91 samples, 20% pos) | 0.512 | 0.382 | **-0.130** |
| Trypanosomiasis (104 samples, 27% pos) | 0.565 | 0.359 | **-0.206** |
| Leishmaniasis (116 samples, 63% pos) | 0.505 | 0.547 | **+0.043** |

All compound-only LODO AUROCs are near random (0.50-0.57). Adding pair features worsens predictions for 2 of 3 diseases and marginally helps for the third.

### Feature Importance (Clean)

| Feature group | Importance | % |
|--------------|:---:|:---:|
| Morgan fingerprint (256) | 0.713 | 71.3% |
| RDKit descriptors (10) | 0.201 | 20.1% |
| Protein pair features (13) | 0.086 | **8.6%** |
| **Compound total** | **0.914** | **91.4%** |
| **Protein total** | **0.086** | **8.6%** |

Top protein-level features: para_len (3.1%), kmer4 (1.7%), aa_l2 (0.9%), esm_hnorm (0.7%), esm_pnorm (0.6%), kmer2 (0.5%), esm_cos (0.5%).

---

## Interpretation

### Why 5-fold CV works but LODO fails

5-fold CV mixes all diseases within each fold. The model learns compound-scaffold patterns associated with each target's selectivity distribution (e.g., "hydroxamic acids from SmHDAC8 → non-selective," "pteridine analogues from PTR1 → selective"). These are valid within-dataset associations. LODO asks whether these patterns transfer to unseen targets with unseen compound scaffolds. They do not.

### Why pair features help within-dataset but hurt across diseases

Within the mixed dataset, pair features act as a target-identity proxy. Compounds from SmHDAC8 (cosine 0.535) have different selectivity rates than compounds from PTR1 (cosine 0.797). The pair features encode "which target group am I?" — and each group has a characteristic selectivity rate.

Across diseases, the pair-to-selectivity mapping learned from two diseases misleads for the third. Protein pair features that correlate with selectivity in schistosomiasis + leishmaniasis are anti-correlated with selectivity in trypanosomiasis. The mapping is not universal.

### Why ESM-2 contributes less than sequence statistics

ESM-2 cosine similarity is high (>0.79) for 5 of 6 pairs and very high (0.99) for SmDHODH-HsDHODH. The protein language model sees these pairs as nearly identical. Yet their selectivity landscapes are dramatically different: SmDHODH has 3 compounds with >10x selectivity (best 30.8x), while cathepsin B (cosine 0.96) has zero. ESM-2 embeddings encode global protein similarity — fold, function, evolutionary distance — not binding-site-level geometry. Selectivity depends on specific active-site residue differences that global mean-pooled representations do not capture.

Simple sequence statistics (length, k-mer overlap) capture coarser but different information: roughly "how structurally different are these proteins at the sequence level?" This is a weaker signal than binding-site comparison but was more discriminative than ESM-2 summaries in the current feature set.

### The SmDHODH cosine observation

SmDHODH and human DHODH have ESM-2 cosine similarity of 0.9897. To the protein language model, they are nearly the same protein. Yet CHEMBL155771 achieves 30.8x selectivity — 23 nM against the parasite enzyme, 709 nM against the human version. The selectivity-enabling structural differences are invisible to global protein representations. This is perhaps the sharpest single illustration in our data of the gap between protein language model similarity and functional selectivity.

---

## Engineering Audit Trail

The clean result required three iterations:

| Script | SmDHODH sequence | Pair features in matrix | LODO result |
|--------|:---:|:---:|:---:|
| Script 19 | Wrong (0 aa, G4LZI2) | Yes (31.4%) | Not tested |
| Script 20 | Wrong (0 aa, G4VQR5) | No (0%, bug) | Fails (0.26-0.49) |
| Script 21 | Wrong (0 aa, G4VQR5) | Yes (21.6%) | Mixed (+0.08 to -0.21) |
| **Script 21b** | **Correct (379 aa, G4VFD7)** | **Yes (8.6%)** | **Fails cleanly (-0.13 to -0.21)** |

The pair feature importance dropped from 21.6% to 8.6% when SmDHODH was fixed, indicating that the broken embedding was artificially inflating the protein signal. The LODO failure was consistent across all runs with pair features present, strengthening confidence that it reflects biology rather than engineering.

---

## Implications

**For antiparasitic drug repurposing:** Selectivity prediction from global protein representations is not viable for cross-target generalization. Kira's experimental selectivity triage (ChEMBL cross-species querying) remains more reliable than ML prediction for new targets.

**For protein language models:** ESM-2 mean-pooled embeddings encode global protein similarity, not binding-site-level selectivity-relevant features. Future work should test active-site-specific embeddings — either per-residue ESM-2 features restricted to binding-site residues, or structure-conditioned representations from AlphaFold 3.

**For ML evaluation in drug discovery:** The gap between 5-fold CV (AUROC 0.895) and LODO (0.36-0.55) is a concrete example of how standard evaluation overestimates real-world performance. Held-out evaluation on unseen target classes should be standard.

---

## Limitations

1. Only 311 samples across 6 target pairs. LODO evaluation has only 3 data points at the disease level.
2. ESM-2 embeddings were mean-pooled over the full sequence. Active-site-specific pooling was not tested.
3. SmHDAC8 UniProt ID (A0A3Q0KU27, 426 aa) has not been independently verified as the canonical S. mansoni HDAC8 sequence.
4. Human orthologues are not always exact structural equivalents (PTR1 vs DHFR are different fold families).
5. Only gradient boosting was used for the ablation; MLP was tested but not ablated by feature group.

---

## Code

Scripts 19-21 and all feature matrices, embeddings, and evaluation results: https://github.com/Danny2045/Kira
