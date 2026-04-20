# Kira v0 to v5 Scientific Summary

## Current identity

Kira is a computational repository for retrospective translational analysis, parasite-vs-human selectivity benchmarking, assay-aware data expansion, and mechanistic hypothesis generation under sparse public medicinal-chemistry data.

It is not yet a production predictor, universal drug-discovery engine, full binding-physics engine, or wet-lab-validated therapeutic platform.

## Version evolution

| Stage | Main contribution | Scientific interpretation |
|---|---|---|
| v0-v1 | Exploratory translational and repurposing analyses | Useful historical context, not a modern benchmark. |
| v2 | More organized selectivity and target-comparison work | Still limited by sparse data and incomplete leakage control. |
| v3 | Retrospective target-pair selectivity benchmark | Too narrow; representation could collapse to target-pair constants. |
| v4 | Compound-conditioned, scaffold-aware benchmark | Real methodological repair. |
| v5 expansion | Assay-aware ChEMBL expansion with corrected target ontology | Expanded evidence substrate under strict rules. |
| v5 exact core | Replicate-collapsed exact-ratio core | Clean but class-degenerate after aggregation. |

## v4 result

v4 uses a 379-feature representation and scaffold-aware cross-validation with groups defined as `target_pair_id::murcko_scaffold`.

| Ablation | Pooled AUROC | Macro pair AUROC | Macro pair Spearman |
|---|---:|---:|---:|
| A0_pair_only | 0.794 | 0.366 | -0.206 |
| A1_compound_only | 0.835 | 0.705 | 0.343 |
| A2_compound_plus_pair | 0.900 | 0.708 | 0.357 |
| A2b_pair_plus_side | 0.898 | 0.697 | 0.348 |
| A2c_pair_plus_compat | 0.896 | 0.686 | 0.345 |
| A3_full_v4 | 0.896 | 0.695 | 0.361 |

The honest interpretation is that compound chemistry currently dominates the predictive signal, while the present coarse compatibility block has not yet shown clear within-pair classification gain.

## v5 result

| Quantity | Count |
|---|---:|
| Curated activity rows | 11,703 |
| Candidate evidence rows | 12,091 |
| Exact matched-ratio candidate rows | 1,227 |
| Exact-core rows | 114 |
| Trainable exact-core rows | 110 |
| Conflicting exact-core rows | 4 |

v5 improved the evidence substrate, but the strict exact core is class-degenerate. The next benchmark should not overclaim v5 modeling performance until B0/B1/B2 is run on a better exact or tiered core.

## Next scientific step

Build `selectivity_v5_tiered_core.py` to aggregate exact plus logically decidable bounded evidence, then run:

- `B0_pair_only`
- `B1_compound_only`
- `B2_compound_plus_pair`

with scaffold-aware grouping.
