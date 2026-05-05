# Parasite Selectivity Benchmark-Repair Report

## Why This Report Exists

Summarize the existing v4 modeling claim, the v5 public-data evidence substrate, and the v6 benchmark-repair ticket campaign so Kira's original empirical work is easier to inspect and act on.

This dossier returns to Kira's strongest empirical public-data substrate: parasite-vs-human selectivity evidence from the existing v4/v5/v6 artifacts.

## Source Artifacts

- `results/selectivity_v4/summary.json`
- `results/selectivity_v4/per_pair_metrics.csv`
- `data/processed/selectivity_v4_summary.json`
- `data/processed/selectivity_v5_expansion_summary.json`
- `data/processed/selectivity_v5_exact_core_summary.json`
- `data/processed/selectivity_v5_tiered_core_summary.json`
- `data/processed/selectivity_v6_campaign_summary.json`
- `data/lab_requests/v6_top_assay_tickets.csv`
- `data/lab_requests/v6_gap_closure_campaign.json`
- `data/reference/selectivity_v5_target_pairs.csv`
- `data/processed/canonical_target_manifest_report.txt`

## How v4, v5, and v6 Connect

v4 is the current modeling claim. v5 expands and cleans the assay-aware public-data evidence substrate. v6 turns gaps in that substrate into benchmark-repair assay tickets.

## v4 Modeling Claim Summary

- v4 remains the current modeling claim.
- Trainable rows: 235.
- Feature count: 379.
- Cross-validation splits: 5.
- Morgan fingerprint bits: 256.
- Interpretation: Compound chemistry carries most current predictive signal: compound-only features outperform pair-only features by macro pair AUROC; compound + pair is the best classifier by macro pair AUROC in the committed v4 summary.

| Ablation | Macro pair AUROC |
|---|---:|
| `A0_pair_only` | 0.366480 |
| `A1_compound_only` | 0.704891 |
| `A2_compound_plus_pair` | 0.708342 |
| `A2b_pair_plus_side` | 0.697192 |
| `A2c_pair_plus_compat` | 0.686001 |
| `A3_full_v4` | 0.695064 |

## v5 Evidence Substrate Summary

| Quantity | Count |
|---|---:|
| Candidate evidence rows | 12091 |
| Curated activity rows | 11703 |
| Exact matched-ratio candidate rows | 1227 |
| Exact-core rows | 114 |
| Trainable exact-core rows | 110 |
| Tiered-core rows | 135 |
| Trainable tiered-core rows | 131 |

Tiered pairs with both classes: LmDHFR, SmHDAC8, TbCathB.

| Pair | Candidate rows |
|---|---:|
| `LmDHFR` | 851 |
| `LmPTR1` | 905 |
| `SmDHODH` | 2631 |
| `SmHDAC8` | 4771 |
| `TbCathB` | 1183 |
| `TbPDEB1` | 1750 |

## v6 Lab-Campaign Summary

- Campaign mode: `benchmark_repair`.
- Benchmark-repair tickets: 100.

| Pair | Tickets |
|---|---:|
| `LmDHFR` | 5 |
| `LmPTR1` | 25 |
| `SmDHODH` | 25 |
| `SmHDAC8` | 15 |
| `TbCathB` | 5 |
| `TbPDEB1` | 25 |

| Missing side | Tickets |
|---|---:|
| `parasite` | 99 |
| `human` | 1 |

## Benchmark Repair Interpretation

| Pair | Evidence state | Exact both classes | Tiered both classes | Missing side | Tickets | Repair priority |
|---|---|---|---|---|---:|---|
| `LmDHFR` | `both_classes_present` | yes | yes | parasite=5 | 5 | `closest-to-ready` |
| `LmPTR1` | `class_degenerate_zero_positive` | no | no | parasite=24, human=1 | 25 | `highest` |
| `SmDHODH` | `class_degenerate_zero_positive` | no | no | parasite=25 | 25 | `highest` |
| `SmHDAC8` | `both_classes_extreme_imbalance` | yes | yes | parasite=15 | 15 | `high` |
| `TbCathB` | `both_classes_present` | yes | yes | parasite=5 | 5 | `closest-to-ready` |
| `TbPDEB1` | `class_degenerate_zero_positive` | no | no | parasite=25 | 25 | `highest` |

### LmDHFR

- Parasite target: Bifunctional dihydrofolate reductase-thymidylate synthase.
- Human comparator: Dihydrofolate reductase (human).
- Exact core: 11 trainable rows (9 positive, 2 negative).
- Tiered core: 11 trainable rows (9 positive, 2 negative).
- Repair note: LmDHFR is closest to benchmark-ready because both classes are present (9 positive, 2 negative); the 5 tickets (parasite=5) mainly improve coverage and matched-evidence robustness.

### LmPTR1

- Parasite target: Pteridine reductase 1.
- Human comparator: Dihydrofolate reductase (human).
- Exact core: 1 trainable row (0 positive, 1 negative).
- Tiered core: 1 trainable row (0 positive, 1 negative).
- Repair note: LmPTR1 has 0 positive and 1 negative trainable tiered rows, so benchmark repair should prioritize 25 missing-side tickets (parasite=24, human=1) that can create matched comparator evidence.

### SmDHODH

- Parasite target: Dihydroorotate dehydrogenase (quinone), mitochondrial.
- Human comparator: Dihydroorotate dehydrogenase (human).
- Exact core: 3 trainable rows (0 positive, 3 negative).
- Tiered core: 3 trainable rows (0 positive, 3 negative).
- Repair note: SmDHODH has 0 positive and 3 negative trainable tiered rows, so benchmark repair should prioritize 25 missing-side tickets (parasite=25) that can create matched comparator evidence.

### SmHDAC8

- Parasite target: Histone deacetylase 8.
- Human comparator: Histone deacetylase 8 (human).
- Exact core: 73 trainable rows (1 positive, 72 negative).
- Tiered core: 87 trainable rows (1 positive, 86 negative).
- Repair note: SmHDAC8 already has both classes but remains highly imbalanced (1 positive, 86 negative); 15 tickets (parasite=15) should add comparator evidence and reduce benchmark skew.

### TbCathB

- Parasite target: Cathepsin B-like cysteine protease.
- Human comparator: Cathepsin B (human).
- Exact core: 6 trainable rows (1 positive, 5 negative).
- Tiered core: 12 trainable rows (7 positive, 5 negative).
- Repair note: TbCathB is closest to benchmark-ready because both classes are present (7 positive, 5 negative); the 5 tickets (parasite=5) mainly improve coverage and matched-evidence robustness.

### TbPDEB1

- Parasite target: Class 1 phosphodiesterase PDEB1.
- Human comparator: Phosphodiesterase 4B (human).
- Exact core: 16 trainable rows (0 positive, 16 negative).
- Tiered core: 17 trainable rows (0 positive, 17 negative).
- Repair note: TbPDEB1 has 0 positive and 17 negative trainable tiered rows, so benchmark repair should prioritize 25 missing-side tickets (parasite=25) that can create matched comparator evidence.

## Next Actions

Closest to benchmark-ready: LmDHFR, TbCathB.
Need parasite-side measurements: LmDHFR, LmPTR1, SmDHODH, SmHDAC8, TbCathB, TbPDEB1.
Need human-side comparator measurements: LmPTR1.

Returned data that would repair the benchmark:

- Return the missing-side comparator assay for the ticketed compound and target pair.
- Use comparable potency units and relation fields so a human-divided-by-parasite ratio can be reconstructed.
- Include returned-data fields: `compound_key`, `pair_id`, `measured_side`, `target_chembl_id`, `assay_type`, `standard_type`, `standard_relation`, `standard_value`, `standard_units`, `replicate_count`, `data_validity_comment`, `assay_chembl_id_or_external_id`, `activity_chembl_id_or_external_id`, `notes`.
- Flag data-validity concerns so repaired rows can be audited before joining any benchmark core.

## Contrast-Driven Inverse-Biology Direction

The report keeps the inverse-biology loop empirical: start from a parasite-vs-human contrast, label the evidence state, expose benchmark weakness, issue missing-measurement tickets, and wait for auditable returned data before making stronger claims.

## Non-Claims

- No drug-discovery claim: this dossier summarizes evidence substrate and repair tickets only.
- No wet-lab validation claim: v6 tickets are measurement requests, not completed experiments.
- No fresh model-performance claim: the report repeats existing v4/v5/v6 artifacts and does not rerun a benchmark.
- No clinical or therapeutic-success claim.
- No claim that parasite disease problems are broadly resolved.

## Validation Commands

```bash
ruff check src/kira/selectivity tests/test_selectivity_benchmark_repair_report.py
pytest -q tests/test_selectivity_benchmark_repair_report.py
ruff check .
pytest -q
```
