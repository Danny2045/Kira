# Selectivity v5 Results

## Summary

v5 adds assay-aware ChEMBL expansion and exact-core aggregation. It is a data-substrate upgrade, not yet a new modeling benchmark.

## Expansion

| Quantity | Count |
|---|---:|
| Curated activity rows | 11,703 |
| Candidate rows | 12,091 |
| Exact matched-ratio rows | 1,227 |
| Lower-bound rows | 72 |
| Upper-bound rows | 69 |
| Interval rows | 4 |
| Unmatched comparable rows | 125 |
| Single-side-only rows | 10,594 |

## Exact core

| Pair | Trainable rows | Positive | Negative | Unique scaffolds |
|---|---:|---:|---:|---:|
| LmDHFR | 11 | 9 | 2 | 10 |
| LmPTR1 | 1 | 0 | 1 | 1 |
| SmDHODH | 3 | 0 | 3 | 3 |
| SmHDAC8 | 73 | 1 | 72 | 37 |
| TbCathB | 6 | 1 | 5 | 4 |
| TbPDEB1 | 16 | 0 | 16 | 15 |

## Interpretation

The v5 exact core is cleaner than the raw candidate evidence, but after replicate collapse it is too class-degenerate for a strong six-pair classifier benchmark. The correct next step is a tiered-evidence core before B0/B1/B2 model evaluation.
