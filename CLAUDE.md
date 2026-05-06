# Kira

Open Causal Discovery for Humanitarian Biology — a contrast-driven inverse-biology
research engine. NTD selectivity is its first concrete domain; AMR and regeneration
are demonstrated cross-domain adapters.

## What this is

A pip-installable Python package (`kira-engine`) with three layers:

1. **Contrast core** — `kira.contrast` — the cross-domain abstraction:
   `intervention -> desired context -> control context -> readout ->
    evidence status -> missing measurement -> experiment ticket`.
2. **Domain adapters** — `kira.amr` (Rwanda AMR scout, AST audit, data-return
   kit), `kira.regeneration` (Levin-flavored CRISPR/regeneration scout).
3. **Selectivity benchmarks** — `kira.experiments.run_selectivity_v3` (original
   pair-level pocket benchmark; preserved for history),
   `run_selectivity_v4` (current modeling claim: compound-conditioned,
   scaffold-aware), `selectivity_v5_*` (data-substrate upgrade, exact + tiered
   cores), `design_v6_lab_campaign` (benchmark-repair tickets), and
   `kira.selectivity.benchmark_repair_report` (V14 dossier that renders
   v4/v5/v6 artifacts without rerunning models).

Approximate-physics utilities (`kira.physics.core`) and pocket comparison
(`kira.causality.binding_site`) are shared with the standalone Physics Auditor
project; both repos currently carry duplicated copies.

## Environment

- `conda activate bio-builder`
- Python 3.11, RDKit, pandas, scikit-learn, JAX (required), PyTorch + ESM-2
  (optional, behind the `[ml]` extra). All data in `data/`.

## Current modeling claim (v4)

| Ablation | Macro pair AUROC |
|---|---:|
| `A0_pair_only` | 0.366 |
| `A1_compound_only` | 0.705 |
| `A2_compound_plus_pair` | 0.708 |
| `A3_full_v4` | 0.695 |

Compound chemistry carries most current predictive signal. Pair-level features
(whether ESM-2 cosine or pocket divergence) add little within-pair beyond
compound features. The v3 LODO improvement of pocket features over ESM-2
(mean 0.519 vs 0.429) is a per-pair-prior effect, not transferable mechanism;
this is documented honestly in `docs/KIRA_V0_TO_V5_SCIENTIFIC_SUMMARY.md` and
the V14 `docs/PARASITE_SELECTIVITY_BENCHMARK_REPAIR_REPORT.md`.

## v5 substrate

| Quantity | Count |
|---|---:|
| Curated activity rows | 11,703 |
| Candidate evidence rows | 12,091 |
| Exact matched-ratio candidate rows | 1,227 |
| Trainable exact-core rows | 110 |
| Trainable tiered-core rows | 131 |

Pairs with both classes in tiered core: LmDHFR, SmHDAC8, TbCathB. The other
three pairs (LmPTR1, SmDHODH, TbPDEB1) are class-degenerate and need v6
benchmark-repair tickets to become trainable.

## v6 lab campaign

100 benchmark-repair tickets, of which 99 are missing-parasite-side measurements
(reflecting the public-data substrate, not a methodology choice).

## Selectivity landscape (descriptive, not predictive)

311 cross-species comparisons across 3 diseases and 6 target pairs:
- 61.7% non-selective overall
- PTR1: best target (median 68.2x selectivity, 34 compounds >10x)
- SmDHODH: top compound CHEMBL155771 at 23 nM, 30.8x selectivity, QED 0.89
- SmHDAC8: 90.4% non-selective (selectivity trap)
- ESM-2 cosine 0.9897 for SmDHODH-HsDHODH yet 30.8x selectivity exists at the
  binding-site level — this is the local-vs-global representation observation,
  not a predictive claim.

## Coding rules

- Run code after writing to verify it works.
- Use existing data files; don't re-query ChEMBL unless needed.
- Commit with descriptive messages (see `.github/PULL_REQUEST_TEMPLATE.md`).
- Document limitations honestly. Use the Non-claims section in PRs.
- Always test on held-out data, not just random CV. For pair-grouped data,
  use `StratifiedGroupKFold` with `target_pair_id::murcko_scaffold` groups
  (the v4 pattern), not vanilla `KFold`.
- Scientific claims must be supported by the numbers in committed artifacts.
- Macro-pair AUROC is the right metric for within-pair model quality;
  pooled AUROC can be inflated by class imbalance across pairs.

## Cross-repo notes

- The JAX physics modules (`kira.physics.core`) and the binding-site
  comparison (`kira.causality.binding_site`) are mirrored in
  `github.com/Danny2045/physics-auditor`. Bug fixes need to be applied
  to both repos until a shared library is extracted.
- The Physics Auditor `ticket_gate` (V10) is a contrast-ticket measurability
  linter, not actual physics; the standalone Physics Auditor project (separate
  repo) is the JAX/LJ structure validator. Names will be disambiguated.
