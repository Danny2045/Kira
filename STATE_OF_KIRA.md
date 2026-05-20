# State of Kira

**Snapshot anchor commit:** `f12fe3d` — *Merge pull request #19 from Danny2045/feat/lodo-rerun-corrected-embeddings*
**Snapshot date:** 2026-05-17
**Repository:** [github.com/Danny2045/Kira](https://github.com/Danny2045/Kira)

---

## How to use this document

This file is the canonical state-of-repository reference for Kira. It exists so that any reader — human collaborator, future Claude window, peer reviewer — can come up to speed on what the repository currently is, what is reproducible from committed state, what claims are bounded by which evidence tier, and what work is outstanding, **without reconstructing that picture from chat history, PR descriptions, or speculative readings of the code**.

Specific guidance:

- **Read this before starting any substantial work** on Kira. The Outstanding Work Register (§9) is the authoritative list of open scientific findings; if you are about to open a PR, check there first to confirm your work is scoped to a real open finding and is not already partially done.
- **This document complements, not replaces, the README and the governance docs.** The [README](README.md) is the marketing-and-onboarding surface. [`docs/SCIENTIFIC_CONSTITUTION.md`](docs/SCIENTIFIC_CONSTITUTION.md), [`docs/CLAIM_INFLATION_AUDIT.md`](docs/CLAIM_INFLATION_AUDIT.md), and [`docs/DISPOSITION_PLAN.md`](docs/DISPOSITION_PLAN.md) are the normative documents that bound what claims may be made and which modules carry which scientific weight. This file is descriptive: it tells you what *is*, not what *ought to be*.
- **This document is updated periodically, not continuously.** Every line carries an implicit "as of commit `f12fe3d`". Branch lists, missing-file tables, finding-status grids, and per-artifact reproducibility lines can drift between snapshots. When you make a substantive structural change to the repository, update this file as part of the same PR; when you only fix code without changing the structure of what's reproducible, you may leave it for the next snapshot. The "Update protocol" section at the end describes the cadence in more detail.
- **Do not introduce new scientific claims or new metrics in this document.** If a number is not already documented in committed state (a result artifact, a provenance sidecar, a governance doc), it does not belong here. This file is a faithful summary; it is not a place to publish.
- **The Outstanding Work Register pairs every open finding with a "Natural next-PR scope" note.** Use those notes as starting points; refine them when you actually open the PR.

---

## 1. Repository overview

Kira is a Rwanda-rooted, contrast-driven inverse-biology system for retrospective translational analysis, pair-level selectivity benchmarking, curated biological comparison, and residue-level mechanistic hypothesis generation against parasitic targets in schistosomiasis, human African trypanosomiasis, and leishmaniasis. The README opens with the scope-limiting paragraph that bounds these claims; the Scientific Constitution and Disposition Plan are the normative documents that enforce the boundary.

The repository contains:

- an **active engine** under `src/kira/` (50 Python files, 12 sub-packages including 4 deprecated empty namespaces),
- a **frozen historical record** under `archived/scripts/` (22 numbered scripts that produced the original — since-withdrawn — preprint analysis between February and April 2026),
- an **active repair toolkit** under `scripts/` (6 files: 5 Python scripts plus a shell wrapper, all added in PRs #17–19 to repair the UniProt mapping cache and rerun affected models),
- a **test suite** of 314 collected tests under `tests/` (309 pass at the snapshot commit; 5 skip because of Finding 9),
- a **data tree** under `data/` (120 committed files, ~13.6 MB) split across `eval/`, `models/`, `processed/`, `publication/`, `lab_requests/`, `docking/`, `validation/`, `reference/`, `reports/`, `leishmania/`, `trypanosoma/`,
- a **governance and documentation tree** under `docs/` (13 files including the Scientific Constitution, Claim Inflation Audit, Disposition Plan, and module-specific design docs; the preprint and six preprint-era narrative documents were withdrawn from main in PR `chore/withdraw-preprint` — git history preserves them),
- and a **results capture** under `results/` (5 files: case-study output, v3 results, v4 per-pair metrics, summary JSON, and a manual provenance markdown).

The full module-level map is in §5 (Methodology map).

---

## 2. Quick-reference dashboard

| Property | Value |
|---|---|
| Snapshot commit | `f12fe3d` |
| Snapshot date | 2026-05-17 |
| Main branch | `main` |
| Tests collected | 314 |
| Tests passing | 309 |
| Tests skipped | 5 (all in `tests/test_target_manifest.py`; see Finding 9) |
| Tests failing | 0 |
| Committed files (tracked) | 279 |
| Committed size | ~14.4 MB |
| Active scripts | 6 (in `scripts/`) |
| Active engine modules | 50 (under `src/kira/`) |
| Archived scripts | 22 (under `archived/scripts/`) |
| Doc files | 20 (under `docs/`) |
| Provenance sidecars | 4 (3 JSON in `data/models/`, 1 Markdown in `results/`) |
| Open findings | 9 (of 12 total — Findings 1, 3, 5-code-portion are closed; 2 partially closed) |
| New findings surfaced by the audit producing this document | 2 (Findings 11 and 12) |
| Constitution forbidden-term violations in active code | 6 |
| Constitution forbidden-term violations in active docs | 4 |
| Modules carrying explicit `Evidence tier: Tier X` docstring declarations | 0 (Finding 7 systemic gap) |

---

## 3. Git history and branch state

### 3.1 Remote configuration

- `origin` → `https://github.com/Danny2045/Kira.git` (fetch + push)
- Refspec: `+refs/heads/*:refs/remotes/origin/*` (standard)

### 3.2 Last 30 commits on `main`

```
f12fe3d Merge pull request #19 from Danny2045/feat/lodo-rerun-corrected-embeddings
023af79 feat: regenerate ESM-2 cache against corrected sequences; partial LODO rerun
de24547 Merge pull request #18 from Danny2045/feat/structured-esm2-cosine-artifact
f4ac2ff feat: add structured ESM-2 DHODH cosine artifact and route case study through it
52adbe7 Merge pull request #17 from Danny2045/fix/uniprot-mapping-repair
7153225 fix: correct UniProt mappings and add verifiable provenance
bd34d00 Merge pull request #16 from Danny2045/docs/readme-rewrite-constitution-compliant
38fbff8 docs: rewrite README and update pyproject Summary for Constitution compliance
c2f05c7 Merge pull request #15 from Danny2045/docs/restore-scientific-governance
f4324be docs: restore scientific governance documents from commit 3604ea4
95b90a0 Merge pull request #14 from Danny2045/feat/parasite-selectivity-benchmark-repair-report
d6e6d29 docs: record parasite selectivity repair report
b12395d feat: add parasite selectivity benchmark-repair report
65e46ed Merge pull request #13 from Danny2045/docs/kira-project-log-and-pr-template
eece856 docs: add Kira project log and PR template
611eac0 Merge pull request #12 from Danny2045/feat/rwanda-amr-ast-data-return-kit
6a7f939 feat: add Rwanda AMR AST data-return kit
237281c Merge pull request #11 from Danny2045/feat/rwanda-amr-ast-completeness-audit
4f7a153 feat: add Rwanda AMR AST completeness audit
bbfaf75 Merge pull request #10 from Danny2045/feat/physics-auditor-ticket-gate
cb2c532 feat: add Physics Auditor ticket gate
efda7f2 Merge pull request #9 from Danny2045/feat/rwanda-amr-contrast-scout
6d5c35b feat: add Rwanda AMR contrast scout
8fe30b5 Merge pull request #8 from Danny2045/feat/regeneration-contrast-scout
47593dd chore: exclude notebooks from ruff
b35ff20 feat: add regeneration contrast scout
f85779a chore: exclude archived scripts from ruff
205b33a Merge pull request #7 from Danny2045/feat/contrast-core-foundation
35c6f7a feat: add contrast-core foundation
7ef6fbe Merge pull request #6 from Danny2045/feat/v6-lab-campaign-designer
```

### 3.3 Merged PRs on main (in order)

PR #2, PR #4, PR #5, PR #6, PR #7, PR #8, PR #9, PR #10, PR #11, PR #12, PR #13, PR #14, PR #15, PR #16, PR #17, PR #18, PR #19. PRs #1 and #3 are gaps in the GitHub PR numbering — they were closed without merge.

### 3.4 Outstanding branches to clean up

The following branches exist at the snapshot commit and are candidates for cleanup. **No branch is deleted by the PR that introduces this document.** Cleanup is left to a deliberate follow-up so that a human can confirm each deletion.

**Local branches that are merged and could be deleted:**

- `docs/readme-rewrite-constitution-compliant` — merged via PR #16
- `docs/restore-scientific-governance` — merged via PR #15
- `feat/lodo-rerun-corrected-embeddings` — merged via PR #19
- `feat/selectivity-v4-minimal` — merged via PR #4
- `feat/structured-esm2-cosine-artifact` — merged via PR #18
- `fix/uniprot-mapping-repair` — merged via PR #17

**Local branch with divergence (needs decision before deletion):**

- `feat/selectivity-v5-expand-data` — merged via PR #5 but the local copy has diverged from `origin/feat/selectivity-v5-expand-data` by 6 ahead and 14 behind. Either rebase or delete after confirming the local-ahead commits carry no unique work.

**Remote-only branches that appear to be post-merge leftovers:**

- `origin/feat/contrast-core-foundation` — merged via PR #7
- `origin/feat/v6-lab-campaign-designer` — merged via PR #6

**Remote-only branch with unknown status (needs human review before deletion):**

- `origin/fixes/technical-pass-1` — present on origin only, not represented in any merge commit on `main`. Provenance unclear; do not delete without confirmation.

---

## 4. File tree inventory

The repository has **279 committed files** totaling **~14.4 MB** at the snapshot commit.

### 4.1 Top-level directory breakdown

| Top-level dir | Files | Committed bytes |
|---|---:|---:|
| `data/` | 120 | 13,594,170 |
| `src/` | 50 | 538,621 |
| `tests/` | 33 | 154,555 |
| `archived/` | 26 | 435,150 |
| `docs/` | 13 | 107,709 |
| `scripts/` | 6 | 25,758 |
| `results/` | 5 | 18,097 |
| `notebooks/` | 3 | 6,559 |
| `examples/` | 2 | 1,252 |
| Root files (`README.md`, `pyproject.toml`, `pixi.toml`, `environment.yml`, `Dockerfile`, `LICENSE`, `CLAUDE.md`, `.gitignore`, `.github/`) | 11 | — |

### 4.2 `data/` sub-breakdown

| Sub-directory | Files | Committed bytes | Purpose |
|---|---:|---:|---|
| `data/docking/ligands/` | 43 | 87,314 | Per-compound PDBQT ligand files used by archived Script 14 (SmTGR docking) |
| `data/docking/receptors/` | 6 | 5,427,361 | SmTGR and HsTrxR1 receptor structures (PDB + PDBQT) |
| `data/eval/` | 7 | 86,618 | Ground-truth evaluation sets (v1 and v2) |
| `data/lab_requests/` | 4 | 778,503 | v6 lab-campaign tickets (gap-closure and potency-discovery modes) |
| `data/leishmania/` | 3 | 364,100 | Leishmaniasis selectivity, activities, targets |
| `data/models/` | 15 | 228,481 | ESM-2 caches, UniProt mappings, model results, provenance sidecars |
| `data/processed/` | 14 | 5,097,157 | v4/v5/v6 selectivity intermediates and summaries, ChEMBL SMILES cache, target-manifest report |
| `data/publication/` | 14 | 55,169 | Cross-disease tables, discovery candidates, platform reports |
| `data/reference/` | 1 | 1,326 | v5 target-pair reference table |
| `data/reports/` | 4 | 30,834 | Historical v1/v2 candidate, novelty, selectivity reports |
| `data/trypanosoma/` | 6 | 206,836 | Trypanosomiasis selectivity, activities, targets, ranking |
| `data/validation/` | 3 | 1,230,471 | PDB structures used by `validate_physics.py` |

### 4.3 Directory map (filtered)

```
.
├── archived/
│   ├── docs/
│   └── scripts/        (22 numbered historical scripts)
├── data/               (see §4.2)
├── docs/               (governance + module design docs)
├── examples/           (Rwanda AMR AST templates)
├── notebooks/          (3 Jupyter notebooks)
├── results/            (case study, v3, v4 outputs + PROVENANCE.md)
├── scripts/            (6 active repair-toolkit scripts)
├── src/
│   └── kira/
│       ├── amr/        (Rwanda AMR scout, audit, data-return)
│       ├── causality/  (binding-site extraction, divergence, energy decomp, selectivity map)
│       ├── contrast/   (domain-agnostic contrast core primitives)
│       ├── data/       (target manifest pipeline)
│       ├── dock/       (empty package — deprecate)
│       ├── eval/       (empty package — deprecate)
│       ├── experiments/(v3/v4/v5/v6 selectivity pipelines, case study, validate_physics)
│       ├── filter/     (empty package — deprecate)
│       ├── graph/      (empty package — deprecate)
│       ├── physics/    (PDB parser, geometry, energy, clashes, ticket gate)
│       ├── regeneration/(CRISPR/regeneration contrast scout)
│       └── selectivity/(benchmark-repair report)
└── tests/              (33 test files, 314 tests, 309 pass / 5 skip)
```

---

## 5. Methodology map

This section documents every module under `src/kira/`, every active script under `scripts/`, and every historical script under `archived/scripts/`, with one-line purpose, IO surface, and test coverage where applicable.

### 5.1 `src/kira/` — active engine

| Module | One-line purpose | Reads | Writes | Tests |
|---|---|---|---|---|
| `src/kira/__init__.py` | Empty package marker. | — | — | — |
| `src/kira/amr/__init__.py` | Re-exports AMR scout, audit, and data-return helpers. | — | — | `tests/test_amr_*.py` |
| `src/kira/amr/audit.py` | AST completeness audit for Rwanda AMR benchmark-readiness. Pure rule logic. | — | — | `test_amr_audit.py` |
| `src/kira/amr/data_return.py` | CSV IO + Markdown rendering around the AST audit. | input CSV path | output CSV + Markdown | `test_amr_data_return.py` |
| `src/kira/amr/scout.py` | Rwanda AMR contrast scout — maps eight seed contrasts into the contrast-core shape. | — | — | `test_amr_scout.py` |
| `src/kira/causality/__init__.py` | Package marker. Disposition: **move to Hypothesis**. | — | — | — |
| `src/kira/causality/binding_site.py` | Extracts pocket residues from a structure plus ligand coordinates or a centroid. | (called with a parsed `Structure`) | — | `test_causality.py`, `test_divergence.py` |
| `src/kira/causality/divergence.py` | Local-vs-global divergence profile over a curated pocket. | — | — | `test_divergence.py` |
| `src/kira/causality/energy_decomp.py` | Per-residue Lennard-Jones decomposition (JAX). | — | — | `test_energy_decomp.py` |
| `src/kira/causality/selectivity_map.py` | Residue-level heuristic attribution from pocket comparison plus residue-energy deltas. | — | — | `test_selectivity_map.py` |
| `src/kira/chemistry.py` | RDKit wrapper for Lipinski-style drug-likeness properties. | — | — | `test_chemistry.py` |
| `src/kira/cli.py` | Typer CLI: `query`, `evaluate`, `selectivity`, `validate`, `info`. | various per command | optional JSON via `--out` | `test_cli_query.py` |
| `src/kira/contrast/schemas.py` | Dataclasses for `ContrastSpec`, `EvidenceRecord`, `EvidenceStatus`, `ExperimentTicket`, `DataReturnSchema`. | — | — | `test_contrast_schemas.py` |
| `src/kira/contrast/tickets.py` | Experiment-ticket helpers and validators. | — | — | `test_contrast_tickets.py` |
| `src/kira/data/target_manifest.py` | Canonical parasite→human target manifest plus dataset-consistency validation pipeline. | `data/processed/schisto_parasite_targets.csv` (MISSING — Finding 9), `data/trypanosoma/tryp_targets.csv`, `data/leishmania/leish_targets.csv`, plus three selectivity CSVs including the missing `data/processed/kira_selectivity_analysis.csv` (Finding 9) | `data/processed/canonical_target_manifest.csv`, `_report.txt`, `_validation.csv` (the report is committed; the CSV outputs are derived) | `test_target_manifest.py` (5 tests skip pending Finding 9 closure) |
| `src/kira/dock/__init__.py` | Empty package. Disposition: **deprecate**. | — | — | — |
| `src/kira/drugs.py` | `CURATED_SMILES` dict for drugs lacking ChEMBL SMILES entries. | — | — | `test_drugs.py` |
| `src/kira/eval/__init__.py` | Empty package. Disposition: **deprecate**. | — | — | — |
| `src/kira/experiments/__init__.py` | `TARGET_PAIRS` registry of seven parasite/human pairs with pocket sequences. Loads ESM-2 cosine from artifact. | `data/models/esm2_dhodh_cosine.json` | — | (indirect via `test_selectivity_v3.py`) |
| `src/kira/experiments/case_study_chembl155771.py` | Curated case study for CHEMBL155771 vs SmDHODH/HsDHODH. Prints narrative to stdout. | `data/models/esm2_dhodh_cosine.json` | stdout only | — |
| `src/kira/experiments/design_v6_lab_campaign.py` | v6 lab-campaign designer; ranks missing-evidence gaps into structured experimental tickets in benchmark-repair and potency-discovery modes. | tiered-core CSV/JSON | `data/lab_requests/v6_top_assay_tickets.csv`, `v6_gap_closure_campaign.json`, `data/processed/selectivity_v6_campaign_summary*.json` | `test_design_v6_lab_campaign.py` |
| `src/kira/experiments/run_selectivity_v3.py` | v3 benchmark: target-pair pocket features vs ESM-2 baseline. | `data/leishmania/leish_selectivity.csv`, `data/trypanosoma/tryp_selectivity*.csv`, `data/publication/discovery_candidates.csv` | `results/selectivity_v3_results.txt` | `test_selectivity_v3.py` |
| `src/kira/experiments/run_selectivity_v4.py` | v4 LogisticRegression benchmark with `StratifiedGroupKFold` grouped on `(target_pair_id, scaffold)`. | `data/processed/selectivity_v4_rows_primary_trainable.csv` (derived — produced by `selectivity_v4_data.py main()`) | `results/selectivity_v4/per_pair_metrics.csv`, `results/selectivity_v4/summary.json` | indirectly via the v4 data/feature tests |
| `src/kira/experiments/selectivity_features.py` | Pocket-feature vector builder (the v3 feature stack). | — | — | `test_selectivity_v3.py` |
| `src/kira/experiments/selectivity_v4_data.py` | Build the v4 trainable table from selectivity CSVs plus ChEMBL SMILES retrieval. | `data/leishmania/leish_selectivity.csv`, `data/trypanosoma/tryp_selectivity_expanded.csv`, `data/publication/discovery_candidates.csv`, `chembl-webresource-client` (network) | `data/processed/selectivity_v4_rows_{all,primary_candidate,primary_trainable}.csv`, `selectivity_v4_summary.json`, `chembl_smiles_cache.json` | `test_selectivity_v4_data.py` |
| `src/kira/experiments/selectivity_v4_features.py` | RDKit Morgan + descriptor + Murcko-scaffold feature blocks for v4. | — | — | `test_selectivity_v4_features.py` |
| `src/kira/experiments/selectivity_v5_exact_core.py` | Build strict exact-matched-ratio core from v5 candidate expansion rows. | `data/processed/selectivity_v5_candidate_rows.csv` | `selectivity_v5_exact_core_rows.csv`, `_exact_core_summary.json` | `test_selectivity_v5_exact_core.py` |
| `src/kira/experiments/selectivity_v5_expand_data.py` | Assay-aware v5 evidence-substrate expansion for curated parasite-vs-human pairs. | input selectivity/activity CSVs | `data/processed/selectivity_v5_candidate_rows.csv`, `_expansion_summary.json` | `test_selectivity_v5_expand_data.py` |
| `src/kira/experiments/selectivity_v5_tiered_core.py` | Build tiered-evidence v5 core (admits one-sided bounded evidence only when the threshold decision is forced). | `data/processed/selectivity_v5_candidate_rows.csv` | `selectivity_v5_tiered_core_rows.csv`, `_tiered_core_summary.json` | `test_selectivity_v5_tiered_core.py` |
| `src/kira/experiments/validate_physics.py` | Downloads real PDB structures and runs the physics checks for plausibility. | network → `data/validation/*.pdb` | stdout | — |
| `src/kira/filter/__init__.py` | Empty package. Disposition: **deprecate**. | — | — | — |
| `src/kira/graph/__init__.py` | Empty package. Disposition: **deprecate**. | — | — | — |
| `src/kira/physics/__init__.py` | Approximate structure-checking toolkit package marker. | — | — | — |
| `src/kira/physics/checks/__init__.py` | Namespace for structure-check routines. | — | — | — |
| `src/kira/physics/checks/clashes.py` | Steric clash detection (van-der-Waals overlap with tolerance). | — | — | `test_physics_energy.py` |
| `src/kira/physics/config.py` | Default thresholds, weights, and check configuration. (Docstring still contains the forbidden term "composite trust score" — Finding 7.) | — | — | — |
| `src/kira/physics/core/__init__.py` | Namespace for core computational modules. | — | — | — |
| `src/kira/physics/core/energy.py` | Lennard-Jones energy kernel (JAX), total / per-atom / per-residue. (Docstring cross-reference to "the causality module" — Finding 7.) | — | — | `test_physics_energy.py` |
| `src/kira/physics/core/geometry.py` | JAX-accelerated distance / angle / dihedral primitives. | — | — | `test_physics_geometry.py` |
| `src/kira/physics/core/parser.py` | PDB-format ATOM/HETATM parser into a `Structure` dataclass. mmCIF is **not** implemented (docstring corrected in earlier work). | PDB text | — | `test_physics_parser.py` |
| `src/kira/physics/core/topology.py` | Bond-graph inference, atom typing, bonded-pair mask generation. | — | — | `test_physics_topology.py` |
| `src/kira/physics/ticket_gate.py` | Deterministic lexical gate for contrast-core experiment tickets. | — | — | `test_physics_ticket_gate.py` |
| `src/kira/regeneration/__init__.py` | Re-exports regeneration/CRISPR contrast scout helpers. | — | — | `test_regeneration_scout.py` |
| `src/kira/regeneration/scout.py` | Six seed contrasts mapping regenerative biology, CRISPR perturbation, and bioelectric examples into the contrast-core abstraction. | — | — | `test_regeneration_scout.py` |
| `src/kira/scoring.py` | Normalize IC50, phase, publication counts into 0-1 scores. | — | — | `test_scoring.py` |
| `src/kira/selectivity/__init__.py` | Re-exports benchmark-repair report helpers. | — | — | — |
| `src/kira/selectivity/benchmark_repair_report.py` | Reader/renderer over committed v4/v5/v6 artifacts that produces a deterministic Markdown dossier. | `data/lab_requests/*`, `data/processed/selectivity_v{4,5,6}_*.json/csv`, `data/reference/selectivity_v5_target_pairs.csv` | Markdown report (path supplied by caller) | `test_selectivity_benchmark_repair_report.py` |
| `src/kira/targets.py` | Curated `TARGET_ESSENTIALITY` floats and `ORTHOLOGUE_MAP` dicts. | — | — | `test_targets.py` |

### 5.2 `scripts/` — active repair toolkit

These six scripts were added in PRs #17–19 to repair the UniProt mapping cache, regenerate the ESM-2 embeddings against the corrected sequences, and partially rerun Script 20 LODO.

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `scripts/verify_uniprot_mappings.py` | Live UniProt verification for the eleven Kira target proteins. Read-only against the canonical files. | network → UniProt REST | `data/models/uniprot_mappings_provenance.json` |
| `scripts/apply_uniprot_corrections.py` | Idempotently promotes verified mappings into the reference files; aborts if any provenance record is not `status: PASS`. | `data/models/uniprot_mappings_provenance.json` | `data/models/uniprot_ids_v2.json`, `data/models/sequences_v2.json` |
| `scripts/generate_esm2_dhodh_cosine.py` | Recomputes the SmDHODH/HsDHODH ESM-2 cosine, Euclidean L2, and 3-mer Jaccard against the corrected sequences. Records a `supersedes` block with the prior incorrect numbers. | `data/models/sequences_v2.json` (G4VFD7, Q02127) | `data/models/esm2_embeddings_dhodh.npz`, `data/models/esm2_dhodh_cosine.json` |
| `scripts/regenerate_esm2_cache_v2.py` | Regenerates the full v2 ESM-2 cache against the corrected sequences. Recipe matches archived Script 20 exactly (model `esm2_t33_650M_UR50D`, layer 33, mean+std pooling over `[1:len(seq)+1]`, max 1022 tokens, float32 dim 1280). | `data/models/uniprot_mappings_provenance.json` | `data/models/esm2_embeddings_v2.npz`, `data/models/esm2_embeddings_v2_provenance.json` |
| `scripts/run_lodo_script20.py` | Wrapper that runs the archived Script 20 LODO in-process with a virtual `__file__` pointing to the repo root, so that the script's relative paths resolve against `data/`. | (delegates to archived `archived/scripts/20_selectivity_model_v2.py`) | `data/models/model_v2_results.json` (partial — only Trypanosomiasis and Leishmaniasis arms; the Schistosomiasis arm cannot be run because its input CSV is not committed — Finding 9) |
| `scripts/reproduce_v4_v5.sh` | Shell wrapper to rerun the v4 and v5 pipelines from a clean repo. | — | (delegates to module entrypoints) |

### 5.3 `archived/scripts/` — frozen historical record

Per the Disposition Plan, these 22 numbered scripts are "move to Evidence". They are preserved verbatim as the historical record of how the original — since-withdrawn — preprint analysis was produced. They are **not** part of the active reproduction path and must not be edited.

| Script | One-line purpose | Output artifact (when re-runnable) |
|---|---|---|
| `01_explore_primekg.py` | PrimeKG exploration (abandoned path). | — |
| `02_query_chembl.py` | Live ChEMBL query for S. mansoni targets and activities. | `data/processed/schisto_filtered_activities.csv` (**missing — Finding 9**) |
| `03_build_eval_set.py` | Build the ground-truth evaluation set. | `data/eval/evaluation_set_v1.csv`, `eval_v2_*.csv` |
| `04_rank_and_evaluate.py` | First composite-ranking algorithm and evaluation. | `data/reports/kira_v1_*_report.txt` |
| `05_structural_similarity.py` | Structural-similarity addition and reranking. | partial `target_pair_features.csv` |
| `06_whole_organism.py` | Whole-organism activity data and V3 ranking. | various reports |
| `07_admet_and_report.py` | ADMET filtering and the final candidate report. | `data/processed/kira_final_ranking_v1.csv` (**missing — Finding 9**) |
| `08_harden_benchmark.py` | Benchmark hardening and methodological self-criticism. | notes only |
| `09_novelty_filter.py` | Novelty filter. | notes only |
| `10_selectivity_analysis.py` | Selectivity vs human orthologues. | `data/processed/kira_selectivity_analysis.csv` (**missing — Finding 9**) |
| `11_selectivity_rerank.py` | Selectivity-adjusted reranking. | `data/reports/kira_v2_selectivity_adjusted_report.txt` |
| `12_publication_analysis.py` | Publication-ready selectivity analysis. | `data/publication/publication_analysis.txt` |
| `13_supply_chain_and_final.py` | Supply-chain reality check and definitive shortlist. | `data/publication/cross_disease_compounds.csv`, `discovery_candidates.csv` |
| `14_docking_smtgr.py` | SmTGR selectivity molecular docking via AutoDock Vina. | `data/processed/smtgr_docking_results.csv` (**missing — Finding 9**) |
| `15_trypanosoma_platform.py` | Trypanosoma platform extension. | `data/trypanosoma/*.csv`, `data/publication/kira_platform_report.txt` |
| `16_tryp_expanded_selectivity.py` | Expanded T. brucei selectivity. | `data/trypanosoma/tryp_selectivity_expanded.csv` |
| `17_cross_disease_analysis.py` | Cross-disease compound analysis. | `data/publication/cross_disease_selectivity*.csv`, `kira_platform_definitive.txt` |
| `18_leishmania_platform.py` | Third-disease (Leishmania) platform extension. | `data/leishmania/*.csv` |
| `19_selectivity_prediction.py` | First ESM-2 selectivity classifier. | `data/models/sequences.json`, `esm2_embeddings.npz`, `embedding_selectivity_correlation.csv` (all v1; superseded by v2) |
| `20_selectivity_model_v2.py` | Selectivity Model v2 with LODO evaluation. | `data/models/model_v2_results.json` (regenerated in PR #19 but partial — see Finding 2) |
| `21_loto_evaluation.py` | Leave-One-Target-pair-Out evaluation. | `data/models/loto_results.json`, `clean_model_results.json` |
| `22_build_target_manifest.py` | Pre-`src/`-package manifest builder. Disposition: **deprecate** (superseded by `src/kira/data/target_manifest.py`). | — |

### 5.4 `tests/` — coverage map

The repository has 33 test files collecting 314 tests. At the snapshot commit, **309 pass and 5 skip** (no failures).

The 5 skipped tests are all in `tests/test_target_manifest.py::TargetManifestPipelineTest` and skip via:

```python
@pytest.mark.skipif(
    not (... / "data" / "processed" / "schisto_parasite_targets.csv").exists(),
    reason="Requires processed data files (data/processed/*.csv) that are not tracked in git.
            Run the pipeline locally to regenerate.",
)
```

This is the test-suite manifestation of Finding 9.

The other test files map to engine modules as follows:

| Test file | Modules under test |
|---|---|
| `test_amr_audit.py` | `kira.amr.audit` |
| `test_amr_data_return.py` | `kira.amr.data_return` |
| `test_amr_scout.py` | `kira.amr.scout`, `kira.contrast` |
| `test_causality.py` | `kira.causality.binding_site`, `kira.physics.core.parser` |
| `test_chemistry.py` | `kira.chemistry` |
| `test_cli_query.py` | `kira.cli` (query command) |
| `test_contrast_schemas.py` | `kira.contrast.schemas` |
| `test_contrast_tickets.py` | `kira.contrast.tickets` |
| `test_design_v6_lab_campaign.py` | `kira.experiments.design_v6_lab_campaign` |
| `test_divergence.py` | `kira.causality.divergence` |
| `test_drugs.py` | `kira.drugs` |
| `test_energy_decomp.py` | `kira.causality.energy_decomp` |
| `test_physics_energy.py` | `kira.physics.core.energy`, `kira.physics.checks.clashes` |
| `test_physics_geometry.py` | `kira.physics.core.geometry` |
| `test_physics_parser.py` | `kira.physics.core.parser` |
| `test_physics_ticket_gate.py` | `kira.physics.ticket_gate` |
| `test_physics_topology.py` | `kira.physics.core.topology` |
| `test_regeneration_scout.py` | `kira.regeneration` |
| `test_scoring.py` | `kira.scoring` |
| `test_selectivity_benchmark_repair_report.py` | `kira.selectivity.benchmark_repair_report` |
| `test_selectivity_map.py` | `kira.causality.selectivity_map` |
| `test_selectivity_v3.py` | `kira.experiments`, `kira.experiments.selectivity_features` |
| `test_selectivity_v4_data.py` | `kira.experiments.selectivity_v4_data` |
| `test_selectivity_v4_features.py` | `kira.experiments.selectivity_v4_features` |
| `test_selectivity_v5_exact_core.py` | `kira.experiments.selectivity_v5_exact_core` |
| `test_selectivity_v5_expand_data.py` | `kira.experiments.selectivity_v5_expand_data` |
| `test_selectivity_v5_tiered_core.py` | `kira.experiments.selectivity_v5_tiered_core` |
| `test_target_manifest.py` | `kira.data.target_manifest` |
| `test_targets.py` | `kira.targets` |

---

## 6. Expected-but-uncommitted file registry

This is the systematic version of the gap that produced Finding 9. Every hardcoded `data/...` path in `src/`, `scripts/`, `archived/scripts/`, and `tests/` is checked against on-disk presence and against `git log --all` history.

### 6.1 Paths referenced in code but **not** in git history

| Path | On disk | In git history | Referenced from |
|---|---|---|---|
| `data/processed/kira_selectivity_analysis.csv` | No | No | **`src/kira/data/target_manifest.py:64` (active)** + 8 archived scripts |
| `data/processed/schisto_filtered_activities.csv` | No | No | 19 references across archived scripts (produced by `archived/scripts/02_query_chembl.py`) |
| `data/processed/smtgr_docking_results.csv` | No | No | 3 archived scripts + listed in `CLAUDE.md` as a key data file ("43 docking results") |
| `data/processed/kira_final_ranking_v1.csv` | No | No | `archived/scripts/07_admet_and_report.py` only |
| `data/processed/schisto_parasite_targets.csv` | No | No | **`src/kira/data/target_manifest.py:46` (active)** |
| `data/processed/canonical_target_manifest.csv` | No | No | **`src/kira/data/target_manifest.py:21` (default OUTPUT)** — derived, regenerable |
| `data/processed/selectivity_v4_rows_primary_trainable.csv` | No | No | `src/kira/experiments/selectivity_v4_data.py:469` (OUTPUT of v4) — derived, regenerable from committed CSVs plus ChEMBL |

### 6.2 Classification of the gaps

- **Active-code dependency, not derivable from committed data — must be addressed (Finding 9):**
  - `kira_selectivity_analysis.csv` (the Schistosomiasis input consumed by `target_manifest.py` validation and by archived Script 20 / Script 21 LODO).
  - `schisto_parasite_targets.csv` (the Schistosomiasis parasite-target table consumed by `target_manifest.py`).
- **Active-code OUTPUT, regenerable by running the module:**
  - `canonical_target_manifest.csv` — output of `target_manifest.py main()`.
  - `selectivity_v4_rows_primary_trainable.csv` — output of `selectivity_v4_data.py main()`.
- **Archived-only OUTPUT, regenerable only by re-running an archived script with its own missing inputs:**
  - `schisto_filtered_activities.csv` — output of archived Script 02 from a live ChEMBL query.
  - `smtgr_docking_results.csv` — output of archived Script 14 from AutoDock Vina (CLAUDE.md cites this as a key data file).
  - `kira_final_ranking_v1.csv` — output of archived Script 07.

### 6.3 Finding 9 framing — broader than originally reported in PR #19

PR #19 (the LODO rerun) reported two missing files (`kira_selectivity_analysis.csv`, `schisto_filtered_activities.csv`) blocking the Schistosomiasis arm of Script 20. The audit producing this document found that the gap is broader:

- A **third** active-code consumer is missing the same family of data: `src/kira/data/target_manifest.py` references both `schisto_parasite_targets.csv` and `kira_selectivity_analysis.csv`. The absence of the parasite-targets file is the structural cause of the five currently-skipped tests in `tests/test_target_manifest.py`.
- A **fourth**, `smtgr_docking_results.csv`, is referenced in `CLAUDE.md` as a key data file ("43 docking results") but is not committed; only the 43 underlying ligand `.pdbqt` files are.
- The active-code OUTPUTs (`canonical_target_manifest.csv`, `selectivity_v4_rows_primary_trainable.csv`) are derived and not technically missing — but the policy question of whether to commit derived intermediates is itself unresolved (see Open Question 6 in §10).

---

## 7. Reproducibility check

This section traces the lineage of every committed result artifact: what produces it, what inputs it requires, whether it carries a provenance sidecar, and whether it can be regenerated from committed state.

### 7.1 Provenance sidecars

Four provenance sidecars exist at the snapshot commit:

- `data/models/uniprot_mappings_provenance.json` — full per-protein verification record (organism, recommended name, expected vs observed length, full canonical sequence, retrieval timestamp, `status: PASS`/`FAIL`). Produced by `scripts/verify_uniprot_mappings.py`. Added in PR #17.
- `data/models/esm2_embeddings_v2_provenance.json` — generation record for the v2 ESM-2 cache. Documents added/removed accession keys vs v1 with per-key explanation. Records model id, layer, pooling, max sequence length, fair-esm version, code git SHA, generation timestamp. Produced by `scripts/regenerate_esm2_cache_v2.py`. Added in PR #19.
- `data/models/model_v2_results_provenance.json` — Script 20 LODO output provenance. Records `diseases_loaded`, `diseases_missing: ["Schistosomiasis"]`, `missing_input_files`, the superseded values from the pre-PR-17 cache, and an explanation of why three of six pair-level protein features were previously zero vectors. Produced by `scripts/run_lodo_script20.py`. Added in PR #19.
- `results/PROVENANCE.md` — manual provenance note dated 2026-04-07, covering the v3 experiment and the case study. **Stale**: not updated through PRs #15–19.

### 7.2 Per-artifact reproducibility table

| Artifact | Reproducible from committed state | Missing inputs |
|---|---|---|
| `data/models/uniprot_ids_v2.json` | Yes | — |
| `data/models/sequences_v2.json` | Yes | — |
| `data/models/uniprot_mappings_provenance.json` | Yes (requires network — UniProt REST) | — |
| `data/models/esm2_embeddings_v2.npz` | Yes (requires ESM-2 weights download, ~2.5 GB) | — |
| `data/models/esm2_embeddings_dhodh.npz` | Yes (requires ESM-2 weights) | — |
| `data/models/esm2_dhodh_cosine.json` | Yes | — |
| `data/models/model_v2_results.json` | **Partial** — only Trypanosomiasis and Leishmaniasis arms recomputed; Schistosomiasis arm absent | `kira_selectivity_analysis.csv`, `schisto_filtered_activities.csv` (Finding 9) |
| `data/models/sequences.json` (v1) | No — stale, contains C1L5Z2 / Q57UX2 / Q8WQ44 / C1LV40 (the pre-PR-17 wrong accessions) | Superseded by `sequences_v2.json`; kept for historical record |
| `data/models/esm2_embeddings.npz` (v1) | No — generated from the wrong sequences | Superseded by `esm2_embeddings_v2.npz`; kept for historical record |
| `data/models/clean_model_results.json` | No | Produced by archived Script 21 against the missing schisto CSV |
| `data/models/loto_results.json` | No | Produced by `run_selectivity_v3.py`; the v3 path is regenerable but the file itself is not provenance-tagged |
| `data/models/target_pair_features.csv` | No | Produced by archived Script 19 against the now-corrected UniProt mapping; not regenerated |
| `data/models/embedding_selectivity_correlation.csv` | No | Same as above |
| `data/processed/selectivity_v5_candidate_rows.csv` | Yes | Regenerable via `selectivity_v5_expand_data` (deterministic) |
| `data/processed/selectivity_v5_exact_core_*` | Yes | Regenerable via `selectivity_v5_exact_core` |
| `data/processed/selectivity_v5_tiered_core_*` | Yes | Regenerable via `selectivity_v5_tiered_core` |
| `data/processed/selectivity_v6_campaign_summary*.json` | Yes | Regenerable via `design_v6_lab_campaign` |
| `data/processed/canonical_target_manifest_report.txt` | No (skips schisto datasets at runtime) | `schisto_parasite_targets.csv`, `kira_selectivity_analysis.csv` (Finding 9) |
| `data/publication/cross_disease_*.csv` | No | Outputs of archived Script 17 with missing schisto inputs |
| `data/publication/discovery_candidates.csv` | No | Output of archived Script 13 |
| `data/publication/kira_platform_report.txt` | No | Archived Script 15 output (Schistosomiasis + Trypanosomiasis only; pre-Leishmania) |
| `data/publication/kira_platform_definitive.txt` | No | Archived Script 17 output (Schistosomiasis + Trypanosomiasis only, despite the "definitive" label — see Finding 4) |
| `data/publication/table*` | No | Outputs of historical publication-analysis scripts |
| `data/lab_requests/v6_*` | Yes | Regenerable via `design_v6_lab_campaign` |
| `results/selectivity_v3_results.txt` | Yes (in principle — the file itself is not provenance-stamped) | — |
| `results/selectivity_v4/per_pair_metrics.csv`, `summary.json` | Yes (in principle) | — |
| `results/case_study_chembl155771.txt` | Yes — `python -m kira.experiments.case_study_chembl155771` regenerates the stdout content, **but the on-disk file is stale**: it shows cosine `0.9897`, distance `6.09`, Jaccard `0.0326` (the pre-PR-17 numbers from C1L5Z2). A live run now emits `0.989732`, `0.979757`, `0.101362` from the corrected `esm2_dhodh_cosine.json` artifact. | — (the file is a one-off manual stdout capture from 2026-04-07) |

### 7.3 Stale cross-references to the pre-PR-18 ESM-2 numbers

The corrected SmDHODH/HsDHODH cosine artifact is at `data/models/esm2_dhodh_cosine.json`:

- New: `cosine_similarity = 0.989732`, `embedding_distance = 0.979757`, `kmer3_jaccard = 0.101362`
- Superseded: `cosine_similarity = 0.9897`, `embedding_distance = 6.09`, `kmer3_jaccard = 0.0326`

Old values still appear in:

- `results/case_study_chembl155771.txt:20-21`

(Two further occurrences in `docs/technical_report.md` and `docs/kira-ml-analysis.md` were retired when those preprint-era documents were withdrawn from main in PR `chore/withdraw-preprint`.)

Note: `scripts/generate_esm2_dhodh_cosine.py:54-55` stores the old values as `PREVIOUS_*` constants for the `supersedes` block — that is correct and intentional.

---

## 8. Constitution compliance

### 8.1 Forbidden-term sweep

The Scientific Constitution at `docs/SCIENTIFIC_CONSTITUTION.md` defines a controlled vocabulary with a list of "Forbidden Overclaim Terms". The list was applied to every committed `.py` and `.md` file (excluding the audit/constitution docs themselves and the project log).

**Active-code violations** (these need to be fixed by the docstring-rewrite PR scoped under Finding 7):

| File:line | Term | Recommended replacement |
|---|---|---|
| `src/kira/causality/__init__.py:1` | `Causality module` | `Mechanistic hypothesis module` (per Disposition Plan, the module name stays `causality` but the docstring narrows) |
| `src/kira/physics/core/energy.py:5` | `causality module` | Cross-reference phrasing |
| `src/kira/experiments/case_study_chembl155771.py:12` | `FULL causality pipeline` | `curated hypothesis workflow` |
| `src/kira/experiments/case_study_chembl155771.py:231` | `predicts selectivity window` | `suggests a possible selectivity window by pocket divergence` |
| `src/kira/experiments/selectivity_features.py:5` | `binding site predicts selectivity` | benchmark-association phrasing |
| `src/kira/physics/config.py:87` | `composite trust score` | `heuristic score` |

**Active-doc violations:**

None remaining in active docs. (Four previously-listed violations in `docs/technical_report.md` — `composite trust score`, `Eight checks applied`, `run through the full pipeline`, `predicts selectivity window` — were retired when that preprint-era document was withdrawn from main in PR `chore/withdraw-preprint`.)

**Stale on-disk artifact:**

- `results/case_study_chembl155771.txt:44` — `predicts selectivity window`. Resolves automatically when the file is regenerated against the updated case-study source (see Finding 8).

**Allowed mentions (preserved by design):**

- `README.md:3, 15, 29` — all three are explicit Constitutional negations (`It is not a … selectivity predictor`, `Kira does not claim to … implement a complete docking engine …`, `not a validated force field`). These must remain.
- `archived/scripts/04_rank_and_evaluate.py:441` and `archived/scripts/19_selectivity_prediction.py:378` — preserved verbatim per Disposition Plan ("move to Evidence"). Do **not** edit.

### 8.2 Evidence-tier declaration audit

The Constitution requires that every module declare its evidence tier in its docstring, using the template:

```
Evidence tier: Tier X
Scientific status: descriptive / retrospective / benchmark proof-of-concept / mechanistic hypothesis / validated experimental
Not a claim of: …
```

Current state:

- The README carries tier declarations for each major package (`Evidence tier: 1`, `3`, `4`, etc.) — this is the only place tiers appear in plain language.
- **Zero `.py` modules** declare `Evidence tier: Tier X` in their docstring. A repository-wide grep for the literal `Evidence tier:` returns only the README and the Constitution itself.

This is a systemic gap (Finding 7). Closing it requires touching roughly 30 module docstrings; the natural unit is per-package (amr, causality, contrast, data, experiments, physics, regeneration, selectivity), and the work is bigger than a single drive-by edit.

### 8.3 TARGET_PAIRS vs verified UniProt IDs

A cross-check of `src/kira/experiments/__init__.py` `TARGET_PAIRS` against `data/models/uniprot_ids_v2.json` surfaces a drift that PR #17 did not propagate. This is Finding 11.

| Pair key in TARGET_PAIRS | parasite_uniprot in TARGET_PAIRS | Verified in `uniprot_ids_v2.json` | Match |
|---|---|---|---|
| SmDHODH | G4VFD7 | G4VFD7 | ✅ |
| SmHDAC8 | A0A3Q0KTZ8 | A5H660 | ❌ |
| SmTGR | Q86LC0 | (not in v2 set — docking-only) | — |
| TbCathB (TARGET_PAIRS key) / TbCatB (verified key) | Q95PM0 | Q6R7Z5 | ❌ + naming inconsistency |
| TbPDEB1 | Q38F42 | Q8WQX9 | ❌ |
| LmPTR1 | Q01782 | Q01782 | ✅ |
| LmDHFR | P07382 | P07382 | ✅ |

Three of the seven entries (four of the six v2-covered pairs) still hold pre-PR-17 wrong UniProt IDs in the `TARGET_PAIRS` dataclass. Grep shows no active code consumes those fields for behavior — the pocket-sequence strings are the load-bearing data, and the parasite/human ID fields are descriptive metadata only — but the dataclass is a public part of the API and is inconsistent with the verified mapping table.

---

## 9. Outstanding work register

The repository carries **twelve** numbered scientific/structural findings. **Findings 11 and 12 are new** — they were surfaced by the audit that produced this document and have not been opened as PRs yet. **Finding 9 is broader than originally framed in PR #19**: three missing files, not two.

Each finding below pairs current status with a "Natural next-PR scope" line. The scope lines are suggestions for the future PR that closes the finding; they are not commitments and the engineer who picks up the work is expected to refine them.

### Finding 1 — UniProt mapping correctness

**Status:** **Closed** by PR #17.

`data/models/uniprot_ids_v2.json` contains eleven verified accessions; `data/models/uniprot_mappings_provenance.json` records organism, recommended name, expected vs observed lengths, full canonical sequence, and `status: PASS` for each. Verification is reproducible via `scripts/verify_uniprot_mappings.py` against the live UniProt REST API.

Note: the closure does **not** extend to the legacy `TARGET_PAIRS` dataclass — see Finding 11.

### Finding 2 — LODO contamination from the mis-mapped ESM-2 cache

**Status:** **Partially closed** by PR #19.

- `data/models/esm2_embeddings_v2.npz` and `data/models/model_v2_results.json` were regenerated against the corrected sequences.
- The new numbers are `gb_cv_auroc = 0.8934`, `mlp_cv_auroc = 0.7688`, `compound_importance = 0.930`, `protein_importance = 0.070`.
- The Schistosomiasis arm of Script 20 was not regenerated because its input CSVs (`kira_selectivity_analysis.csv`, `schisto_filtered_activities.csv`) are not committed — see Finding 9.

**Natural next-PR scope:** Once Finding 9 is closed, rerun `scripts/run_lodo_script20.py` and either replace `data/models/model_v2_results.json` with the full three-disease results or commit an additional `_schisto.json` slice with a new provenance sidecar.

### Finding 3 — Case-study cosine artifact

**Status:** **Closed** by PR #18.

- `data/models/esm2_dhodh_cosine.json` exists with the corrected numbers (`cosine_similarity = 0.989732`, `embedding_distance = 0.979757`, `kmer3_jaccard = 0.101362`) and a `supersedes` block recording the prior values.
- Both `src/kira/experiments/case_study_chembl155771.py:58` and `src/kira/experiments/__init__.py:30` read the artifact dynamically at import time.

Residual: the on-disk `results/case_study_chembl155771.txt` is stale (it is a 2026-04-07 stdout capture, not a dynamically regenerated file). Re-running the case study regenerates the printed output but does not write back. This is the residual that Finding 8 addresses.

### Finding 4 — Platform-report reconciliation

**Status:** **Open.**

`data/publication/kira_platform_report.txt` (archived Script 15 output — "1,035 unique compounds, Schistosomiasis + Trypanosomiasis") and `data/publication/kira_platform_definitive.txt` (archived Script 17 output — "1,210 unique compounds, two-disease coverage") still coexist in the publication directory with overlapping but inconsistent scope and numbers. Neither file's lineage was retouched after PR #11 (Leishmaniasis added). Both predate the three-disease framing documented in `CLAUDE.md` ("3 diseases, 2,699 compounds").

**Natural next-PR scope:** Pick one canonical platform summary, deprecate the other (or move to `archived/data/publication/`), and add a provenance sidecar that lists committed inputs and disease coverage. If the three-disease summary cannot be regenerated until Finding 9 is closed, document that dependency explicitly.

### Finding 5 — CLI consistency

**Status:** **Code fixed; audit doc not updated.**

`src/kira/cli.py` `query()` (lines 39–86) and `selectivity()` (lines 204+) both implement the corrected behavior. `query` iterates `TARGET_ESSENTIALITY` correctly (treating values as float essentiality scores) and uses `ORTHOLOGUE_MAP.get(tgt, {})`; `selectivity` takes PDB paths and ligand centroid coordinates, exactly as described in CLAIM_INFLATION_AUDIT.md's "honest narrower wording".

The audit document at `docs/CLAIM_INFLATION_AUDIT.md:18-20` still lists these rows as open violations.

**Natural next-PR scope:** A documentation-only update marking the two CLAIM_INFLATION_AUDIT.md rows as resolved, with a reference to the implementing commit/PR.

### Finding 6 — Pf-to-Sm structural alignment

**Status:** **Open.**

`src/kira/causality/binding_site.py:199-202` and `selectivity_map.py:182-185` still default to a sequential alignment with the documented caveat *"should come from a proper structural alignment"*. The same limitation applies to distant homologs (e.g. LmPTR1 vs HsDHFR, which are different enzyme families).

**Natural next-PR scope:** Introduce a hand-curated alignment table for each non-trivially-aligned pair, **or** explicitly restrict the pipeline to pairs where sequential alignment is defensible and surface a runtime error for the rest.

### Finding 7 — Docstring rewrites for Constitution compliance

**Status:** **Open (substantially).**

Many phrasing fixes have already shipped (e.g. the `causality` modules carry "mechanistic hypothesis" framing in many places; the README is fully Constitution-compliant since PR #16). Remaining cleanup is itemized in §8.1 above — six source-file lines and four docs-file lines.

The bigger systemic task: **add `Evidence tier: Tier X` declarations to every module docstring** (currently zero modules have one in code).

**Natural next-PR scope:** A single "constitutional cleanup" PR that fixes the ten forbidden-term lines and adds per-package tier declarations following the README template. Optionally introduce a small `ruff` or pytest gate that enforces the tier-declaration presence going forward.

### Finding 8 — Case-study structured output

**Status:** **Partially closed.**

The cosine/distance/Jaccard numbers now live in a structured JSON (`esm2_dhodh_cosine.json`) read at runtime — the most fragile piece is fixed. The case study itself still writes its full output via stdout `print()` calls only; `results/case_study_chembl155771.txt` is a manual capture from 2026-04-07 and is now stale (see §7.3).

**Natural next-PR scope:** Have `case_study_chembl155771.py` write a structured JSON sidecar (compound metadata, ESM-2 metrics from the artifact, per-position analysis array, pocket-feature numbers) and a regenerated text summary, both with `code_git_sha` and `computed_at_utc` fields.

### Finding 9 — Missing intermediate data (broader than first reported in PR #19)

**Status:** **Open.**

The original PR #19 framing identified two missing files (`kira_selectivity_analysis.csv`, `schisto_filtered_activities.csv`) blocking the Schistosomiasis arm of Script 20 LODO. The audit producing this document found the gap is **broader**:

- **Three** missing files are active-code dependencies, not two:
  - `data/processed/kira_selectivity_analysis.csv` (referenced by active `target_manifest.py:64` plus 8 archived scripts).
  - `data/processed/schisto_filtered_activities.csv` (referenced by 19 archived scripts).
  - `data/processed/schisto_parasite_targets.csv` (referenced by active `target_manifest.py:46`).
- One more missing file is cited in `CLAUDE.md` as a key data file but only the upstream PDBQT ligands are committed:
  - `data/processed/smtgr_docking_results.csv` (output of archived Script 14 docking).
- One historical-only output is also absent (lowest priority):
  - `data/processed/kira_final_ranking_v1.csv` (output of archived Script 07).

The structural cause of the **five currently-skipped tests** in `tests/test_target_manifest.py` is the absence of `schisto_parasite_targets.csv`.

**Natural next-PR scope:** Re-run archived Script 02 (`schisto_filtered_activities.csv`) and Script 10 (`kira_selectivity_analysis.csv`) against the current ChEMBL snapshot, commit both with a provenance sidecar recording the ChEMBL release version, the script version, and a SHA-256 of each output. The schisto target table is small and could be hand-curated from the existing `data/trypanosoma/tryp_targets.csv` and `data/leishmania/leish_targets.csv` patterns. After closing this finding, also re-enable the five skipped tests.

### Finding 10 — Archived-script silent-regeneration fragility

**Status:** **Open.**

The pattern surfaced by PR #19: an archived script (Script 20) reads from a `MODEL_DIR` derived from `os.path.dirname(__file__)`, runs to completion against whatever happens to be in the cache, and silently emits results whose protein-feature block was a zero vector — because three of six pair embeddings were missing — without any complaint or warning. The `run_lodo_script20.py` wrapper repaired the path problem but cannot repair the silent-fallback problem, which lives inside Script 20.

The Disposition Plan blocks rewriting archived scripts (they are preserved verbatim as historical evidence). The mitigation must live outside them.

**Natural next-PR scope:** Add a "preflight" helper in `scripts/` that, before any wrapper runs an archived ML script, verifies (against `uniprot_mappings_provenance.json` and `esm2_embeddings_v2_provenance.json`) that every required protein key resolves to a non-zero embedding, and aborts with a loud error otherwise. Apply the preflight inside `run_lodo_script20.py` and any future wrappers.

### Finding 11 — TARGET_PAIRS UniProt drift (NEW, surfaced by this audit)

**Status:** **Open.** Not previously opened as a PR.

`src/kira/experiments/__init__.py` `TARGET_PAIRS` carries pre-PR-17 (wrong) UniProt IDs for SmHDAC8, TbCathB, and TbPDEB1, and uses an inconsistent key (`TbCathB` in TARGET_PAIRS versus the verified `TbCatB` in `uniprot_ids_v2.json`). The full discrepancy table is in §8.3.

No active code reads those `parasite_uniprot` / `human_uniprot` fields for behavior — the pocket-sequence strings are the load-bearing data — but the dataclass is a public part of the API and is exactly the kind of "two sources of truth" that PR #17's provenance pattern was meant to eliminate. The naming inconsistency is also a small but real footgun.

**Natural next-PR scope:** Replace the literal uniprot strings in `TARGET_PAIRS` with a runtime lookup against `data/models/uniprot_ids_v2.json` (or a generated registry derived from it), reconcile the `TbCathB`/`TbCatB` spelling (decide which form is canonical and propagate), and add a test that asserts equality between the two sources for every key present in both.

### Finding 12 — Project log lapsed past PR #14 (NEW, surfaced by this audit)

**Status:** **Open.** Not previously opened as a PR.

`docs/KIRA_PROJECT_LOG.md` is the durable per-PR record described in its own preamble as "the permanent repo record for meaningful Kira branches… GitHub PRs capture the review event. This log captures why each branch mattered scientifically after the PR has merged." The log currently ends at PR #14. PRs #15, #16, #17, #18, and #19 — five meaningful PRs covering governance restoration, README rewrite, UniProt mapping repair, structured cosine artifact, and the partial LODO rerun — have no entries.

**Natural next-PR scope:** Backfill log entries for PRs #15–19 from each PR's body and the docstring/provenance changes those PRs made, then add a release-discipline rule (a checklist entry in the PR template at `.github/PULL_REQUEST_TEMPLATE.md` if one exists, otherwise create one) so the log is updated as part of every future "meaningful" PR.

### Status summary

| Finding | Status | Origin |
|---|---|---|
| 1. UniProt mappings | **Closed** (PR #17) | Pre-audit |
| 2. LODO contamination | **Partially closed** (PR #19) | Pre-audit |
| 3. Case-study cosine artifact | **Closed** (PR #18) | Pre-audit |
| 4. Platform-report reconciliation | Open | Pre-audit |
| 5. CLI consistency | Code fixed; audit doc still flags | Pre-audit |
| 6. Pf-to-Sm structural alignment | Open | Pre-audit |
| 7. Docstring rewrites + tier declarations | Open (substantially) | Pre-audit |
| 8. Case-study structured output | Partially closed | Pre-audit |
| 9. Missing intermediate data (broader scope) | Open — three missing files, not two | Pre-audit, **scope expanded by this audit** |
| 10. Archived-script silent-regeneration fragility | Open | Pre-audit (opened by PR #19) |
| 11. TARGET_PAIRS UniProt drift | Open | **New — surfaced by this audit** |
| 12. Project log lapsed past PR #14 | Open | **New — surfaced by this audit** |

---

## 10. Open questions and known unknowns

These are the things the audit producing this document surfaced that do not have an answer in committed state. They are decisions for human review, not actions.

1. **`origin/fixes/technical-pass-1`** — present on origin only, not in any merge commit on `main`. Is this an abandoned branch, a fix queued for a future PR, or someone else's work? Needs a human check before deciding whether to delete the remote ref.
2. **Stale local feature branches** — the six post-merge local branches listed in §3.4 and the locally-diverged `feat/selectivity-v5-expand-data`. Should a single dedicated cleanup PR delete them, or is the user comfortable doing it interactively? No deletion happens automatically.
3. **Schistosomiasis input substrate** — is the right path forward to (a) re-query ChEMBL against the current release and version-tag the regenerated snapshot, or (b) freeze a snapshot from a specific ChEMBL release and version-tag it? PR #19's provenance language implies (a) but does not commit to it. (The original preprint, withdrawn from main in PR `chore/withdraw-preprint`, is no longer the reproduction target; the substrate stands on its own provenance.)
4. **Platform-report canonicalization** — is `kira_platform_definitive.txt` actually superseded by `kira_platform_report.txt` (or vice versa), or are both meant to coexist as Tryp-only and Schisto+Tryp-only snapshots? Naming suggests the former; content suggests the latter.
5. **TARGET_PAIRS scope** — does SmTGR (which is not in `uniprot_ids_v2.json`) belong in the verified-mapping set, or is it intentionally docking-only? If the former, `verify_uniprot_mappings.py` should be extended to include it.
6. **Should derived files like `selectivity_v4_rows_primary_trainable.csv` be committed?** They are deterministic outputs of a committed pipeline against committed inputs. Not committing them keeps the repo small; committing them removes "did you remember to run the prep step?" as a footgun. Whichever way: pick one and document it.
7. **Disposition of empty packages** — `src/kira/dock`, `src/kira/eval`, `src/kira/filter`, `src/kira/graph` are marked `deprecate` in the Disposition Plan but still ship as importable packages. Is the plan to leave them as namespace placeholders, or to actually remove them?
8. **`results/PROVENANCE.md` cadence** — the manual provenance file was last updated 2026-04-07 and none of PRs #15–19 touched it. Is it dead, or is it supposed to be maintained alongside the per-artifact JSON sidecars introduced in PRs #17–19?
9. **Tier-declaration enforcement** — once tier strings are added to module docstrings as part of Finding 7 closure, should a `ruff` plugin or simple pytest fixture enforce their presence going forward? Worth deciding before writing all thirty docstrings.
10. **Notebook reproducibility** — `notebooks/` has three Jupyter files (`00_v0_to_v5_evolution.ipynb`, `01_reproduce_v4_benchmark.ipynb`, `02_v5_expansion_exact_core.ipynb`). None were touched by the audit producing this document. Are they current with the codebase, or themselves stale?

---

## 11. Update protocol

This document is intended to be refreshed periodically, not continuously.

- **When to update inline (same PR):** when your PR changes the structure of what is committed or what is reproducible — adding a new active-code dependency on a data file, adding a new provenance sidecar, changing the test count by more than a handful, opening or closing a numbered finding, merging or deleting a branch listed in §3.4.
- **When to defer to the next snapshot:** small documentation tweaks, bug fixes in already-tested code paths, formatting and lint cleanups.
- **When to produce a new snapshot:** at least once per major milestone (e.g. after closing two or more findings), or whenever the divergence between this document and reality becomes uncomfortable enough to slow down a new contributor. Re-run the per-file inventories (§4), the expected-but-uncommitted registry (§6), and the per-artifact reproducibility table (§7). Bump the snapshot commit hash and date at the top.
- **What this document does *not* track:** ephemeral task lists, in-progress branches, chat-level discussion. Those belong in the PR description, the project log (`docs/KIRA_PROJECT_LOG.md`), or in scratchpad notes outside the repository.
- **Authority:** when this document and the code disagree, the **code is authoritative** and this document needs to be updated. When this document and a chat-level summary disagree, **this document is authoritative**.
