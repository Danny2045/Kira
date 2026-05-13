# Kira Claim Inflation Audit

This document records places where repository language overstates what the current code actually does.

## README

| File path | Quoted phrase | Why it is overstated | Honest narrower wording | Severity |
|---|---|---|---|---|
| `README.md` | `Open Causal Discovery` | The active package does not perform causal identification or intervention-based inference. The active `causality` package computes pocket divergence and residue attributions from curated inputs and approximate energies. | `Open computational selectivity analysis and mechanistic hypothesis generation` | High |
| `README.md` | `A unified computational engine for drug selectivity analysis, physics validation, and mechanistic explanation` | The active package is not a fully unified engine for all repository claims. Major functionality, especially docking and much of the historical translational pipeline, lives only in archived scripts. | `A package combining retrospective selectivity analysis, approximate structure checks, and residue-level hypothesis tools` | Medium |
| `README.md` | `Binding-site divergence predicts cross-disease drug selectivity` | `run_selectivity_v3.py` uses one constant pocket feature vector per target pair for all compounds in that pair. The code supports a pair-level correlation/proof-of-concept benchmark, not a strong compound-level predictive claim. | `Curated binding-site divergence features correlate with target-pair selectivity patterns in a small cross-disease benchmark` | High |
| `README.md` | `selectivity prediction` | The active benchmark is not true compound-specific prediction in the usual sense, because all compounds from one target pair share the same 15 pocket features. | `target-pair-level selectivity benchmark` | High |
| `README.md` | `Physics validation engine` | The active validator computes clashes, LJ totals, and backbone dihedral counts. It does not implement the full documented structure-validation stack. | `approximate structure-checking utilities` | Medium |
| `README.md` | `Validates AI-generated protein structures ... with 8 physics checks` | The active `kira validate` path reports only steric clashes, Lennard-Jones totals, and backbone dihedrals. The README’s 8-check list is not implemented in the active validator. | `checks structures with a limited active set of geometry and steric heuristics` | High |
| `README.md` | `Produces a Trust Report with accept/relax/discard recommendation` | The active code does not produce a multi-check trust report. In `cli.validate()`, `global_score` is effectively the clash subscore. | `produces a clash-dominated heuristic report with accept/short_md/discard thresholds` | High |
| `README.md` | `Mechanistic explanation layer — Explains WHY selectivity exists` | The active explanation path uses curated pockets plus protein-only LJ decomposition, not ligand-conditioned binding calculations or experimental mechanism validation. | `Mechanistic hypothesis layer that proposes residue-level explanations for selectivity` | High |
| `README.md` | `attributes selectivity to specific pocket differences` | The active attribution is built from pocket identity plus residue energy deltas from protein-only decomposition, not true protein-ligand binding attribution. | `highlights pocket differences associated with the model’s selectivity heuristic` | High |
| `README.md` | `kira validate structure.pdb      # Physics validation with Trust Report` | CLI help overpromises relative to actual outputs. No full trust report with the documented eight checks is produced. | `kira validate structure.pdb      # Approximate clash/LJ/backbone summary` | High |
| `README.md` | `kira query --target SmTGR` | The active `query()` implementation is internally inconsistent and likely broken because it treats float essentiality scores as dicts and asks for nonexistent orthologue keys. | `query command currently needs repair before being presented as stable CLI functionality` | High |
| `README.md` | `kira selectivity --parasite SmDHODH --human HsDHODH` | The active CLI expects PDB file paths, not target names, and requires a ligand centroid, not a compound or docked complex. | `kira selectivity --parasite parasite.pdb --human human.pdb --ligand-x ... --ligand-y ... --ligand-z ...` | High |
| `README.md` | `checks/           # 8 physics checks + composite scorer` | The active package does not expose eight implemented checks in the validator, and the composite score is not a weighted eight-term score in practice. | `checks/           # steric clash logic and related validation scaffolding` | High |
| `README.md` | `selectivity_map.py # Selectivity attribution` | The file does not attribute ligand-conditioned binding energetics. It attributes residue-energy differences from protein-only decompositions merged with pocket divergence. | `selectivity_map.py # residue-level heuristic attribution` | High |

## Docs

| File path | Quoted phrase | Why it is overstated | Honest narrower wording | Severity |
|---|---|---|---|---|
| `docs/technical_report.md` | `Physics validation engine ... checks them for physical validity` | The active code performs limited heuristic checks, not a broad physical-validity assessment. | `approximate structure-quality checker` | Medium |
| `docs/technical_report.md` | `Produces a composite trust score with an accept/relax/discard recommendation` | In the active validator, `composite = clash_result.subscore`; there is no active weighted combination over the documented eight checks. | `produces a clash-based heuristic score and recommendation thresholds` | High |
| `docs/technical_report.md` | `Mechanistic explanation layer ... addresses the question "why is this compound selective?"` | The active code never loads a ligand structure in the selectivity path and does not compute protein-ligand interaction energy. | `addresses a weaker question: which curated pocket differences and residue-energy heuristics are associated with observed selectivity?` | High |
| `docs/technical_report.md` | `mean leave-one-disease-out (LODO) AUROC of 0.519 for selectivity prediction` | The benchmark is not clean compound-level prediction; all compounds from a pair share the same pocket vector. | `mean LODO AUROC of 0.519 for a target-pair-feature selectivity benchmark` | High |
| `docs/technical_report.md` | `Compound-level selectivity data was assembled` | The labels are compound-level, but the active v3 feature matrix is not compound-specific. The wording can mislead about what is actually predicted. | `compound labels were assembled, then paired with target-pair-level pocket features` | Medium |
| `docs/technical_report.md` | `CHEMBL155771 ... was run through the full pipeline` | `case_study_chembl155771.py` is not a full structure-to-docking-to-attribution pipeline; it uses hard-coded compound metadata and curated pocket sequences. | `CHEMBL155771 was analyzed through a curated case-study script using stored metadata and pocket definitions` | High |
| `docs/technical_report.md` | `predicts no selectivity` / `predicts selectivity window` | The case study does not run a trained predictor there; it narrates interpretation from cosine similarity and pocket identity. | `is consistent with a no-selectivity intuition` / `suggests a possible selectivity window` | Medium |
| `docs/kira-technical-breakdown.md` | `Binding energy comparison → computational selectivity estimate` | The active package has no docking module and no binding-energy engine. The historical docking script produces docking scores, not rigorous binding energies. | `docking-score comparison → rough computational selectivity hypothesis` | High |
| `docs/kira-technical-breakdown.md` | `Docking validated by significant correlation with experimental IC50` | Correlation in one retrospective docking setup is not validation in the broad scientific sense, especially with non-equivalent sites and rigid docking caveats. | `docking setup showed a modest retrospective correlation with experimental IC50 in this case study` | Medium |
| `docs/kira-technical-breakdown.md` | `the non-selectivity finding is genuine` | The docking result may be suggestive, but the caveats listed in the same paragraph mean `genuine` is too strong. | `the docking result is suggestive but limited by rigid docking and site non-equivalence` | High |

## Docstrings

| File path | Quoted phrase | Why it is overstated | Honest narrower wording | Severity |
|---|---|---|---|---|
| `src/kira/physics/core/parser.py` | `Reads PDB and minimal mmCIF files` | The implementation only parses PDB-like ATOM/HETATM records. There is no mmCIF parser. | `Reads PDB-format ATOM/HETATM records into a Structure dataclass` | High |
| `src/kira/causality/__init__.py` | `Causality module` | The module performs attribution and divergence analysis, not causal inference. | `Mechanistic hypothesis module` | High |
| `src/kira/causality/binding_site.py` | `Given a protein structure and a ligand, extracts binding pocket residues` | In active use, the `ligand` can be just a single centroid point from CLI coordinates, not an explicit ligand structure. | `Given a protein structure and ligand coordinates or a centroid, extracts nearby residues` | Medium |
| `src/kira/causality/binding_site.py` | `This is the mechanistic explanation of selectivity: which positions differ and how that affects binding` | The function only compares residue identities at aligned pocket positions. It does not model how those changes affect ligand binding. | `This identifies which pocket positions differ between two curated sites` | High |
| `src/kira/causality/energy_decomp.py` | `given a compound docked into a parasite target and its human orthologue` | The function receives no ligand or docking pose. It compares two residue-energy decompositions. | `given two pocket residue-energy decompositions` | High |
| `src/kira/causality/selectivity_map.py` | `For a compound docked into both a target and its human ortholog` | No ligand pose is passed into `build_selectivity_map()` or the active CLI path. | `For two homologous pocket decompositions` | High |
| `src/kira/causality/selectivity_map.py` | `binding energy difference` | The underlying deltas are not ligand-conditioned binding energies. | `residue-energy difference` | High |
| `src/kira/causality/selectivity_map.py` | `mechanistic explanation of selectivity` | The module generates heuristic attributions from pocket comparison plus protein-only LJ decomposition. | `mechanistic hypothesis for selectivity` | High |
| `src/kira/experiments/case_study_chembl155771.py` | `End-to-end selectivity explanation` | The script is not end-to-end in the structural sense; it uses hard-coded compound metadata and curated pocket sequences, not live docking or ligand-conditioned energy calculation. | `Curated selectivity case study` | High |
| `src/kira/experiments/case_study_chembl155771.py` | `FULL causality pipeline` | The script does not invoke the active pocket extraction, energy decomposition, or selectivity-map pipeline from structures and ligand pose. | `curated explanatory workflow based on stored metadata and pocket sequences` | High |
| `src/kira/experiments/case_study_chembl155771.py` | `predicts NO selectivity` / `predicts selectivity window` | The script does not run predictive models in those steps; it prints interpretation. | `suggests no selectivity by global similarity` / `suggests a possible selectivity window by pocket divergence` | Medium |
| `src/kira/experiments/run_selectivity_v3.py` | `Load compound-level selectivity data` | The labels are compound-level, but the active features are pair-level constants reused for all compounds in the pair. | `Load compound labels for a target-pair-feature benchmark` | Medium |
| `src/kira/experiments/run_selectivity_v3.py` | `Selectivity Prediction v3` | The benchmark is weaker than the title suggests because it does not use compound-specific features in the current version. | `Selectivity Benchmark v3 using target-pair pocket features` | High |
| `src/kira/experiments/run_selectivity_v3.py` | `Binding-site divergence captures transferable selectivity signal` | The code supports a small retrospective benchmark with target-pair-level features. | `binding-site divergence features show modest pair-level transfer signal in this small benchmark` | Medium |
| `src/kira/experiments/selectivity_features.py` | `binding site predicts selectivity better than global protein similarity` | The code computes feature vectors; prediction performance comes from the later benchmark, and even there the evidence is limited by pair-level feature reuse. | `binding-site feature vectors are used in a small benchmark against global similarity baselines` | Medium |
| `src/kira/experiments/validate_physics.py` | `Validate physics engine against real PDB structures` | The script checks plausibility on a few structures, not full validation of a physics engine. | `Check whether the structure-analysis code behaves plausibly on selected real PDB structures` | Medium |
| `src/kira/experiments/validate_physics.py` | `Physics engine validated on real protein structures` | Passing a few sanity checks on a few proteins is not broad validation. | `physics checks behaved plausibly on this small test set` | High |
| `src/kira/targets.py` | `Best-validated target` / `Well-validated` / `Selectivity window exists` | These comments are strong biological claims embedded as priors, not results established by the active code. | `high-priority curated target` / `commonly studied target` / `possible selectivity window has been proposed` | Medium |

## CLI

| File path | Quoted phrase | Why it is overstated | Honest narrower wording | Severity |
|---|---|---|---|---|
| `src/kira/cli.py` | `Kira Engine CLI — unified entry point for causal drug discovery` | The CLI does not implement causal discovery; it wraps heuristic scoring, structure checks, and attribution utilities. | `Kira CLI — entry point for computational selectivity and structure-analysis utilities` | High |
| `src/kira/cli.py` | `A unified scientific engine for drug repurposing, physics validation, and selectivity analysis` | The active CLI surface does not cover the full historical repurposing pipeline and exposes reduced validation functionality. | `a CLI exposing parts of the repository’s scoring, structure-checking, and selectivity-analysis workflow` | Medium |
| `src/kira/cli.py` | `binding energy differential` | `selectivity()` does not compute protein-ligand binding energy. It computes residue-energy differences from protein-only decomposition after pocket selection by centroid. | `residue-energy differential within centroid-defined pockets` | High |

## Summary

The most severe inflation patterns are:

- `causal` language applied to hypothesis-generation modules
- `binding` and `binding energy` language applied to non-ligand-conditioned code paths
- `prediction` language applied to a pair-level benchmark
- `validation`, `Trust Report`, and `8 checks` language applied to a reduced active validator
- `full pipeline` and `end-to-end` language applied to curated narrative workflows
