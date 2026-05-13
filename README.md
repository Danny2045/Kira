# Kira: Open Selectivity Engine for Neglected Tropical Diseases

Kira is a computational repository for retrospective translational analysis, pair-level selectivity benchmarking, curated biological comparison, and residue-level mechanistic hypothesis generation against parasitic targets in schistosomiasis, human African trypanosomiasis, and leishmaniasis. It is an evidence substrate, a benchmark engine, and an experiment-design compiler. It is not a validated drug-discovery engine, a compound-conditioned selectivity predictor, a binding-affinity predictor, or a source of wet-lab validated mechanism.

**Evidence tier policy.** This repository documents work at Constitution evidence tiers one through four (descriptive, retrospective translational analysis, benchmark correlation, and mechanistic hypothesis generation). No module in this repository qualifies for tier five (validated causal or experimental evidence). Every public-facing claim in this README, in the docs, in the docstrings, and in the CLI is bounded by the Scientific Constitution at `docs/SCIENTIFIC_CONSTITUTION.md`, the Claim Inflation Audit at `docs/CLAIM_INFLATION_AUDIT.md`, and the Disposition Plan at `docs/DISPOSITION_PLAN.md`.

## What Kira does

Kira performs systematic selectivity analysis on chemical compound activities against parasite and human enzyme targets. It pulls activity measurements from ChEMBL, computes selectivity ratios for compounds with both parasite and human measurements, classifies the resulting evidence by assay matchedness and statistical decidability, and produces structured benchmark substrates from the resulting evidence. It generates residue-level mechanistic hypotheses about why specific compounds are selective for parasite enzymes over their human counterparts. It compiles missing-evidence gaps into ranked experimental tickets that can be returned to wet-lab collaborators. It provides domain-general schemas that extend the contrast pattern to antimicrobial resistance surveillance and regenerative biology.

The empirical work covers three diseases, seven target pairs, and approximately twenty-seven hundred unique compounds drawn from ChEMBL. The headline finding, supported by the data, is that the majority of antiparasitic compounds tested against conserved targets are non-selective when assessed against their human counterparts. Across seven target pairs, weighted non-selectivity is approximately sixty-two percent. The most selectivity-favorable target in the dataset is leishmania pteridine reductase 1 with a median sixty-eight-fold selectivity over human dihydrofolate reductase. The single most-selective compound identified is CHEMBL155771 at thirty-fold selectivity for *Schistosoma mansoni* dihydroorotate dehydrogenase over the human ortholog.

## What Kira does not claim

Kira does not claim to discover causal mechanisms, validate drug candidates experimentally, predict compound-conditioned selectivity from explicit protein-ligand physics, compute true binding free energies, implement a complete docking engine in the active package, perform the full eight-check structure validation workflow described in early documentation, parse mmCIF files, or generate clinical-grade recommendations. These boundaries are formalized in the Scientific Constitution and enforced by the Disposition Plan's per-module assignments.

## Components

The active codebase is organized into the following components.

**Selectivity pipeline.** The `src/kira/experiments/` package contains the v3, v4, v5, and v6 selectivity pipelines. The v4 benchmark is a compound-conditioned scaffold-aware retrospective benchmark using StratifiedGroupKFold grouped on target pair and Murcko scaffold. The v5 evidence-substrate expansion pipeline classifies ChEMBL records into evidence-status categories including exact matched ratios, intervals, lower and upper bounds, single-side observations, and unmatched comparables. The v5 exact-core and tiered-core builders distill the candidate substrate into strict and logically-decidable trainable rows. The v6 lab-campaign designer ranks missing-evidence gaps into structured experimental tickets in benchmark-repair and potency-discovery modes. Evidence tier: 3 (benchmark proof-of-concept).

**Contrast core.** The `src/kira/contrast/` package defines domain-agnostic primitives for biological contrasts: comparison of a desired biological context against a control or failure context under a measurable readout. It provides typed schemas for contrast specifications, evidence records, evidence statuses, experiment tickets, and data-return schemas, with no chemistry, model, or assay-runtime dependencies. The contrast core is the substrate that the AMR, regeneration, and physics-auditor adapters build on. Evidence tier: 1 (infrastructure).

**Domain scouts.** Three scout adapters extend the contrast pattern to specific biological domains. The `src/kira/amr/` package implements a Rwanda antimicrobial resistance contrast scout with eight seed contrasts covering AST completeness, sentinel surveillance, stewardship, outbreak clustering, genomic confirmation, antibiotic consumption, AST quality control, and benchmark-complete records, together with an AST completeness audit and a collaborator-facing CSV data-return kit. The `src/kira/regeneration/` package implements a regeneration and CRISPR contrast scout with six seed contrasts covering bioelectric modulation, planarian target-morphology perturbation, assembloid CRISPR, three-dimensional organoid CRISPR, partial reprogramming, and morphogen organoid phenotyping. The `src/kira/selectivity/benchmark_repair_report.py` module renders a deterministic dossier over the v4, v5, and v6 artifacts for the original parasite selectivity domain. Evidence tier: 1 (domain adapters over infrastructure), with explicit non-claims of therapeutic discovery, wet-lab validation, or model-performance results.

**Physics-auditor ticket gate.** The `src/kira/physics/ticket_gate.py` module is a deterministic lexical gate that checks whether an experiment ticket is observable, falsifiable, bounded in units and timing, and grounded in an operational scale. It does not evaluate whether an intervention works. Evidence tier: 1 (infrastructure).

**Approximate structure-checking toolkit.** The `src/kira/physics/` package contains a PDB parser, topology inference, geometry kernels, Lennard-Jones energy computation, and steric clash detection. The toolkit is approximate, not a validated force field, and exposes only the currently implemented checks rather than the full eight-check suite described in early planning documents. Evidence tier: 1 (descriptive structure analysis).

**Mechanistic hypothesis module.** The `src/kira/causality/` package, despite its name, performs mechanistic hypothesis generation rather than causal inference. It computes binding-site extraction from ligand coordinates, per-residue Lennard-Jones decomposition, residue-level selectivity attribution between curated pocket pairs, and ESM-2 protein language model divergence between full sequences and pocket-restricted residues. The outputs are interpretive hypotheses suitable for guiding follow-up experimentation, not validated mechanisms. Evidence tier: 4 (mechanistic hypothesis generation).

**Historical pipeline.** The `archived/scripts/` directory contains the twenty-one numbered scripts that produced the original preprint analysis from February through April 2026. These scripts are preserved as the historical evidence record per the Disposition Plan. They are not part of the active reproduction path. Evidence tier: 2 (retrospective translational analysis), historical.

## Repository layout

```
src/kira/
  __init__.py
  chemistry.py            # standard medicinal-chemistry heuristics
  cli.py                  # command-line interface
  contrast/               # domain-agnostic contrast schemas
  amr/                    # Rwanda AMR scout, audit, data-return kit
  regeneration/           # regeneration and CRISPR scout
  selectivity/            # parasite selectivity benchmark-repair report
  experiments/            # v3, v4, v5, v6 selectivity pipeline
  causality/              # mechanistic hypothesis module
  physics/                # approximate structure-checking + ticket gate
  data/                   # target manifest validation

archived/
  scripts/                # historical 21-script pipeline (Disposition: Evidence)

data/
  publication/            # current 3-disease results
  models/                 # ESM-2 embeddings, model results
  lab_requests/           # v6 campaign tickets
  reference/              # target-pair configurations

docs/
  SCIENTIFIC_CONSTITUTION.md      # binding governance
  CLAIM_INFLATION_AUDIT.md        # overclaim catalog and replacements
  DISPOSITION_PLAN.md             # per-module Evidence/Engine/Hypothesis assignments
  CONTRAST_ENGINE.md              # contrast core reference
  PHYSICS_AUDITOR_TICKET_GATE.md  # ticket gate reference
  REGENERATION_CONTRAST_SCOUT.md  # regeneration scout reference
  RWANDA_AMR_*.md                 # AMR scout, audit, data-return references
  PARASITE_SELECTIVITY_*.md       # selectivity dossier reference
  KIRA_PROJECT_LOG.md             # branch and PR history
  KIRA_SCIENTIFIC_OPERATING_MODEL.md  # process documentation

results/
  selectivity_v3_results.txt
  selectivity_v4/         # v4 ablation summaries
  case_study_chembl155771.txt
```

## Quick start

The repository can be installed in editable mode and exercised against committed local artifacts without network access for most workflows.

Install:

```bash
pip install -e .
```

Run the full test suite:

```bash
pytest -q
```

Reproduce the v4 ablation study:

```bash
python -m kira.experiments.run_selectivity_v4
```

Generate the parasite selectivity benchmark-repair report:

```bash
python -m kira.selectivity.benchmark_repair_report
```

Inspect contrast scouts:

```bash
python -m kira.amr.scout
python -m kira.regeneration.scout
```

## Related work

This repository is one half of a two-repository platform. The companion repository, Physics Auditor at `github.com/Danny2045/physics-auditor`, provides protein structure provenance verification, approximate physics validation, and per-residue selectivity attribution as an independently-installable package. Kira and Physics Auditor share scientific lineage but maintain separate dependencies and release surfaces. Future integration between the two repositories operates through the structure-audit request and result interface described in `docs/PHYSICS_AUDITOR_TICKET_GATE.md`.

## Scientific governance

All public claims in this repository are bounded by three governance documents committed to `docs/`:

The **Scientific Constitution** at `docs/SCIENTIFIC_CONSTITUTION.md` defines what Kira is and is not allowed to claim, lists forbidden overclaim terms with approved replacements, and specifies the evidence-tier policy that every public-facing claim must declare.

The **Claim Inflation Audit** at `docs/CLAIM_INFLATION_AUDIT.md` catalogs specific overclaim language in earlier versions of the README, docs, docstrings, and CLI, with severity ratings and approved narrower wording.

The **Disposition Plan** at `docs/DISPOSITION_PLAN.md` assigns each module to one of five categories — Evidence, Engine, Hypothesis, Split, or Deprecate — with scientific reasoning for each assignment.

These documents are binding governance. They were committed on April 8, 2026, reverted the same day, and restored on May 13, 2026, after independent audits rediscovered most of their findings. Subsequent work on the repository is expected to comply with them or to propose amendments through a separate PR.

## License

MIT License. See `LICENSE`.

## Citation

If you use Kira in academic work, please cite the repository directly until peer-reviewed publication is available. The historical preprint at `docs/kira-final-preprint-v2.docx` is preserved as a historical artifact and is being revised; it does not represent the current state of the repository's claims.

## Contact

Daniel Ngabonziza. Franklin, Tennessee. Independent researcher.
