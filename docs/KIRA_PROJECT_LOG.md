# Kira Project Log

This is the permanent repo record for meaningful Kira branches. GitHub PRs
capture the review event. This log captures why each branch mattered
scientifically after the PR has merged.

## Current Strategic Frame

Kira is a Rwanda-rooted, contrast-driven inverse-biology system.

Its job is not to decorate biology with architecture. Its job is to turn messy
evidence into bounded contrasts, expose missing measurements, create experiment
or data tickets, audit returned data, and repair benchmark substrates without
overstating claims.

## Core Kira Loop

```text
messy evidence
-> contrast
-> evidence status
-> benchmark core
-> missing measurement
-> experiment/data ticket
-> Physics Auditor gate
-> returned-data audit
-> benchmark-ready rows
-> next repair ticket
```

Every meaningful branch should improve at least one step in this loop and leave
two records:

- a GitHub record: PR title/body, scientific purpose, changed files,
  validation, non-claims, risk, and next step
- a repo record: project log or operating-model update that survives after the
  PR is merged

## PR #5: v5 Assay-Aware Expansion and Evidence Audit

Branch/PR name: PR #5, v5 assay-aware expansion and evidence audit.

What changed: v5 expanded the parasite selectivity evidence substrate with
assay-aware ChEMBL processing, evidence auditing, and exact/tiered benchmark
cores.

Important numbers:

| Quantity | Count |
|---|---:|
| v5 candidate evidence rows | 12,091 |
| v5 curated activity rows | 11,703 |
| v5 exact matched-ratio candidate rows | 1,227 |
| v5 exact-core rows | 114 |
| v5 trainable exact-core rows | 110 |
| v5 tiered-core rows | 135 |
| v5 trainable tiered-core rows | 131 |

Tiered pairs with both classes: LmDHFR, SmHDAC8, TbCathB.

Scientific meaning: PR #5 made the public-data substrate cleaner and more
auditable. It showed where selectivity evidence is real enough to structure,
but also where class balance and matched evidence remain too weak for broad
claims.

Validation if known: v5 result artifacts were documented in
`docs/V5_RESULTS.md`.

Non-claims: v5 did not prove a new model-performance result. It did not claim
drug discovery or wet-lab validation.

Next consequence: v5 created the need for a repair-oriented lab campaign rather
than premature benchmark theater.

## PR #6: v6 Lab-Campaign Designer

Branch/PR name: PR #6, v6 lab-campaign designer.

What changed: v6 generated benchmark-repair tickets from the selectivity
substrate.

Important numbers:

| Target pair | benchmark_repair tickets |
|---|---:|
| LmDHFR | 5 |
| LmPTR1 | 25 |
| SmDHODH | 25 |
| SmHDAC8 | 15 |
| TbCathB | 5 |
| TbPDEB1 | 25 |
| Total | 100 |

Missing side distribution:

| Missing side | Tickets |
|---|---:|
| parasite | 99 |
| human | 1 |

Scientific meaning: PR #6 translated benchmark weakness into concrete repair
work. The ticket distribution showed that most repair need was parasite-side
measurement, not human-side data.

Validation if known: validation details are not recorded in this log.

Non-claims: v6 tickets are not wet-lab validated. They are repair priorities,
not biological proof.

Next consequence: the project needed durable contrast and ticket objects so
future scouts and repair branches could speak a shared language.

## PR #7: Contrast-Core Foundation

Branch/PR name: PR #7, contrast-core foundation.

What changed: PR #7 added the shared contrast-core object model and helper
functions.

Core objects:

- `ContrastSpec`
- `EvidenceRecord`
- `EvidenceStatus`
- `ExperimentTicket`
- `DataReturnSchema`

Core helper functions:

- `validate_contrast_spec`
- `validate_experiment_ticket`
- `ticket_to_dict`
- `ticket_from_dict`
- `make_data_return_schema`

Scientific meaning: PR #7 gave Kira stable internal contracts for contrasts,
evidence status, experiment tickets, and returned-data schemas. This made later
domain scouts less ad hoc.

Validation if known: validation details are not recorded in this log.

Non-claims: the contrast core is infrastructure. It does not prove a biological
result, validate a ticket, or establish benchmark performance.

Next consequence: the shared contracts made frontier scouts possible while
keeping them tied to explicit evidence and ticket schemas.

## PR #8: Regeneration / CRISPR / Bioelectric Contrast Scout

Branch/PR name: PR #8, regeneration / CRISPR / bioelectric contrast scout.

What changed: PR #8 added a frontier biology scout for regeneration, CRISPR,
bioelectric state, morphology, and reprogramming contrasts.

Scientific meaning: PR #8 tested whether the contrast-core language could hold
frontier biology without turning uncertainty into hype. It made regeneration
work expressible as contrast tickets with evidence status and data-return
needs.

Validation after merge:

```bash
ruff check .
pytest -q
# 273 passed, 5 skipped
```

Non-claims: Kira did not solve regeneration, did not solve CRISPR, and did not
wet-lab validate any frontier ticket.

Next consequence: the scout needed stronger measurement discipline so frontier
tickets could be audited before downstream use.

## PR #9: Rwanda AMR Contrast Scout

Branch/PR name: PR #9, Rwanda AMR contrast scout.

What changed: PR #9 added the first Rwanda-facing AMR evidence,
surveillance, stewardship, and benchmark-repair scout.

Scientific meaning: PR #9 grounded Kira in a practical Rwanda-relevant AMR
measurement problem. It connected contrast infrastructure to AST completeness,
surveillance data, and stewardship boundaries without giving clinical advice.

Validation after merge:

```bash
ruff check .
pytest -q
# 279 passed, 5 skipped
```

Non-claims: PR #9 made no clinical recommendations, gave no prescribing advice,
and did not claim that Kira solves AMR.

Next consequence: AMR and regeneration tickets both needed a gate that checks
whether a ticket is measurable, observable, bounded, falsifiable, and
data-return ready.

## PR #10: Physics Auditor Ticket Gate

Branch/PR name: PR #10, Physics Auditor ticket gate.

What changed: PR #10 added a ticket gate for auditing measurement readiness
across scout tickets.

Important ticket-gate summary:

| Status | Tickets |
|---|---:|
| pass | 4 |
| needs_measurement_detail | 10 |
| blocked | 0 |
| total audited scout tickets | 14 |

By domain:

| Domain | pass | needs_measurement_detail | blocked | total |
|---|---:|---:|---:|---:|
| rwanda_amr | 4 | 4 | 0 | 8 |
| regeneration_crispr | 0 | 6 | 0 | 6 |

Scientific meaning: Kira can now audit whether tickets are measurable,
observable, bounded, falsifiable, and data-return ready before treating them as
repair work.

Validation after merge:

```bash
ruff check .
pytest -q
# 283 passed, 5 skipped
```

Non-claims: the gate does not evaluate whether interventions work. It audits
ticket measurability and claim boundaries.

Next consequence: the next AMR step was to apply the measurement discipline to
returned AST records and benchmark readiness.

## PR #11: Rwanda AMR AST Completeness Audit

Branch/PR name: PR #11, Rwanda AMR AST completeness audit.

What changed: PR #11 added row-level audit logic for isolate-level and
aggregate/facility-period AMR records.

Scientific meaning: PR #11 turned AMR scout intent into a concrete
data-quality and benchmark-readiness audit. It checks whether rows contain the
organism, specimen, antibiotic, AST method/result or count denominator,
breakpoint version, and QC fields needed for benchmark construction.

Validation:

```bash
ruff check .
pytest -q
# 292 passed, 5 skipped
```

Non-claims: PR #11 made no clinical recommendations, gave no prescribing
advice, inferred no facility performance, and made no public-health outcome
claim.

Next consequence: the audit needed a collaborator-facing data-return kit so a
template, example CSV, report, and repair actions could form a returned-data
loop.

## PR #12: Rwanda AMR AST Data-Return Kit

Branch/PR name: PR #12, Rwanda AMR AST data-return kit.

What changed: PR #12 added Kira's first collaborator-facing returned-data loop.

Files:

- `src/kira/amr/data_return.py`
- `tests/test_amr_data_return.py`
- `docs/RWANDA_AMR_AST_DATA_RETURN_KIT.md`
- `examples/rwanda_amr_ast_template.csv`
- `examples/rwanda_amr_ast_example.csv`

Scientific meaning: PR #12 made the AMR completeness audit usable as a
returned-data workflow. A collaborator can start from a CSV template, return
AST rows, and receive a deterministic benchmark-readiness report with missing
field repair actions.

Validation after merge:

```bash
ruff check .
pytest -q
# 301 passed, 5 skipped
```

Non-claims: PR #12 used no real patient data and no real facility names. It
gave no clinical advice, did not claim Kira solves AMR, and made no
partnership or equivalence claims with Biohub, Arc, Ginkgo, IGI, Levin,
Doudna, or any institution or scientist.

Next consequence: the project can return to its strongest empirical public-data
substrate and generate a parasite selectivity benchmark-repair dossier.

## PR #14 — Parasite Selectivity Benchmark-Repair Report

Branch: `feat/parasite-selectivity-benchmark-repair-report`.

Changed files:

- `src/kira/selectivity/__init__.py`
- `src/kira/selectivity/benchmark_repair_report.py`
- `tests/test_selectivity_benchmark_repair_report.py`
- `docs/PARASITE_SELECTIVITY_BENCHMARK_REPAIR_REPORT.md`

Scientific meaning:

- returns to Kira's strongest empirical public-data substrate
- connects the v4 modeling claim, v5 evidence substrate, exact/tiered cores,
  and v6 benchmark-repair tickets
- produces a deterministic benchmark-repair dossier from committed local
  artifacts
- identifies all six v6 target pairs: LmDHFR, LmPTR1, SmDHODH, SmHDAC8,
  TbCathB, TbPDEB1
- records v5/v6 counts from local committed result files
- identifies v5 tiered pairs with both classes: LmDHFR, SmHDAC8, TbCathB
- records the v6 missing-side distribution: parasite 99, human 1

v5 counts:

| Quantity | Count |
|---|---:|
| candidate evidence rows | 12,091 |
| curated activity rows | 11,703 |
| exact matched-ratio candidate rows | 1,227 |
| exact-core rows | 114 |
| trainable exact-core rows | 110 |
| tiered-core rows | 135 |
| trainable tiered-core rows | 131 |

v6 benchmark-repair tickets:

| Target pair | Tickets |
|---|---:|
| LmDHFR | 5 |
| LmPTR1 | 25 |
| SmDHODH | 25 |
| SmHDAC8 | 15 |
| TbCathB | 5 |
| TbPDEB1 | 25 |
| Total | 100 |

v6 missing-side distribution:

| Missing side | Tickets |
|---|---:|
| parasite | 99 |
| human | 1 |

Validation:

```bash
ruff check src/kira/selectivity tests/test_selectivity_benchmark_repair_report.py
pytest -q tests/test_selectivity_benchmark_repair_report.py
# 8 passed

ruff check .
# All checks passed

pytest -q
# 309 passed, 5 skipped
```

Non-claims:

- no drug-discovery claim
- no wet-lab validation claim
- no clinical claim
- no new model-performance claim beyond the committed v4 summary
- no claim that parasite diseases are solved
- no claim that any v6 ticket has already been executed in a wet lab

Next consequence: use the report to build a collaborator-facing selectivity
benchmark-repair data-return kit.

Future branch: `feat/selectivity-benchmark-repair-data-return-kit`.
