# Kira Scientific Operating Model

Kira must not become architecture theater. A branch is useful only when it moves
the project toward data, measurement, benchmark repair, experiment or data
tickets, returned-data audit, or stricter claim discipline.

This document defines how Kira work should be built, reviewed, and bounded.

## Operating Roles

Kira uses one writer and many critics.

- Human scientific owner: Daniel
- Builder: Codex
- Scientific reviewer/architect: ChatGPT

Daniel owns the scientific direction and final judgment. Codex implements
bounded changes in the repo. ChatGPT is used as a scientific reviewer and
architecture critic. The goal is disciplined pressure on the work, not a swarm
of agents producing disconnected artifacts.

## Hard Gates

Every meaningful branch should pass through hard gates:

- Git: clean branch intent, reviewable diff, no hidden local-state dependency
- Ruff: style and static checks
- pytest: focused tests plus full project tests when applicable
- CI: remote repeatability
- claim review: explicit non-claims and bounded conclusions

The default full-project validation commands are:

```bash
ruff check .
pytest -q
```

Branches that touch a narrow subsystem should also run focused tests and list
the exact commands in the PR body.

## Branch Standard

Each branch should answer:

- What data, measurement, benchmark, ticket, audit, or claim boundary changed?
- Which files changed?
- Which commands validate the change?
- What does the change not claim?
- What risk remains?
- What is the next repair step?

Documentation-only branches should say so and avoid scientific code changes.
Scientific-code branches should update docs when the change alters the project
record, operating model, benchmark substrate, or claim boundary.

## Agent Council as Review Discipline

The agent council is a review discipline, not a swarm. These lenses can be used
to criticize a branch before it is treated as serious scientific progress:

- Orchestrator: checks branch scope, sequencing, and whether the work moves the
  Kira loop forward.
- Lattner Architecture Agent: checks whether abstractions are simple, typed
  where useful, maintainable, and integrated with existing repo patterns.
- Carmack Reproducibility Agent: checks exact commands, clean diffs,
  deterministic artifacts, and absence of local-state ambiguity.
- Sutskever/Karpathy Benchmark Agent: checks benchmark truth, leakage
  avoidance, class balance, splits, metric meaning, and whether a result is
  actually supported.
- Assay/Ontology Agent: checks assay comparability, units, target naming,
  organism boundaries, evidence status, and ontology drift.
- Lab Translation Agent: checks whether tickets can become real measurements,
  protocols, returned-data schemas, or repair actions.
- Physics Auditor Agent: checks observability, falsifiability, bounded units,
  timing, scale, and returned-data readiness.
- Hype/Claim Boundary Agent: removes unsupported clinical, therapeutic,
  public-health, model-performance, wet-lab, partnership, or equivalence
  claims.
- Rwanda Relevance Agent: checks whether Rwanda-facing work is practical,
  respectful, locally relevant, and careful about public-health authority.
- Frontier Biology Judgment Gate: checks whether CRISPR, regeneration,
  bioelectric, morphology, or reprogramming language stays tied to measurable
  state changes and avoids frontier-biology overclaiming.

These names are internal review lenses. They do not imply endorsement,
partnership, affiliation, or equivalence with any named person, institution, or
lab.

## Carmack-Doudna-Ilya-Levin Principle

The Carmack-Doudna-Ilya-Levin principle is an operating shorthand for four
review pressures. It is not a claim that John Carmack, Jennifer Doudna, Ilya
Sutskever, Andrej Karpathy, Michael Levin, or any associated institution has
endorsed, reviewed, partnered with, or validated Kira.

Carmack lens: reproducibility. A Kira result should have exact commands,
reviewable diffs, deterministic artifacts when possible, and no hidden local
state. If another reviewer cannot rerun the branch from the repo, the result is
not yet operationally real.

Doudna lens: perturbation precision and safety boundaries. CRISPR or
perturbation language must identify the measured causal intervention, the
readout, the control or failure mode, and the safety boundary. Kira should not
turn gene-editing vocabulary into broad therapeutic claims.

Ilya/Sutskever/Karpathy lens: benchmark truth. Model and benchmark claims must
respect leakage avoidance, train/test separation, class balance, assay
comparability, calibration of metrics, and result support. Kira should not
perform metric theater or imply model progress from dataset scaffolding.

Levin lens: state transitions and pattern restoration. Regeneration,
bioelectric, and morphology work must specify state, transition, readout,
timing, scale, and failure criteria. Kira should not claim regeneration from a
ticket, scaffold, or unvalidated contrast.

These lenses are useful because they keep Kira's ambition tied to measurement.
They are inspiration lenses and operating principles only.

## Claim Discipline

Default non-claims for Kira branches:

- Kira did not discover drugs unless a branch explicitly supports that claim
  with appropriate evidence.
- Kira did not solve AMR.
- Kira did not solve regeneration.
- Kira does not provide clinical recommendations.
- Kira does not provide prescribing advice.
- Kira did not wet-lab validate tickets unless wet-lab evidence is explicitly
  present and reviewed.
- Kira did not prove a new model-performance result unless the benchmark was
  rerun and the result is explicitly supported.
- Kira does not claim endorsement, partnership, affiliation, or equivalence with
  named scientists, labs, companies, institutions, or public-health bodies
  unless formal evidence is present.

The claim boundary is part of the science. A branch that improves claims from
"impressive-sounding" to "actually supported" has moved Kira forward.

## Definition of Progress

Progress means at least one of these changed in a durable way:

- messier evidence became auditable evidence
- a contrast became explicit
- evidence status became clearer
- a benchmark core became cleaner
- a missing measurement became visible
- an experiment or data ticket became more measurable
- a Physics Auditor gate caught ambiguity
- returned data became audit-ready
- benchmark-ready rows were produced
- the next repair ticket became sharper

Anything else should be treated as narrative until it connects back to data,
measurement, benchmark repair, or claim discipline.
