# Kira Contrast Engine

Kira is becoming a contrast-driven discovery engine. The shared abstraction is a
measurable comparison between a desired biological context and a control, safety,
or failure context for a specific intervention and readout.

This is an infrastructure step. Kira v4 remains the current modeling claim.
Kira v5 and v6 established the first complete proof domain: neglected tropical
disease parasite-vs-human selectivity. The contrast core generalizes the shape
of that evidence problem without changing the existing NTD selectivity pipeline
and without making a new disease or predictive-model claim.

## Core Contract

Every contrast domain must define:

- `contrast_id`: stable identifier for the comparison.
- `domain`: adapter or evidence domain, such as `ntd_selectivity`.
- `intervention_type` and `intervention_id`: the thing being tested.
- `desired_context`: the biological context where the intervention is intended
  to have the desired effect.
- `control_context`: the comparator, safety, or failure context.
- `readout_type` and `readout_units`: the measurable output.
- `evidence_status`: whether evidence is proposed, one-sided, paired, bounded,
  conflicting, or otherwise explicitly described.
- `known_side` and `missing_side`: which side of the contrast is measured or
  absent.
- `label`: benchmark or workflow label, not a scientific conclusion by itself.
- `uncertainty_note`: scope limits and unresolved assumptions.
- `provenance`: source of the record.
- `experiment_question`: executable question for a lab, data, or benchmark task.
- `expected_benchmark_impact`: why this returned data would improve an evidence
  table or benchmark.

The contrast core is deliberately lightweight. It is a schema and ticket layer,
not a chemistry toolkit, docking layer, model runtime, or wet-lab validation
system.

## First Proof Domain

NTD selectivity is the first hard proof domain. In that setting, the desired
context can be parasite target activity, the control context can be the matched
human target comparator, and the readout can be a potency or selectivity ratio.
The v5/v6 work demonstrated how missing comparator evidence can be represented
and turned into executable lab-campaign tickets.

The contrast core preserves that pattern while keeping the existing v4, v5, and
v6 code intact.

## Future Adapters

Future adapters may include oncology, regeneration, CRISPR, biologics, and
biomanufacturing. Those adapters must map their domain records into the same
contrast contract: contrast, readout, evidence status, benchmark label, and
executable experiment ticket.

This document does not claim that Kira has solved oncology, regeneration,
CRISPR, biologics, or biomanufacturing. Those domains would require their own
evidence definitions, validation datasets, benchmark labels, and executable
experiments before any scientific claim could be evaluated.
