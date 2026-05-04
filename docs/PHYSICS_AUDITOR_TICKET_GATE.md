# Physics Auditor Ticket Gate

## Why This Gate Exists

Kira now emits contrast-core experiment and data tickets from the regeneration
scout and Rwanda AMR scout, while v6 provides earlier benchmark-repair ticket
machinery. Those tickets can
look scientific even when they are still too vague to measure. The Physics
Auditor ticket gate prevents that failure mode by checking whether a ticket is
observable, falsifiable, bounded in units and timing, and grounded in an
operational scale.

The gate does not evaluate whether an intervention works. It evaluates whether
the ticket says enough about the measurement for a future data return to repair,
reject, or downgrade the contrast.

## What It Checks

The public module is `kira.physics.ticket_gate`.

Core entry points:

- `audit_ticket(ticket)` returns one `TicketGateReport`.
- `audit_tickets(tickets)` returns one report per ticket.
- `report_to_dict(report)` returns a stable JSON-ready payload.
- `summarize_gate_reports(reports)` counts reports by status and domain.

Gate statuses:

- `pass`: required fields and measurement boundaries are present.
- `needs_measurement_detail`: the ticket is usable but lacks a sharper unit,
  timing, falsifiability, or domain-specific detail.
- `blocked`: required fields, observability, data-return schema, falsifiability,
  or non-claim boundaries fail.

Checks:

- Required contrast fields: verifies the contrast-core ticket fields needed to
  identify the contrast, context, readout, evidence status, provenance, question,
  and benchmark impact.
- Units: requires non-empty `readout_units` and warns on generic unit language.
- Scale: infers a measurement scale such as molecular, cellular, organoid,
  tissue, organism, facility, surveillance, or population.
- Timescale: looks for a time window, reporting period, assay duration,
  perturbation timing, or bounded timing language.
- Observability: requires a returned data kind such as an AST record,
  resistance-rate table, morphology score, marker panel, bioelectric map,
  image/phenotype table, facility-month data, sequence/marker call, or QC record.
- Falsifiability: checks that the ticket states how measurement can change the
  contrast, not just decorate it.
- Data-return schema: verifies that the contrast-core required fields are in the
  schema. For current `ExperimentTicket` objects, the gate can synthesize the
  schema with `make_data_return_schema()` because the core dataclass does not
  store that object as a field.
- Domain-specific checks: adds AMR and regeneration measurement boundaries.
- Non-claims: blocks broad claim language that is outside an audit report.

## Rwanda AMR Tickets

For `rwanda_amr` tickets, the gate checks AST, resistance-rate, and surveillance
records for operational fields such as denominator or counts, organism or
pathogen boundary, antibiotic or antimicrobial field, specimen or source
provenance, reporting period, AST linkage, and breakpoint or QC fields where
those are applicable.

This is surveillance and benchmark infrastructure. The gate does not provide
patient-level treatment direction or public-health outcome assertions.

## Regeneration Tickets

For `regeneration_crispr` tickets, the gate checks morphology, bioelectric,
organoid, CRISPR, and reprogramming tickets for dose or perturbation timing,
readout scale, control or failure arm, replicate or blinded scoring language,
and safety or toxicity boundary where applicable.

This is measurement discipline for contrast tickets. The gate does not assert a
therapy, a completed lab result, or domain completion.

## Validation Commands

Run the focused checks:

```bash
ruff check src/kira/physics tests/test_physics_ticket_gate.py
pytest -q tests/test_physics_ticket_gate.py
```

Run the related contrast and scout tests:

```bash
pytest -q tests/test_contrast_schemas.py tests/test_contrast_tickets.py tests/test_regeneration_scout.py tests/test_amr_scout.py tests/test_physics_ticket_gate.py
```

Run the full project checks:

```bash
ruff check .
pytest -q
```

## How This Moves Kira Toward Measurable Science

The scout branches created useful scaffolding: they turn scientific domains into
shared contrast records and experiment tickets. The ticket gate adds a stricter
next layer. It asks whether each ticket has a measurable readout, units, scale,
timing, observable return data, falsifiability rule, and domain-specific data
boundary.

That keeps future Kira branches from treating labels, narratives, or incomplete
schemas as evidence. A ticket can still be useful while needing measurement
detail, but the report makes that limitation explicit before downstream code or
documentation makes it look stronger than it is.
