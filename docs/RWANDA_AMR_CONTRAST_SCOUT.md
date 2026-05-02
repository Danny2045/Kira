# Rwanda AMR contrast scout

This branch adds a narrow Rwanda AMR scout adapter that turns surveillance,
stewardship, infection-control, antibiotic-use, and lab-data gaps into Kira's
contrast-core shape:

```text
intervention
-> desired AMR context
-> control / failure / missing-data context
-> measurable readout
-> evidence status
-> missing measurement
-> executable experiment or data ticket
```

## Scope

The scout covers eight fixed seed contrasts:

1. AST completeness for priority bacterial isolate records.
2. Sentinel surveillance sampling completeness.
3. Stewardship evidence gap for AST-linked antibiotic-use audit records.
4. Infection-control outbreak-cluster contrast.
5. Genomic resistance marker confirmation against phenotypic AST.
6. Antibiotic consumption paired with resistance-rate measurement.
7. AST quality-control lab-data ticket.
8. Benchmark-complete pathogen-antibiotic pair record.

The adapter is deliberately small. It uses Rwanda as the first provenance and
traction setting, but the output shape is the same world-scale contrast
infrastructure needed for comparable AMR rows across places, facilities,
periods, pathogens, antibiotics, and evidence states.

## Non-claims

The scout does not provide individual-patient treatment advice. It does not say
Kira can fix AMR, does not assert that any lab validation has been completed,
does not infer facility performance, and does not turn surveillance gaps into
intervention success claims.

AMR is treated here as an evidence, surveillance, stewardship, and
benchmark-repair problem. The strongest object this adapter emits is an
executable data ticket for a missing or incomplete contrast row.

## Evidence-status discipline

The seed set intentionally includes multiple evidence states:

| Status | Meaning in this scout |
|---|---|
| `paired_complete` | Desired and control sides are represented enough to form a first benchmark row. |
| `single_side_only` | One useful side exists, but a paired comparator or denominator is missing. |
| `proposed` | The seed is only an executable ticket shape. |
| `bounded` | The contrast is useful only when time, place, denominator, QC, or case-definition boundaries are explicit. |
| `conflicting` | Marker, AST, or surveillance signals require an adjudicating paired measurement. |

Every seed defines a contrast, readout type, evidence status, benchmark label,
missing measurement, and executable experiment/data ticket.

## Code surface

New package:

```text
src/kira/amr/__init__.py
src/kira/amr/scout.py
tests/test_amr_scout.py
docs/RWANDA_AMR_CONTRAST_SCOUT.md
```

Public helpers:

```python
from kira.amr import (
    amr_scout_seeds,
    build_amr_scout,
    build_contrast_specs,
    build_evidence_records,
    build_experiment_tickets,
    summarize_evidence_statuses,
    tickets_as_dicts,
    validate_amr_seed,
)
```

## Expected validation commands

```bash
ruff check src/kira/amr tests/test_amr_scout.py
pytest -q tests/test_amr_scout.py
pytest -q tests/test_contrast_schemas.py tests/test_contrast_tickets.py tests/test_amr_scout.py
ruff check .
pytest -q
```

## Source anchors

These anchors motivate seed shapes. They are not claims that Kira has produced
new AMR findings.

- Rwanda Biomedical Centre, *2nd National Action Plan on AMR (2025-2029)*:
  `https://www.rbc.gov.rw/fileadmin/user_upload/strategy/2nd_NATIONAL_ACTION_PLAN_ON_AMR__NAP_2025-2029_.pdf`
- WHO publication page for *Rwanda: National action plan on antimicrobial resistance 2020-2024*:
  `https://www.who.int/publications/m/item/rwanda-national-action-plan-on-antimicrobial-resistance-2020-2024`
- Fleming Fund Rwanda country page:
  `https://www.flemingfund.org/countries/rwanda/`
- WHO GLASS-AMR routine data surveillance:
  `https://www.who.int/initiatives/glass/glass-routine-data-surveillance`
- WHO GLASS-AMU module:
  `https://www.who.int/initiatives/glass/glass-amc-module`
- WHO GLASS methodology for national antimicrobial consumption surveillance:
  `https://www.who.int/publications/i/item/9789240012639`
- WHO GLASS manual for AMR surveillance in common bacteria:
  `https://www.who.int/publications/i/item/9789240076600`
- CLSI M100 antimicrobial susceptibility testing standards page:
  `https://clsi.org/shop/standards/m100/`

## Rwanda-first traction, world-scale infrastructure

Rwanda is a useful first scout setting because its AMR planning and One Health
coordination make the missing-measurement problem concrete: AST completeness,
sampling denominators, site coverage, antibiotic-use linkage, IPC boundaries,
genotype-phenotype pairing, and pathogen-antibiotic benchmark rows.

The contrast-core framing keeps that local traction portable. A completed
Rwanda AMR row has the same schema elements a global benchmark needs:
provenance, desired side, control side, readout, evidence state, missing side,
and a ticket that states exactly what data must come back.
