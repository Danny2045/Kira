# Rwanda AMR AST completeness audit

This audit moves the Rwanda AMR work from scout tickets toward concrete
data-quality checks. It inspects isolate-level or facility-period AMR records
and decides whether each row has enough non-empty AST metadata to support an AMR
contrast or benchmark row.

## Why AST completeness matters

AMR benchmark rows are only useful when the tested organism, specimen source,
antibiotic, AST method, result or count denominator, breakpoint version, and QC
status are present. Missing AST metadata can make a row look comparable while
leaving the measurement boundary unclear.

The audit therefore checks completeness before any downstream benchmark use. It
does not interpret patient care, rank facilities, assert lab completion, or turn
surveillance gaps into outcome claims.

## Required isolate-level fields

An isolate-level row is benchmark-ready only when all of these fields are
present and non-empty:

| Field | Meaning |
|---|---|
| `facility_id` | Reporting facility identifier. |
| `reporting_period` | Bounded reporting period, such as facility-month. |
| `organism` | Organism identification. |
| `specimen_source` | Specimen source or specimen type. |
| `antibiotic` | Tested antibiotic. |
| `ast_method` | AST method used for the isolate-antibiotic result. |
| `ast_result` | Isolate-antibiotic result, such as S/I/R or a reported category. |
| `breakpoint_version` | Breakpoint table or version used for interpretation. |
| `qc_status` | AST quality-control status for the row. |

Useful optional fields include `isolate_id`, `patient_age_group`, `ward`,
`specimen_date`, `mic_value`, `sir_category`, `resistant_count`,
`susceptible_count`, and `tested_count`.

## Isolate-level vs aggregate rows

The audit supports dictionaries, dataclass instances, objects with `to_dict()`,
and simple objects with attributes.

Rows with an explicit aggregate `record_type` or facility-period count fields
are treated as aggregate rows when they do not carry a per-isolate `ast_result`.
Aggregate rows use count completeness instead of isolate-result completeness.

Aggregate rows require:

| Field | Meaning |
|---|---|
| `facility_id` | Reporting facility identifier. |
| `reporting_period` | Bounded reporting period. |
| `organism` | Organism or pathogen boundary for the count. |
| `specimen_source` | Specimen source or specimen type boundary. |
| `antibiotic` | Tested antibiotic. |
| `tested_count` | Count of tested isolates in the facility-period row. |
| `resistant_count` or `susceptible_count` | At least one result-side count. |
| `ast_method` | AST method used for the count row. |
| `breakpoint_version` | Breakpoint table or version. |
| `qc_status` | AST quality-control status for the count row. |

## Benchmark-ready rule

The rule is intentionally narrow:

```text
benchmark-ready = all required fields for the detected row type are present
                  and non-empty
```

Blank strings, `None`, NaN floats, and empty containers count as missing. Numeric
zero is present, so `resistant_count = 0` or `susceptible_count = 0` can satisfy
aggregate count completeness.

The report computes:

- `total_records`
- `complete_records`
- `incomplete_records`
- `completeness_percent`
- `missing_field_counts`
- `benchmark_ready_count`
- `non_benchmark_ready_count`
- `benchmark_ready_rows`
- `incomplete_rows`
- repair actions grouped by missing field

## Code surface

Public helpers live in `kira.amr.audit` and are re-exported from `kira.amr`:

```python
from kira.amr import (
    AmrAuditRecord,
    AmrCompletenessReport,
    MissingFieldFinding,
    audit_ast_completeness,
    audit_record,
    benchmark_ready_records,
    data_repair_tickets,
    report_to_dict,
    summarize_ast_completeness,
)
```

`report_to_dict(report)` returns a deterministic JSON-ready payload. Repair
actions are grouped by missing field and point back to the AST completeness scout
ticket.

## Provenance and the mixed-file refusal boundary

Synthetic-vs-real provenance is classified from the `synthetic_data_notice`
column (the single source of truth). The mixed-provenance / inconsistent-notice
refusal guarantee — a file with some rows marked synthetic and others unmarked,
or with differing notice strings, raises `ValueError` instead of being processed
— holds at the **CSV entry points**: `audit_amr_csv()` and `make_amr_csv_report()`
(and `write_amr_csv_report()`). Use these whenever the input is a CSV.

`audit_ast_completeness(load_amr_csv(...))` is a deliberate lower-level
composition: `load_amr_csv()` is a pure CSV reader that returns rows only and
**intentionally does not classify provenance or refuse mixed files**, and
`audit_ast_completeness()` consumes already-in-memory records. This primitive
path is for callers that have established provenance by other means; it does not
carry the synthetic banner or the mixed-file guarantee. For CSV-sourced reports,
go through `make_amr_csv_report()` / `write_amr_csv_report()` so the report is
provenance-stamped and mixed files are refused.

## Non-claims

This audit is only a data-completeness and benchmark-readiness check.

It does not:

- give patient-care or medication-use guidance
- claim Kira has fixed AMR
- assert completed wet lab validation
- infer facility performance
- make public-health outcome claims

## Rwanda AMR scout and Physics Auditor connection

The audit is tied to the Rwanda AMR scout concept:

```text
rwanda-amr-ast-completeness-priority-isolates
```

Reports include this value as `source_ticket_id`. They also include a narrow
Physics Auditor ticket-gate summary for that AST-completeness scout ticket:
`ticket_gate_status` and `ticket_gate_check_sections`.

This is not a global AMR-ticket gate. It connects only the AST completeness audit
to the existing scout ticket so downstream work can trace why the audit exists
and which measurement concept it repairs.


## Why this matters beyond AMR

This audit is Kira's first concrete returned-data readiness pattern. It is small
by design: AMR AST completeness is a practical Rwanda-relevant measurement
problem. But the same pattern generalizes to future virtual-biology models,
CRISPR perturbation screens, regeneration readouts, autonomous labs, and
biofoundry-style experimental systems.

Kira's role is not to replace predictive models, wet labs, autonomous labs, or
public-health authorities. Kira defines the contrast, checks whether the ticket
is measurable through the Physics Auditor gate, and specifies the data that must
return before a claim becomes benchmark-ready.

This document does not claim partnership or equivalence with Biohub, Arc,
Ginkgo, IGI, Levin, Doudna, or any institution or scientist.

## Validation commands

```bash
ruff check src/kira/amr tests/test_amr_audit.py
pytest -q tests/test_amr_audit.py
pytest -q tests/test_amr_scout.py tests/test_physics_ticket_gate.py tests/test_amr_audit.py
ruff check .
pytest -q
```
