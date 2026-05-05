# Rwanda AMR AST data-return kit

This kit turns the merged Rwanda AMR AST completeness audit into a practical
returned-data loop. A collaborator can start from a CSV template, fill synthetic
placeholder rows with local isolate-level or facility-period AST data, run the
audit, and receive a benchmark-readiness report with concrete missing-data
repair actions.

The files are:

- `src/kira/amr/data_return.py`
- `examples/rwanda_amr_ast_template.csv`
- `examples/rwanda_amr_ast_example.csv`
- `tests/test_amr_data_return.py`

The example and template rows are synthetic. They contain no real patient data
and use generic fake facility identifiers.

## Why this kit exists

The Rwanda AMR scout identified AST completeness for priority bacterial isolate
records as the first data-return target:

```text
rwanda-amr-ast-completeness-priority-isolates
```

The Physics Auditor ticket gate then checks whether that ticket is measurable:
readout, units, scale, timing, observable returned data, falsifiability, and
domain-specific AMR boundaries. The AST completeness audit applies the narrow
row-level rule that decides whether a returned row has enough non-empty AST
metadata to become a benchmark-ready row.

This kit is the collaborator-facing layer around that audit. It does not
rewrite the audit logic. It handles CSV loading, template writing, report
rendering, and report file output.

## Required isolate-level columns

Each isolate-level AST row must include non-empty values for:

| Column | Meaning |
|---|---|
| `facility_id` | Generic facility identifier or local facility code. |
| `reporting_period` | Bounded period such as facility-month. |
| `organism` | Organism identification. |
| `specimen_source` | Specimen source or specimen type. |
| `antibiotic` | Tested antibiotic. |
| `ast_method` | AST method used for the isolate-antibiotic result. |
| `ast_result` | Isolate-antibiotic result category or reported result. |
| `breakpoint_version` | Breakpoint table or version used for interpretation. |
| `qc_status` | AST quality-control status for the row. |

Optional columns, such as `isolate_id`, can be included. The audit ignores
extra columns.

## Aggregate count guidance

Facility-period count rows should use `record_type=aggregate` and include:

| Column | Meaning |
|---|---|
| `facility_id` | Generic facility identifier or local facility code. |
| `reporting_period` | Bounded period such as facility-month. |
| `organism` | Organism boundary for the count. |
| `specimen_source` | Specimen source or specimen type boundary. |
| `antibiotic` | Tested antibiotic. |
| `tested_count` | Count of tested isolates in the period. |
| `resistant_count` or `susceptible_count` | At least one result-side count. |
| `ast_method` | AST method used for the count row. |
| `breakpoint_version` | Breakpoint table or version. |
| `qc_status` | AST quality-control status for the count row. |

For aggregate rows, `tested_count=0`, `resistant_count=0`, or
`susceptible_count=0` are present values. Blank cells are missing.

## How to fill the template

Start from:

```text
examples/rwanda_amr_ast_template.csv
```

The template includes two synthetic placeholder rows:

- `record_type=isolate_level` for isolate-antibiotic AST rows.
- `record_type=aggregate` for facility-period count rows.

Replace the placeholder cells with local data before running the audit. Keep
generic facility codes if the file is shared outside the source data system.
Do not add patient identifiers. The audit only needs facility, period, organism,
specimen, antibiotic, AST method/result or counts, breakpoint version, and QC
status.

## Python usage

```python
from kira.amr import (
    audit_amr_csv,
    benchmark_ready_records,
    load_amr_csv,
    make_markdown_report,
    report_to_dict,
    write_data_return_template,
    write_markdown_report,
)

write_data_return_template("rwanda_amr_ast_template.csv")

records = load_amr_csv("examples/rwanda_amr_ast_example.csv")
report = audit_amr_csv("examples/rwanda_amr_ast_example.csv")

ready_records = benchmark_ready_records(records)
payload = report_to_dict(report)
markdown = make_markdown_report(report)
write_markdown_report(report, "rwanda_amr_ast_readiness_report.md")
```

`load_amr_csv()` accepts a filesystem path or an already-open text stream. It
uses only Python standard library CSV handling.

## What the report means

`make_markdown_report(report)` produces a deterministic markdown report with:

- Source ticket id.
- Ticket gate status and check-section statuses.
- Total, complete, and incomplete record counts.
- Completeness percent.
- Missing field counts.
- Benchmark-ready row indices.
- Incomplete row indices and missing fields.
- Repair actions grouped by missing field.
- Non-claims.

The benchmark-ready rule is intentionally narrow:

```text
benchmark-ready = all required fields for the detected row type are present
                  and non-empty
```

The report says whether rows are ready for benchmark construction under that
schema. It does not say that a pathogen-antibiotic signal is meaningful, that a
facility performed better or worse than another, or that a data gap changed a
health outcome.

## Non-claims

This kit is only for data completeness and benchmark-readiness.

It does not:

- provide patient-care or medication-use direction
- assert that AMR is fixed
- assert completed external validation
- order facilities or sites by performance
- attribute health outcomes to a measurement gap or repair action
- assert collaboration with any named external institution

## Validation commands

Run the focused checks:

```bash
ruff check src/kira/amr tests/test_amr_data_return.py
pytest -q tests/test_amr_data_return.py
```

Run the related AMR and ticket-gate checks:

```bash
pytest -q tests/test_amr_audit.py tests/test_amr_scout.py tests/test_physics_ticket_gate.py tests/test_amr_data_return.py
```

Run the full project checks:

```bash
ruff check .
pytest -q
```

## Why this is Kira's first returned-data readiness loop

The earlier AMR scout created the missing-measurement ticket. The Physics
Auditor gate checked whether that ticket was measurable. The AST completeness
audit turned the ticket into a row-level readiness rule. This kit closes the
loop by giving a collaborator an input template, a loader, an audit call, a
deterministic readiness report, and grouped repair actions.

That is the first complete Kira pattern where a scout ticket can ask for data,
returned rows can be checked, and missing fields can be converted into concrete
repair work before any downstream benchmark claim is made.
