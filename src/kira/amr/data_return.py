"""CSV data-return helpers for the Rwanda AMR AST completeness audit.

This module handles collaborator-facing CSV IO and markdown rendering around
the existing AST completeness audit. It does not change the benchmark-ready
rule; all completeness decisions come from :mod:`kira.amr.audit`.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping
from os import PathLike
from pathlib import Path
from typing import Any, TextIO

from kira.amr.audit import (
    ISOLATE_REQUIRED_FIELDS,
    audit_ast_completeness,
    data_repair_tickets,
    report_to_dict,
    summarize_ast_completeness,
)

PathOrTextStream = str | PathLike[str] | TextIO

AGGREGATE_CSV_COLUMNS = (
    "facility_id",
    "reporting_period",
    "organism",
    "specimen_source",
    "antibiotic",
    "tested_count",
    "resistant_count",
    "susceptible_count",
    "ast_method",
    "breakpoint_version",
    "qc_status",
)

TEMPLATE_COLUMNS = (
    "synthetic_data_notice",
    "record_type",
    "facility_id",
    "reporting_period",
    "organism",
    "specimen_source",
    "antibiotic",
    "ast_method",
    "ast_result",
    "breakpoint_version",
    "qc_status",
    "tested_count",
    "resistant_count",
    "susceptible_count",
    "isolate_id",
)


def required_isolate_columns() -> tuple[str, ...]:
    """Return required CSV columns for isolate-level AST rows."""

    return tuple(ISOLATE_REQUIRED_FIELDS)


def required_aggregate_columns() -> tuple[str, ...]:
    """Return CSV columns for aggregate facility-period count rows.

    Aggregate rows require `tested_count` plus at least one result-side count:
    `resistant_count` or `susceptible_count`.
    """

    return AGGREGATE_CSV_COLUMNS


def load_amr_csv(path_or_file: PathOrTextStream) -> tuple[dict[str, str], ...]:
    """Load AMR AST CSV rows from a filesystem path or text stream."""

    if _is_text_stream(path_or_file):
        return _read_csv_rows(path_or_file)

    with Path(path_or_file).open("r", encoding="utf-8-sig", newline="") as handle:
        return _read_csv_rows(handle)


def audit_amr_csv(path_or_file: PathOrTextStream):
    """Load an AMR AST CSV and run the existing completeness audit."""

    return audit_ast_completeness(load_amr_csv(path_or_file))


def make_markdown_report(report: Any) -> str:
    """Render a deterministic benchmark-readiness report in markdown."""

    payload = report_to_dict(report)
    lines = [
        "# Rwanda AMR AST Data-Return Readiness Report",
        "",
        "## Summary",
        "",
        summarize_ast_completeness(report),
        "",
        "| Field | Value |",
        "|---|---:|",
        f"| Source ticket id | `{payload['source_ticket_id']}` |",
        f"| Ticket gate status | `{payload['ticket_gate_status']}` |",
        f"| Total records | {payload['total_records']} |",
        f"| Complete records | {payload['complete_records']} |",
        f"| Incomplete records | {payload['incomplete_records']} |",
        f"| Completeness percent | {payload['completeness_percent']:.2f}% |",
        "",
    ]

    lines.extend(_ticket_gate_section(payload))
    lines.extend(_missing_field_counts_section(payload))
    lines.extend(_row_index_section("Benchmark-Ready Rows", payload["benchmark_ready_rows"]))
    lines.extend(_incomplete_rows_section(payload))
    lines.extend(_repair_actions_section(report))
    lines.extend(_non_claims_section())
    return "\n".join(lines).rstrip() + "\n"


def make_data_return_template_rows() -> tuple[dict[str, str], ...]:
    """Return synthetic placeholder rows for the data-return CSV template."""

    empty_row = {column: "" for column in TEMPLATE_COLUMNS}
    isolate_row = {
        **empty_row,
        "synthetic_data_notice": "SYNTHETIC TEMPLATE ROW - replace before use",
        "record_type": "isolate_level",
    }
    aggregate_row = {
        **empty_row,
        "synthetic_data_notice": "SYNTHETIC TEMPLATE ROW - replace before use",
        "record_type": "aggregate",
    }
    return (isolate_row, aggregate_row)


def write_data_return_template(path: str | PathLike[str]) -> Path:
    """Write the synthetic AMR AST data-return CSV template."""

    output_path = Path(path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TEMPLATE_COLUMNS)
        writer.writeheader()
        writer.writerows(make_data_return_template_rows())
    return output_path


def write_markdown_report(report: Any, path: str | PathLike[str]) -> Path:
    """Write a markdown report for an AMR completeness audit."""

    output_path = Path(path)
    output_path.write_text(make_markdown_report(report), encoding="utf-8")
    return output_path


def _is_text_stream(value: PathOrTextStream) -> bool:
    return callable(getattr(value, "read", None))


def _read_csv_rows(handle: TextIO) -> tuple[dict[str, str], ...]:
    reader = csv.DictReader(handle)
    if reader.fieldnames is None:
        return ()
    reader.fieldnames = [_normalize_header(field) for field in reader.fieldnames]
    return tuple(_normalize_row(row) for row in reader)


def _normalize_header(field: str | None) -> str:
    if field is None:
        return ""
    return field.strip()


def _normalize_row(row: Mapping[str | None, str | None]) -> dict[str, str]:
    normalized: dict[str, str] = {}
    for key, value in row.items():
        if key is None:
            continue
        normalized[_normalize_header(key)] = "" if value is None else value.strip()
    return normalized


def _ticket_gate_section(payload: Mapping[str, Any]) -> list[str]:
    sections = payload["ticket_gate_check_sections"]
    lines = ["## Ticket Gate Sections", ""]
    if not sections:
        lines.extend(["No ticket-gate section details were returned.", ""])
        return lines

    lines.extend(["| Section | Status |", "|---|---|"])
    for section, status in sorted(sections.items()):
        lines.append(f"| {section} | `{status}` |")
    lines.append("")
    return lines


def _missing_field_counts_section(payload: Mapping[str, Any]) -> list[str]:
    counts = payload["missing_field_counts"]
    lines = ["## Missing Field Counts", ""]
    if not counts:
        lines.extend(["No missing required fields.", ""])
        return lines

    lines.extend(["| Missing field | Count |", "|---|---:|"])
    for field_name, count in sorted(counts.items()):
        lines.append(f"| `{field_name}` | {count} |")
    lines.append("")
    return lines


def _row_index_section(title: str, rows: list[int]) -> list[str]:
    lines = [f"## {title}", ""]
    if rows:
        lines.append(", ".join(str(row) for row in rows))
    else:
        lines.append("None.")
    lines.append("")
    return lines


def _incomplete_rows_section(payload: Mapping[str, Any]) -> list[str]:
    lines = ["## Incomplete Rows", ""]
    incomplete_rows = set(payload["incomplete_rows"])
    if not incomplete_rows:
        lines.extend(["None.", ""])
        return lines

    lines.extend(["| Row index | Record type | Missing fields |", "|---:|---|---|"])
    for record in payload["record_reports"]:
        if record["row_index"] in incomplete_rows:
            missing = ", ".join(f"`{field}`" for field in record["missing_fields"])
            lines.append(
                f"| {record['row_index']} | {record['record_type']} | {missing} |"
            )
    lines.append("")
    return lines


def _repair_actions_section(report: Any) -> list[str]:
    tickets = data_repair_tickets(report)
    lines = ["## Repair Actions By Missing Field", ""]
    if not tickets:
        lines.extend(["No missing-data repair actions.", ""])
        return lines

    for ticket in tickets:
        rows = ", ".join(str(row) for row in ticket["row_indices"])
        record_types = ", ".join(ticket["record_types"])
        lines.extend(
            [
                f"### `{ticket['field']}`",
                "",
                f"- Missing count: {ticket['missing_count']}",
                f"- Row indices: {rows}",
                f"- Record types: {record_types}",
                f"- Action: {ticket['action']}",
                "",
            ]
        )
    return lines


def _non_claims_section() -> list[str]:
    return [
        "## Non-Claims",
        "",
        "- Data-completeness and benchmark-readiness check only.",
        "- No patient-care or medication-use direction.",
        "- No assertion that AMR is fixed.",
        "- No assertion that external lab validation is complete.",
        "- No facility or site performance ordering.",
        "- No health-outcome attribution.",
        "- No named-institution collaboration assertion.",
        "",
    ]


__all__ = [
    "AGGREGATE_CSV_COLUMNS",
    "TEMPLATE_COLUMNS",
    "audit_amr_csv",
    "load_amr_csv",
    "make_data_return_template_rows",
    "make_markdown_report",
    "required_aggregate_columns",
    "required_isolate_columns",
    "write_data_return_template",
    "write_markdown_report",
]
