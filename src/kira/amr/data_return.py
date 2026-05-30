"""CSV data-return helpers for the Rwanda AMR AST completeness audit.

This module handles collaborator-facing CSV IO and markdown rendering around
the existing AST completeness audit. It does not change the benchmark-ready
rule; all completeness decisions come from :mod:`kira.amr.audit`.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
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

SYNTHETIC_NOTICE_COLUMN = "synthetic_data_notice"

DATA_STATUS_SYNTHETIC = "SYNTHETIC"
DATA_STATUS_REAL = "REAL"
DATA_STATUS_EMPTY = "EMPTY"

SYNTHETIC_BANNER = "# ⚠ SYNTHETIC DATA — NOT A REAL AMR FINDING — tooling demonstration only"
_EMPTY_DATA_NOTE = (
    "> No data audited: the input contained no AMR rows; "
    "this report implies no real AMR finding."
)


@dataclass(frozen=True, slots=True)
class AmrProvenance:
    """Synthetic-vs-real provenance for one AMR input, tied to its exact bytes.

    ``content_sha256`` is the SHA-256 of the raw CSV bytes, so the same input
    always yields a byte-identical report. No report-generation timestamp is
    stored: that belongs in a log line or filename, not the reproducible artifact.
    """

    data_status: str
    source: str
    n_rows: int
    content_sha256: str


def classify_data_status(rows: Sequence[Mapping[str, Any]]) -> str:
    """Classify rows as SYNTHETIC, REAL, or EMPTY from `synthetic_data_notice`.

    This is the single source of truth for provenance status. It looks only at
    the ``synthetic_data_notice`` column VALUES (never at facility IDs or other
    fields). Rows are partitioned into populated (notice non-empty after strip)
    versus empty:

    - all rows populated with one notice value -> ``SYNTHETIC``
    - all rows empty -> ``REAL``
    - no rows -> ``EMPTY``
    - some populated and some empty -> raises ``ValueError`` (mixed-provenance)
    - all populated but with differing notice values -> raises ``ValueError``
      (inconsistent synthetic-status column)
    """

    populated = [row for row in rows if _notice_value(row)]
    empty = [row for row in rows if not _notice_value(row)]

    if not rows:
        return DATA_STATUS_EMPTY
    if populated and empty:
        raise ValueError(
            f"mixed-provenance file: {len(populated)} rows marked synthetic, "
            f"{len(empty)} rows unmarked — separate synthetic and real rows into "
            "different files before auditing"
        )
    if populated:
        distinct = sorted({_notice_value(row) for row in populated})
        if len(distinct) > 1:
            values = ", ".join(repr(value) for value in distinct)
            raise ValueError(
                f"inconsistent synthetic_data_notice column: rows carry differing "
                f"notice values [{values}] — the synthetic-status column must be "
                "uniform; split differing values into separate files before auditing"
            )
        return DATA_STATUS_SYNTHETIC
    return DATA_STATUS_REAL


def _notice_value(row: Mapping[str, Any]) -> str:
    return str(row.get(SYNTHETIC_NOTICE_COLUMN, "") or "").strip()


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
    """Load AMR AST CSV rows from a filesystem path or text stream.

    This low-level reader returns rows only and does not classify provenance.
    Use :func:`load_amr_csv_with_provenance` (or :func:`audit_amr_csv` /
    :func:`make_amr_csv_report`) when provenance status and the mixed-file
    refusal are required.
    """

    _source, text, _raw = _read_source_text(path_or_file)
    return _parse_csv_rows(text)


def load_amr_csv_with_provenance(
    path_or_file: PathOrTextStream,
) -> tuple[tuple[dict[str, str], ...], AmrProvenance]:
    """Load AMR AST CSV rows and their provenance from a path or text stream.

    Reads the input exactly once (so streams are safe) and records the source
    identifier and a SHA-256 of the raw CSV bytes. Classifies the rows via
    :func:`classify_data_status`, which raises ``ValueError`` on a
    mixed-provenance or inconsistent-notice file.
    """

    source, text, raw = _read_source_text(path_or_file)
    rows = _parse_csv_rows(text)
    data_status = classify_data_status(rows)
    provenance = AmrProvenance(
        data_status=data_status,
        source=source,
        n_rows=len(rows),
        content_sha256=hashlib.sha256(raw).hexdigest(),
    )
    return rows, provenance


def audit_amr_csv(path_or_file: PathOrTextStream):
    """Load an AMR AST CSV and run the existing completeness audit.

    Routes through :func:`load_amr_csv_with_provenance`, so a mixed-provenance
    or inconsistent-notice file raises ``ValueError`` rather than being audited.
    """

    rows, _provenance = load_amr_csv_with_provenance(path_or_file)
    return audit_ast_completeness(rows)


def make_amr_csv_report(path_or_file: PathOrTextStream) -> str:
    """Load an AMR AST CSV and render its provenance-stamped markdown report.

    Raises ``ValueError`` on a mixed-provenance or inconsistent-notice file. The
    report leads with a synthetic-data banner for SYNTHETIC inputs and always
    carries a Data Provenance block tied to the input by SHA-256.
    """

    rows, provenance = load_amr_csv_with_provenance(path_or_file)
    report = audit_ast_completeness(rows)
    return make_markdown_report(report, provenance)


def make_markdown_report(report: Any, provenance: AmrProvenance | None = None) -> str:
    """Render a deterministic benchmark-readiness report in markdown.

    When ``provenance`` is supplied, a Data Provenance block is inserted below
    the title and (for SYNTHETIC inputs only) a synthetic-data banner is placed
    above it. These additions are the only difference: every numeric audit
    section from ``## Summary`` onward is byte-identical with or without
    provenance.

    For CSV-sourced data, use :func:`make_amr_csv_report` (or
    :func:`write_amr_csv_report`) instead: it classifies the file's provenance,
    stamps it into the report, and refuses mixed-provenance files. Calling this
    function bare (``provenance=None``) produces an UNMARKED report and is only
    appropriate for already-in-memory, non-CSV records whose provenance you have
    established by other means.
    """

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

    if provenance is not None:
        lines = _inject_provenance(lines, provenance)
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
    """Write a markdown report for an AMR completeness audit.

    For CSV-sourced data use :func:`write_amr_csv_report` instead: this helper
    renders the bare, UNMARKED report (no provenance block or synthetic banner)
    and is only appropriate for already-in-memory, non-CSV records.
    """

    output_path = Path(path)
    output_path.write_text(make_markdown_report(report), encoding="utf-8")
    return output_path


def write_amr_csv_report(path_or_file: PathOrTextStream, path: str | PathLike[str]) -> Path:
    """Write a provenance-stamped markdown report rendered from an AMR AST CSV.

    Raises ``ValueError`` on a mixed-provenance or inconsistent-notice input.
    """

    output_path = Path(path)
    output_path.write_text(make_amr_csv_report(path_or_file), encoding="utf-8")
    return output_path


def _inject_provenance(body_lines: list[str], provenance: AmrProvenance) -> list[str]:
    """Splice the banner (SYNTHETIC only) and Data Provenance block into a report.

    The title stays first (after any banner); the provenance block goes between
    the title and ``## Summary`` so every numeric section below is untouched.
    """

    title, rest = body_lines[0], body_lines[1:]
    banner = (
        [SYNTHETIC_BANNER, ""] if provenance.data_status == DATA_STATUS_SYNTHETIC else []
    )
    block = ["", *_provenance_block_lines(provenance)]
    if provenance.data_status == DATA_STATUS_EMPTY:
        block.extend(["", _EMPTY_DATA_NOTE])
    return [*banner, title, *block, *rest]


def _provenance_block_lines(provenance: AmrProvenance) -> list[str]:
    return [
        "## Data Provenance",
        "",
        "| Field | Value |",
        "|---|---|",
        f"| Data status | `{provenance.data_status}` |",
        f"| Source | `{provenance.source}` |",
        f"| Rows audited | {provenance.n_rows} |",
        f"| Content sha256 | `{provenance.content_sha256}` |",
    ]


def _is_text_stream(value: PathOrTextStream) -> bool:
    return callable(getattr(value, "read", None))


def _read_source_text(path_or_file: PathOrTextStream) -> tuple[str, str, bytes]:
    """Return ``(source, text, raw_bytes)`` from a path or stream in one read.

    ``raw_bytes`` is what gets hashed for provenance: the file's bytes for a
    path, or the stream text encoded as UTF-8 for a stream.
    """

    if _is_text_stream(path_or_file):
        text = path_or_file.read()
        return "<stream>", text, text.encode("utf-8")

    path = Path(path_or_file)
    raw_bytes = path.read_bytes()
    return str(path), raw_bytes.decode("utf-8-sig"), raw_bytes


def _parse_csv_rows(text: str) -> tuple[dict[str, str], ...]:
    return _read_csv_rows(StringIO(text))


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
    "DATA_STATUS_EMPTY",
    "DATA_STATUS_REAL",
    "DATA_STATUS_SYNTHETIC",
    "SYNTHETIC_BANNER",
    "SYNTHETIC_NOTICE_COLUMN",
    "TEMPLATE_COLUMNS",
    "AmrProvenance",
    "audit_amr_csv",
    "classify_data_status",
    "load_amr_csv",
    "load_amr_csv_with_provenance",
    "make_amr_csv_report",
    "make_data_return_template_rows",
    "make_markdown_report",
    "required_aggregate_columns",
    "required_isolate_columns",
    "write_amr_csv_report",
    "write_data_return_template",
    "write_markdown_report",
]
