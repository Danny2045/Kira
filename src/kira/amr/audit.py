"""AST completeness audit for Rwanda AMR benchmark-readiness.

The audit checks whether isolate-level or facility-period AMR records contain
the fields needed for a contrast or benchmark row. It is a data-quality check
only: it does not provide patient-care guidance, infer facility performance,
assert completed lab validation, or make public-health outcome claims.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, is_dataclass
from dataclasses import field as dataclass_field
from math import isnan
from typing import Any

from kira.amr.scout import build_experiment_tickets
from kira.physics.ticket_gate import audit_ticket
from kira.physics.ticket_gate import report_to_dict as gate_report_to_dict

AST_COMPLETENESS_SCOUT_TICKET_ID = "rwanda-amr-ast-completeness-priority-isolates"

ISOLATE_REQUIRED_FIELDS = (
    "facility_id",
    "reporting_period",
    "organism",
    "specimen_source",
    "antibiotic",
    "ast_method",
    "ast_result",
    "breakpoint_version",
    "qc_status",
)

AGGREGATE_BASE_REQUIRED_FIELDS = (
    "facility_id",
    "reporting_period",
    "organism",
    "specimen_source",
    "antibiotic",
    "tested_count",
    "ast_method",
    "breakpoint_version",
    "qc_status",
)

COUNT_RESULT_FIELD = "resistant_count_or_susceptible_count"
AGGREGATE_REQUIRED_FIELDS = (*AGGREGATE_BASE_REQUIRED_FIELDS, COUNT_RESULT_FIELD)

_COUNT_FIELDS = ("tested_count", "resistant_count", "susceptible_count")
_AGGREGATE_RECORD_TYPES = {
    "aggregate",
    "facility_period",
    "facility-period",
    "facility_period_count",
    "facility-period-count",
}
_ISOLATE_RECORD_TYPES = {"isolate", "isolate_level", "isolate-level"}


@dataclass(frozen=True, slots=True)
class MissingFieldFinding:
    """One missing field in one AMR audit row."""

    row_index: int
    record_type: str
    field: str
    source_ticket_id: str
    action: str


@dataclass(frozen=True, slots=True)
class AmrAuditRecord:
    """Normalized readiness result for one input record."""

    row_index: int
    record_type: str
    benchmark_ready: bool
    required_fields: tuple[str, ...]
    missing_fields: tuple[str, ...]
    source_ticket_id: str = AST_COMPLETENESS_SCOUT_TICKET_ID


@dataclass(frozen=True, slots=True)
class AmrCompletenessReport:
    """Aggregated AST completeness audit report."""

    total_records: int
    complete_records: int
    incomplete_records: int
    completeness_percent: float
    missing_field_counts: dict[str, int]
    benchmark_ready_count: int
    non_benchmark_ready_count: int
    benchmark_ready_rows: tuple[int, ...]
    incomplete_rows: tuple[int, ...]
    record_reports: tuple[AmrAuditRecord, ...]
    findings: tuple[MissingFieldFinding, ...]
    source_ticket_id: str = AST_COMPLETENESS_SCOUT_TICKET_ID
    ticket_gate_status: str = ""
    ticket_gate_check_sections: dict[str, str] = dataclass_field(default_factory=dict)


def audit_record(record: Mapping[str, Any] | Any, row_index: int = 0) -> AmrAuditRecord:
    """Audit one AMR record for AST completeness."""

    values = _record_to_mapping(record)
    record_type = _infer_record_type(values)
    required_fields = _required_fields(record_type)
    missing_fields = _missing_fields(values, record_type)
    return AmrAuditRecord(
        row_index=row_index,
        record_type=record_type,
        benchmark_ready=not missing_fields,
        required_fields=required_fields,
        missing_fields=missing_fields,
    )


def audit_ast_completeness(records: Iterable[Mapping[str, Any] | Any]) -> AmrCompletenessReport:
    """Audit a deterministic sequence of AMR records for benchmark-readiness."""

    record_reports = tuple(
        audit_record(record, row_index=row_index) for row_index, record in enumerate(records)
    )
    missing_field_counts = dict(
        sorted(Counter(field for report in record_reports for field in report.missing_fields).items())
    )
    findings = tuple(
        MissingFieldFinding(
            row_index=report.row_index,
            record_type=report.record_type,
            field=missing_field,
            source_ticket_id=AST_COMPLETENESS_SCOUT_TICKET_ID,
            action=_repair_action(missing_field),
        )
        for report in record_reports
        for missing_field in report.missing_fields
    )
    benchmark_ready_rows = tuple(
        report.row_index for report in record_reports if report.benchmark_ready
    )
    incomplete_rows = tuple(report.row_index for report in record_reports if not report.benchmark_ready)
    complete_records = len(benchmark_ready_rows)
    total_records = len(record_reports)
    ticket_gate_status, ticket_gate_sections = _ast_ticket_gate_summary()

    return AmrCompletenessReport(
        total_records=total_records,
        complete_records=complete_records,
        incomplete_records=total_records - complete_records,
        completeness_percent=_percent(complete_records, total_records),
        missing_field_counts=missing_field_counts,
        benchmark_ready_count=complete_records,
        non_benchmark_ready_count=total_records - complete_records,
        benchmark_ready_rows=benchmark_ready_rows,
        incomplete_rows=incomplete_rows,
        record_reports=record_reports,
        findings=findings,
        ticket_gate_status=ticket_gate_status,
        ticket_gate_check_sections=ticket_gate_sections,
    )


def report_to_dict(report: AmrCompletenessReport) -> dict[str, Any]:
    """Return a stable JSON-ready dictionary for an AST completeness report."""

    return {
        "source_ticket_id": report.source_ticket_id,
        "ticket_gate_status": report.ticket_gate_status,
        "ticket_gate_check_sections": dict(sorted(report.ticket_gate_check_sections.items())),
        "total_records": report.total_records,
        "complete_records": report.complete_records,
        "incomplete_records": report.incomplete_records,
        "completeness_percent": report.completeness_percent,
        "missing_field_counts": dict(sorted(report.missing_field_counts.items())),
        "benchmark_ready_count": report.benchmark_ready_count,
        "non_benchmark_ready_count": report.non_benchmark_ready_count,
        "benchmark_ready_rows": list(report.benchmark_ready_rows),
        "incomplete_rows": list(report.incomplete_rows),
        "record_reports": [
            {
                "row_index": record.row_index,
                "record_type": record.record_type,
                "benchmark_ready": record.benchmark_ready,
                "required_fields": list(record.required_fields),
                "missing_fields": list(record.missing_fields),
                "source_ticket_id": record.source_ticket_id,
            }
            for record in report.record_reports
        ],
        "findings": [
            {
                "row_index": finding.row_index,
                "record_type": finding.record_type,
                "field": finding.field,
                "source_ticket_id": finding.source_ticket_id,
                "action": finding.action,
            }
            for finding in report.findings
        ],
        "repair_tickets": list(data_repair_tickets(report)),
    }


def summarize_ast_completeness(report: AmrCompletenessReport) -> str:
    """Return a compact human-readable AST completeness summary."""

    if report.missing_field_counts:
        missing = ", ".join(
            f"{field}={count}" for field, count in sorted(report.missing_field_counts.items())
        )
    else:
        missing = "no missing required fields"
    return (
        f"AST completeness audit for {report.source_ticket_id}: "
        f"{report.benchmark_ready_count}/{report.total_records} records benchmark-ready "
        f"({report.completeness_percent:.2f}%). Missing required data: {missing}."
    )


def benchmark_ready_records(records: Iterable[Mapping[str, Any] | Any]) -> tuple[Mapping[str, Any] | Any, ...]:
    """Return the original records that pass the AST benchmark-ready rule."""

    record_tuple = tuple(records)
    ready_rows = set(audit_ast_completeness(record_tuple).benchmark_ready_rows)
    return tuple(record for row_index, record in enumerate(record_tuple) if row_index in ready_rows)


def data_repair_tickets(report: AmrCompletenessReport) -> tuple[dict[str, Any], ...]:
    """Group missing-data repair actions by field."""

    findings_by_field: dict[str, list[MissingFieldFinding]] = defaultdict(list)
    for finding in report.findings:
        findings_by_field[finding.field].append(finding)

    tickets: list[dict[str, Any]] = []
    for field_name in sorted(findings_by_field):
        field_findings = tuple(findings_by_field[field_name])
        tickets.append(
            {
                "ticket_id": f"{report.source_ticket_id}:missing-{field_name}",
                "source_ticket_id": report.source_ticket_id,
                "field": field_name,
                "missing_count": len(field_findings),
                "row_indices": sorted(finding.row_index for finding in field_findings),
                "record_types": sorted({finding.record_type for finding in field_findings}),
                "action": _repair_action(field_name),
            }
        )
    return tuple(tickets)


def _record_to_mapping(record: Mapping[str, Any] | Any) -> dict[str, Any]:
    if isinstance(record, Mapping):
        return {str(key): value for key, value in record.items()}
    if is_dataclass(record) and not isinstance(record, type):
        return {str(key): value for key, value in asdict(record).items()}
    to_dict = getattr(record, "to_dict", None)
    if callable(to_dict):
        payload = to_dict()
        if isinstance(payload, Mapping):
            return {str(key): value for key, value in payload.items()}
    if hasattr(record, "__dict__"):
        return {str(key): value for key, value in vars(record).items()}
    raise TypeError("AMR audit records must be mappings, dataclass instances, or to_dict objects")


def _infer_record_type(values: Mapping[str, Any]) -> str:
    explicit_type = str(values.get("record_type") or values.get("row_type") or "").strip().lower()
    if explicit_type in _AGGREGATE_RECORD_TYPES:
        return "aggregate"
    if explicit_type in _ISOLATE_RECORD_TYPES:
        return "isolate_level"
    if _is_present(values.get("isolate_id")):
        return "isolate_level"
    if any(field in values for field in _COUNT_FIELDS) and not _is_present(values.get("ast_result")):
        return "aggregate"
    return "isolate_level"


def _required_fields(record_type: str) -> tuple[str, ...]:
    if record_type == "aggregate":
        return AGGREGATE_REQUIRED_FIELDS
    return ISOLATE_REQUIRED_FIELDS


def _missing_fields(values: Mapping[str, Any], record_type: str) -> tuple[str, ...]:
    if record_type == "aggregate":
        missing = [
            field for field in AGGREGATE_BASE_REQUIRED_FIELDS if not _is_present(values.get(field))
        ]
        if not (
            _is_present(values.get("resistant_count"))
            or _is_present(values.get("susceptible_count"))
        ):
            missing.append(COUNT_RESULT_FIELD)
        return tuple(missing)
    return tuple(field for field in ISOLATE_REQUIRED_FIELDS if not _is_present(values.get(field)))


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, float) and isnan(value):
        return False
    if isinstance(value, Mapping):
        return bool(value)
    if isinstance(value, (list, tuple, set, frozenset)):
        return bool(value)
    return True


def _percent(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return round((numerator / denominator) * 100.0, 2)


def _ast_ticket_gate_summary() -> tuple[str, dict[str, str]]:
    for ticket in build_experiment_tickets():
        if ticket.contrast_id == AST_COMPLETENESS_SCOUT_TICKET_ID:
            payload = gate_report_to_dict(audit_ticket(ticket))
            return str(payload["status"]), {
                str(section): str(status)
                for section, status in payload["check_sections"].items()
            }
    return "unknown", {}


def _repair_action(field_name: str) -> str:
    actions = {
        "facility_id": "Backfill the reporting facility identifier for each listed row.",
        "reporting_period": "Backfill a bounded reporting period for each listed row.",
        "organism": "Backfill the organism identification for each listed row.",
        "specimen_source": "Backfill the specimen source or specimen type for each listed row.",
        "antibiotic": "Backfill the tested antibiotic name for each listed row.",
        "ast_method": "Backfill the AST method for each listed row.",
        "ast_result": "Backfill the isolate-antibiotic AST result for each listed row.",
        "breakpoint_version": "Backfill the breakpoint table or version for each listed row.",
        "qc_status": "Backfill the AST quality-control status for each listed row.",
        "tested_count": "Backfill the tested isolate count for each listed facility-period row.",
        COUNT_RESULT_FIELD: (
            "Backfill at least one result count, resistant_count or susceptible_count, "
            "for each listed facility-period row."
        ),
    }
    return actions.get(field_name, f"Backfill non-empty {field_name} values for each listed row.")


__all__ = [
    "AGGREGATE_REQUIRED_FIELDS",
    "AST_COMPLETENESS_SCOUT_TICKET_ID",
    "ISOLATE_REQUIRED_FIELDS",
    "AmrAuditRecord",
    "AmrCompletenessReport",
    "MissingFieldFinding",
    "audit_ast_completeness",
    "audit_record",
    "benchmark_ready_records",
    "data_repair_tickets",
    "report_to_dict",
    "summarize_ast_completeness",
]
