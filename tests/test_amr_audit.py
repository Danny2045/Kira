from __future__ import annotations

import json
from dataclasses import dataclass

from kira.amr import (
    AST_COMPLETENESS_SCOUT_TICKET_ID,
    audit_ast_completeness,
    audit_record,
    benchmark_ready_records,
    data_repair_tickets,
    report_to_dict,
    summarize_ast_completeness,
)

FORBIDDEN_OUTPUT_TERMS = (
    "clinical recommendation",
    "prescribe",
    "solved AMR",
    "wet-lab validated",
    "clinical cure",
    "partnership with Biohub",
    "partnership with Arc",
    "partnership with Ginkgo",
)


def _complete_isolate(**overrides):
    record = {
        "facility_id": "rwa-facility-001",
        "reporting_period": "2026-04",
        "organism": "Escherichia coli",
        "specimen_source": "urine",
        "antibiotic": "ceftriaxone",
        "ast_method": "disk diffusion",
        "ast_result": "resistant",
        "breakpoint_version": "CLSI M100 2026",
        "qc_status": "pass",
        "isolate_id": "iso-001",
    }
    record.update(overrides)
    return record


def _complete_aggregate(**overrides):
    record = {
        "record_type": "facility_period",
        "facility_id": "rwa-facility-001",
        "reporting_period": "2026-04",
        "organism": "Klebsiella pneumoniae",
        "specimen_source": "blood",
        "antibiotic": "meropenem",
        "tested_count": 12,
        "resistant_count": 3,
        "ast_method": "broth microdilution",
        "breakpoint_version": "CLSI M100 2026",
        "qc_status": "pass",
    }
    record.update(overrides)
    return record


def test_complete_isolate_level_record_is_benchmark_ready() -> None:
    record = _complete_isolate()
    audited_record = audit_record(record)
    report = audit_ast_completeness([record])

    assert audited_record.benchmark_ready is True
    assert audited_record.record_type == "isolate_level"
    assert audited_record.missing_fields == ()
    assert report.total_records == 1
    assert report.complete_records == 1
    assert report.incomplete_records == 0
    assert report.completeness_percent == 100.0
    assert report.benchmark_ready_count == 1
    assert report.non_benchmark_ready_count == 0
    assert report.benchmark_ready_rows == (0,)
    assert report.incomplete_rows == ()
    assert report.source_ticket_id == AST_COMPLETENESS_SCOUT_TICKET_ID
    assert report.ticket_gate_status
    assert benchmark_ready_records([record]) == (record,)


def test_missing_fields_are_counted_correctly() -> None:
    records = [
        _complete_isolate(ast_result=None, qc_status=""),
        _complete_isolate(organism=None, ast_result=""),
        _complete_isolate(),
    ]

    report = audit_ast_completeness(records)

    assert report.total_records == 3
    assert report.complete_records == 1
    assert report.incomplete_records == 2
    assert report.completeness_percent == 33.33
    assert report.missing_field_counts == {
        "ast_result": 2,
        "organism": 1,
        "qc_status": 1,
    }
    assert report.benchmark_ready_rows == (2,)
    assert report.incomplete_rows == (0, 1)


def test_blank_strings_count_as_missing() -> None:
    audited_record = audit_record(_complete_isolate(organism="   ", ast_result="\t"))

    assert audited_record.benchmark_ready is False
    assert audited_record.missing_fields == ("organism", "ast_result")


def test_dataclass_like_records_are_supported() -> None:
    @dataclass(frozen=True)
    class AstRow:
        facility_id: str
        reporting_period: str
        organism: str
        specimen_source: str
        antibiotic: str
        ast_method: str
        ast_result: str
        breakpoint_version: str
        qc_status: str

    row = AstRow(
        facility_id="rwa-facility-002",
        reporting_period="2026-04",
        organism="Staphylococcus aureus",
        specimen_source="wound",
        antibiotic="oxacillin",
        ast_method="disk diffusion",
        ast_result="susceptible",
        breakpoint_version="CLSI M100 2026",
        qc_status="pass",
    )

    assert audit_record(row).benchmark_ready is True


def test_aggregate_count_records_can_be_benchmark_ready() -> None:
    resistant_count_record = _complete_aggregate(resistant_count=0)
    susceptible_count_record = _complete_aggregate(resistant_count=None, susceptible_count=9)
    report = audit_ast_completeness([resistant_count_record, susceptible_count_record])

    assert all(record.record_type == "aggregate" for record in report.record_reports)
    assert report.complete_records == 2
    assert report.benchmark_ready_count == 2
    assert report.missing_field_counts == {}
    assert benchmark_ready_records([resistant_count_record, susceptible_count_record]) == (
        resistant_count_record,
        susceptible_count_record,
    )


def test_incomplete_aggregate_records_are_not_benchmark_ready() -> None:
    record = _complete_aggregate(tested_count="", resistant_count=None)
    audited_record = audit_record(record)
    report = audit_ast_completeness([record])

    assert audited_record.record_type == "aggregate"
    assert audited_record.benchmark_ready is False
    assert audited_record.missing_fields == (
        "tested_count",
        "resistant_count_or_susceptible_count",
    )
    assert report.benchmark_ready_count == 0
    assert report.non_benchmark_ready_count == 1
    assert report.missing_field_counts == {
        "resistant_count_or_susceptible_count": 1,
        "tested_count": 1,
    }


def test_report_to_dict_is_deterministic_and_json_serializable() -> None:
    report = audit_ast_completeness(
        [
            _complete_isolate(),
            _complete_isolate(ast_method=" "),
            _complete_aggregate(resistant_count=None, susceptible_count=7),
        ]
    )
    payload = report_to_dict(report)

    assert payload == report_to_dict(report)
    json.dumps(payload, sort_keys=True)
    assert set(payload) == {
        "source_ticket_id",
        "ticket_gate_status",
        "ticket_gate_check_sections",
        "total_records",
        "complete_records",
        "incomplete_records",
        "completeness_percent",
        "missing_field_counts",
        "benchmark_ready_count",
        "non_benchmark_ready_count",
        "benchmark_ready_rows",
        "incomplete_rows",
        "record_reports",
        "findings",
        "repair_tickets",
    }
    assert payload["source_ticket_id"] == AST_COMPLETENESS_SCOUT_TICKET_ID
    assert payload["record_reports"][1]["missing_fields"] == ["ast_method"]
    assert payload["repair_tickets"][0]["field"] == "ast_method"


def test_data_repair_tickets_return_concrete_missing_data_actions() -> None:
    report = audit_ast_completeness([_complete_isolate(ast_result="", qc_status="")])
    tickets = data_repair_tickets(report)
    tickets_by_field = {ticket["field"]: ticket for ticket in tickets}

    assert set(tickets_by_field) == {"ast_result", "qc_status"}
    assert tickets_by_field["ast_result"]["source_ticket_id"] == AST_COMPLETENESS_SCOUT_TICKET_ID
    assert tickets_by_field["ast_result"]["missing_count"] == 1
    assert tickets_by_field["ast_result"]["row_indices"] == [0]
    assert tickets_by_field["ast_result"]["action"].startswith("Backfill")
    assert "AST result" in tickets_by_field["ast_result"]["action"]
    json.dumps(tickets, sort_keys=True)


def test_audit_outputs_do_not_contain_forbidden_claims() -> None:
    report = audit_ast_completeness(
        [
            _complete_isolate(ast_result=""),
            _complete_aggregate(tested_count="", resistant_count=None),
        ]
    )
    output_text = json.dumps(
        {
            "report": report_to_dict(report),
            "summary": summarize_ast_completeness(report),
            "repair_tickets": data_repair_tickets(report),
        },
        sort_keys=True,
    ).lower()

    for claim in FORBIDDEN_OUTPUT_TERMS:
        assert claim.lower() not in output_text
