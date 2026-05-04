from __future__ import annotations

import json

from kira.amr import build_experiment_tickets as build_amr_tickets
from kira.physics.ticket_gate import (
    REQUIRED_CHECK_SECTIONS,
    TicketGateStatus,
    audit_ticket,
    audit_tickets,
    report_to_dict,
    summarize_gate_reports,
)
from kira.regeneration import build_experiment_tickets as build_regeneration_tickets

FORBIDDEN_REPORT_TERMS = (
    "clinical recommendation",
    "prescribe",
    "solved AMR",
    "solved regeneration",
    "wet-lab validated",
    "clinical cure",
)


def _scout_tickets():
    return (*build_amr_tickets(), *build_regeneration_tickets())


def _bad_ticket() -> dict[str, str]:
    return {
        "contrast_id": "bad-ticket-gate-fixture",
        "domain": "rwanda_amr",
        "intervention_type": "platform",
        "intervention_id": "platform::bad",
        "desired_context": "better state",
        "control_context": "baseline state",
        "readout_type": "vibes",
        "readout_units": "",
        "evidence_status": "proposed",
        "known_side": "none",
        "missing_side": "both",
        "label": "bad_fixture",
        "uncertainty_note": "solved AMR",
        "provenance": "unit-test fixture",
        "experiment_question": "See if it works.",
        "expected_benchmark_impact": "Big impact.",
    }


def test_ticket_gate_audits_current_scout_tickets_and_report_shape() -> None:
    tickets = _scout_tickets()
    reports = audit_tickets(tickets)

    assert len(tickets) == 14
    assert len(reports) == 14
    assert {report.domain for report in reports} == {"rwanda_amr", "regeneration_crispr"}
    assert any(report.status == TicketGateStatus.NEEDS_MEASUREMENT_DETAIL for report in reports)
    assert all(report.status != TicketGateStatus.BLOCKED for report in reports)

    for report in reports:
        assert report.contrast_id
        assert report.domain
        assert report.status in set(TicketGateStatus)
        assert isinstance(report.findings, tuple)
        assert set(report.check_sections) == set(REQUIRED_CHECK_SECTIONS)
        assert set(report_to_dict(report)["check_sections"]) == set(REQUIRED_CHECK_SECTIONS)
        assert report.inferred_scale in {
            "molecular",
            "cellular",
            "organoid",
            "tissue",
            "organism",
            "facility",
            "surveillance",
            "population",
        }
        assert report.observable_data_kinds


def test_ticket_gate_summary_counts_status_and_domain() -> None:
    reports = audit_tickets(_scout_tickets())
    summary = summarize_gate_reports(reports)

    assert summary["total"] == 14
    assert summary["by_status"]["pass"] + summary["by_status"]["needs_measurement_detail"] == 14
    assert summary["by_status"]["blocked"] == 0
    assert summary["by_domain"]["rwanda_amr"]["total"] == 8
    assert summary["by_domain"]["regeneration_crispr"]["total"] == 6


def test_ticket_gate_blocks_bad_vague_ticket_dict() -> None:
    report = audit_ticket(_bad_ticket())
    payload = report_to_dict(report)

    assert report.status == TicketGateStatus.BLOCKED
    assert payload["status"] == "blocked"
    assert payload["contrast_id"] == "bad-ticket-gate-fixture"
    assert payload["check_sections"]["required_fields"] == "blocked"
    assert payload["check_sections"]["observability"] == "blocked"
    assert payload["check_sections"]["falsifiability"] == "blocked"
    assert payload["check_sections"]["non_claims"] == "blocked"


def test_report_to_dict_is_stable_and_claim_redacted() -> None:
    reports = audit_tickets(_scout_tickets())
    payloads = [report_to_dict(report) for report in reports]
    bad_payload = report_to_dict(audit_ticket(_bad_ticket()))

    assert payloads == [report_to_dict(report) for report in reports]
    json.dumps(payloads + [bad_payload], sort_keys=True)

    first = payloads[0]
    assert set(first) == {
        "contrast_id",
        "domain",
        "status",
        "check_sections",
        "inferred_scale",
        "observable_data_kinds",
        "findings",
    }

    report_text = json.dumps(payloads + [bad_payload], sort_keys=True).lower()
    for claim in FORBIDDEN_REPORT_TERMS:
        assert claim.lower() not in report_text
