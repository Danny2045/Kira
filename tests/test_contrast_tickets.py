from __future__ import annotations

import pytest

from kira.contrast.schemas import EvidenceStatus, make_data_return_schema
from kira.contrast.tickets import (
    ExperimentTicket,
    ticket_from_dict,
    ticket_to_dict,
    validate_experiment_ticket,
)


def _ticket_payload(**overrides):
    payload = {
        "contrast_id": "ntd-selectivity::compound-a::target-pair-1",
        "domain": "ntd_selectivity",
        "intervention_type": "small_molecule",
        "intervention_id": "compound-a",
        "desired_context": "parasite target inhibition",
        "control_context": "human target safety comparator",
        "readout_type": "IC50 ratio human_div_parasite",
        "readout_units": "fold",
        "evidence_status": EvidenceStatus.SINGLE_SIDE_ONLY,
        "known_side": "parasite",
        "missing_side": "human",
        "label": "benchmark_repair_candidate",
        "uncertainty_note": "Human comparator measurement is missing.",
        "provenance": "v6 lab-campaign ticket fixture",
        "experiment_question": (
            "Measure compound-a against the missing human comparator using an IC50 assay."
        ),
        "expected_benchmark_impact": (
            "Completes one parasite-vs-human contrast needed for selectivity evidence."
        ),
    }
    payload.update(overrides)
    return payload


def test_ticket_round_trips_to_and_from_dict() -> None:
    ticket = ExperimentTicket(**_ticket_payload())

    payload = ticket_to_dict(ticket)
    restored = ticket_from_dict(payload)

    assert payload["evidence_status"] == "single_side_only"
    assert restored == ticket


def test_missing_required_ticket_fields_fail_validation() -> None:
    payload = _ticket_payload()
    del payload["experiment_question"]

    with pytest.raises(ValueError, match="experiment_question"):
        validate_experiment_ticket(payload)


def test_ticket_data_return_schema_contains_required_fields() -> None:
    ticket = ExperimentTicket(**_ticket_payload())

    schema = make_data_return_schema(ticket)

    assert schema.experiment_question == ticket.experiment_question
    assert schema.expected_benchmark_impact == ticket.expected_benchmark_impact
    assert "contrast_id" in schema.required_fields
    assert "expected_benchmark_impact" in schema.required_fields
