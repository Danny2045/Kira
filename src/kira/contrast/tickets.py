"""Experiment ticket helpers for the contrast core."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
from typing import Any

from kira.contrast.schemas import (
    DATA_RETURN_REQUIRED_FIELDS,
    ExperimentTicket,
    evidence_status_to_string,
    normalize_evidence_status,
)

EXPERIMENT_TICKET_REQUIRED_FIELDS = DATA_RETURN_REQUIRED_FIELDS


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _missing_required_fields(values: Mapping[str, Any]) -> list[str]:
    return [field for field in EXPERIMENT_TICKET_REQUIRED_FIELDS if not _is_present(values.get(field))]


def _as_mapping(ticket: ExperimentTicket | Mapping[str, Any]) -> Mapping[str, Any]:
    if isinstance(ticket, ExperimentTicket):
        return asdict(ticket)
    if isinstance(ticket, Mapping):
        return ticket
    raise TypeError("experiment ticket must be an ExperimentTicket or mapping")


def validate_experiment_ticket(ticket: ExperimentTicket | Mapping[str, Any]) -> ExperimentTicket:
    """Validate required ticket fields and return a normalized ExperimentTicket."""

    values = _as_mapping(ticket)
    missing = _missing_required_fields(values)
    if missing:
        raise ValueError(f"experiment ticket is missing required fields: {missing}")

    payload = {
        field: values[field]
        for field in EXPERIMENT_TICKET_REQUIRED_FIELDS
    }
    payload["evidence_status"] = normalize_evidence_status(payload["evidence_status"])
    return ExperimentTicket(**payload)


def ticket_to_dict(ticket: ExperimentTicket | Mapping[str, Any]) -> dict[str, Any]:
    """Return a JSON-ready dictionary for an experiment ticket."""

    validated = validate_experiment_ticket(ticket)
    payload = asdict(validated)
    payload["evidence_status"] = evidence_status_to_string(validated.evidence_status)
    return payload


def ticket_from_dict(payload: Mapping[str, Any]) -> ExperimentTicket:
    """Build an ExperimentTicket from a serialized dictionary."""

    return validate_experiment_ticket(payload)
