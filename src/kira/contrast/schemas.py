"""Lightweight schemas for domain-agnostic biological contrasts.

A contrast compares a desired biological context against a control, safety, or
failure context under a measurable readout. These structures intentionally avoid
chemistry, model, and assay-runtime dependencies so that future adapters can map
their own domain records onto the same small contract.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Any


class EvidenceStatus(str, Enum):
    """Coarse evidence states shared by contrast domains."""

    PROPOSED = "proposed"
    SINGLE_SIDE_ONLY = "single_side_only"
    PAIRED_COMPLETE = "paired_complete"
    BOUNDED = "bounded"
    CONFLICTING = "conflicting"
    UNKNOWN = "unknown"


CONTRAST_SPEC_REQUIRED_FIELDS = (
    "contrast_id",
    "domain",
    "intervention_type",
    "intervention_id",
    "desired_context",
    "control_context",
    "readout_type",
    "readout_units",
    "evidence_status",
    "known_side",
    "missing_side",
    "label",
    "uncertainty_note",
    "provenance",
)

EXPERIMENT_CONTEXT_FIELDS = (
    "experiment_question",
    "expected_benchmark_impact",
)

DATA_RETURN_REQUIRED_FIELDS = (
    *CONTRAST_SPEC_REQUIRED_FIELDS,
    *EXPERIMENT_CONTEXT_FIELDS,
)


@dataclass(frozen=True)
class ContrastSpec:
    """Domain-neutral definition of a measurable contrast."""

    contrast_id: str
    domain: str
    intervention_type: str
    intervention_id: str
    desired_context: str
    control_context: str
    readout_type: str
    readout_units: str
    evidence_status: EvidenceStatus | str
    known_side: str
    missing_side: str
    label: str
    uncertainty_note: str
    provenance: str


@dataclass(frozen=True)
class EvidenceRecord:
    """A recorded or imported evidence row for a contrast."""

    contrast_id: str
    domain: str
    intervention_type: str
    intervention_id: str
    desired_context: str
    control_context: str
    readout_type: str
    readout_units: str
    evidence_status: EvidenceStatus | str
    known_side: str
    missing_side: str
    label: str
    uncertainty_note: str
    provenance: str


@dataclass(frozen=True)
class DataReturnSchema:
    """Expected fields for data returned by an executable experiment ticket."""

    contrast_id: str
    domain: str
    intervention_type: str
    intervention_id: str
    desired_context: str
    control_context: str
    readout_type: str
    readout_units: str
    evidence_status: EvidenceStatus | str
    known_side: str
    missing_side: str
    label: str
    uncertainty_note: str
    provenance: str
    experiment_question: str
    expected_benchmark_impact: str
    required_fields: tuple[str, ...] = DATA_RETURN_REQUIRED_FIELDS


@dataclass(frozen=True)
class ExperimentTicket:
    """Executable question that can return evidence for one contrast."""

    contrast_id: str
    domain: str
    intervention_type: str
    intervention_id: str
    desired_context: str
    control_context: str
    readout_type: str
    readout_units: str
    evidence_status: EvidenceStatus | str
    known_side: str
    missing_side: str
    label: str
    uncertainty_note: str
    provenance: str
    experiment_question: str
    expected_benchmark_impact: str


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _missing_required_fields(values: Mapping[str, Any], required_fields: tuple[str, ...]) -> list[str]:
    return [field for field in required_fields if not _is_present(values.get(field))]


def _status_value(status: EvidenceStatus | str) -> str:
    if isinstance(status, EvidenceStatus):
        return status.value
    return str(status).strip()


def normalize_evidence_status(status: EvidenceStatus | str) -> EvidenceStatus | str:
    """Return a known EvidenceStatus enum when possible, otherwise a clean string."""

    if isinstance(status, EvidenceStatus):
        return status
    if not _is_present(status):
        raise ValueError("evidence_status is required")
    value = str(status).strip()
    try:
        return EvidenceStatus(value)
    except ValueError:
        return value


def _as_mapping(spec: ContrastSpec | Mapping[str, Any]) -> Mapping[str, Any]:
    if isinstance(spec, ContrastSpec):
        return {
            field: getattr(spec, field)
            for field in CONTRAST_SPEC_REQUIRED_FIELDS
        }
    if isinstance(spec, Mapping):
        return spec
    raise TypeError("contrast spec must be a ContrastSpec or mapping")


def validate_contrast_spec(spec: ContrastSpec | Mapping[str, Any]) -> ContrastSpec:
    """Validate required contrast fields and return a normalized ContrastSpec."""

    values = _as_mapping(spec)
    missing = _missing_required_fields(values, CONTRAST_SPEC_REQUIRED_FIELDS)
    if missing:
        raise ValueError(f"contrast spec is missing required fields: {missing}")

    payload = {
        field: values[field]
        for field in CONTRAST_SPEC_REQUIRED_FIELDS
    }
    payload["evidence_status"] = normalize_evidence_status(payload["evidence_status"])
    return ContrastSpec(**payload)


def _read_source_field(source: Any, field: str, default: str = "") -> Any:
    if isinstance(source, Mapping):
        return source.get(field, default)
    return getattr(source, field, default)


def make_data_return_schema(source: Any) -> DataReturnSchema:
    """Create the expected data-return schema from a contrast spec or ticket-like object."""

    values = {
        field: _read_source_field(source, field)
        for field in DATA_RETURN_REQUIRED_FIELDS
    }
    missing = _missing_required_fields(values, CONTRAST_SPEC_REQUIRED_FIELDS)
    if missing:
        raise ValueError(f"data return schema source is missing required fields: {missing}")

    values["evidence_status"] = normalize_evidence_status(values["evidence_status"])
    values["experiment_question"] = str(values.get("experiment_question") or "")
    values["expected_benchmark_impact"] = str(values.get("expected_benchmark_impact") or "")
    return DataReturnSchema(**values)


def evidence_status_to_string(status: EvidenceStatus | str) -> str:
    """Serialize evidence status values for JSON-ready ticket dictionaries."""

    return _status_value(status)
