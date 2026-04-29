"""Domain-agnostic contrast core primitives for Kira."""

from kira.contrast.schemas import (
    CONTRAST_SPEC_REQUIRED_FIELDS,
    DATA_RETURN_REQUIRED_FIELDS,
    ContrastSpec,
    DataReturnSchema,
    EvidenceRecord,
    EvidenceStatus,
    ExperimentTicket,
    make_data_return_schema,
    validate_contrast_spec,
)
from kira.contrast.tickets import (
    EXPERIMENT_TICKET_REQUIRED_FIELDS,
    ticket_from_dict,
    ticket_to_dict,
    validate_experiment_ticket,
)

__all__ = [
    "CONTRAST_SPEC_REQUIRED_FIELDS",
    "DATA_RETURN_REQUIRED_FIELDS",
    "EXPERIMENT_TICKET_REQUIRED_FIELDS",
    "ContrastSpec",
    "DataReturnSchema",
    "EvidenceRecord",
    "EvidenceStatus",
    "ExperimentTicket",
    "make_data_return_schema",
    "ticket_from_dict",
    "ticket_to_dict",
    "validate_contrast_spec",
    "validate_experiment_ticket",
]
