"""Physics-auditor gate for contrast-core experiment tickets.

The gate is intentionally lexical and deterministic. It does not score biology,
predict outcomes, or validate a wet-lab result. It checks whether a ticket is
measurable enough to deserve scientific-looking treatment downstream.
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from typing import Any

from kira.contrast.schemas import (
    DATA_RETURN_REQUIRED_FIELDS,
    DataReturnSchema,
    ExperimentTicket,
    evidence_status_to_string,
    make_data_return_schema,
)


class TicketGateStatus(str, Enum):
    """Coarse status assigned by the physics ticket gate."""

    PASS = "pass"
    NEEDS_MEASUREMENT_DETAIL = "needs_measurement_detail"
    BLOCKED = "blocked"


@dataclass(frozen=True, slots=True)
class TicketGateFinding:
    """One deterministic ticket-gate finding."""

    check: str
    severity: str
    code: str
    message: str
    fields: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class TicketGateReport:
    """Structured audit result for one experiment ticket."""

    contrast_id: str
    domain: str
    status: TicketGateStatus
    findings: tuple[TicketGateFinding, ...]
    check_sections: dict[str, TicketGateStatus]
    inferred_scale: str
    observable_data_kinds: tuple[str, ...]


REQUIRED_CHECK_SECTIONS = (
    "required_fields",
    "units",
    "scale",
    "timescale",
    "observability",
    "falsifiability",
    "data_return_schema",
    "domain_specific",
    "non_claims",
)

_TEXT_FIELDS = (
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
    "experiment_question",
    "expected_benchmark_impact",
    "missing_measurement",
    "measurement_gap",
    "rationale",
    "model_system",
    "recommended_assay_type",
    "recommended_standard_type",
    "recommended_units",
)

_FORBIDDEN_CLAIM_TERMS = (
    "clinical recommendation",
    "prescribe",
    "solved amr",
    "solved regeneration",
    "wet-lab validated",
    "clinical cure",
)

_CONCRETE_UNIT_TERMS = (
    "%",
    "percent",
    "mic",
    "nm",
    "ddd",
    "dot",
    "count",
    "rate",
    "days",
    "month",
    "detected",
    "not_detected",
    "fold",
    "ratio",
)

_VAGUE_UNIT_TERMS = (
    "assay-specific",
    "assay specific",
    "unitless",
    "qualitative",
    "semiquantitative",
    "tbd",
    "to be determined",
    "category",
    "score",
    "signal",
)

_TIMING_TERMS = (
    "time",
    "timing",
    "time-bounded",
    "time bounded",
    "time-window",
    "window",
    "period",
    "same-period",
    "reporting-period",
    "reporting period",
    "month",
    "day",
    "days",
    "duration",
    "course",
    "schedule",
    "after",
    "before",
    "transient",
    "repeatability",
    "bounded",
)

_OBSERVABLE_DATA_KINDS = (
    (
        "AST record",
        (
            "ast result",
            "ast record",
            "ast panel",
            "phenotypic ast",
            "ast method",
            "ast availability",
            "isolate-antibiotic result",
        ),
    ),
    (
        "resistance-rate table",
        (
            "resistance rate",
            "percent resistant",
            "resistant count",
            "susceptible count",
            "ast denominator",
            "pathogen-antibiotic",
        ),
    ),
    (
        "morphology score",
        (
            "morphology score",
            "morphology atlas",
            "head-shape labels",
            "wound-repair phenotype",
        ),
    ),
    (
        "marker panel",
        (
            "marker panel",
            "marker-panel",
            "marker call",
            "marker-call",
            "lineage-retention",
            "epigenetic-age",
        ),
    ),
    (
        "bioelectric map",
        (
            "membrane-voltage map",
            "voltage-pattern",
            "bioelectric pattern",
            "membrane-voltage",
        ),
    ),
    (
        "image/phenotype table",
        (
            "image-based",
            "phenotype",
            "migration index",
            "organoid phenotype",
            "functional-rescue",
            "functional rescue",
        ),
    ),
    (
        "facility-month data",
        (
            "facility-month",
            "facility month",
            "facility aggregate",
            "facility-level",
            "site-by-specimen",
            "ward exposure",
        ),
    ),
    (
        "sequence/marker call",
        (
            "sequence-quality",
            "sequence quality",
            "marker-call",
            "marker call",
            "genomic resistance marker",
            "linked isolate id",
        ),
    ),
    (
        "QC record",
        (
            "quality-control",
            "quality control",
            "breakpoint",
            "qc",
            "repeatability",
        ),
    ),
)

_READOUT_SCALE_HINTS = {
    "ast result": "organism",
    "resistance rate": "surveillance",
    "antibiotic consumption": "facility",
    "treatment outcome": "facility",
    "infection-control metric": "facility",
    "genomic resistance marker": "molecular",
    "surveillance completeness": "surveillance",
    "bioelectric pattern": "tissue",
    "functional rescue": "tissue",
    "marker panel delta": "cellular",
    "morphology score": "organism",
    "organoid phenotype": "organoid",
    "wound closure": "tissue",
    "ic50 ratio human_div_parasite": "molecular",
}

_SCALE_TERMS = (
    ("molecular", ("marker", "sequence", "mic", "compound", "target", "ic50", "guide")),
    ("cellular", ("cell", "single-cell", "cell-state", "lineage", "interneuron")),
    ("organoid", ("organoid", "assembloid")),
    ("tissue", ("tissue", "wound", "tail", "bioelectric", "morphology")),
    ("organism", ("organism", "pathogen", "bacterial", "isolate", "xenopus", "planarian")),
    ("facility", ("facility", "ward", "hospital", "site", "ipc")),
    ("surveillance", ("surveillance", "sentinel", "sampling", "reporting", "glass")),
    ("population", ("population", "public health", "one health")),
)


def audit_ticket(ticket: ExperimentTicket | Mapping[str, Any] | Any) -> TicketGateReport:
    """Audit one contrast-core experiment ticket and return a structured report."""

    values = _as_ticket_mapping(ticket)
    contrast_id = _display_value(values, "contrast_id", fallback="unknown")
    domain = _display_value(values, "domain", fallback="unknown")
    text = _ticket_text(values)
    lower_text = text.lower()

    sections: dict[str, list[TicketGateFinding]] = {name: [] for name in REQUIRED_CHECK_SECTIONS}
    inferred_scale = infer_measurement_scale(values)
    observable_data_kinds = infer_observable_data_kinds(values)

    sections["required_fields"].extend(_audit_required_fields(values))
    sections["units"].extend(_audit_units(values))
    sections["scale"].extend(_audit_scale(inferred_scale))
    sections["timescale"].extend(_audit_timescale(lower_text))
    sections["observability"].extend(_audit_observability(observable_data_kinds))
    sections["falsifiability"].extend(_audit_falsifiability(lower_text))
    sections["data_return_schema"].extend(_audit_data_return_schema(values))
    sections["domain_specific"].extend(_audit_domain_specific(values, lower_text))
    sections["non_claims"].extend(_audit_non_claims(lower_text))

    findings: list[TicketGateFinding] = []
    check_sections: dict[str, TicketGateStatus] = {}
    for section in REQUIRED_CHECK_SECTIONS:
        section_findings = tuple(sections[section])
        findings.extend(section_findings)
        check_sections[section] = _status_from_findings(section_findings)

    return TicketGateReport(
        contrast_id=contrast_id,
        domain=domain,
        status=_status_from_findings(tuple(findings)),
        findings=tuple(findings),
        check_sections=check_sections,
        inferred_scale=inferred_scale,
        observable_data_kinds=observable_data_kinds,
    )


def audit_tickets(tickets: Iterable[ExperimentTicket | Mapping[str, Any]]) -> tuple[TicketGateReport, ...]:
    """Audit a deterministic sequence of tickets."""

    return tuple(audit_ticket(ticket) for ticket in tickets)


def report_to_dict(report: TicketGateReport) -> dict[str, Any]:
    """Return a stable JSON-ready dictionary for a ticket-gate report."""

    return {
        "contrast_id": report.contrast_id,
        "domain": report.domain,
        "status": _status_value(report.status),
        "check_sections": {
            section: _status_value(report.check_sections[section])
            for section in REQUIRED_CHECK_SECTIONS
        },
        "inferred_scale": report.inferred_scale,
        "observable_data_kinds": list(report.observable_data_kinds),
        "findings": [
            {
                "check": finding.check,
                "severity": finding.severity,
                "code": finding.code,
                "message": finding.message,
                "fields": list(finding.fields),
            }
            for finding in report.findings
        ],
    }


def summarize_gate_reports(reports: Iterable[TicketGateReport]) -> dict[str, Any]:
    """Summarize gate reports by status and domain."""

    report_tuple = tuple(reports)
    by_status = {status.value: 0 for status in TicketGateStatus}
    by_status.update(Counter(_status_value(report.status) for report in report_tuple))

    domains = sorted({report.domain for report in report_tuple})
    by_domain: dict[str, dict[str, int]] = {}
    for domain in domains:
        domain_reports = [report for report in report_tuple if report.domain == domain]
        counts = {status.value: 0 for status in TicketGateStatus}
        counts.update(Counter(_status_value(report.status) for report in domain_reports))
        counts["total"] = len(domain_reports)
        by_domain[domain] = counts

    return {
        "total": len(report_tuple),
        "by_status": by_status,
        "by_domain": by_domain,
    }


def infer_measurement_scale(ticket: ExperimentTicket | Mapping[str, Any] | Any) -> str:
    """Infer the most likely biological or operational measurement scale."""

    values = _as_ticket_mapping(ticket)
    readout_type = _display_value(values, "readout_type").lower()
    if readout_type in _READOUT_SCALE_HINTS:
        return _READOUT_SCALE_HINTS[readout_type]

    lower_text = _ticket_text(values).lower()
    for scale, terms in _SCALE_TERMS:
        if _contains_any(lower_text, terms):
            return scale
    return "unknown"


def infer_observable_data_kinds(ticket: ExperimentTicket | Mapping[str, Any] | Any) -> tuple[str, ...]:
    """Return observable data kinds detected in the ticket text."""

    lower_text = _ticket_text(_as_ticket_mapping(ticket)).lower()
    kinds = [
        data_kind
        for data_kind, terms in _OBSERVABLE_DATA_KINDS
        if _contains_any(lower_text, terms)
    ]
    return tuple(kinds)


def _audit_required_fields(values: Mapping[str, Any]) -> list[TicketGateFinding]:
    missing = tuple(field for field in DATA_RETURN_REQUIRED_FIELDS if not _is_present(values.get(field)))
    if not missing:
        return []
    return [
        TicketGateFinding(
            check="required_fields",
            severity="error",
            code="missing_required_fields",
            message="Ticket is missing contrast-core fields required for audit.",
            fields=missing,
        )
    ]


def _audit_units(values: Mapping[str, Any]) -> list[TicketGateFinding]:
    units = values.get("readout_units")
    if not _is_present(units):
        return [
            TicketGateFinding(
                check="units",
                severity="error",
                code="missing_readout_units",
                message="Ticket readout units must be present and non-empty.",
                fields=("readout_units",),
            )
        ]

    unit_text = str(units).strip().lower()
    if "assay-specific" in unit_text or (
        _contains_any(unit_text, _VAGUE_UNIT_TERMS)
        and not _contains_any(unit_text, _CONCRETE_UNIT_TERMS)
    ):
        return [
            TicketGateFinding(
                check="units",
                severity="warning",
                code="vague_readout_units",
                message="Ticket readout units are generic and need a bounded measurement definition.",
                fields=("readout_units",),
            )
        ]
    return []


def _audit_scale(inferred_scale: str) -> list[TicketGateFinding]:
    if inferred_scale != "unknown":
        return []
    return [
        TicketGateFinding(
            check="scale",
            severity="warning",
            code="missing_measurement_scale",
            message="Ticket does not expose a biological or operational measurement scale.",
            fields=("domain", "readout_type"),
        )
    ]


def _audit_timescale(lower_text: str) -> list[TicketGateFinding]:
    if _contains_any(lower_text, _TIMING_TERMS):
        return []
    return [
        TicketGateFinding(
            check="timescale",
            severity="warning",
            code="missing_bounded_timing",
            message="Ticket needs a time window, reporting period, or bounded assay timing.",
            fields=("experiment_question", "expected_benchmark_impact"),
        )
    ]


def _audit_observability(data_kinds: tuple[str, ...]) -> list[TicketGateFinding]:
    if data_kinds:
        return []
    return [
        TicketGateFinding(
            check="observability",
            severity="error",
            code="missing_observable_data_kind",
            message="Ticket does not identify a returned data kind that can prove the readout.",
            fields=("readout_type", "expected_benchmark_impact"),
        )
    ]


def _audit_falsifiability(lower_text: str) -> list[TicketGateFinding]:
    measurement_terms = (
        "measure",
        "measurement",
        "returning",
        "returns",
        "paired",
        "compare",
        "control",
        "failure",
        "failed",
        "missing",
        "bounded",
        "adjudicat",
    )
    explicit_outcome_terms = (
        "repair",
        "reject",
        "downgrade",
        "adjudicat",
    )
    if not _contains_any(lower_text, measurement_terms):
        return [
            TicketGateFinding(
                check="falsifiability",
                severity="error",
                code="not_falsifiable",
                message="Ticket lacks a concrete measurement path that could change the contrast.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        ]

    if not _contains_any(lower_text, explicit_outcome_terms):
        return [
            TicketGateFinding(
                check="falsifiability",
                severity="warning",
                code="implicit_falsifiability_rule",
                message="Ticket should state the measurement rule for repair, rejection, or downgrade.",
                fields=("expected_benchmark_impact",),
            )
        ]
    return []


def _audit_data_return_schema(values: Mapping[str, Any]) -> list[TicketGateFinding]:
    if _missing_core_fields(values):
        return [
            TicketGateFinding(
                check="data_return_schema",
                severity="error",
                code="schema_source_missing_core_fields",
                message="Data-return schema cannot be verified until core ticket fields are present.",
                fields=tuple(_missing_core_fields(values)),
            )
        ]

    schema = values.get("data_return_schema")
    if not _is_present(schema):
        try:
            schema = make_data_return_schema(values)
        except (TypeError, ValueError) as exc:
            return [
                TicketGateFinding(
                    check="data_return_schema",
                    severity="error",
                    code="schema_construction_failed",
                    message=f"Contrast helper could not construct a data-return schema: {type(exc).__name__}.",
                    fields=("data_return_schema",),
                )
            ]

    schema_fields = _schema_required_fields(schema)
    missing_schema_fields = tuple(
        field for field in DATA_RETURN_REQUIRED_FIELDS if field not in schema_fields
    )
    if missing_schema_fields:
        return [
            TicketGateFinding(
                check="data_return_schema",
                severity="error",
                code="schema_missing_required_fields",
                message="Data-return schema omits fields required by the contrast core.",
                fields=missing_schema_fields,
            )
        ]

    schema_payload = _schema_to_json_ready(schema)
    try:
        json.dumps(schema_payload, sort_keys=True)
    except TypeError as exc:
        return [
            TicketGateFinding(
                check="data_return_schema",
                severity="error",
                code="schema_not_serializable",
                message=f"Data-return schema is not JSON-serializable: {type(exc).__name__}.",
                fields=("data_return_schema",),
            )
        ]
    return []


def _audit_domain_specific(values: Mapping[str, Any], lower_text: str) -> list[TicketGateFinding]:
    domain = _display_value(values, "domain").lower()
    if domain == "rwanda_amr":
        return _audit_rwanda_amr(lower_text)
    if domain == "regeneration_crispr":
        return _audit_regeneration(values, lower_text)
    return []


def _audit_non_claims(lower_text: str) -> list[TicketGateFinding]:
    if not _contains_any(lower_text, _FORBIDDEN_CLAIM_TERMS):
        return []
    return [
        TicketGateFinding(
            check="non_claims",
            severity="error",
            code="forbidden_broad_claim_language",
            message="Ticket contains broad claim language outside the gate scope.",
            fields=("uncertainty_note", "experiment_question", "expected_benchmark_impact"),
        )
    ]


def _audit_rwanda_amr(lower_text: str) -> list[TicketGateFinding]:
    findings: list[TicketGateFinding] = []
    readout_is_rate_or_surveillance = _contains_any(
        lower_text,
        ("resistance rate", "surveillance completeness", "antibiotic consumption"),
    )
    readout_is_ast_or_marker = _contains_any(
        lower_text,
        ("ast result", "phenotypic ast", "genomic resistance marker"),
    )

    if readout_is_rate_or_surveillance and not _contains_any(
        lower_text,
        ("denominator", "count", "counts", "percent", "ddd", "dot", "completeness"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_denominator",
                message="AMR surveillance tickets need a denominator or count field.",
                fields=("expected_benchmark_impact",),
            )
        )
    if (readout_is_rate_or_surveillance or readout_is_ast_or_marker) and not _contains_any(
        lower_text,
        ("organism", "pathogen", "bacterial", "isolate", "species"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_organism_or_pathogen",
                message="AMR tickets need an organism, pathogen, species, or isolate boundary.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        )
    if (readout_is_rate_or_surveillance or readout_is_ast_or_marker) and not _contains_any(
        lower_text,
        ("antibiotic", "antimicrobial"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_antibiotic",
                message="AMR tickets need an antibiotic or antimicrobial field where applicable.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        )
    if readout_is_rate_or_surveillance and not _contains_any(
        lower_text,
        ("specimen", "source", "site", "ward", "facility", "isolate", "accession"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_specimen_or_source",
                message="AMR surveillance tickets need specimen, source, site, or isolate provenance.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        )
    if readout_is_rate_or_surveillance and not _contains_any(
        lower_text,
        ("reporting period", "reporting-period", "facility-month", "same-period", "period"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_reporting_period",
                message="AMR surveillance tickets need a reporting period boundary.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        )
    if readout_is_ast_or_marker and "ast" not in lower_text:
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_ast_result",
                message="AMR AST or marker tickets need an AST result linkage.",
                fields=("readout_type", "expected_benchmark_impact"),
            )
        )
    if "ast result" in lower_text and not _contains_any(
        lower_text,
        ("breakpoint", "quality-control", "quality control", "qc", "repeatability"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="amr_missing_breakpoint_or_qc",
                message="AMR AST tickets need breakpoint or QC boundary fields where applicable.",
                fields=("expected_benchmark_impact",),
            )
        )
    return findings


def _audit_regeneration(values: Mapping[str, Any], lower_text: str) -> list[TicketGateFinding]:
    findings: list[TicketGateFinding] = []
    if not _contains_any(
        lower_text,
        ("dose", "timing", "time-window", "time window", "schedule", "time course", "transient"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="regen_missing_dose_or_timing",
                message="Regeneration tickets need dose, perturbation timing, or assay-duration detail.",
                fields=("experiment_question", "expected_benchmark_impact"),
            )
        )

    if _display_value(values, "readout_units").strip().lower() == "assay-specific":
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="regen_generic_readout_units",
                message="Regeneration tickets need a readout-specific scale or scoring unit.",
                fields=("readout_units",),
            )
        )

    if not _contains_any(lower_text, ("control", "failure", "failed", "malformed", "toxic")):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="regen_missing_control_or_failure_arm",
                message="Regeneration tickets need a control or failure arm.",
                fields=("control_context",),
            )
        )

    if not _contains_any(
        lower_text,
        ("replicate", "blinded", "positive/negative", "negative repair controls", "guide-level"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="regen_missing_replicate_or_blinded_scoring",
                message="Regeneration tickets need replicate support or blinded scoring language.",
                fields=("expected_benchmark_impact",),
            )
        )

    safety_scope_text = _selected_text(
        values,
        (
            "intervention_type",
            "intervention_id",
            "desired_context",
            "control_context",
            "readout_type",
            "experiment_question",
            "expected_benchmark_impact",
            "label",
        ),
    ).lower()
    safety_applicable = _contains_any(
        safety_scope_text,
        ("crispr", "organoid", "reprogramming", "small-molecule", "drug", "osk", "oskm"),
    )
    if safety_applicable and not _contains_any(
        lower_text,
        ("safety", "toxicity", "toxic", "tumor", "identity-loss", "identity loss", "non-diseased"),
    ):
        findings.append(
            TicketGateFinding(
                check="domain_specific",
                severity="warning",
                code="regen_missing_safety_boundary",
                message="Regeneration perturbation tickets need a safety or toxicity boundary.",
                fields=("control_context", "expected_benchmark_impact"),
            )
        )
    return findings


def _as_ticket_mapping(ticket: ExperimentTicket | Mapping[str, Any] | Any) -> Mapping[str, Any]:
    if isinstance(ticket, ExperimentTicket):
        return asdict(ticket)
    if isinstance(ticket, Mapping):
        return dict(ticket)
    if is_dataclass(ticket) and not isinstance(ticket, type):
        return asdict(ticket)
    to_dict = getattr(ticket, "to_dict", None)
    if callable(to_dict):
        payload = to_dict()
        if isinstance(payload, Mapping):
            return dict(payload)
    raise TypeError("ticket must be an ExperimentTicket, mapping, dataclass instance, or to_dict object")


def _display_value(values: Mapping[str, Any], field: str, fallback: str = "") -> str:
    value = values.get(field, fallback)
    if not _is_present(value):
        return fallback
    return _scalar_to_text(value)


def _ticket_text(values: Mapping[str, Any]) -> str:
    return " ".join(_scalar_to_text(values[field]) for field in _TEXT_FIELDS if field in values)


def _selected_text(values: Mapping[str, Any], fields: Iterable[str]) -> str:
    return " ".join(_scalar_to_text(values[field]) for field in fields if field in values)


def _scalar_to_text(value: Any) -> str:
    if hasattr(value, "value") and isinstance(getattr(value, "value"), str):
        return str(value.value)
    if value is None:
        return ""
    return str(value)


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _contains_any(text: str, terms: Iterable[str]) -> bool:
    return any(term in text for term in terms)


def _missing_core_fields(values: Mapping[str, Any]) -> list[str]:
    return [field for field in DATA_RETURN_REQUIRED_FIELDS if not _is_present(values.get(field))]


def _schema_required_fields(schema: Any) -> tuple[str, ...]:
    if isinstance(schema, DataReturnSchema):
        return tuple(str(field) for field in schema.required_fields)
    if isinstance(schema, str):
        try:
            parsed = json.loads(schema)
        except json.JSONDecodeError:
            return ()
        return _schema_required_fields(parsed)
    if is_dataclass(schema) and not isinstance(schema, type):
        return _schema_required_fields(asdict(schema))
    if isinstance(schema, Mapping):
        required_fields = schema.get("required_fields")
        if isinstance(required_fields, (list, tuple, set)):
            return tuple(str(field) for field in required_fields)
        return tuple(str(field) for field in schema.keys())
    return ()


def _schema_to_json_ready(schema: Any) -> Any:
    if isinstance(schema, DataReturnSchema):
        payload = asdict(schema)
        payload["evidence_status"] = evidence_status_to_string(schema.evidence_status)
        return _schema_to_json_ready(payload)
    if is_dataclass(schema) and not isinstance(schema, type):
        return _schema_to_json_ready(asdict(schema))
    if isinstance(schema, Mapping):
        return {str(key): _schema_to_json_ready(value) for key, value in schema.items()}
    if isinstance(schema, tuple):
        return [_schema_to_json_ready(value) for value in schema]
    if isinstance(schema, list):
        return [_schema_to_json_ready(value) for value in schema]
    if hasattr(schema, "value") and isinstance(getattr(schema, "value"), str):
        return schema.value
    return schema


def _status_from_findings(findings: tuple[TicketGateFinding, ...]) -> TicketGateStatus:
    if any(finding.severity == "error" for finding in findings):
        return TicketGateStatus.BLOCKED
    if any(finding.severity == "warning" for finding in findings):
        return TicketGateStatus.NEEDS_MEASUREMENT_DETAIL
    return TicketGateStatus.PASS


def _status_value(status: TicketGateStatus | str) -> str:
    if isinstance(status, TicketGateStatus):
        return status.value
    return str(status)


__all__ = [
    "REQUIRED_CHECK_SECTIONS",
    "TicketGateFinding",
    "TicketGateReport",
    "TicketGateStatus",
    "audit_ticket",
    "audit_tickets",
    "infer_measurement_scale",
    "infer_observable_data_kinds",
    "report_to_dict",
    "summarize_gate_reports",
]
