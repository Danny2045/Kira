"""Narrow regeneration/CRISPR scout built on contrast-core primitives.

This module does not make a therapeutic or discovery claim. It defines a small,
inspectable seed set that maps regenerative biology, CRISPR perturbation, and
bioelectric pattern-control examples into the project-wide contrast abstraction:

    intervention -> desired context -> control context -> readout -> evidence
    status -> missing measurement -> experiment ticket

The adapter intentionally keeps source records compact. Its job is to make the
new domain testable by the contrast core, not to act as a literature database.
"""

from __future__ import annotations

import inspect
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Any, TypeVar

try:  # Prefer the public package surface added by contrast-core.
    from kira.contrast import (
        ContrastSpec,
        DataReturnSchema,
        EvidenceRecord,
        EvidenceStatus,
        ExperimentTicket,
        make_data_return_schema,
        ticket_to_dict,
        validate_contrast_spec,
        validate_experiment_ticket,
    )
except ImportError:  # pragma: no cover - compatibility for narrower exports.
    from kira.contrast.schemas import (  # type: ignore[no-redef]
        ContrastSpec,
        DataReturnSchema,
        EvidenceRecord,
        EvidenceStatus,
        validate_contrast_spec,
    )
    from kira.contrast.tickets import (  # type: ignore[no-redef]
        ExperimentTicket,
        ticket_to_dict,
        validate_experiment_ticket,
    )

    try:
        from kira.contrast.schemas import make_data_return_schema  # type: ignore[no-redef]
    except ImportError:
        from kira.contrast.tickets import make_data_return_schema  # type: ignore[no-redef]

REGENERATION_SCOUT_VERSION = "0.1.0"

_ALLOWED_INTERVENTION_FAMILIES = frozenset(
    {
        "bioelectric perturbation",
        "crispr perturbation",
        "partial reprogramming",
        "small molecule",
        "morphogen/control signal",
    }
)

_ALLOWED_READOUT_TYPES = frozenset(
    {
        "bioelectric pattern",
        "functional rescue",
        "marker panel delta",
        "morphology score",
        "organoid phenotype",
        "wound closure",
    }
)

_STATUS_VALUES = frozenset(
    {
        "proposed",
        "single-side-only",
        "paired",
        "bounded",
        "conflicting",
    }
)

_T = TypeVar("_T")


@dataclass(frozen=True, slots=True)
class RegenerationScoutSeed:
    """Compact adapter input before conversion into contrast-core objects."""

    contrast_id: str
    intervention_family: str
    intervention: str
    desired_context: str
    control_context: str
    readout_type: str
    evidence_status: str
    missing_measurement: str
    benchmark_label: str
    priority: int
    rationale: str
    source_label: str
    source_url: str
    model_system: str
    notes: str = ""


@dataclass(frozen=True, slots=True)
class RegenerationScoutRecord:
    """Generated scout bundle for one seed."""

    seed: RegenerationScoutSeed
    contrast_spec: ContrastSpec
    evidence_record: EvidenceRecord
    experiment_ticket: ExperimentTicket


def regeneration_scout_seeds() -> tuple[RegenerationScoutSeed, ...]:
    """Return the fixed narrow seed set for the first scout branch.

    The examples are chosen to cover the requested frontier without broad
    platform claims: Michael Levin-style bioelectric control, CRISPR/organoid
    perturbation, Yamanaka-factor partial reprogramming, and morphology-first
    regeneration readouts.
    """

    seeds = (
        RegenerationScoutSeed(
            contrast_id="regen-bioelectric-xenopus-tail-vmem",
            intervention_family="bioelectric perturbation",
            intervention=(
                "membrane-voltage / proton-pump modulation after Xenopus "
                "tail amputation"
            ),
            desired_context=(
                "regenerative tail outgrowth with pattern-restored morphology"
            ),
            control_context=(
                "regeneration-refractory tail wound or malformed repair state"
            ),
            readout_type="bioelectric pattern",
            evidence_status="paired",
            missing_measurement=(
                "blinded paired morphology score plus membrane-voltage map "
                "time course under a standardized perturbation schedule"
            ),
            benchmark_label="bioelectric_pattern_restoration",
            priority=5,
            rationale=(
                "Canonical bioelectric regeneration contrast: the missing "
                "branch-worthy measurement is a reproducible desired-vs-control "
                "voltage-pattern and morphology panel."
            ),
            source_label="Michael Levin lab / Xenopus bioelectric regeneration",
            source_url="https://doi.org/10.1242/dev.02812",
            model_system="Xenopus laevis tail-amputation model",
            notes="Keep as contrast-scout seed, not as a therapeutic claim.",
        ),
        RegenerationScoutSeed(
            contrast_id="regen-bioelectric-planarian-target-morphology",
            intervention_family="bioelectric perturbation",
            intervention=(
                "gap-junction or endogenous bioelectric-gradient perturbation "
                "in planarian regeneration"
            ),
            desired_context="species-typical restored head/body morphology",
            control_context="ectopic, species-shifted, or malformed regenerated anatomy",
            readout_type="morphology score",
            evidence_status="bounded",
            missing_measurement=(
                "dose-bounded morphology atlas with replicate-level head-shape "
                "labels and recovery/failure annotations"
            ),
            benchmark_label="target_morphology_recovery",
            priority=4,
            rationale=(
                "Pattern-control examples are useful only when bounded by "
                "explicit dose, timing, and morphology-return measurements."
            ),
            source_label="Michael Levin lab / planarian target morphology",
            source_url="https://doi.org/10.3390/ijms161126061",
            model_system="Girardia / planarian regeneration model",
            notes="Useful scout record for morphology-as-readout discipline.",
        ),
        RegenerationScoutSeed(
            contrast_id="regen-crispr-assembloid-interneuron-migration",
            intervention_family="crispr perturbation",
            intervention=(
                "pooled CRISPR loss-of-function perturbation in human "
                "forebrain assembloids"
            ),
            desired_context="corrected or preserved interneuron generation and migration",
            control_context="failed interneuron migration or malformed circuit assembly",
            readout_type="organoid phenotype",
            evidence_status="paired",
            missing_measurement=(
                "paired perturbation-vs-control migration index with marker-panel "
                "delta and guide-level replicate support"
            ),
            benchmark_label="organoid_crispr_migration_rescue",
            priority=5,
            rationale=(
                "CRISPR/assembloid systems give Kira a clean perturbation-to-"
                "phenotype contrast shape without claiming clinical repair."
            ),
            source_label="Human forebrain assembloid CRISPR screen",
            source_url="https://doi.org/10.1038/s41586-023-06564-w",
            model_system="hiPSC-derived forebrain organoids / assembloids",
            notes="Treat as a development/phenotype scout, not regeneration proof.",
        ),
        RegenerationScoutSeed(
            contrast_id="regen-crispr-3d-organoid-recovery",
            intervention_family="crispr perturbation",
            intervention=(
                "CRISPR knockout, CRISPRi, CRISPRa, or single-cell CRISPR "
                "perturbation in primary human 3D organoids"
            ),
            desired_context="cell-state recovery or preserved organoid architecture",
            control_context="drug-injured, toxic, or failed-recovery organoid state",
            readout_type="marker panel delta",
            evidence_status="single-side-only",
            missing_measurement=(
                "paired non-diseased organoid safety arm with single-cell "
                "state recovery and toxicity readouts"
            ),
            benchmark_label="organoid_recovery_safety_pair",
            priority=4,
            rationale=(
                "Human 3D organoid perturbation screens are powerful, but this "
                "scout should demand an explicit desired-state plus safety/failure arm."
            ),
            source_label="Primary human 3D organoid CRISPR perturbation screens",
            source_url="https://doi.org/10.1038/s41467-025-62818-3",
            model_system="primary human 3D organoids",
            notes="Scout asks for missing paired safety side before any stronger claim.",
        ),
        RegenerationScoutSeed(
            contrast_id="regen-yamanaka-partial-reprogramming-window",
            intervention_family="partial reprogramming",
            intervention=(
                "transient OSK/OSKM or non-genetic partial-reprogramming signal"
            ),
            desired_context="repaired or rejuvenated cell state with lineage retained",
            control_context="loss of identity, pluripotency drift, tumor risk, or toxicity",
            readout_type="functional rescue",
            evidence_status="conflicting",
            missing_measurement=(
                "dose/time-window experiment returning lineage-retention, "
                "tumor-risk, epigenetic-age, and functional-rescue measurements"
            ),
            benchmark_label="partial_reprogramming_safety_window",
            priority=5,
            rationale=(
                "Yamanaka-factor work belongs in Kira as a bounded contrast: "
                "beneficial repair signals must be paired against identity-loss "
                "and safety failure states."
            ),
            source_label="Yamanaka-factor partial reprogramming review literature",
            source_url="https://doi.org/10.1038/s41467-024-46020-5",
            model_system="cell, mouse, and tissue-specific partial-reprogramming models",
            notes="Explicitly retain as a safety-window ticket, not rejuvenation hype.",
        ),
        RegenerationScoutSeed(
            contrast_id="regen-organoid-morphology-phenoscreen",
            intervention_family="morphogen/control signal",
            intervention=(
                "morphogen, pathway, or small-molecule control signal in a "
                "regeneration-oriented organoid assay"
            ),
            desired_context="regenerative organoid morphology or wound-repair phenotype",
            control_context="failed organoid formation, scar-like state, or malformed repair",
            readout_type="morphology score",
            evidence_status="proposed",
            missing_measurement=(
                "image-based morphology benchmark with positive/negative repair "
                "controls and marker-panel confirmation"
            ),
            benchmark_label="morphology_first_regeneration_screen",
            priority=3,
            rationale=(
                "This seed steers frontier systems into Kira's measurable-ticket "
                "discipline before expanding to broader regenerative assays."
            ),
            source_label="Image-based organoid regeneration phenotyping literature",
            source_url="https://doi.org/10.1038/s41586-020-2776-9",
            model_system="intestinal or tissue-specific organoid regeneration assay",
            notes="Starter seed for future lab-campaign repair tickets.",
        ),
    )
    for seed in seeds:
        validate_regeneration_seed(seed)
    return seeds


def validate_regeneration_seed(seed: RegenerationScoutSeed) -> None:
    """Validate local scout constraints before using contrast-core validators."""

    required_strings = {
        "contrast_id": seed.contrast_id,
        "intervention_family": seed.intervention_family,
        "intervention": seed.intervention,
        "desired_context": seed.desired_context,
        "control_context": seed.control_context,
        "readout_type": seed.readout_type,
        "evidence_status": seed.evidence_status,
        "missing_measurement": seed.missing_measurement,
        "benchmark_label": seed.benchmark_label,
        "rationale": seed.rationale,
        "source_label": seed.source_label,
        "source_url": seed.source_url,
        "model_system": seed.model_system,
    }
    blanks = [name for name, value in required_strings.items() if not value.strip()]
    if blanks:
        raise ValueError(f"blank regeneration scout fields: {', '.join(blanks)}")

    if seed.intervention_family not in _ALLOWED_INTERVENTION_FAMILIES:
        raise ValueError(
            f"unsupported intervention family {seed.intervention_family!r}; "
            f"expected one of {sorted(_ALLOWED_INTERVENTION_FAMILIES)!r}"
        )
    if seed.readout_type not in _ALLOWED_READOUT_TYPES:
        raise ValueError(
            f"unsupported readout type {seed.readout_type!r}; "
            f"expected one of {sorted(_ALLOWED_READOUT_TYPES)!r}"
        )
    if _normalize_status(seed.evidence_status) not in _STATUS_VALUES:
        raise ValueError(
            f"unsupported evidence status {seed.evidence_status!r}; "
            f"expected one of {sorted(_STATUS_VALUES)!r}"
        )
    if not 1 <= seed.priority <= 5:
        raise ValueError("priority must be on the 1..5 scout scale")
    if not seed.source_url.startswith(("https://doi.org/", "https://")):
        raise ValueError("source_url must be an https URL or DOI resolver URL")


def build_contrast_specs(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> tuple[ContrastSpec, ...]:
    """Convert regeneration scout seeds into validated ContrastSpec objects."""

    specs: list[ContrastSpec] = []
    for seed in _ordered_seeds(seeds):
        values = _contrast_values(seed)
        spec = _construct(ContrastSpec, values)
        validate_contrast_spec(spec)
        specs.append(spec)
    return tuple(specs)


def build_evidence_records(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> tuple[EvidenceRecord, ...]:
    """Convert scout seeds into EvidenceRecord objects."""

    records: list[EvidenceRecord] = []
    for seed in _ordered_seeds(seeds):
        values = _evidence_values(seed)
        records.append(_construct(EvidenceRecord, values))
    return tuple(records)


def build_experiment_tickets(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> tuple[ExperimentTicket, ...]:
    """Convert scout seeds into validated ExperimentTicket objects."""

    tickets: list[ExperimentTicket] = []
    for seed in _ordered_seeds(seeds):
        values = _ticket_values(seed)
        ticket = _construct(ExperimentTicket, values)
        validate_experiment_ticket(ticket)
        tickets.append(ticket)
    return tuple(tickets)


def build_regeneration_scout(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> tuple[RegenerationScoutRecord, ...]:
    """Build the full seed/spec/evidence/ticket scout bundle."""

    ordered = _ordered_seeds(seeds)
    specs = build_contrast_specs(ordered)
    evidence_records = build_evidence_records(ordered)
    tickets = build_experiment_tickets(ordered)
    return tuple(
        RegenerationScoutRecord(
            seed=seed,
            contrast_spec=spec,
            evidence_record=evidence_record,
            experiment_ticket=ticket,
        )
        for seed, spec, evidence_record, ticket in zip(
            ordered, specs, evidence_records, tickets, strict=True
        )
    )


def tickets_as_dicts(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> tuple[dict[str, Any], ...]:
    """Return experiment tickets in the contrast-core dict representation."""

    return tuple(ticket_to_dict(ticket) for ticket in build_experiment_tickets(seeds))


def summarize_evidence_statuses(
    seeds: Iterable[RegenerationScoutSeed] | None = None,
) -> dict[str, int]:
    """Count scout seeds by evidence status."""

    return dict(Counter(seed.evidence_status for seed in _ordered_seeds(seeds)))


def _ordered_seeds(
    seeds: Iterable[RegenerationScoutSeed] | None,
) -> tuple[RegenerationScoutSeed, ...]:
    raw_seeds = regeneration_scout_seeds() if seeds is None else tuple(seeds)
    for seed in raw_seeds:
        validate_regeneration_seed(seed)
    return tuple(sorted(raw_seeds, key=lambda seed: (-seed.priority, seed.contrast_id)))


def _contrast_values(seed: RegenerationScoutSeed) -> dict[str, Any]:
    status = _normalize_status(seed.evidence_status)
    known_side, missing_side = _known_and_missing_sides(status)
    uncertainty_note = _uncertainty_note(seed)
    intervention_type = _normalize_field(seed.intervention_family).replace("/", "_")

    return {
        "id": seed.contrast_id,
        "contrast_id": seed.contrast_id,
        "name": seed.contrast_id,
        "domain": "regeneration_crispr",
        "intervention_type": intervention_type,
        "intervention_id": f"{intervention_type}::{seed.contrast_id}",
        "intervention": seed.intervention,
        "intervention_family": seed.intervention_family,
        "desired_context": seed.desired_context,
        "control_context": seed.control_context,
        "readout": seed.readout_type,
        "readout_type": seed.readout_type,
        "readout_units": "assay-specific",
        "evidence_status": _coerce_evidence_status(seed.evidence_status),
        "known_side": known_side,
        "missing_side": missing_side,
        "label": seed.benchmark_label,
        "uncertainty_note": uncertainty_note,
        "provenance": f"{seed.source_label}; {seed.source_url}",
        "missing_measurement": seed.missing_measurement,
        "benchmark_label": seed.benchmark_label,
        "model_system": seed.model_system,
        "source_label": seed.source_label,
        "source_url": seed.source_url,
        "notes": seed.notes,
    }


def _evidence_values(seed: RegenerationScoutSeed) -> dict[str, Any]:
    return {
        **_contrast_values(seed),
        "id": f"evidence-{seed.contrast_id}",
        "record_id": f"evidence-{seed.contrast_id}",
        "evidence_id": f"evidence-{seed.contrast_id}",
        "status": _coerce_evidence_status(seed.evidence_status),
        "evidence_status": _coerce_evidence_status(seed.evidence_status),
        "evidence_source": seed.source_label,
        "source": seed.source_label,
        "source_url": seed.source_url,
        "measurement_status": seed.evidence_status,
        "missing_measurement": seed.missing_measurement,
        "has_desired_side": _normalize_status(seed.evidence_status)
        in {"paired", "bounded", "conflicting"},
        "has_control_side": _normalize_status(seed.evidence_status)
        in {"paired", "bounded"},
    }


def _ticket_values(seed: RegenerationScoutSeed) -> dict[str, Any]:
    schema = _make_regeneration_data_return_schema(seed)
    return {
        **_contrast_values(seed),
        "id": f"ticket-{seed.contrast_id}",
        "ticket_id": f"ticket-{seed.contrast_id}",
        "experiment_id": f"ticket-{seed.contrast_id}",
        "title": f"Regeneration contrast scout: {seed.benchmark_label}",
        "missing_experiment": seed.missing_measurement,
        "missing_measurement": seed.missing_measurement,
        "measurement_gap": seed.missing_measurement,
        "priority": seed.priority,
        "priority_score": seed.priority,
        "rank": seed.priority,
        "rationale": seed.rationale,
        "why": seed.rationale,
        "experiment_question": _experiment_question(seed),
        "expected_benchmark_impact": _expected_benchmark_impact(seed),
        "data_return_schema": schema,
        "return_schema": schema,
        "schema": schema,
        "source_label": seed.source_label,
        "source_url": seed.source_url,
    }


def _make_regeneration_data_return_schema(seed: RegenerationScoutSeed) -> DataReturnSchema:
    values = {
        **_contrast_values(seed),
        "experiment_question": _experiment_question(seed),
        "expected_benchmark_impact": _expected_benchmark_impact(seed),
    }
    return make_data_return_schema(values)


_STATUS_ALIASES = {
    "paired": "paired-complete",
    "single-side-only": "single-side-only",
}


def _known_and_missing_sides(status: str) -> tuple[str, str]:
    if status == "paired":
        return "desired_and_control", "none"
    if status == "bounded":
        return "desired_and_control_bounded", "dose_timing_or_safety_boundary"
    if status == "single-side-only":
        return "desired", "control_or_safety"
    if status == "conflicting":
        return "conflicting_evidence", "matched_resolution_measurement"
    if status == "proposed":
        return "none", "desired_and_control"
    return "unknown", "unknown"


def _experiment_question(seed: RegenerationScoutSeed) -> str:
    return (
        f"Measure {seed.intervention} in {seed.model_system} to compare "
        f"{seed.desired_context} against {seed.control_context} using "
        f"{seed.readout_type}."
    )


def _expected_benchmark_impact(seed: RegenerationScoutSeed) -> str:
    return (
        f"Returns the missing measurement for {seed.benchmark_label}: "
        f"{seed.missing_measurement}"
    )


def _uncertainty_note(seed: RegenerationScoutSeed) -> str:
    suffix = seed.notes or "No therapeutic, efficacy, or wet-lab validation claim."
    return f"Scout seed only; {seed.rationale} {suffix}"


def _coerce_evidence_status(value: str) -> EvidenceStatus:
    normalized = _STATUS_ALIASES.get(_normalize_status(value), _normalize_status(value))
    members = _evidence_status_members()
    try:
        return members[normalized]
    except KeyError as exc:
        raise ValueError(f"unknown evidence status: {value!r}") from exc


def _evidence_status_members() -> dict[str, EvidenceStatus]:
    members: dict[str, EvidenceStatus] = {}
    if hasattr(EvidenceStatus, "__members__"):
        for member in EvidenceStatus:  # type: ignore[union-attr]
            _add_status_member(members, member)
    for attr in dir(EvidenceStatus):
        if attr.startswith("_"):
            continue
        maybe_member = getattr(EvidenceStatus, attr)
        if isinstance(maybe_member, (str, Enum)):
            _add_status_member(members, maybe_member)
    return members


def _add_status_member(members: dict[str, EvidenceStatus], member: Any) -> None:
    raw_values = [member]
    if isinstance(member, Enum):
        raw_values.extend([member.name, member.value])
    for raw_value in raw_values:
        members[_normalize_status(str(raw_value))] = member


def _normalize_status(value: str) -> str:
    return value.strip().lower().replace("_", "-")


def _construct(factory: Callable[..., _T], values: Mapping[str, Any]) -> _T:
    return _call_with_supported_kwargs(factory, values)


def _call_with_supported_kwargs(
    factory: Callable[..., _T], values: Mapping[str, Any]
) -> _T:
    signature = inspect.signature(factory)
    parameters = tuple(signature.parameters.values())
    if any(parameter.kind == inspect.Parameter.VAR_KEYWORD for parameter in parameters):
        return factory(**dict(values))

    normalized_values = {_normalize_field(name): value for name, value in values.items()}
    kwargs: dict[str, Any] = {}
    missing: list[str] = []
    for parameter in parameters:
        if parameter.kind in {
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        }:
            continue
        value = _value_for_parameter(parameter.name, normalized_values)
        if value is _MISSING:
            if parameter.default is inspect.Parameter.empty:
                missing.append(parameter.name)
            continue
        kwargs[parameter.name] = value

    if missing:
        raise TypeError(
            f"cannot construct {factory!r}; unsupported required fields: "
            f"{', '.join(missing)}"
        )
    return factory(**kwargs)


_MISSING = object()

_FIELD_ALIASES = {
    "assay_id": ("ticket_id", "experiment_id", "id"),
    "benchmark": ("benchmark_label",),
    "benchmark_label": ("label", "benchmark"),
    "contrast_key": ("contrast_id", "id"),
    "data_schema": ("data_return_schema", "return_schema", "schema"),
    "desired_biology": ("desired_context",),
    "desired_state": ("desired_context",),
    "evidence_id": ("record_id", "id"),
    "failure_context": ("control_context",),
    "label": ("benchmark_label",),
    "measurement_gap": ("missing_measurement", "missing_experiment"),
    "missing_experiment": ("missing_measurement", "measurement_gap"),
    "name": ("contrast_id", "id"),
    "negative_context": ("control_context",),
    "positive_context": ("desired_context",),
    "priority_score": ("priority", "rank"),
    "readout": ("readout_type",),
    "return_schema": ("data_return_schema", "schema"),
    "safety_context": ("control_context",),
    "spec_id": ("contrast_id", "id"),
    "status": ("evidence_status",),
    "target_context": ("desired_context",),
    "ticket_key": ("ticket_id", "id"),
    "uid": ("contrast_id", "ticket_id", "record_id", "id"),
}


def _value_for_parameter(
    parameter_name: str, normalized_values: Mapping[str, Any]
) -> Any:
    normalized_name = _normalize_field(parameter_name)
    candidate_names = [normalized_name]
    candidate_names.extend(_FIELD_ALIASES.get(normalized_name, ()))
    if "contrast" in normalized_name and "id" in normalized_name:
        candidate_names.extend(["contrast_id", "id"])
    if "ticket" in normalized_name and "id" in normalized_name:
        candidate_names.extend(["ticket_id", "experiment_id", "id"])
    if "record" in normalized_name and "id" in normalized_name:
        candidate_names.extend(["record_id", "evidence_id", "id"])
    if "source" in normalized_name and "url" not in normalized_name:
        candidate_names.extend(["source_label", "source"])
    if "url" in normalized_name:
        candidate_names.append("source_url")

    for candidate_name in candidate_names:
        normalized_candidate = _normalize_field(candidate_name)
        if normalized_candidate in normalized_values:
            return normalized_values[normalized_candidate]
    return _MISSING


def _normalize_field(value: str) -> str:
    return value.strip().lower().replace("-", "_").replace(" ", "_")
