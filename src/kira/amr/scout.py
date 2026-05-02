"""Rwanda AMR contrast scout built on contrast-core primitives.

This module is a narrow evidence adapter. It maps a fixed Rwanda-relevant AMR
seed set into Kira's shared contrast shape:

    intervention -> desired context -> control context -> readout -> evidence
    status -> missing measurement -> experiment ticket

The scout is for surveillance, stewardship evidence, infection-control
benchmarks, and benchmark-repair tickets. It does not provide individual-patient
treatment advice, assert lab completion, or claim that AMR is solved.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Any

from kira.contrast import (
    ContrastSpec,
    EvidenceRecord,
    EvidenceStatus,
    ExperimentTicket,
    make_data_return_schema,
    ticket_to_dict,
    validate_contrast_spec,
    validate_experiment_ticket,
)

AMR_SCOUT_VERSION = "0.1.0"

_ALLOWED_INTERVENTION_FAMILIES = frozenset(
    {
        "stewardship action",
        "diagnostic/AST measurement",
        "infection-control measure",
        "surveillance sampling",
        "antibiotic exposure policy",
        "lab validation assay",
    }
)

_ALLOWED_DESIRED_CONTEXTS = frozenset(
    {
        "susceptible infection",
        "reduced resistance signal",
        "appropriate antibiotic use",
        "controlled transmission",
        "benchmark-complete AMR evidence",
    }
)

_ALLOWED_CONTROL_CONTEXTS = frozenset(
    {
        "resistant infection",
        "inappropriate antibiotic use",
        "treatment failure",
        "outbreak cluster",
        "missing AST",
        "unpaired surveillance evidence",
    }
)

_ALLOWED_READOUT_TYPES = frozenset(
    {
        "AST result",
        "resistance rate",
        "antibiotic consumption",
        "treatment outcome",
        "infection-control metric",
        "genomic resistance marker",
        "surveillance completeness",
    }
)

_STATUS_VALUES = frozenset(
    {
        "proposed",
        "single_side_only",
        "paired_complete",
        "bounded",
        "conflicting",
    }
)

_FORBIDDEN_CLAIMS = frozenset(
    {
        "clinical recommendation",
        "clinical recommendations",
        "prescribe",
        "solved amr",
        "solve amr",
        "wet-lab validated",
        "wet lab validated",
        "cure",
    }
)

_READOUT_UNITS = {
    "AST result": "susceptible/intermediate/resistant category or MIC",
    "resistance rate": "percent resistant",
    "antibiotic consumption": "DDD, DOT, or facility denominator",
    "treatment outcome": "aggregate outcome category",
    "infection-control metric": "incidence rate, cluster count, or days since last case",
    "genomic resistance marker": "marker detected/not_detected with linked isolate id",
    "surveillance completeness": "percent complete",
}


@dataclass(frozen=True, slots=True)
class AmrScoutSeed:
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
    rwanda_context: str
    notes: str = ""


@dataclass(frozen=True, slots=True)
class AmrScoutRecord:
    """Generated scout bundle for one AMR seed."""

    seed: AmrScoutSeed
    contrast_spec: ContrastSpec
    evidence_record: EvidenceRecord
    experiment_ticket: ExperimentTicket


def amr_scout_seeds() -> tuple[AmrScoutSeed, ...]:
    """Return the fixed narrow Rwanda AMR scout seed set."""

    seeds = (
        AmrScoutSeed(
            contrast_id="rwanda-amr-ast-completeness-priority-isolates",
            intervention_family="diagnostic/AST measurement",
            intervention=(
                "attach organism identification and AST panel completeness flags "
                "to priority bacterial isolate records"
            ),
            desired_context="benchmark-complete AMR evidence",
            control_context="missing AST",
            readout_type="surveillance completeness",
            evidence_status="single_side_only",
            missing_measurement=(
                "facility-month denominator of priority bacterial isolates with "
                "species identification, specimen source, AST panel, and AST result"
            ),
            benchmark_label="priority_isolate_ast_completeness",
            priority=5,
            rationale=(
                "Completeness is the first repair target before resistance-rate "
                "comparisons are treated as benchmark evidence."
            ),
            source_label="WHO GLASS-AMR routine data surveillance",
            source_url="https://www.who.int/initiatives/glass/glass-routine-data-surveillance",
            rwanda_context=(
                "Rwanda sentinel and reference laboratory stream for common "
                "human bacterial pathogens"
            ),
            notes="Evidence and benchmark-repair ticket only.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-sentinel-sampling-gap",
            intervention_family="surveillance sampling",
            intervention=(
                "map sentinel site, specimen type, and reporting-period coverage "
                "before comparing AMR rates"
            ),
            desired_context="benchmark-complete AMR evidence",
            control_context="unpaired surveillance evidence",
            readout_type="surveillance completeness",
            evidence_status="proposed",
            missing_measurement=(
                "site-by-specimen sampling frame with isolate counts, zero-reporting "
                "flags, and reporting-period coverage"
            ),
            benchmark_label="sentinel_sampling_completeness",
            priority=4,
            rationale=(
                "A Rwanda-first scout should separate true resistance signals from "
                "unpaired or missing surveillance denominators."
            ),
            source_label="Rwanda NAP 2.0 One Health AMR surveillance objective",
            source_url=(
                "https://www.rbc.gov.rw/fileadmin/user_upload/strategy/"
                "2nd_NATIONAL_ACTION_PLAN_ON_AMR__NAP_2025-2029_.pdf"
            ),
            rwanda_context=(
                "human-health, animal-health, agriculture, and environment "
                "coordination channels referenced by Rwanda NAP 2.0"
            ),
            notes="Sampling-frame gap ticket; no incidence claim.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-stewardship-ast-linked-use-audit",
            intervention_family="stewardship action",
            intervention=(
                "facility antibiotic-use audit linked to diagnostic result "
                "availability and documented indication category"
            ),
            desired_context="appropriate antibiotic use",
            control_context="inappropriate antibiotic use",
            readout_type="treatment outcome",
            evidence_status="single_side_only",
            missing_measurement=(
                "paired facility aggregate linking antibiotic-use category, AST "
                "availability, indication group, and outcome completeness"
            ),
            benchmark_label="ast_linked_stewardship_evidence_gap",
            priority=5,
            rationale=(
                "Stewardship belongs in this scout as an evidence-gap contrast, "
                "not as individual treatment instruction."
            ),
            source_label="WHO GLASS antimicrobial use surveillance module",
            source_url="https://www.who.int/initiatives/glass/glass-amc-module",
            rwanda_context=(
                "Rwanda hospital stewardship and antimicrobial-use monitoring "
                "data streams"
            ),
            notes="Aggregate-use evidence ticket only.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-ipc-outbreak-cluster-contrast",
            intervention_family="infection-control measure",
            intervention=(
                "infection-control line list, isolate linkage, and ward-time "
                "exposure review for a suspected resistant-organism cluster"
            ),
            desired_context="controlled transmission",
            control_context="outbreak cluster",
            readout_type="infection-control metric",
            evidence_status="bounded",
            missing_measurement=(
                "time-bounded cluster curve with isolate linkage, ward exposure "
                "window, and days since last linked case"
            ),
            benchmark_label="ipc_cluster_control_contrast",
            priority=4,
            rationale=(
                "Outbreak contrasts are useful only with explicit time, ward, "
                "case-definition, and isolate-linkage boundaries."
            ),
            source_label="Rwanda NAP 2.0 infection prevention and control objective",
            source_url=(
                "https://www.rbc.gov.rw/fileadmin/user_upload/strategy/"
                "2nd_NATIONAL_ACTION_PLAN_ON_AMR__NAP_2025-2029_.pdf"
            ),
            rwanda_context="Rwanda facility IPC and laboratory reporting interface",
            notes="Bounded transmission-control benchmark ticket.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-genomic-marker-phenotype-confirmation",
            intervention_family="lab validation assay",
            intervention=(
                "confirm candidate resistance marker calls against phenotypic AST "
                "for the same isolate accession"
            ),
            desired_context="benchmark-complete AMR evidence",
            control_context="resistant infection",
            readout_type="genomic resistance marker",
            evidence_status="conflicting",
            missing_measurement=(
                "paired marker-call, sequence-quality, organism identification, "
                "antibiotic, and phenotypic AST record for each isolate"
            ),
            benchmark_label="genotype_phenotype_marker_confirmation",
            priority=5,
            rationale=(
                "Genomic markers should repair AMR benchmarks only when tied to "
                "the matching isolate phenotype and quality fields."
            ),
            source_label="WHO GLASS AMR surveillance methods",
            source_url="https://www.who.int/publications/i/item/9789240076600",
            rwanda_context=(
                "Rwanda reference-laboratory isolate records that may later be "
                "linked to sequence-derived marker calls"
            ),
            notes="Marker-confirmation ticket; no standalone resistance claim.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-consumption-resistance-paired-measure",
            intervention_family="antibiotic exposure policy",
            intervention=(
                "pair hospital antibiotic consumption summary with facility-level "
                "pathogen-antibiotic resistance rate"
            ),
            desired_context="reduced resistance signal",
            control_context="resistant infection",
            readout_type="antibiotic consumption",
            evidence_status="paired_complete",
            missing_measurement=(
                "same-period DDD or DOT denominator paired with pathogen-antibiotic "
                "resistance rate and AST completeness denominator"
            ),
            benchmark_label="consumption_resistance_pair",
            priority=4,
            rationale=(
                "Consumption and resistance become useful contrast infrastructure "
                "only when paired by facility, period, pathogen, and antibiotic."
            ),
            source_label="WHO GLASS methodology for antimicrobial consumption",
            source_url="https://www.who.int/publications/i/item/9789240012639",
            rwanda_context=(
                "Rwanda facility-level antimicrobial consumption and AMR "
                "surveillance records"
            ),
            notes="Paired measurement scaffold, not policy advice.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-ast-qc-lab-assay-ticket",
            intervention_family="lab validation assay",
            intervention=(
                "record AST quality-control strain, method, breakpoint version, "
                "and repeatability flag for priority isolate testing"
            ),
            desired_context="benchmark-complete AMR evidence",
            control_context="missing AST",
            readout_type="AST result",
            evidence_status="proposed",
            missing_measurement=(
                "AST method, breakpoint table version, quality-control outcome, "
                "repeatability flag, and isolate-antibiotic result"
            ),
            benchmark_label="ast_qc_assay_validation_ticket",
            priority=3,
            rationale=(
                "The scout can issue an executable lab-data ticket without "
                "asserting that validation has already happened."
            ),
            source_label="CLSI antimicrobial susceptibility testing standards",
            source_url="https://clsi.org/shop/standards/m100/",
            rwanda_context=(
                "Rwanda bacteriology laboratory AST quality-control and reporting "
                "workflow"
            ),
            notes="Executable data-ticket shape only.",
        ),
        AmrScoutSeed(
            contrast_id="rwanda-amr-pathogen-antibiotic-benchmark-pair",
            intervention_family="diagnostic/AST measurement",
            intervention=(
                "compile one benchmark-ready pathogen-antibiotic pair with AST "
                "denominator, susceptible count, and resistant count"
            ),
            desired_context="susceptible infection",
            control_context="resistant infection",
            readout_type="resistance rate",
            evidence_status="paired_complete",
            missing_measurement=(
                "pathogen, antibiotic, specimen type, reporting period, AST "
                "denominator, susceptible count, and resistant count"
            ),
            benchmark_label="pathogen_antibiotic_pair_record",
            priority=5,
            rationale=(
                "A complete pathogen-antibiotic pair is the smallest world-scale "
                "AMR contrast row that still keeps Rwanda-specific provenance."
            ),
            source_label="WHO GLASS priority bacterial pathogen surveillance",
            source_url="https://www.who.int/initiatives/glass/glass-routine-data-surveillance",
            rwanda_context=(
                "Rwanda GLASS-aligned bacterial pathogen and antibiotic reporting "
                "record"
            ),
            notes="Benchmark row, not a treatment statement.",
        ),
    )
    for seed in seeds:
        validate_amr_seed(seed)
    return seeds


def validate_amr_seed(seed: AmrScoutSeed) -> None:
    """Validate local AMR scout constraints before contrast-core conversion."""

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
        "rwanda_context": seed.rwanda_context,
    }
    blanks = [name for name, value in required_strings.items() if not value.strip()]
    if blanks:
        raise ValueError(f"blank AMR scout fields: {', '.join(blanks)}")

    if seed.intervention_family not in _ALLOWED_INTERVENTION_FAMILIES:
        raise ValueError(
            f"unsupported intervention family {seed.intervention_family!r}; "
            f"expected one of {sorted(_ALLOWED_INTERVENTION_FAMILIES)!r}"
        )
    if seed.desired_context not in _ALLOWED_DESIRED_CONTEXTS:
        raise ValueError(
            f"unsupported desired context {seed.desired_context!r}; "
            f"expected one of {sorted(_ALLOWED_DESIRED_CONTEXTS)!r}"
        )
    if seed.control_context not in _ALLOWED_CONTROL_CONTEXTS:
        raise ValueError(
            f"unsupported control context {seed.control_context!r}; "
            f"expected one of {sorted(_ALLOWED_CONTROL_CONTEXTS)!r}"
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
    if not seed.source_url.startswith("https://"):
        raise ValueError("source_url must be an https URL")

    seed_text = repr(seed).lower()
    claims = sorted(claim for claim in _FORBIDDEN_CLAIMS if claim in seed_text)
    if claims:
        raise ValueError(f"AMR scout seed contains forbidden claim text: {claims}")


def build_contrast_specs(
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> tuple[ContrastSpec, ...]:
    """Convert AMR scout seeds into validated ContrastSpec objects."""

    specs: list[ContrastSpec] = []
    for seed in _ordered_seeds(seeds):
        spec = validate_contrast_spec(_contrast_values(seed))
        specs.append(spec)
    return tuple(specs)


def build_evidence_records(
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> tuple[EvidenceRecord, ...]:
    """Convert AMR scout seeds into EvidenceRecord objects."""

    records: list[EvidenceRecord] = []
    for seed in _ordered_seeds(seeds):
        values = validate_contrast_spec(_contrast_values(seed))
        records.append(EvidenceRecord(**values.__dict__))
    return tuple(records)


def build_experiment_tickets(
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> tuple[ExperimentTicket, ...]:
    """Convert AMR scout seeds into validated ExperimentTicket objects."""

    tickets: list[ExperimentTicket] = []
    for seed in _ordered_seeds(seeds):
        ticket = validate_experiment_ticket(_ticket_values(seed))
        tickets.append(ticket)
    return tuple(tickets)


def build_amr_scout(
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> tuple[AmrScoutRecord, ...]:
    """Build the full seed/spec/evidence/ticket AMR scout bundle."""

    ordered = _ordered_seeds(seeds)
    specs = build_contrast_specs(ordered)
    evidence_records = build_evidence_records(ordered)
    tickets = build_experiment_tickets(ordered)
    return tuple(
        AmrScoutRecord(
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
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> tuple[dict[str, Any], ...]:
    """Return AMR experiment tickets in the contrast-core dict representation."""

    return tuple(ticket_to_dict(ticket) for ticket in build_experiment_tickets(seeds))


def summarize_evidence_statuses(
    seeds: Iterable[AmrScoutSeed] | None = None,
) -> dict[str, int]:
    """Count scout seeds by evidence status."""

    return dict(Counter(seed.evidence_status for seed in _ordered_seeds(seeds)))


def _ordered_seeds(seeds: Iterable[AmrScoutSeed] | None) -> tuple[AmrScoutSeed, ...]:
    raw_seeds = amr_scout_seeds() if seeds is None else tuple(seeds)
    for seed in raw_seeds:
        validate_amr_seed(seed)
    return tuple(sorted(raw_seeds, key=lambda seed: (-seed.priority, seed.contrast_id)))


def _contrast_values(seed: AmrScoutSeed) -> dict[str, Any]:
    status = _normalize_status(seed.evidence_status)
    known_side, missing_side = _known_and_missing_sides(status)
    intervention_type = _normalize_field(seed.intervention_family)
    return {
        "contrast_id": seed.contrast_id,
        "domain": "rwanda_amr",
        "intervention_type": intervention_type,
        "intervention_id": f"{intervention_type}::{seed.contrast_id}",
        "desired_context": seed.desired_context,
        "control_context": seed.control_context,
        "readout_type": seed.readout_type,
        "readout_units": _READOUT_UNITS[seed.readout_type],
        "evidence_status": _coerce_evidence_status(seed.evidence_status),
        "known_side": known_side,
        "missing_side": missing_side,
        "label": seed.benchmark_label,
        "uncertainty_note": _uncertainty_note(seed),
        "provenance": f"{seed.source_label}; {seed.source_url}",
    }


def _ticket_values(seed: AmrScoutSeed) -> dict[str, Any]:
    values = {
        **_contrast_values(seed),
        "experiment_question": _experiment_question(seed),
        "expected_benchmark_impact": _expected_benchmark_impact(seed),
    }
    values["data_return_schema"] = make_data_return_schema(values)
    return values


def _known_and_missing_sides(status: str) -> tuple[str, str]:
    if status == "paired_complete":
        return "desired_and_control", "none"
    if status == "bounded":
        return "desired_and_control_bounded", "time_place_denominator_or_qc_boundary"
    if status == "single_side_only":
        return "one_side", "paired_comparator"
    if status == "conflicting":
        return "conflicting_evidence", "adjudicating_paired_measurement"
    if status == "proposed":
        return "none", "desired_and_control"
    return "unknown", "unknown"


def _experiment_question(seed: AmrScoutSeed) -> str:
    return (
        f"Execute the Rwanda AMR scout ticket for {seed.benchmark_label} in "
        f"{seed.rwanda_context}: compare {seed.desired_context} with "
        f"{seed.control_context} under {seed.intervention}, returning "
        f"{seed.readout_type}."
    )


def _expected_benchmark_impact(seed: AmrScoutSeed) -> str:
    return (
        f"Repairs or completes the AMR contrast row for {seed.benchmark_label} "
        f"by returning: {seed.missing_measurement}"
    )


def _uncertainty_note(seed: AmrScoutSeed) -> str:
    suffix = seed.notes or "Evidence-scout record only."
    return (
        f"Scout seed only; {seed.rationale} {suffix} This is surveillance, "
        "stewardship, or benchmark evidence rather than individual-patient "
        "treatment advice."
    )


_STATUS_ALIASES = {
    "single-side-only": "single_side_only",
    "paired-complete": "paired_complete",
}


def _coerce_evidence_status(value: str) -> EvidenceStatus:
    normalized = _normalize_status(value)
    try:
        return EvidenceStatus(normalized)
    except ValueError as exc:
        raise ValueError(f"unknown evidence status: {value!r}") from exc


def _normalize_status(value: str) -> str:
    normalized = value.strip().lower().replace("-", "_")
    return _STATUS_ALIASES.get(normalized, normalized)


def _normalize_field(value: str) -> str:
    return value.strip().lower().replace("/", "_").replace("-", "_").replace(" ", "_")


def _status_values() -> tuple[str, ...]:
    values: list[str] = []
    for status in EvidenceStatus:
        if isinstance(status, Enum):
            values.append(str(status.value))
    return tuple(values)


assert _STATUS_VALUES.issubset(set(_status_values()))
