from __future__ import annotations

import pytest

from kira.amr import (
    amr_scout_seeds,
    build_amr_scout,
    build_contrast_specs,
    build_evidence_records,
    build_experiment_tickets,
    summarize_evidence_statuses,
    tickets_as_dicts,
    validate_amr_seed,
)
from kira.amr.scout import AmrScoutSeed
from kira.contrast import (
    ContrastSpec,
    EvidenceRecord,
    ExperimentTicket,
    ticket_from_dict,
    ticket_to_dict,
)
from kira.contrast.schemas import EvidenceStatus

EXPECTED_INTERVENTION_FAMILIES = {
    "stewardship action",
    "diagnostic/AST measurement",
    "infection-control measure",
    "surveillance sampling",
    "antibiotic exposure policy",
    "lab validation assay",
}

EXPECTED_READOUT_TYPES = {
    "AST result",
    "resistance rate",
    "antibiotic consumption",
    "treatment outcome",
    "infection-control metric",
    "genomic resistance marker",
    "surveillance completeness",
}


def test_amr_scout_seeds_are_fixed_narrow_and_valid() -> None:
    seeds = amr_scout_seeds()

    assert len(seeds) == 8
    assert {seed.intervention_family for seed in seeds} == EXPECTED_INTERVENTION_FAMILIES
    assert {seed.readout_type for seed in seeds} == EXPECTED_READOUT_TYPES
    assert all(seed.contrast_id.startswith("rwanda-amr-") for seed in seeds)
    assert all(seed.benchmark_label for seed in seeds)
    assert all(seed.missing_measurement for seed in seeds)

    for seed in seeds:
        validate_amr_seed(seed)


def test_amr_scout_status_distribution() -> None:
    status_counts = summarize_evidence_statuses()

    assert status_counts == {
        "paired_complete": 2,
        "single_side_only": 2,
        "proposed": 2,
        "bounded": 1,
        "conflicting": 1,
    }


def test_amr_scout_rejects_broad_or_unmeasured_seed() -> None:
    bad_seed = AmrScoutSeed(
        contrast_id="bad-amr-hype",
        intervention_family="platform claim",
        intervention="",
        desired_context="AMR solved",
        control_context="",
        readout_type="press release",
        evidence_status="validated",
        missing_measurement="",
        benchmark_label="",
        priority=9,
        rationale="",
        source_label="",
        source_url="http://example.invalid",
        rwanda_context="",
    )

    with pytest.raises(ValueError):
        validate_amr_seed(bad_seed)


def test_amr_scout_builds_contrast_core_objects() -> None:
    specs = build_contrast_specs()
    evidence_records = build_evidence_records()
    tickets = build_experiment_tickets()
    scout_records = build_amr_scout()

    assert len(specs) == 8
    assert len(evidence_records) == 8
    assert len(tickets) == 8
    assert len(scout_records) == 8
    assert all(isinstance(spec, ContrastSpec) for spec in specs)
    assert all(isinstance(record, EvidenceRecord) for record in evidence_records)
    assert all(isinstance(ticket, ExperimentTicket) for ticket in tickets)
    assert {spec.domain for spec in specs} == {"rwanda_amr"}
    assert all(isinstance(spec.evidence_status, EvidenceStatus) for spec in specs)


def test_amr_scout_tickets_are_ranked_and_round_trip() -> None:
    scout_records = build_amr_scout()
    tickets = build_experiment_tickets()
    ticket_dicts = [ticket_to_dict(ticket) for ticket in tickets]
    restored = [ticket_from_dict(ticket_dict) for ticket_dict in ticket_dicts]

    assert [record.seed.priority for record in scout_records] == sorted(
        [record.seed.priority for record in scout_records], reverse=True
    )
    assert [record.seed.benchmark_label for record in scout_records[:4]] == [
        "priority_isolate_ast_completeness",
        "genotype_phenotype_marker_confirmation",
        "pathogen_antibiotic_pair_record",
        "ast_linked_stewardship_evidence_gap",
    ]
    assert ticket_dicts == list(tickets_as_dicts())
    assert restored == list(tickets)
    assert {ticket_dict["evidence_status"] for ticket_dict in ticket_dicts} == {
        "proposed",
        "single_side_only",
        "paired_complete",
        "bounded",
        "conflicting",
    }


def test_amr_scout_does_not_encode_clinical_or_grand_claims() -> None:
    forbidden_claims = {
        "clinical recommendation",
        "prescribe",
        "solved amr",
        "wet-lab validated",
        "cure",
    }
    scout_text = repr(amr_scout_seeds()).lower()
    ticket_text = repr(tickets_as_dicts()).lower()

    for claim in forbidden_claims:
        assert claim not in scout_text
        assert claim not in ticket_text
