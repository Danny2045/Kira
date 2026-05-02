from __future__ import annotations

import pytest

from kira.contrast import ticket_from_dict, ticket_to_dict
from kira.regeneration import (
    build_contrast_specs,
    build_evidence_records,
    build_experiment_tickets,
    build_regeneration_scout,
    regeneration_scout_seeds,
    summarize_evidence_statuses,
    tickets_as_dicts,
    validate_regeneration_seed,
)
from kira.regeneration.scout import RegenerationScoutSeed


def test_regeneration_scout_seeds_are_narrow_and_valid() -> None:
    seeds = regeneration_scout_seeds()

    assert len(seeds) == 6
    assert {seed.intervention_family for seed in seeds} == {
        "bioelectric perturbation",
        "crispr perturbation",
        "partial reprogramming",
        "morphogen/control signal",
    }
    assert any("Michael Levin" in seed.source_label for seed in seeds)
    assert any("Yamanaka" in seed.source_label for seed in seeds)
    assert all(seed.missing_measurement for seed in seeds)

    for seed in seeds:
        validate_regeneration_seed(seed)


def test_regeneration_scout_status_distribution() -> None:
    status_counts = summarize_evidence_statuses()

    assert status_counts == {
        "paired": 2,
        "bounded": 1,
        "single-side-only": 1,
        "conflicting": 1,
        "proposed": 1,
    }


def test_regeneration_scout_rejects_broad_or_unmeasured_seed() -> None:
    bad_seed = RegenerationScoutSeed(
        contrast_id="bad-regeneration-hype",
        intervention_family="magic platform",
        intervention="",
        desired_context="regeneration solved",
        control_context="",
        readout_type="press release",
        evidence_status="validated",
        missing_measurement="",
        benchmark_label="",
        priority=9,
        rationale="",
        source_label="",
        source_url="http://example.invalid",
        model_system="",
    )

    with pytest.raises(ValueError):
        validate_regeneration_seed(bad_seed)


def test_regeneration_scout_builds_contrast_core_objects() -> None:
    specs = build_contrast_specs()
    evidence_records = build_evidence_records()
    tickets = build_experiment_tickets()
    scout_records = build_regeneration_scout()

    assert len(specs) == 6
    assert len(evidence_records) == 6
    assert len(tickets) == 6
    assert len(scout_records) == 6


def test_regeneration_scout_tickets_are_ranked_and_round_trip() -> None:
    tickets = build_experiment_tickets()
    ticket_dicts = [ticket_to_dict(ticket) for ticket in tickets]
    restored = [ticket_from_dict(ticket_dict) for ticket_dict in ticket_dicts]

    assert ticket_dicts == list(tickets_as_dicts())
    assert len(restored) == len(tickets)
    assert "partial_reprogramming_safety_window" in repr(ticket_dicts)
    assert "organoid_crispr_migration_rescue" in repr(ticket_dicts)


def test_regeneration_scout_does_not_encode_grand_claims() -> None:
    forbidden_claims = {
        "discovered drugs",
        "wet-lab validated",
        "solved regeneration",
        "solved crispr",
        "clinical cure",
    }
    scout_text = repr(regeneration_scout_seeds()).lower()

    for claim in forbidden_claims:
        assert claim not in scout_text
