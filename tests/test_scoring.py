"""Tests for kira.scoring — the composite ranking algorithm's signal functions."""

import math

import pytest

from kira.scoring import (
    WEIGHTS_V1,
    WEIGHTS_V2,
    WEIGHTS_V3,
    classify_novelty,
    classify_selectivity,
    compute_confidence_score,
    compute_drug_stage_score,
    compute_potency_score,
    compute_target_score,
)

# ---------------------------------------------------------------------------
# Weight dictionaries
# ---------------------------------------------------------------------------

class TestWeights:
    """Verify weight dictionaries sum to 1.0 and have expected keys."""

    def test_v1_weights_sum_to_one(self):
        assert math.isclose(sum(WEIGHTS_V1.values()), 1.0, rel_tol=1e-9)

    def test_v2_weights_sum_to_one(self):
        assert math.isclose(sum(WEIGHTS_V2.values()), 1.0, rel_tol=1e-9)

    def test_v3_weights_sum_to_one(self):
        assert math.isclose(sum(WEIGHTS_V3.values()), 1.0, rel_tol=1e-9)

    def test_v1_has_expected_keys(self):
        expected = {"potency", "target", "confidence", "drug_stage", "multitarget"}
        assert set(WEIGHTS_V1.keys()) == expected

    def test_v2_adds_similarity(self):
        assert "similarity" in WEIGHTS_V2

    def test_v3_adds_whole_org(self):
        assert "whole_org" in WEIGHTS_V3


# ---------------------------------------------------------------------------
# Potency score
# ---------------------------------------------------------------------------

class TestPotencyScore:
    """compute_potency_score converts IC50 (nM) -> 0-1 via pIC50."""

    def test_1nm_gives_max_score(self):
        score = compute_potency_score(1.0)
        assert score == 1.0

    def test_100000nm_gives_zero(self):
        score = compute_potency_score(100_000.0)
        assert score == pytest.approx(0.0, abs=1e-6)

    def test_1000nm_gives_midrange(self):
        # pIC50 = 6.0 -> (6-4)/(9-4) = 0.4
        score = compute_potency_score(1000.0)
        assert score == pytest.approx(0.4, abs=1e-6)

    def test_10nm_gives_high_score(self):
        # pIC50 = 8.0 -> (8-4)/5 = 0.8
        score = compute_potency_score(10.0)
        assert score == pytest.approx(0.8, abs=1e-6)

    def test_nan_returns_zero(self):
        assert compute_potency_score(float("nan")) == 0.0

    def test_zero_returns_zero(self):
        assert compute_potency_score(0.0) == 0.0

    def test_negative_returns_zero(self):
        assert compute_potency_score(-100.0) == 0.0

    def test_very_potent_clamped_at_one(self):
        # IC50 = 0.01 nM -> pIC50 = 11 -> (11-4)/5 = 1.4 -> clamped to 1.0
        assert compute_potency_score(0.01) == 1.0

    def test_very_weak_clamped_at_zero(self):
        # IC50 = 1e8 nM -> pIC50 = 1 -> (1-4)/5 = -0.6 -> clamped to 0.0
        assert compute_potency_score(1e8) == 0.0

    def test_monotonically_decreasing(self):
        """More potent (lower IC50) should give higher scores."""
        scores = [compute_potency_score(ic50) for ic50 in [1, 10, 100, 1000, 10000]]
        for i in range(len(scores) - 1):
            assert scores[i] > scores[i + 1]


# ---------------------------------------------------------------------------
# Target essentiality score
# ---------------------------------------------------------------------------

class TestTargetScore:
    def test_known_target(self):
        assert compute_target_score("Thioredoxin glutathione reductase") == 1.0

    def test_hdac8(self):
        assert compute_target_score("Histone deacetylase 8") == 0.85

    def test_unknown_target_gets_default(self):
        score = compute_target_score("Some unknown target")
        assert score == 0.3  # DEFAULT_ESSENTIALITY

    def test_negative_control(self):
        assert compute_target_score("None (negative control)") == 0.0


# ---------------------------------------------------------------------------
# Confidence score
# ---------------------------------------------------------------------------

class TestConfidenceScore:
    def test_zero_measurements(self):
        assert compute_confidence_score(0) == 0.1

    def test_one_measurement(self):
        score = compute_confidence_score(1)
        assert 0.1 <= score <= 0.3

    def test_ten_measurements_near_one(self):
        score = compute_confidence_score(10)
        assert score >= 0.75

    def test_hundred_measurements_capped(self):
        score = compute_confidence_score(100)
        assert score == 1.0

    def test_negative_measurements(self):
        assert compute_confidence_score(-1) == 0.1

    def test_monotonically_increasing(self):
        """More measurements should give higher confidence."""
        scores = [compute_confidence_score(n) for n in [1, 2, 5, 10, 50]]
        for i in range(len(scores) - 1):
            assert scores[i] <= scores[i + 1]


# ---------------------------------------------------------------------------
# Drug stage score
# ---------------------------------------------------------------------------

class TestDrugStageScore:
    def test_approved_drug(self):
        assert compute_drug_stage_score(4.0) == 1.0

    def test_phase_three(self):
        assert compute_drug_stage_score(3.0) == 0.7

    def test_phase_two(self):
        assert compute_drug_stage_score(2.0) == 0.5

    def test_phase_one(self):
        assert compute_drug_stage_score(1.0) == 0.3

    def test_preclinical(self):
        assert compute_drug_stage_score(0.0) == 0.1

    def test_nan_returns_default(self):
        assert compute_drug_stage_score(float("nan")) == 0.1

    def test_unknown_phase_returns_default(self):
        assert compute_drug_stage_score(99.0) == 0.1


# ---------------------------------------------------------------------------
# Selectivity classification
# ---------------------------------------------------------------------------

class TestClassifySelectivity:
    def test_selective(self):
        assert classify_selectivity(15.0) == "SELECTIVE"

    def test_moderate(self):
        assert classify_selectivity(5.0) == "MODERATE"

    def test_poor(self):
        assert classify_selectivity(1.5) == "POOR"

    def test_counter_selective(self):
        assert classify_selectivity(0.5) == "COUNTER-SELECTIVE"

    def test_boundary_selective(self):
        assert classify_selectivity(10.0) == "SELECTIVE"

    def test_boundary_moderate(self):
        assert classify_selectivity(3.0) == "MODERATE"

    def test_boundary_poor(self):
        assert classify_selectivity(1.0) == "POOR"


# ---------------------------------------------------------------------------
# Novelty classification
# ---------------------------------------------------------------------------

class TestClassifyNovelty:
    def test_novel(self):
        assert classify_novelty(0) == "NOVEL"

    def test_emerging(self):
        assert classify_novelty(1) == "EMERGING"
        assert classify_novelty(3) == "EMERGING"

    def test_known(self):
        assert classify_novelty(4) == "KNOWN"
        assert classify_novelty(10) == "KNOWN"

    def test_well_known(self):
        assert classify_novelty(11) == "WELL-KNOWN"
        assert classify_novelty(100) == "WELL-KNOWN"

    def test_error_on_negative(self):
        assert classify_novelty(-1) == "ERROR"
