"""Tests for kira.chemistry — ADMET filtering and drug-likeness."""


import pytest

from kira.chemistry import (
    classify_drug_likeness,
    compute_admet_multiplier,
    count_lipinski_violations,
)


class TestCountLipinskiViolations:
    """Count violations of Lipinski's Rule of Five."""

    def test_no_violations(self):
        # All values within limits
        assert count_lipinski_violations(mw=300, logp=2.0, hbd=2, hba=5) == 0

    def test_mw_violation(self):
        assert count_lipinski_violations(mw=600, logp=2.0, hbd=2, hba=5) == 1

    def test_logp_violation(self):
        assert count_lipinski_violations(mw=300, logp=6.0, hbd=2, hba=5) == 1

    def test_hbd_violation(self):
        assert count_lipinski_violations(mw=300, logp=2.0, hbd=6, hba=5) == 1

    def test_hba_violation(self):
        assert count_lipinski_violations(mw=300, logp=2.0, hbd=2, hba=12) == 1

    def test_all_violations(self):
        assert count_lipinski_violations(mw=600, logp=6.0, hbd=6, hba=12) == 4

    def test_boundary_mw(self):
        # Exactly at boundary — should NOT violate
        assert count_lipinski_violations(mw=500, logp=2.0, hbd=2, hba=5) == 0

    def test_boundary_logp(self):
        assert count_lipinski_violations(mw=300, logp=5.0, hbd=2, hba=5) == 0

    def test_boundary_hbd(self):
        assert count_lipinski_violations(mw=300, logp=2.0, hbd=5, hba=5) == 0

    def test_boundary_hba(self):
        assert count_lipinski_violations(mw=300, logp=2.0, hbd=2, hba=10) == 0


class TestClassifyDrugLikeness:
    def test_excellent(self):
        assert classify_drug_likeness(0) == "excellent"

    def test_acceptable(self):
        assert classify_drug_likeness(1) == "acceptable"

    def test_marginal(self):
        assert classify_drug_likeness(2) == "marginal"

    def test_poor(self):
        assert classify_drug_likeness(3) == "poor"
        assert classify_drug_likeness(4) == "poor"


class TestComputeAdmetMultiplier:
    def test_no_violations(self):
        mult = compute_admet_multiplier(0, 100.0)
        assert mult == 1.0

    def test_one_violation(self):
        mult = compute_admet_multiplier(1, 100.0)
        assert mult == pytest.approx(0.85)

    def test_two_violations(self):
        mult = compute_admet_multiplier(2, 100.0)
        assert mult == pytest.approx(0.60)

    def test_three_violations(self):
        mult = compute_admet_multiplier(3, 100.0)
        assert mult == pytest.approx(0.35)

    def test_high_tpsa_penalty(self):
        # 0 violations but TPSA > 140 -> 1.0 * 0.75 = 0.75
        mult = compute_admet_multiplier(0, 160.0)
        assert mult == pytest.approx(0.75)

    def test_violations_plus_high_tpsa(self):
        # 1 violation + TPSA > 140 -> 0.85 * 0.75 = 0.6375
        mult = compute_admet_multiplier(1, 160.0)
        assert mult == pytest.approx(0.6375)

    def test_none_violations_returns_default(self):
        mult = compute_admet_multiplier(None, 100.0)
        assert mult == pytest.approx(0.8)

    def test_nan_violations_returns_default(self):
        mult = compute_admet_multiplier(float("nan"), 100.0)
        assert mult == pytest.approx(0.8)

    def test_none_tpsa_no_penalty(self):
        mult = compute_admet_multiplier(0, None)
        assert mult == 1.0

    def test_multiplier_between_zero_and_one(self):
        """Multiplier should always be in [0, 1]."""
        for v in range(5):
            for tpsa in [50.0, 100.0, 150.0, 200.0]:
                mult = compute_admet_multiplier(v, tpsa)
                assert 0.0 <= mult <= 1.0
