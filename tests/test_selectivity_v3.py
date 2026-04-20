"""Tests for selectivity v3 feature extraction."""

import numpy as np
import pytest

from kira.experiments import TARGET_PAIRS, get_pairs_for_disease, get_target_pair
from kira.experiments.selectivity_features import (
    PocketFeatures,
    compute_pocket_features,
)


class TestTargetPairs:
    """Test target pair definitions."""

    def test_all_pairs_have_sequences(self):
        for pair_id, pair in TARGET_PAIRS.items():
            assert len(pair.parasite_pocket_sequence) > 0, f"{pair_id} missing parasite pocket seq"
            assert len(pair.human_pocket_sequence) > 0, f"{pair_id} missing human pocket seq"

    def test_pocket_sequences_same_length(self):
        for pair_id, pair in TARGET_PAIRS.items():
            assert len(pair.parasite_pocket_sequence) == len(pair.human_pocket_sequence), \
                f"{pair_id}: pocket sequences differ in length"

    def test_three_diseases_covered(self):
        diseases = {p.disease for p in TARGET_PAIRS.values()}
        assert "Schistosomiasis" in diseases
        assert "Trypanosomiasis" in diseases
        assert "Leishmaniasis" in diseases

    def test_get_pair_by_name(self):
        pair = get_target_pair("SmDHODH")
        assert pair.disease == "Schistosomiasis"

    def test_get_pairs_for_disease(self):
        leish = get_pairs_for_disease("Leishmaniasis")
        assert len(leish) >= 2

    def test_unknown_pair_raises(self):
        with pytest.raises(KeyError):
            get_target_pair("NotARealTarget")


class TestPocketFeatures:
    """Test pocket feature computation."""

    def test_identical_pockets_zero_divergence(self):
        """Identical sequences should have zero divergence."""
        pf = compute_pocket_features("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY")

        assert pf.pocket_identity == 1.0
        assert pf.pocket_divergence == 0.0
        assert pf.mean_physicochemical_distance == 0.0
        assert pf.n_charge_changes == 0
        assert pf.n_hydrophobicity_flips == 0
        assert pf.n_volume_changes == 0

    def test_completely_different_pockets(self):
        """Maximally different sequences should have high divergence."""
        pf = compute_pocket_features("AAAA", "RRRR")

        assert pf.pocket_identity == 0.0
        assert pf.pocket_divergence == 1.0
        assert pf.mean_physicochemical_distance > 0.0

    def test_charge_change_detection(self):
        """Should detect charge changes (D/E → K/R or reverse)."""
        # D (negative) → K (positive) = charge change
        pf = compute_pocket_features("D", "K")
        assert pf.n_charge_changes == 1

    def test_hydrophobicity_flip_detection(self):
        """Should detect polar/nonpolar flips."""
        # I (hydrophobic, 4.5) → D (hydrophilic, -3.5)
        pf = compute_pocket_features("I", "D")
        assert pf.n_hydrophobicity_flips == 1

    def test_volume_change_detection(self):
        """Should detect large volume changes."""
        # G (60.1 A³) → W (227.8 A³) = 167.7 A³ difference
        pf = compute_pocket_features("G", "W")
        assert pf.n_volume_changes == 1

    def test_feature_array_shape(self):
        """Feature vector should have 15 dimensions."""
        pf = compute_pocket_features("ACDEF", "ACDEF")
        arr = pf.to_array()
        assert arr.shape == (15,)
        assert arr.dtype == np.float32

    def test_feature_names_match_array(self):
        """Feature names should match array length."""
        names = PocketFeatures.feature_names()
        assert len(names) == 15

    def test_ptr1_high_divergence(self):
        """PTR1 (different enzyme family from DHFR) should show high divergence."""
        pair = get_target_pair("LmPTR1")
        pf = compute_pocket_features(
            pair.parasite_pocket_sequence,
            pair.human_pocket_sequence,
        )
        # Different enzyme families → low identity
        assert pf.pocket_identity < 0.8

    def test_hdac8_low_divergence(self):
        """SmHDAC8 (conserved active site) should show low divergence."""
        pair = get_target_pair("SmHDAC8")
        pf = compute_pocket_features(
            pair.parasite_pocket_sequence,
            pair.human_pocket_sequence,
        )
        # Conserved active site → high identity
        assert pf.pocket_identity > 0.8

    def test_divergence_correlates_with_selectivity(self):
        """Target pairs with higher selectivity should have higher pocket divergence.

        PTR1 (75.6% selective) should be more divergent than HDAC8 (1.4% selective).
        This is the core hypothesis of the experiment.
        """
        ptr1 = get_target_pair("LmPTR1")
        hdac8 = get_target_pair("SmHDAC8")

        pf_ptr1 = compute_pocket_features(
            ptr1.parasite_pocket_sequence, ptr1.human_pocket_sequence,
        )
        pf_hdac8 = compute_pocket_features(
            hdac8.parasite_pocket_sequence, hdac8.human_pocket_sequence,
        )

        # PTR1 (high selectivity) should have HIGHER divergence
        # than HDAC8 (low selectivity)
        assert pf_ptr1.pocket_divergence > pf_hdac8.pocket_divergence, \
            f"PTR1 divergence ({pf_ptr1.pocket_divergence:.3f}) should be > " \
            f"HDAC8 ({pf_hdac8.pocket_divergence:.3f})"

    def test_empty_sequences(self):
        """Should handle empty sequences gracefully."""
        pf = compute_pocket_features("", "")
        assert pf.pocket_size == 0
