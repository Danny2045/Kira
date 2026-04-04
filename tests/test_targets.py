"""Tests for kira.targets — target essentiality scores and orthologue mappings."""

from kira.targets import (
    DEFAULT_ESSENTIALITY,
    ORTHOLOGUE_MAP,
    TARGET_ESSENTIALITY,
    compute_target_score,
)


class TestTargetEssentiality:
    """Verify the TARGET_ESSENTIALITY dictionary is well-formed."""

    def test_all_scores_between_zero_and_one(self):
        for target, score in TARGET_ESSENTIALITY.items():
            assert 0.0 <= score <= 1.0, f"{target} has score {score} outside [0, 1]"

    def test_smtgr_is_highest(self):
        """SmTGR is the best-validated target — should score 1.0."""
        assert TARGET_ESSENTIALITY["Thioredoxin glutathione reductase"] == 1.0

    def test_negative_control_is_zero(self):
        assert TARGET_ESSENTIALITY["None (negative control)"] == 0.0

    def test_default_essentiality_is_reasonable(self):
        assert 0.0 < DEFAULT_ESSENTIALITY < 1.0

    def test_at_least_ten_targets(self):
        assert len(TARGET_ESSENTIALITY) >= 10


class TestComputeTargetScore:
    def test_known_target_returns_essentiality(self):
        score = compute_target_score("Thioredoxin glutathione reductase")
        assert score == 1.0

    def test_unknown_target_returns_default(self):
        score = compute_target_score("Nonexistent target")
        assert score == DEFAULT_ESSENTIALITY


class TestOrthologueMap:
    """Verify the orthologue mapping structure."""

    def test_has_hdac8(self):
        assert "Histone deacetylase 8" in ORTHOLOGUE_MAP

    def test_has_smtgr(self):
        assert "Thioredoxin glutathione reductase" in ORTHOLOGUE_MAP

    def test_has_dhodh(self):
        assert "Dihydroorotate dehydrogenase (quinone), mitochondrial" in ORTHOLOGUE_MAP

    def test_each_entry_has_required_keys(self):
        required = {"parasite_id", "human_name", "human_id", "notes"}
        for target, info in ORTHOLOGUE_MAP.items():
            assert required.issubset(set(info.keys())), (
                f"{target} missing keys: {required - set(info.keys())}"
            )

    def test_chembl_ids_start_with_chembl(self):
        for target, info in ORTHOLOGUE_MAP.items():
            assert info["parasite_id"].startswith("CHEMBL"), (
                f"{target} parasite_id invalid"
            )
            assert info["human_id"].startswith("CHEMBL"), (
                f"{target} human_id invalid"
            )
