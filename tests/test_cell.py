"""Tests for the virtual cell context layer.

Tests expression profiles, perturbation response modeling,
selectivity-expression integration, and on/off-target comparison.
"""

from __future__ import annotations

import pytest

from kira.cell.expression import (
    EXPRESSION_PROFILES,
    ExpressionLevel,
    GeneExpressionProfile,
    SelectivityExpressionContext,
    TissueExpression,
    compute_selectivity_expression_context,
    get_expression_profile,
)
from kira.cell.perturbation import (
    CellStateChange,
    Perturbation,
    PerturbationType,
    compare_on_off_target,
    predict_perturbation_response,
)


# ---------------------------------------------------------------------------
# Expression profile tests
# ---------------------------------------------------------------------------


class TestExpressionLevel:
    """Tests for expression level classification."""

    def test_not_detected(self) -> None:
        assert ExpressionLevel.from_tpm(0.0) == ExpressionLevel.NOT_DETECTED
        assert ExpressionLevel.from_tpm(0.5) == ExpressionLevel.NOT_DETECTED

    def test_low(self) -> None:
        assert ExpressionLevel.from_tpm(1.0) == ExpressionLevel.LOW
        assert ExpressionLevel.from_tpm(5.0) == ExpressionLevel.LOW

    def test_medium(self) -> None:
        assert ExpressionLevel.from_tpm(10.0) == ExpressionLevel.MEDIUM
        assert ExpressionLevel.from_tpm(50.0) == ExpressionLevel.MEDIUM

    def test_high(self) -> None:
        assert ExpressionLevel.from_tpm(100.0) == ExpressionLevel.HIGH
        assert ExpressionLevel.from_tpm(500.0) == ExpressionLevel.HIGH


class TestGeneExpressionProfile:
    """Tests for gene expression profiles."""

    def test_dhodh_profile_exists(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert profile.gene_name == "DHODH"
        assert profile.organism == "Homo sapiens"

    def test_hdac8_profile_exists(self) -> None:
        profile = get_expression_profile("HDAC8")
        assert profile is not None
        assert profile.gene_name == "HDAC8"

    def test_txnrd1_profile_exists(self) -> None:
        profile = get_expression_profile("TXNRD1")
        assert profile is not None
        assert profile.gene_name == "TXNRD1"

    def test_case_insensitive_lookup(self) -> None:
        assert get_expression_profile("dhodh") is not None
        assert get_expression_profile("Dhodh") is not None

    def test_unknown_gene_returns_none(self) -> None:
        assert get_expression_profile("FAKE_GENE_XYZ") is None

    def test_n_tissues_expressed(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert profile.n_tissues_expressed > 5

    def test_max_tpm_positive(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert profile.max_tpm > 0

    def test_tissue_specificity_between_0_and_1(self) -> None:
        for name, profile in EXPRESSION_PROFILES.items():
            spec = profile.tissue_specificity
            assert 0 <= spec <= 1, f"{name} specificity out of range: {spec}"

    def test_ubiquitous_genes_low_specificity(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert profile.ubiquitous
        # Ubiquitous genes should have lower specificity (< 0.7)
        assert profile.tissue_specificity < 0.7

    def test_get_expression_in_tissue(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        liver = profile.get_expression_in("liver")
        assert liver is not None
        assert liver.tpm > 0

    def test_get_expression_in_unknown_tissue(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert profile.get_expression_in("mars_colony") is None

    def test_summary_not_empty(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert len(profile.summary()) > 20

    def test_to_dict_serialization(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        d = profile.to_dict()
        assert d["gene_name"] == "DHODH"
        assert "tissue_expression" in d
        assert len(d["tissue_expression"]) > 0

    def test_primary_tissues_set(self) -> None:
        profile = get_expression_profile("DHODH")
        assert profile is not None
        assert len(profile.primary_tissues) > 0

    def test_empty_profile_specificity(self) -> None:
        profile = GeneExpressionProfile(
            gene_name="EMPTY",
            organism="test",
        )
        assert profile.n_tissues_expressed == 0
        assert profile.max_tpm == 0.0


# ---------------------------------------------------------------------------
# Selectivity-expression integration tests
# ---------------------------------------------------------------------------


class TestSelectivityExpressionContext:
    """Tests for selectivity-expression integration."""

    def test_high_selectivity_low_risk(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmDHODH", "DHODH", selectivity_ratio=30.0,
        )
        assert ctx is not None
        assert ctx.selectivity_ratio == 30.0
        # All tissue risks should be low with 30x selectivity
        assert all(r < 0.5 for r in ctx.tissue_risk.values())

    def test_low_selectivity_high_risk(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmDHODH", "DHODH", selectivity_ratio=1.0,
        )
        assert ctx is not None
        # Some tissues should have high risk with 1x selectivity
        assert max(ctx.tissue_risk.values()) > 0.5

    def test_safest_and_riskiest_tissues(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmDHODH", "DHODH", selectivity_ratio=10.0,
        )
        assert ctx is not None
        assert len(ctx.safest_delivery_tissues) > 0
        assert len(ctx.riskiest_tissues) > 0
        # Riskiest should be where expression is highest
        assert ctx.riskiest_tissues[0] in [
            te.tissue
            for te in ctx.human_profile.tissue_expression
        ]

    def test_unknown_orthologue_returns_none(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmFAKE", "FAKEGENE", selectivity_ratio=10.0,
        )
        assert ctx is None

    def test_summary_not_empty(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmDHODH", "DHODH", selectivity_ratio=30.0,
        )
        assert ctx is not None
        assert len(ctx.summary()) > 20

    def test_to_dict_serialization(self) -> None:
        ctx = compute_selectivity_expression_context(
            "SmDHODH", "DHODH", selectivity_ratio=30.0,
        )
        assert ctx is not None
        d = ctx.to_dict()
        assert "tissue_risk" in d
        assert "safest_delivery_tissues" in d


# ---------------------------------------------------------------------------
# Perturbation response tests
# ---------------------------------------------------------------------------


class TestPerturbationResponse:
    """Tests for perturbation response modeling."""

    def test_dhodh_inhibition_kills_parasite(self) -> None:
        pert = Perturbation(
            name="test",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmDHODH",
            concentration_nm=1000.0,
            selectivity_ratio=1.0,
        )
        response = predict_perturbation_response(
            pert, "schistosomulum", "Schistosoma mansoni",
        )
        # High concentration should reduce viability
        assert response.viability < 0.8
        assert len(response.pathway_effects) > 0

    def test_selective_compound_spares_human(self) -> None:
        pert = Perturbation(
            name="selective",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="HsDHODH",
            concentration_nm=1000.0,
            selectivity_ratio=30.0,  # 30x selective for parasite
        )
        response = predict_perturbation_response(
            pert, "hepatocyte", "Homo sapiens",
        )
        # With 30x selectivity, human cells should be mostly spared
        assert response.viability > 0.5

    def test_zero_concentration_no_effect(self) -> None:
        pert = Perturbation(
            name="zero",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmDHODH",
            concentration_nm=0.0,
        )
        response = predict_perturbation_response(pert)
        assert response.predicted_state_change == CellStateChange.UNCHANGED
        assert response.viability > 0.9

    def test_pathway_effects_have_correct_fields(self) -> None:
        pert = Perturbation(
            name="fields",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmDHODH",
            concentration_nm=500.0,
        )
        response = predict_perturbation_response(pert)
        for pe in response.pathway_effects:
            assert pe.pathway_name
            assert pe.direction in {"up", "down", "mixed"}
            assert 0 <= pe.magnitude <= 1
            assert 0 <= pe.confidence <= 1
            assert pe.mechanism

    def test_tgr_inhibition_pathway_effects(self) -> None:
        pert = Perturbation(
            name="tgr",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmTGR",
            concentration_nm=1000.0,
        )
        response = predict_perturbation_response(
            pert, "schistosomulum", "Schistosoma mansoni",
        )
        pathway_names = [pe.pathway_name for pe in response.pathway_effects]
        assert "thioredoxin_glutathione_redox" in pathway_names
        assert "reactive_oxygen_species" in pathway_names

    def test_hdac8_inhibition_pathways(self) -> None:
        pert = Perturbation(
            name="hdac8",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmHDAC8",
            concentration_nm=500.0,
        )
        response = predict_perturbation_response(pert)
        pathway_names = [pe.pathway_name for pe in response.pathway_effects]
        assert "chromatin_remodeling" in pathway_names

    def test_to_dict_serialization(self) -> None:
        pert = Perturbation(
            name="dict_test",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmDHODH",
            concentration_nm=100.0,
        )
        response = predict_perturbation_response(pert)
        d = response.to_dict()
        assert "perturbation" in d
        assert "pathway_effects" in d
        assert "viability" in d

    def test_summary_string(self) -> None:
        pert = Perturbation(
            name="summary",
            perturbation_type=PerturbationType.CHEMICAL,
            target_gene="SmDHODH",
            concentration_nm=100.0,
        )
        response = predict_perturbation_response(pert)
        assert len(response.summary) > 20


# ---------------------------------------------------------------------------
# On/off-target comparison tests
# ---------------------------------------------------------------------------


class TestOnOffTargetComparison:
    """Tests for comparing on-target vs off-target effects."""

    def test_selective_compound_good_window(self) -> None:
        result = compare_on_off_target(
            compound_name="CHEMBL155771",
            compound_id="CHEMBL155771",
            parasite_target="SmDHODH",
            human_orthologue="HsDHODH",
            parasite_ic50_nm=50.0,
            selectivity_ratio=30.0,
        )
        assert result["selectivity_ratio"] == 30.0
        # On-target should have lower viability than off-target
        assert result["on_target_viability"] < result["off_target_viability"]

    def test_non_selective_compound_poor_window(self) -> None:
        result = compare_on_off_target(
            compound_name="NonSelective",
            compound_id="TEST001",
            parasite_target="SmDHODH",
            human_orthologue="HsDHODH",
            parasite_ic50_nm=50.0,
            selectivity_ratio=1.0,
        )
        assert result["selectivity_at_cell_level"] in {"POOR", "MODERATE"}

    def test_therapeutic_concentration(self) -> None:
        result = compare_on_off_target(
            compound_name="test",
            compound_id="TEST",
            parasite_target="SmDHODH",
            human_orthologue="HsDHODH",
            parasite_ic50_nm=100.0,
            selectivity_ratio=10.0,
        )
        # Therapeutic concentration should be 10x IC50
        assert result["therapeutic_concentration_nm"] == 1000.0

    def test_result_has_all_keys(self) -> None:
        result = compare_on_off_target(
            compound_name="keys",
            compound_id="TEST",
            parasite_target="SmDHODH",
            human_orthologue="HsDHODH",
            parasite_ic50_nm=50.0,
            selectivity_ratio=10.0,
        )
        expected_keys = {
            "compound", "compound_id", "therapeutic_concentration_nm",
            "selectivity_ratio", "on_target", "off_target",
            "therapeutic_window", "on_target_viability",
            "off_target_viability", "selectivity_at_cell_level",
        }
        assert expected_keys.issubset(set(result.keys()))
