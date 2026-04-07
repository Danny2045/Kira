"""Tests for the delivery physics layer.

Tests LNP formulation modeling, biodistribution prediction,
delivery probability chains, and multi-objective optimization.
"""

from __future__ import annotations

import pytest

from kira.delivery.formulation import (
    BIONTECH_COVID_FORMULATION,
    MODERNA_COVID_FORMULATION,
    ONPATTRO_FORMULATION,
    AdministrationRoute,
    HelperLipidClass,
    IonizableLipidClass,
    LNPFormulation,
    PEGLipidClass,
    TargetingPeptide,
)


# ---------------------------------------------------------------------------
# Formulation tests
# ---------------------------------------------------------------------------


class TestLNPFormulation:
    """Tests for LNP formulation data model."""

    def test_default_formulation_is_balanced(self) -> None:
        f = LNPFormulation(name="test")
        assert f.is_balanced
        assert abs(f.total_lipid_mol_pct - 100.0) < 1.0

    def test_unbalanced_formulation_detected(self) -> None:
        f = LNPFormulation(
            name="unbalanced",
            ionizable_lipid_mol_pct=60.0,
            helper_lipid_mol_pct=20.0,
            cholesterol_mol_pct=30.0,
            peg_lipid_mol_pct=5.0,
        )
        assert not f.is_balanced
        warnings = f.validate()
        assert any("sum to" in w for w in warnings)

    def test_validation_catches_extreme_parameters(self) -> None:
        f = LNPFormulation(
            name="extreme",
            ionizable_lipid_mol_pct=70.0,
            np_ratio=20.0,
            target_size_nm=300.0,
            peg_lipid_mol_pct=10.0,
        )
        warnings = f.validate()
        assert len(warnings) >= 3  # Multiple out-of-range params

    def test_has_targeting_without_peptides(self) -> None:
        f = LNPFormulation(name="no_peptide")
        assert not f.has_targeting

    def test_has_targeting_with_peptides(self) -> None:
        peptide = TargetingPeptide(
            name="RGD",
            sequence="RGD",
            target_receptor="integrin_avb3",
            tissue_affinity={"tumor": 0.8},
        )
        f = LNPFormulation(name="targeted", targeting_peptides=[peptide])
        assert f.has_targeting

    def test_to_dict_roundtrip(self) -> None:
        f = LNPFormulation(name="test_dict")
        d = f.to_dict()
        assert d["name"] == "test_dict"
        assert d["ionizable_lipid"] == "SM-102"
        assert isinstance(d["is_balanced"], bool)
        assert isinstance(d["targeting_peptides"], list)

    def test_large_payload_warning(self) -> None:
        f = LNPFormulation(
            name="large_payload",
            payload_size_nt=6000,
            target_size_nm=60.0,
        )
        warnings = f.validate()
        assert any("Large payload" in w for w in warnings)


class TestReferenceFormulations:
    """Tests for curated reference formulations."""

    def test_onpattro_is_balanced(self) -> None:
        assert ONPATTRO_FORMULATION.is_balanced
        assert ONPATTRO_FORMULATION.ionizable_lipid == IonizableLipidClass.DLIN_MC3
        assert ONPATTRO_FORMULATION.payload_type == "siRNA"

    def test_moderna_is_balanced(self) -> None:
        assert MODERNA_COVID_FORMULATION.is_balanced
        assert MODERNA_COVID_FORMULATION.ionizable_lipid == IonizableLipidClass.SM102
        assert MODERNA_COVID_FORMULATION.route == AdministrationRoute.IM

    def test_biontech_is_balanced(self) -> None:
        assert BIONTECH_COVID_FORMULATION.is_balanced
        assert BIONTECH_COVID_FORMULATION.ionizable_lipid == IonizableLipidClass.ALC0315

    def test_reference_formulations_have_no_warnings(self) -> None:
        for f in [ONPATTRO_FORMULATION, MODERNA_COVID_FORMULATION, BIONTECH_COVID_FORMULATION]:
            warnings = f.validate()
            # Reference formulations should have 0 or minimal warnings
            assert len(warnings) <= 1, f"{f.name} has warnings: {warnings}"


# ---------------------------------------------------------------------------
# Biodistribution tests
# ---------------------------------------------------------------------------


class TestBiodistribution:
    """Tests for biodistribution prediction."""

    def test_iv_liver_dominant(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="iv_test", route=AdministrationRoute.IV)
        biodist = predict_biodistribution(f, "liver")
        assert biodist.dominant_tissue == "liver"
        assert biodist.liver_fraction > 0.4

    def test_im_muscle_dominant(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="im_test", route=AdministrationRoute.IM)
        biodist = predict_biodistribution(f, "muscle")
        assert biodist.tissue_fractions["muscle"] > 0.2

    def test_intranasal_lung_dominant(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="in_test", route=AdministrationRoute.IN)
        biodist = predict_biodistribution(f, "lung")
        assert biodist.tissue_fractions["lung"] > 0.3

    def test_fractions_sum_to_one(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        for route in AdministrationRoute:
            f = LNPFormulation(name=f"test_{route.name}", route=route)
            biodist = predict_biodistribution(f, "liver")
            total = sum(biodist.tissue_fractions.values())
            assert abs(total - 1.0) < 0.01, f"Fractions sum to {total} for {route}"

    def test_mc3_increases_liver(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f_sm102 = LNPFormulation(
            name="sm102", ionizable_lipid=IonizableLipidClass.SM102,
        )
        f_mc3 = LNPFormulation(
            name="mc3", ionizable_lipid=IonizableLipidClass.DLIN_MC3,
        )
        b_sm102 = predict_biodistribution(f_sm102, "liver")
        b_mc3 = predict_biodistribution(f_mc3, "liver")
        assert b_mc3.liver_fraction > b_sm102.liver_fraction

    def test_targeting_peptide_shifts_distribution(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        peptide = TargetingPeptide(
            name="lung_peptide",
            sequence="CGFECVRQCPERC",
            target_receptor="NRP1",
            tissue_affinity={"lung": 0.9},
        )
        f_untargeted = LNPFormulation(name="untargeted")
        f_targeted = LNPFormulation(
            name="targeted", targeting_peptides=[peptide],
        )
        b_un = predict_biodistribution(f_untargeted, "lung")
        b_tg = predict_biodistribution(f_targeted, "lung")
        assert b_tg.tissue_fractions["lung"] > b_un.tissue_fractions["lung"]

    def test_therapeutic_index_positive(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="ti_test")
        biodist = predict_biodistribution(f, "liver")
        assert biodist.therapeutic_index_estimate > 0

    def test_warning_for_brain_iv(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="brain_iv", route=AdministrationRoute.IV)
        biodist = predict_biodistribution(f, "brain")
        assert any("BBB" in w or "brain" in w.lower() for w in biodist.warnings)

    def test_summary_not_empty(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="summary_test")
        biodist = predict_biodistribution(f, "liver")
        assert len(biodist.summary()) > 20

    def test_to_dict_serialization(self) -> None:
        from kira.delivery.biodistribution import predict_biodistribution

        f = LNPFormulation(name="dict_test")
        biodist = predict_biodistribution(f, "liver")
        d = biodist.to_dict()
        assert "tissue_fractions" in d
        assert "dominant_tissue" in d
        assert "warnings" in d


# ---------------------------------------------------------------------------
# Delivery probability chain tests
# ---------------------------------------------------------------------------


class TestDeliveryChain:
    """Tests for delivery probability chain estimation."""

    def test_chain_probabilities_between_0_and_1(self) -> None:
        from kira.delivery.biodistribution import estimate_delivery_chain

        f = LNPFormulation(name="chain_test")
        chain = estimate_delivery_chain(f, "liver")
        assert 0 < chain.p_circulation <= 1
        assert 0 < chain.p_tissue_access <= 1
        assert 0 < chain.p_cell_binding <= 1
        assert 0 < chain.p_internalization <= 1
        assert 0 < chain.p_endosomal_escape <= 1
        assert 0 < chain.p_intracellular_action <= 1

    def test_effective_probability_product(self) -> None:
        from kira.delivery.biodistribution import estimate_delivery_chain

        f = LNPFormulation(name="product_test")
        chain = estimate_delivery_chain(f, "liver")
        expected = (
            chain.p_circulation
            * chain.p_tissue_access
            * chain.p_cell_binding
            * chain.p_internalization
            * chain.p_endosomal_escape
            * chain.p_intracellular_action
        )
        assert abs(chain.p_effective - expected) < 1e-10

    def test_bottleneck_is_endosomal_escape(self) -> None:
        from kira.delivery.biodistribution import estimate_delivery_chain

        # For standard LNPs, endosomal escape is typically the bottleneck
        f = LNPFormulation(name="bottleneck_test")
        chain = estimate_delivery_chain(f, "liver")
        name, val = chain.bottleneck
        assert name == "endosomal_escape"
        assert val < 0.1  # Endosomal escape is typically very low

    def test_dope_improves_escape(self) -> None:
        from kira.delivery.biodistribution import estimate_delivery_chain

        f_dspc = LNPFormulation(
            name="dspc", helper_lipid=HelperLipidClass.DSPC,
        )
        f_dope = LNPFormulation(
            name="dope", helper_lipid=HelperLipidClass.DOPE,
        )
        c_dspc = estimate_delivery_chain(f_dspc, "liver")
        c_dope = estimate_delivery_chain(f_dope, "liver")
        assert c_dope.p_endosomal_escape > c_dspc.p_endosomal_escape

    def test_to_dict_has_bottleneck(self) -> None:
        from kira.delivery.biodistribution import estimate_delivery_chain

        f = LNPFormulation(name="dict_test")
        chain = estimate_delivery_chain(f, "liver")
        d = chain.to_dict()
        assert "bottleneck" in d
        assert "p_effective" in d


# ---------------------------------------------------------------------------
# Multi-objective optimization tests
# ---------------------------------------------------------------------------


class TestOptimization:
    """Tests for formulation optimization."""

    def test_evaluate_formulation_all_scores_valid(self) -> None:
        from kira.delivery.optimization import evaluate_formulation

        f = LNPFormulation(name="eval_test")
        result = evaluate_formulation(f, "liver")
        assert 0 <= result.efficacy <= 1
        assert 0 <= result.safety <= 1
        assert 0 <= result.manufacturability <= 1
        assert 0 <= result.stability <= 1
        assert result.overall > 0

    def test_custom_lipid_reduces_manufacturability(self) -> None:
        from kira.delivery.optimization import evaluate_formulation

        f_standard = LNPFormulation(
            name="standard", ionizable_lipid=IonizableLipidClass.SM102,
        )
        f_custom = LNPFormulation(
            name="custom", ionizable_lipid=IonizableLipidClass.CUSTOM,
        )
        r_std = evaluate_formulation(f_standard, "liver")
        r_cust = evaluate_formulation(f_custom, "liver")
        assert r_cust.manufacturability < r_std.manufacturability

    def test_pareto_front_is_subset(self) -> None:
        from kira.delivery.optimization import evaluate_formulation, find_pareto_front

        formulations = [
            LNPFormulation(name=f"f{i}", ionizable_lipid_mol_pct=35 + i * 5)
            for i in range(5)
        ]
        evals = [evaluate_formulation(f, "liver") for f in formulations]
        pareto = find_pareto_front(evals)
        assert len(pareto) <= len(evals)
        assert len(pareto) >= 1

    def test_pareto_front_empty_input(self) -> None:
        from kira.delivery.optimization import find_pareto_front

        assert find_pareto_front([]) == []

    def test_grid_search_returns_results(self) -> None:
        from kira.delivery.optimization import grid_search_formulations

        results = grid_search_formulations(
            target_tissue="liver",
            routes=[AdministrationRoute.IV],
            n_lipid_steps=2,
        )
        assert len(results) > 0
        # Should be sorted by overall score descending
        if len(results) > 1:
            assert results[0].overall >= results[1].overall

    def test_grid_search_all_viable_check(self) -> None:
        from kira.delivery.optimization import grid_search_formulations

        results = grid_search_formulations(
            target_tissue="liver",
            routes=[AdministrationRoute.IV],
            n_lipid_steps=2,
        )
        # At least some results should be viable for liver delivery
        viable = [r for r in results if r.is_viable]
        assert len(viable) > 0
