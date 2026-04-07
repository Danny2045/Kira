"""Tests for the autonomous lab interface layer.

Tests experiment protocol design, discovery campaign management,
and the design-build-test-learn closed loop.
"""

from __future__ import annotations

import pytest

from kira.lab.protocol import (
    AssayType,
    CompoundSpec,
    ConcentrationSeries,
    ControlSpec,
    ExperimentProtocol,
    PlateFormat,
    ReadoutType,
    TargetSpec,
    design_cytotoxicity_counter_screen,
    design_dose_response,
    design_selectivity_assay,
)
from kira.lab.loop import (
    CompoundStatus,
    CycleIteration,
    CyclePhase,
    DiscoveryCampaign,
    ExperimentResult,
    SelectivityResult,
)


# ---------------------------------------------------------------------------
# ConcentrationSeries tests
# ---------------------------------------------------------------------------


class TestConcentrationSeries:
    """Tests for concentration series generation."""

    def test_default_series(self) -> None:
        cs = ConcentrationSeries()
        assert cs.n_points == 8
        assert cs.n_replicates == 3
        assert cs.total_wells == 24

    def test_concentrations_descending(self) -> None:
        cs = ConcentrationSeries(
            top_concentration_um=100.0,
            dilution_factor=3.0,
            n_points=8,
        )
        concs = cs.concentrations_um
        assert len(concs) == 8
        assert concs[0] == 100.0
        # Each subsequent concentration should be smaller
        for i in range(1, len(concs)):
            assert concs[i] < concs[i - 1]

    def test_nm_conversion(self) -> None:
        cs = ConcentrationSeries(top_concentration_um=100.0)
        um = cs.concentrations_um
        nm = cs.concentrations_nm
        for u, n in zip(um, nm):
            assert abs(n - u * 1000) < 0.01

    def test_dilution_factor_applied(self) -> None:
        cs = ConcentrationSeries(
            top_concentration_um=100.0,
            dilution_factor=10.0,
            n_points=3,
        )
        concs = cs.concentrations_um
        assert abs(concs[0] - 100.0) < 0.01
        assert abs(concs[1] - 10.0) < 0.01
        assert abs(concs[2] - 1.0) < 0.01

    def test_total_wells_calculation(self) -> None:
        cs = ConcentrationSeries(n_points=10, n_replicates=4)
        assert cs.total_wells == 40

    def test_to_dict(self) -> None:
        cs = ConcentrationSeries()
        d = cs.to_dict()
        assert "concentrations_um" in d
        assert "total_wells" in d


# ---------------------------------------------------------------------------
# Protocol tests
# ---------------------------------------------------------------------------


class TestExperimentProtocol:
    """Tests for experiment protocol generation."""

    def test_protocol_id_generated(self) -> None:
        proto = ExperimentProtocol(name="test")
        assert proto.protocol_id.startswith("KIRA-")
        assert len(proto.protocol_id) > 5

    def test_created_at_timestamp(self) -> None:
        proto = ExperimentProtocol(name="test")
        assert proto.created_at.endswith("Z")

    def test_quality_criteria_defaults(self) -> None:
        proto = ExperimentProtocol(name="test")
        assert proto.quality_criteria["z_prime_min"] == 0.5
        assert proto.quality_criteria["cv_max_pct"] == 20.0

    def test_total_wells_calculation(self) -> None:
        compounds = [CompoundSpec(compound_id=f"C{i}", name=f"c{i}") for i in range(5)]
        targets = [TargetSpec(target_name="T1")]
        proto = ExperimentProtocol(
            name="test",
            compounds=compounds,
            targets=targets,
            concentration_series=ConcentrationSeries(n_points=8, n_replicates=3),
            controls=ControlSpec(n_control_wells=8),
        )
        # 5 compounds * 1 target * 8 points * 3 reps = 120 data wells
        # + 8 control wells * 2 * 1 target = 16 control wells
        assert proto.total_data_wells == 120
        assert proto.total_control_wells == 16
        assert proto.total_wells == 136

    def test_n_plates_calculation(self) -> None:
        compounds = [CompoundSpec(compound_id=f"C{i}", name=f"c{i}") for i in range(50)]
        proto = ExperimentProtocol(
            name="many_compounds",
            compounds=compounds,
            targets=[TargetSpec(target_name="T1")],
            plate_format=PlateFormat.PLATE_384,
        )
        assert proto.n_plates >= 1

    def test_estimated_cost_positive(self) -> None:
        proto = ExperimentProtocol(
            name="cost",
            compounds=[CompoundSpec(compound_id="C1", name="c1")],
            targets=[TargetSpec(target_name="T1")],
        )
        assert proto.estimated_cost_usd > 0

    def test_summary_string(self) -> None:
        proto = ExperimentProtocol(
            name="summary_test",
            compounds=[CompoundSpec(compound_id="C1", name="c1")],
            targets=[TargetSpec(target_name="T1")],
        )
        summary = proto.summary()
        assert "KIRA-" in summary
        assert "summary_test" in summary

    def test_to_dict_serialization(self) -> None:
        proto = ExperimentProtocol(
            name="dict_test",
            compounds=[CompoundSpec(compound_id="C1", name="c1", smiles="CCO")],
            targets=[TargetSpec(target_name="T1", organism="test")],
        )
        d = proto.to_dict()
        assert d["name"] == "dict_test"
        assert len(d["compounds"]) == 1
        assert d["compounds"][0]["smiles"] == "CCO"


# ---------------------------------------------------------------------------
# Protocol generator tests
# ---------------------------------------------------------------------------


class TestProtocolGenerators:
    """Tests for protocol generator functions."""

    def test_selectivity_assay(self) -> None:
        compounds = [
            CompoundSpec(compound_id="CHEMBL155771", name="CHEMBL155771"),
            CompoundSpec(compound_id="CHEMBL234567", name="Compound 2"),
        ]
        proto = design_selectivity_assay(
            compounds=compounds,
            parasite_target=TargetSpec(
                target_name="SmDHODH",
                organism="Schistosoma mansoni",
            ),
            human_target=TargetSpec(
                target_name="HsDHODH",
                organism="Homo sapiens",
            ),
        )
        assert proto.assay_type == AssayType.SELECTIVITY_PANEL
        assert proto.n_compounds == 2
        assert proto.n_targets == 2  # Both parasite and human
        assert "Selectivity" in proto.name

    def test_dose_response(self) -> None:
        compounds = [CompoundSpec(compound_id="C1", name="c1")]
        proto = design_dose_response(
            compounds=compounds,
            target=TargetSpec(target_name="SmDHODH"),
        )
        assert proto.assay_type == AssayType.DOSE_RESPONSE
        assert proto.n_compounds == 1
        assert proto.n_targets == 1

    def test_cytotoxicity_counter_screen(self) -> None:
        compounds = [CompoundSpec(compound_id="C1", name="c1")]
        proto = design_cytotoxicity_counter_screen(
            compounds=compounds,
            cell_type="HEK293",
            incubation_hours=48.0,
        )
        assert proto.assay_type == AssayType.CYTOTOXICITY
        assert proto.incubation_hours == 48.0
        assert "staurosporine" in proto.controls.positive_control

    def test_selectivity_assay_custom_concentrations(self) -> None:
        compounds = [CompoundSpec(compound_id="C1", name="c1")]
        proto = design_selectivity_assay(
            compounds=compounds,
            parasite_target=TargetSpec(target_name="SmTGR"),
            human_target=TargetSpec(target_name="HsTXNRD1"),
            top_concentration_um=50.0,
            n_dose_points=12,
            n_replicates=4,
        )
        assert proto.concentration_series.top_concentration_um == 50.0
        assert proto.concentration_series.n_points == 12
        assert proto.concentration_series.n_replicates == 4


# ---------------------------------------------------------------------------
# Discovery campaign tests
# ---------------------------------------------------------------------------


class TestDiscoveryCampaign:
    """Tests for discovery campaign management."""

    def test_create_campaign(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="KIRA-NTD-001",
            name="SmDHODH selectivity screen",
            target_pair=("SmDHODH", "HsDHODH"),
            disease="schistosomiasis",
        )
        assert campaign.current_iteration == 0
        assert campaign.total_compounds_tested == 0

    def test_start_iteration(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
        )
        iteration = campaign.start_iteration()
        assert iteration.iteration == 1
        assert iteration.phase == CyclePhase.DESIGN
        assert campaign.current_iteration == 1

    def test_record_results_hits(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
            activity_threshold_nm=1000.0,
        )
        campaign.start_iteration()

        results = [
            ExperimentResult(
                compound_id="C1",
                target_name="SmDHODH",
                ic50_nm=50.0,
                is_active=True,
            ),
            ExperimentResult(
                compound_id="C2",
                target_name="SmDHODH",
                ic50_nm=5000.0,
                is_active=True,
            ),
            ExperimentResult(
                compound_id="C3",
                target_name="SmDHODH",
                is_active=False,
            ),
        ]

        iteration = campaign.record_results(results)
        assert iteration.n_hits == 1  # Only C1 below threshold
        assert campaign.active_compounds["C1"] == CompoundStatus.HIT
        assert campaign.active_compounds["C3"] == CompoundStatus.ELIMINATED

    def test_selectivity_promotes_to_lead(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
            selectivity_threshold=10.0,
        )
        campaign.start_iteration()

        results = [
            ExperimentResult(
                compound_id="C1",
                target_name="SmDHODH",
                ic50_nm=50.0,
                is_active=True,
            ),
        ]
        sel_results = [
            SelectivityResult(
                compound_id="C1",
                parasite_target="SmDHODH",
                human_target="HsDHODH",
                parasite_ic50_nm=50.0,
                human_ic50_nm=1500.0,  # 30x selective
            ),
        ]

        campaign.record_results(results, sel_results)
        assert campaign.active_compounds["C1"] == CompoundStatus.LEAD
        assert campaign.total_hits == 1  # LEAD counts as hit

    def test_recommendations_generated(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
        )
        campaign.start_iteration()

        results = [
            ExperimentResult(
                compound_id="C1",
                target_name="SmDHODH",
                ic50_nm=50.0,
                is_active=True,
            ),
        ]
        iteration = campaign.record_results(results)
        assert len(iteration.next_actions) > 0

    def test_multiple_iterations(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
        )

        # Iteration 1
        campaign.start_iteration()
        campaign.record_results([
            ExperimentResult(
                compound_id="C1", target_name="SmDHODH",
                ic50_nm=50.0, is_active=True,
            ),
        ])

        # Iteration 2
        campaign.start_iteration()
        campaign.record_results([
            ExperimentResult(
                compound_id="C2", target_name="SmDHODH",
                ic50_nm=25.0, is_active=True,
            ),
        ])

        assert campaign.current_iteration == 2
        assert campaign.total_compounds_tested == 2

    def test_selectivity_result_auto_calculation(self) -> None:
        sel = SelectivityResult(
            compound_id="C1",
            parasite_target="SmDHODH",
            human_target="HsDHODH",
            parasite_ic50_nm=50.0,
            human_ic50_nm=1500.0,
            selectivity_threshold=10.0,
        )
        assert sel.selectivity_ratio is not None
        assert abs(sel.selectivity_ratio - 30.0) < 0.01
        assert sel.is_selective

    def test_selectivity_result_not_selective(self) -> None:
        sel = SelectivityResult(
            compound_id="C2",
            parasite_target="SmDHODH",
            human_target="HsDHODH",
            parasite_ic50_nm=50.0,
            human_ic50_nm=100.0,  # Only 2x
            selectivity_threshold=10.0,
        )
        assert sel.selectivity_ratio is not None
        assert abs(sel.selectivity_ratio - 2.0) < 0.01
        assert not sel.is_selective

    def test_campaign_summary(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="KIRA-001",
            name="SmDHODH screen",
            target_pair=("SmDHODH", "HsDHODH"),
        )
        summary = campaign.summary()
        assert "KIRA-001" in summary
        assert "SmDHODH" in summary

    def test_campaign_to_dict(self) -> None:
        campaign = DiscoveryCampaign(
            campaign_id="test",
            name="test",
            target_pair=("SmDHODH", "HsDHODH"),
        )
        d = campaign.to_dict()
        assert d["campaign_id"] == "test"
        assert d["target_pair"] == ["SmDHODH", "HsDHODH"]

    def test_hit_rate_zero_when_no_tests(self) -> None:
        iteration = CycleIteration(iteration=1)
        assert iteration.hit_rate == 0.0
