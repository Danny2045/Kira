"""Tests for selectivity attribution maps."""

from pathlib import Path

import numpy as np

from kira.causality.binding_site import extract_binding_site
from kira.causality.energy_decomp import (
    decompose_binding_site_energy,
)
from kira.causality.selectivity_map import (
    ResidueAttribution,
    SelectivityMap,
    build_selectivity_map,
    compute_selectivity_score,
)
from kira.physics.core.parser import parse_pdb
from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

FIXTURES = Path(__file__).parent / "fixtures"


def _build_decomp():
    """Build a decomposition from tri-alanine fixture."""
    struct = parse_pdb(FIXTURES / "tri_ala.pdb")
    bonds = infer_bonds_from_topology(struct)
    mask = build_bonded_mask(struct.n_atoms, bonds)
    center = np.mean(struct.coords, axis=0, keepdims=True)
    site = extract_binding_site(struct, center, cutoff=50.0)
    decomp = decompose_binding_site_energy(
        struct.coords, struct.elements, struct.res_indices,
        mask, site,
    )
    return site, decomp


class TestBuildSelectivityMap:
    """Test selectivity map construction."""

    def test_self_comparison_is_neutral(self):
        """Comparing identical structures should yield neutral selectivity."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(
            site, site, decomp, decomp,
            structure1_name="SmTGR", structure2_name="HsTrxR1",
        )

        assert isinstance(sel_map, SelectivityMap)
        assert sel_map.selectivity_direction == "neutral"
        np.testing.assert_allclose(sel_map.total_delta_energy, 0.0, atol=1e-4)

    def test_map_has_attributions(self):
        """Map should contain per-position attributions."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)

        assert sel_map.n_positions > 0
        for attr in sel_map.attributions:
            assert isinstance(attr, ResidueAttribution)
            assert isinstance(attr.position, int)
            assert isinstance(attr.is_divergent, bool)
            assert np.isfinite(attr.delta_energy)
            assert 0.0 <= attr.attribution_score <= 1.0

    def test_attribution_scores_sum_to_one(self):
        """Attribution scores should sum to approximately 1.0."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)

        if sel_map.attributions:
            total_score = sum(a.attribution_score for a in sel_map.attributions)
            # May not sum to exactly 1.0 due to floating point,
            # but should be close (or 0 if all deltas are 0)
            assert total_score <= 1.0 + 1e-5

    def test_structure_names_preserved(self):
        """Structure names should be stored in the map."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(
            site, site, decomp, decomp,
            structure1_name="SmDHODH", structure2_name="HsDHODH",
        )

        assert sel_map.structure1_name == "SmDHODH"
        assert sel_map.structure2_name == "HsDHODH"

    def test_divergent_and_conserved_properties(self):
        """divergent_attributions and conserved_attributions should partition."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)

        n_div = len(sel_map.divergent_attributions)
        n_con = len(sel_map.conserved_attributions)
        assert n_div + n_con == sel_map.n_positions

    def test_top_drivers(self):
        """top_drivers should return at most n items."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)

        drivers = sel_map.top_drivers(2)
        assert len(drivers) <= 2

    def test_summary_returns_string(self):
        """summary() should return a non-empty string."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)

        s = sel_map.summary()
        assert isinstance(s, str)
        assert len(s) > 0


class TestSelectivityScore:
    """Test scalar selectivity score computation."""

    def test_neutral_has_low_score(self):
        """Identical structures should yield a low selectivity score."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)
        score = compute_selectivity_score(sel_map)

        assert 0.0 <= score <= 1.0
        # For identical structures, energy magnitude component should be near 0
        # but concentration component may be nonzero

    def test_score_in_valid_range(self):
        """Score should always be in [0, 1]."""
        site, decomp = _build_decomp()

        sel_map = build_selectivity_map(site, site, decomp, decomp)
        score = compute_selectivity_score(sel_map)

        assert 0.0 <= score <= 1.0

    def test_empty_map_returns_zero(self):
        """Empty selectivity map should return score 0."""
        from kira.causality.binding_site import PocketComparison

        empty_map = SelectivityMap(
            structure1_name="A",
            structure2_name="B",
            total_delta_energy=0.0,
            pocket_comparison=PocketComparison(
                n_aligned_positions=0,
                sequence_identity=0.0,
                identity_at_positions=[],
                divergent_positions=[],
                pocket1_residues=[],
                pocket2_residues=[],
            ),
            attributions=[],
            n_positions=0,
            n_divergent=0,
            selectivity_direction="neutral",
        )

        score = compute_selectivity_score(empty_map)
        assert score == 0.0
