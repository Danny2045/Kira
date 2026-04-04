"""Tests for local vs global divergence analysis."""

from pathlib import Path

import numpy as np

from kira.causality.binding_site import extract_binding_site
from kira.causality.divergence import (
    AA_PROPERTIES,
    DivergenceProfile,
    PositionDivergence,
    _physicochemical_distance,
    compute_divergence_profile,
    rank_divergent_positions,
)
from kira.physics.core.parser import parse_pdb

FIXTURES = Path(__file__).parent / "fixtures"


def _make_site():
    """Build a binding site from tri-alanine fixture."""
    struct = parse_pdb(FIXTURES / "tri_ala.pdb")
    center = np.mean(struct.coords, axis=0, keepdims=True)
    return extract_binding_site(struct, center, cutoff=50.0)


class TestPhysicochemicalDistance:
    """Test amino acid property distance computation."""

    def test_identical_residues_zero_distance(self):
        """Same residue should have zero distance."""
        d = _physicochemical_distance("ALA", "ALA")
        assert d == 0.0

    def test_charge_change_has_large_distance(self):
        """ASP(-) -> LYS(+) should have large distance (charge flip)."""
        d = _physicochemical_distance("ASP", "LYS")
        assert d > 0.5  # Significant distance due to charge change

    def test_conservative_substitution_small_distance(self):
        """ILE -> LEU should have small distance (similar properties)."""
        d = _physicochemical_distance("ILE", "LEU")
        assert d < 0.2  # Very similar residues

    def test_distance_is_symmetric(self):
        """Distance should be symmetric."""
        d1 = _physicochemical_distance("PHE", "ALA")
        d2 = _physicochemical_distance("ALA", "PHE")
        np.testing.assert_allclose(d1, d2)

    def test_unknown_residue_uses_defaults(self):
        """Unknown residues should use default properties."""
        d = _physicochemical_distance("XYZ", "XYZ")
        assert d == 0.0  # Same defaults -> zero distance


class TestComputeDivergenceProfile:
    """Test divergence profile computation."""

    def test_identical_pockets_full_identity(self):
        """Same pocket compared to itself should have 100% local identity."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        assert isinstance(profile, DivergenceProfile)
        assert profile.local_sequence_identity == 1.0

    def test_identical_pockets_zero_changes(self):
        """Identical pockets should have no charge/hydro/volume changes."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        assert profile.n_charge_changes == 0
        assert profile.n_hydrophobicity_flips == 0
        assert profile.n_volume_changes == 0

    def test_per_position_divergence_count(self):
        """Should have one entry per aligned position."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        assert len(profile.per_position_divergence) == site.n_residues

    def test_per_position_fields(self):
        """Each PositionDivergence should have expected fields."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        for pos in profile.per_position_divergence:
            assert isinstance(pos, PositionDivergence)
            assert isinstance(pos.position, int)
            assert isinstance(pos.is_identical, bool)
            assert np.isfinite(pos.hydrophobicity_delta)
            assert np.isfinite(pos.charge_delta)
            assert np.isfinite(pos.volume_delta)
            assert np.isfinite(pos.physicochemical_distance)

    def test_with_global_sequences(self):
        """Should compute global identity when sequences provided."""
        site = _make_site()

        # Provide identical global sequences
        profile = compute_divergence_profile(
            site, site,
            global_seq1="AAAGGG",
            global_seq2="AAAGGG",
        )

        assert profile.global_sequence_identity == 1.0
        assert profile.identity_gap == 0.0  # Both 100%

    def test_divergent_global_sequences(self):
        """Should detect lower global identity with different sequences."""
        site = _make_site()

        profile = compute_divergence_profile(
            site, site,
            global_seq1="AAAAAA",
            global_seq2="AADDDD",
        )

        assert profile.global_sequence_identity < 1.0

    def test_is_locally_divergent_property(self):
        """is_locally_divergent should reflect identity gap."""
        site = _make_site()

        # Identical pockets, divergent global -> identity_gap > 0
        # means pocket is more conserved than global, so NOT locally divergent
        profile = compute_divergence_profile(
            site, site,
            global_seq1="AAAAAA",
            global_seq2="DDDDDD",
        )

        # Global identity is 0%, local is 100% -> gap = 0 - 1 = -1.0
        # Negative gap means pocket is MORE conserved -> not locally divergent
        assert not profile.is_locally_divergent

    def test_summary_returns_string(self):
        """summary() should return a descriptive string."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        s = profile.summary()
        assert isinstance(s, str)
        assert "identity" in s.lower()

    def test_local_physicochemical_distance_zero_for_identical(self):
        """Identical pockets should have zero local physicochemical distance."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        assert profile.local_physicochemical_distance == 0.0


class TestRankDivergentPositions:
    """Test position ranking by divergence."""

    def test_identical_pockets_no_ranked_positions(self):
        """Identical pockets should have no divergent positions to rank."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        ranked = rank_divergent_positions(profile)
        assert len(ranked) == 0  # All positions are identical

    def test_filter_by_min_distance(self):
        """Should filter out positions below min_distance."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        # With very high min_distance, nothing should pass
        ranked = rank_divergent_positions(profile, min_distance=100.0)
        assert len(ranked) == 0

    def test_sorted_by_distance_descending(self):
        """Ranked positions should be sorted by distance descending."""
        site = _make_site()
        profile = compute_divergence_profile(site, site)

        ranked = rank_divergent_positions(profile, min_distance=0.0)
        if len(ranked) >= 2:
            distances = [p.physicochemical_distance for p in ranked]
            assert distances == sorted(distances, reverse=True)


class TestAAProperties:
    """Test amino acid property data."""

    def test_all_standard_residues_present(self):
        """All 20 standard amino acids should have properties."""
        standard = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL",
        ]
        for aa in standard:
            assert aa in AA_PROPERTIES, f"Missing properties for {aa}"

    def test_property_keys_consistent(self):
        """All amino acids should have the same property keys."""
        expected_keys = {"hydrophobicity", "charge", "volume", "hbd", "hba"}
        for aa, props in AA_PROPERTIES.items():
            assert set(props.keys()) == expected_keys, f"Wrong keys for {aa}"

    def test_asp_glu_negative_charge(self):
        """ASP and GLU should have negative charge."""
        assert AA_PROPERTIES["ASP"]["charge"] < 0
        assert AA_PROPERTIES["GLU"]["charge"] < 0

    def test_lys_arg_positive_charge(self):
        """LYS and ARG should have positive charge."""
        assert AA_PROPERTIES["LYS"]["charge"] > 0
        assert AA_PROPERTIES["ARG"]["charge"] > 0
