"""Tests for per-residue energy decomposition at binding sites."""

from pathlib import Path

import numpy as np

from kira.causality.binding_site import extract_binding_site
from kira.causality.energy_decomp import (
    BindingSiteDecomposition,
    ResidueEnergyContribution,
    compare_residue_energies,
    decompose_binding_site_energy,
)
from kira.physics.core.parser import parse_pdb
from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

FIXTURES = Path(__file__).parent / "fixtures"


class TestDecomposeBindingSiteEnergy:
    """Test per-residue energy decomposition."""

    def _setup_structure(self):
        """Parse tri-alanine and build topology."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)
        return struct, mask

    def test_decomposition_returns_correct_type(self):
        """Should return BindingSiteDecomposition."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        assert isinstance(decomp, BindingSiteDecomposition)

    def test_decomposition_has_residues(self):
        """Decomposition should contain contributions for pocket residues."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        assert decomp.n_residues > 0
        assert decomp.n_favorable + decomp.n_unfavorable == decomp.n_residues

    def test_residue_contribution_fields(self):
        """Each ResidueEnergyContribution should have expected fields."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        for contrib in decomp.residue_contributions:
            assert isinstance(contrib, ResidueEnergyContribution)
            assert isinstance(contrib.residue_index, int)
            assert isinstance(contrib.res_name, str)
            assert np.isfinite(contrib.total_energy)
            assert contrib.n_atoms > 0
            assert isinstance(contrib.favorable, bool)
            assert len(contrib.per_atom_energies) == contrib.n_atoms

    def test_total_energy_is_sum_of_residues(self):
        """Total energy should equal sum of residue contributions."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        residue_sum = sum(r.total_energy for r in decomp.residue_contributions)
        np.testing.assert_allclose(decomp.total_energy, residue_sum, rtol=1e-5)

    def test_sorted_by_absolute_energy(self):
        """Contributions should be sorted by |energy| descending."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        energies = [abs(r.total_energy) for r in decomp.residue_contributions]
        assert energies == sorted(energies, reverse=True)

    def test_empty_site_returns_empty_decomposition(self):
        """Empty binding site should return zero-energy decomposition."""
        struct, mask = self._setup_structure()
        # Place ligand far away so no residues are captured
        ligand_coords = np.array([[1000.0, 1000.0, 1000.0]], dtype=np.float32)
        site = extract_binding_site(struct, ligand_coords, cutoff=0.1)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        assert decomp.total_energy == 0.0
        assert decomp.n_residues == 0
        assert decomp.n_favorable == 0
        assert decomp.n_unfavorable == 0
        assert len(decomp.hotspot_indices) == 0

    def test_mean_residue_energy(self):
        """mean_residue_energy property should work."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        if decomp.n_residues > 0:
            expected = decomp.total_energy / decomp.n_residues
            np.testing.assert_allclose(decomp.mean_residue_energy, expected)

    def test_top_favorable_and_unfavorable(self):
        """top_favorable and top_unfavorable should return correct subsets."""
        struct, mask = self._setup_structure()
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        decomp = decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

        for r in decomp.top_favorable(5):
            assert r.favorable is True
        for r in decomp.top_unfavorable(5):
            assert r.favorable is False


class TestCompareResidueEnergies:
    """Test cross-structure energy comparison."""

    def _make_decomp(self):
        """Build a decomposition from tri-alanine."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)
        return decompose_binding_site_energy(
            struct.coords, struct.elements, struct.res_indices,
            mask, site,
        )

    def test_self_comparison_zero_delta(self):
        """Comparing a decomposition to itself should yield zero deltas."""
        decomp = self._make_decomp()
        comparisons = compare_residue_energies(decomp, decomp)

        for c in comparisons:
            np.testing.assert_allclose(c["delta_energy"], 0.0, atol=1e-5)
            assert c["same_residue"] is True

    def test_comparison_sorted_by_delta(self):
        """Comparisons should be sorted by |delta_energy| descending."""
        decomp = self._make_decomp()
        comparisons = compare_residue_energies(decomp, decomp)

        deltas = [abs(c["delta_energy"]) for c in comparisons]
        assert deltas == sorted(deltas, reverse=True)

    def test_comparison_has_expected_keys(self):
        """Each comparison dict should have all expected keys."""
        decomp = self._make_decomp()
        comparisons = compare_residue_energies(decomp, decomp)

        expected_keys = {
            "pos", "res1_name", "res2_name", "res1_index", "res2_index",
            "energy1", "energy2", "delta_energy", "same_residue",
        }
        for c in comparisons:
            assert set(c.keys()) == expected_keys
