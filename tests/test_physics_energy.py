"""Tests for LJ energy computation and steric clash detection."""

from pathlib import Path

import jax.numpy as jnp
import numpy as np

from kira.physics.checks.clashes import check_clashes
from kira.physics.core.energy import (
    get_lj_params_arrays,
    run_lj_analysis,
)
from kira.physics.core.geometry import compute_distance_matrix
from kira.physics.core.parser import parse_pdb
from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

FIXTURES = Path(__file__).parent / "fixtures"


class TestLJEnergy:
    """Test Lennard-Jones energy computation."""

    def test_lj_params_shape(self):
        """LJ parameter arrays should match atom count."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        sigma, epsilon = get_lj_params_arrays(struct.elements)

        assert sigma.shape == (15,)
        assert epsilon.shape == (15,)
        assert np.all(sigma > 0)
        assert np.all(epsilon > 0)

    def test_lj_energy_finite(self):
        """Total LJ energy should be finite."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = run_lj_analysis(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        assert np.isfinite(result["total_energy"])

    def test_per_residue_sums_to_total(self):
        """Per-residue energies should approximately sum to total."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = run_lj_analysis(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        per_res_sum = np.sum(result["per_residue_energy"])
        # Per-residue uses per-atom (half of row sums), total uses upper triangle
        # They should be close
        np.testing.assert_allclose(
            per_res_sum, result["total_energy"], rtol=0.01
        )

    def test_per_atom_energy_shape(self):
        """Per-atom energy should have correct shape."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = run_lj_analysis(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        assert result["per_atom_energy"].shape == (15,)
        assert result["per_residue_energy"].shape == (3,)


class TestClashDetection:
    """Test steric clash detection."""

    def test_no_clashes_in_good_structure(self):
        """Well-built tri-alanine should have few or no clashes."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = check_clashes(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        # Good structure should have low clashscore
        assert result.clashscore < 20.0, f"Clashscore {result.clashscore} too high for good structure"
        assert result.subscore > 0.3

    def test_clashes_in_bad_structure(self):
        """Clashing fixture should have detectable clashes."""
        struct = parse_pdb(FIXTURES / "clashing.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = check_clashes(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        # Should detect at least one clash (CB atoms at 0.87 Å apart)
        assert result.n_clashes > 0, "Expected clashes in clashing fixture"

    def test_clash_result_fields(self):
        """ClashResult should have all expected fields."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        result = check_clashes(
            dmat, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues,
        )

        assert hasattr(result, "n_clashes")
        assert hasattr(result, "n_severe_clashes")
        assert hasattr(result, "clashscore")
        assert hasattr(result, "worst_overlap")
        assert hasattr(result, "subscore")
        assert hasattr(result, "per_residue_clashes")
        assert result.per_residue_clashes.shape == (3,)
        assert 0.0 <= result.subscore <= 1.0
