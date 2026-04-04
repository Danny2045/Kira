"""Tests for JAX geometry kernels."""

from pathlib import Path

import jax.numpy as jnp
import numpy as np

from kira.physics.core.geometry import (
    compute_bond_angles,
    compute_dihedral_angles,
    compute_distance_matrix,
    extract_backbone_dihedrals,
)
from kira.physics.core.parser import parse_pdb

FIXTURES = Path(__file__).parent / "fixtures"


class TestDistanceMatrix:
    """Test pairwise distance computation."""

    def test_simple_distances(self):
        """Known distances for simple coordinates."""
        coords = jnp.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        dmat = compute_distance_matrix(coords)

        assert dmat.shape == (3, 3)
        # d(0,1) = 1.0
        np.testing.assert_allclose(float(dmat[0, 1]), 1.0, atol=1e-3)
        # d(0,2) = 1.0
        np.testing.assert_allclose(float(dmat[0, 2]), 1.0, atol=1e-3)
        # d(1,2) = sqrt(2)
        np.testing.assert_allclose(float(dmat[1, 2]), np.sqrt(2), atol=1e-3)

    def test_symmetric(self):
        """Distance matrix should be symmetric."""
        coords = jnp.array(np.random.randn(10, 3).astype(np.float32))
        dmat = compute_distance_matrix(coords)

        np.testing.assert_allclose(
            np.asarray(dmat), np.asarray(dmat.T), atol=1e-5
        )

    def test_diagonal_near_zero(self):
        """Diagonal should be near zero (epsilon for stability)."""
        coords = jnp.array(np.random.randn(5, 3).astype(np.float32))
        dmat = compute_distance_matrix(coords)

        for i in range(5):
            assert float(dmat[i, i]) < 0.001

    def test_from_structure(self):
        """Distance matrix from parsed structure should have correct shape."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        coords = jnp.array(struct.coords)
        dmat = compute_distance_matrix(coords)

        assert dmat.shape == (15, 15)


class TestDihedralAngles:
    """Test dihedral angle computation."""

    def test_trans_dihedral(self):
        """Trans configuration should give ~180 degrees."""
        # Four atoms in trans: 0-1-2-3 along x-axis with z-offsets
        p0 = jnp.array([[0.0, 0.0, 1.0]])
        p1 = jnp.array([[1.0, 0.0, 0.0]])
        p2 = jnp.array([[2.0, 0.0, 0.0]])
        p3 = jnp.array([[3.0, 0.0, -1.0]])

        angle = compute_dihedral_angles(p0, p1, p2, p3)
        angle_deg = float(jnp.degrees(angle[0]))

        # Trans should be close to ±180
        assert abs(abs(angle_deg) - 180.0) < 5.0, f"Expected ~180, got {angle_deg}"

    def test_cis_dihedral(self):
        """Cis configuration should give ~0 degrees."""
        p0 = jnp.array([[0.0, 0.0, 1.0]])
        p1 = jnp.array([[1.0, 0.0, 0.0]])
        p2 = jnp.array([[2.0, 0.0, 0.0]])
        p3 = jnp.array([[3.0, 0.0, 1.0]])

        angle = compute_dihedral_angles(p0, p1, p2, p3)
        angle_deg = float(jnp.degrees(angle[0]))

        assert abs(angle_deg) < 5.0, f"Expected ~0, got {angle_deg}"

    def test_batch_dihedrals(self):
        """Vectorized computation of multiple dihedrals."""
        n = 5
        p0 = jnp.array(np.random.randn(n, 3).astype(np.float32))
        p1 = jnp.array(np.random.randn(n, 3).astype(np.float32))
        p2 = jnp.array(np.random.randn(n, 3).astype(np.float32))
        p3 = jnp.array(np.random.randn(n, 3).astype(np.float32))

        angles = compute_dihedral_angles(p0, p1, p2, p3)
        assert angles.shape == (5,)
        # All angles should be in [-pi, pi]
        assert jnp.all(angles >= -jnp.pi - 0.01)
        assert jnp.all(angles <= jnp.pi + 0.01)


class TestBondAngles:
    """Test bond angle computation."""

    def test_right_angle(self):
        """90-degree angle."""
        p0 = jnp.array([[1.0, 0.0, 0.0]])
        p1 = jnp.array([[0.0, 0.0, 0.0]])  # vertex
        p2 = jnp.array([[0.0, 1.0, 0.0]])

        angle = compute_bond_angles(p0, p1, p2)
        angle_deg = float(jnp.degrees(angle[0]))

        np.testing.assert_allclose(angle_deg, 90.0, atol=0.1)

    def test_straight_angle(self):
        """180-degree (linear) angle."""
        p0 = jnp.array([[-1.0, 0.0, 0.0]])
        p1 = jnp.array([[0.0, 0.0, 0.0]])
        p2 = jnp.array([[1.0, 0.0, 0.0]])

        angle = compute_bond_angles(p0, p1, p2)
        angle_deg = float(jnp.degrees(angle[0]))

        np.testing.assert_allclose(angle_deg, 180.0, atol=0.1)


class TestBackboneDihedrals:
    """Test backbone dihedral extraction from structures."""

    def test_tri_ala_dihedrals(self):
        """Tri-alanine should have phi, psi, and omega angles."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        dihedrals = extract_backbone_dihedrals(
            struct.coords, struct.atom_names, struct.res_indices,
            struct.is_protein_mask, struct.chain_ids_array,
        )

        # 3 residues: phi defined for residues 2,3; psi for 1,2; omega for 2,3
        assert "phi" in dihedrals
        assert "psi" in dihedrals
        assert "omega" in dihedrals

        assert len(dihedrals["phi"]["res_indices"]) == 2  # residues 2 and 3
        assert len(dihedrals["psi"]["res_indices"]) == 2  # residues 1 and 2
        assert len(dihedrals["omega"]["res_indices"]) == 2

        # Omega angles should be roughly trans (~180 degrees)
        omega_angles = np.asarray(dihedrals["omega"]["angles"])
        for angle in omega_angles:
            assert abs(abs(float(angle)) - 180.0) < 60.0, \
                f"Omega angle {float(angle)} too far from trans"
