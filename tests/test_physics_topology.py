"""Tests for topology builder — bond inference and mask generation."""

from pathlib import Path

import numpy as np

from kira.physics.core.parser import parse_pdb
from kira.physics.core.topology import (
    build_bonded_mask,
    get_vdw_radii_array,
    infer_bonds_from_topology,
)

FIXTURES = Path(__file__).parent / "fixtures"


class TestBondInference:
    """Test covalent bond inference."""

    def test_tri_ala_bond_count(self):
        """Tri-alanine should have known number of bonds."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)

        # 3 residues × (N-CA, CA-C, C-O, CA-CB) = 12 intra-residue bonds
        # + 2 peptide bonds (C1-N2, C2-N3) = 14 total
        # But GLY has no CB, ALA has 4 per residue → 3×4=12 + 2 = 14
        assert len(bonds) >= 12  # At minimum backbone + CB bonds
        assert len(bonds) <= 16  # Upper bound with some tolerance

    def test_bonds_are_sorted_pairs(self):
        """All bonds should be (i, j) with i < j."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)

        for i, j in bonds:
            assert i < j, f"Bond ({i}, {j}) not sorted"

    def test_peptide_bonds_present(self):
        """Peptide bonds between residues should be detected."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)

        # C of residue 1 (atom index 2) to N of residue 2 (atom index 5)
        # C of residue 2 (atom index 7) to N of residue 3 (atom index 10)
        # But actual indices depend on ordering — check by atom names
        bond_set = set(bonds)

        # Find C and N atoms at residue boundaries
        c_atoms = [i for i, a in enumerate(struct.atoms) if a.name == "C"]
        n_atoms = [i for i, a in enumerate(struct.atoms) if a.name == "N"]

        # Should have peptide bond between C[0]-N[1] and C[1]-N[2]
        peptide_found = 0
        for c_idx in c_atoms:
            for n_idx in n_atoms:
                pair = tuple(sorted((c_idx, n_idx)))
                if pair in bond_set:
                    # Verify these are from adjacent residues
                    c_res = struct.res_indices[c_idx]
                    n_res = struct.res_indices[n_idx]
                    if n_res == c_res + 1:
                        peptide_found += 1

        assert peptide_found == 2, f"Expected 2 peptide bonds, found {peptide_found}"


class TestBondedMask:
    """Test non-bonded pair mask generation."""

    def test_mask_shape(self):
        """Mask should be (N, N) boolean."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        assert mask.shape == (15, 15)
        assert mask.dtype == bool

    def test_diagonal_is_false(self):
        """Self-interactions should be excluded."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        assert not np.any(np.diag(mask))

    def test_bonded_pairs_excluded(self):
        """Directly bonded atoms should be masked out."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        for i, j in bonds:
            assert not mask[i, j], f"Bonded pair ({i}, {j}) should be masked"
            assert not mask[j, i], f"Bonded pair ({j}, {i}) should be masked"

    def test_mask_is_symmetric(self):
        """Non-bonded mask should be symmetric."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        np.testing.assert_array_equal(mask, mask.T)


class TestVdWRadii:
    """Test van der Waals radii lookup."""

    def test_radii_shape(self):
        """Should return one radius per atom."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        radii = get_vdw_radii_array(struct)

        assert radii.shape == (15,)
        assert radii.dtype == np.float32

    def test_known_radii(self):
        """Carbon and nitrogen should have standard radii."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        radii = get_vdw_radii_array(struct)

        # First atom is N (radius ~1.55)
        assert abs(radii[0] - 1.55) < 0.01
        # Second atom is C (radius ~1.70)
        assert abs(radii[1] - 1.70) < 0.01
