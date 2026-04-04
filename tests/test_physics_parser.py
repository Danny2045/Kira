"""Tests for the PDB parser."""

from pathlib import Path

import numpy as np
import pytest

from kira.physics.core.parser import parse_pdb, parse_pdb_string

FIXTURES = Path(__file__).parent / "fixtures"


class TestParsePDB:
    """Test PDB file parsing."""

    def test_parse_tri_ala(self):
        """Parse tri-alanine fixture and verify basic properties."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        assert struct.name == "tri_ala"
        assert struct.n_atoms == 15
        assert struct.n_residues == 3
        assert struct.n_chains == 1

    def test_coordinates_shape(self):
        """Verify coordinate array has correct shape."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        assert struct.coords.shape == (15, 3)
        assert struct.coords.dtype == np.float32

    def test_first_atom_coordinates(self):
        """Verify first atom (N of ALA 1) has correct coordinates."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        np.testing.assert_allclose(
            struct.coords[0], [1.458, 0.0, 0.0], atol=1e-3
        )

    def test_element_parsing(self):
        """Verify elements are correctly parsed."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        assert struct.elements[0] == "N"   # First atom is nitrogen
        assert struct.elements[1] == "C"   # CA is carbon
        assert struct.elements[3] == "O"   # O is oxygen

    def test_protein_mask(self):
        """All atoms in tri-ala should be protein."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        assert np.all(struct.is_protein_mask)

    def test_backbone_mask(self):
        """N, CA, C, O should be backbone; CB should not."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        # Atom 0=N, 1=CA, 2=C, 3=O are backbone; 4=CB is not
        assert struct.is_backbone_mask[0]  # N
        assert struct.is_backbone_mask[1]  # CA
        assert struct.is_backbone_mask[2]  # C
        assert struct.is_backbone_mask[3]  # O
        assert not struct.is_backbone_mask[4]  # CB

    def test_residue_indices(self):
        """Verify residue index mapping."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        # First 5 atoms belong to residue 0
        assert struct.res_indices[0] == 0
        assert struct.res_indices[4] == 0
        # Next 5 to residue 1
        assert struct.res_indices[5] == 1
        assert struct.res_indices[9] == 1
        # Last 5 to residue 2
        assert struct.res_indices[10] == 2
        assert struct.res_indices[14] == 2

    def test_chain_sequence(self):
        """Verify chain sequence extraction."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        chain = struct.chains["A"]
        assert chain.sequence == "AAA"
        assert chain.is_protein

    def test_residue_atom_lookup(self):
        """Verify residue-level atom access."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        rid = ("A", 1, "")
        residue = struct.residues[rid]
        assert residue.res_name == "ALA"
        assert "CA" in residue.atoms
        assert "CB" in residue.atoms

        ca_coord = residue.get_coord("CA")
        assert ca_coord is not None
        np.testing.assert_allclose(ca_coord, [2.009, 1.420, 0.0], atol=1e-3)

    def test_heavy_atom_coords(self):
        """Heavy atom coords should exclude hydrogens (none in fixture)."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        heavy = struct.heavy_atom_coords()
        # No hydrogens in fixture, so same count
        assert heavy.shape == (15, 3)

    def test_file_not_found(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_pdb("/nonexistent/path.pdb")

    def test_parse_string(self):
        """Parse PDB from string."""
        pdb_text = (FIXTURES / "tri_ala.pdb").read_text()
        struct = parse_pdb_string(pdb_text, name="from_string")

        assert struct.name == "from_string"
        assert struct.n_atoms == 15

    def test_clashing_fixture_parses(self):
        """Verify clashing fixture parses correctly."""
        struct = parse_pdb(FIXTURES / "clashing.pdb")

        assert struct.n_atoms == 10
        assert struct.n_residues == 2
