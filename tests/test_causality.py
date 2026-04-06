"""Tests for causality binding site extraction and comparison."""

from pathlib import Path

import numpy as np

from kira.causality.binding_site import (
    BindingSite,
    PocketComparison,
    compare_binding_sites,
    extract_binding_site,
)
from kira.physics.core.parser import parse_pdb

FIXTURES = Path(__file__).parent / "fixtures"


class TestBindingSiteExtraction:
    """Test binding pocket extraction."""

    def test_extract_near_origin(self):
        """Ligand near first residue should capture nearby residues."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        # Place a fake ligand near residue 1
        ligand_coords = np.array([[2.0, 1.4, 0.0]], dtype=np.float32)

        site = extract_binding_site(struct, ligand_coords, cutoff=5.0)

        assert isinstance(site, BindingSite)
        assert site.n_residues > 0
        assert site.n_atoms > 0
        assert site.coords.shape[1] == 3
        assert site.centroid.shape == (3,)

    def test_large_cutoff_gets_all(self):
        """Very large cutoff should capture all residues."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        # Ligand at center of structure
        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        assert site.n_residues == 3

    def test_zero_cutoff_gets_none(self):
        """Cutoff of 0 should get no residues."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        ligand_coords = np.array([[100.0, 100.0, 100.0]], dtype=np.float32)
        site = extract_binding_site(struct, ligand_coords, cutoff=0.1)

        assert site.n_residues == 0

    def test_sequence_property(self):
        """Binding site should have correct sequence."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")

        center = np.mean(struct.coords, axis=0, keepdims=True)
        site = extract_binding_site(struct, center, cutoff=50.0)

        assert site.sequence == "AAA"


class TestPocketComparison:
    """Test binding pocket comparison."""

    def test_identical_pockets(self):
        """Two identical pockets should have 100% identity."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        center = np.mean(struct.coords, axis=0, keepdims=True)

        site1 = extract_binding_site(struct, center, cutoff=50.0)
        site2 = extract_binding_site(struct, center, cutoff=50.0)

        comp = compare_binding_sites(site1, site2)

        assert isinstance(comp, PocketComparison)
        assert comp.sequence_identity == 1.0
        assert comp.n_divergent == 0
        assert comp.divergence_fraction == 0.0

    def test_comparison_fields(self):
        """PocketComparison should have all expected fields."""
        struct = parse_pdb(FIXTURES / "tri_ala.pdb")
        center = np.mean(struct.coords, axis=0, keepdims=True)

        site = extract_binding_site(struct, center, cutoff=50.0)
        comp = compare_binding_sites(site, site)

        assert hasattr(comp, "n_aligned_positions")
        assert hasattr(comp, "sequence_identity")
        assert hasattr(comp, "divergent_positions")
        assert hasattr(comp, "pocket1_residues")
        assert hasattr(comp, "pocket2_residues")
        assert comp.n_aligned_positions == 3
