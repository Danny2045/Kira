"""Binding site extraction and comparison.

Given a protein structure and ligand coordinates or a representative
centroid, extracts nearby pocket residues. Given two homologous proteins,
aligns their pockets and computes local divergence metrics associated with
selectivity.

This is the module that catches what ESM-2 global embeddings miss:
99.97% cosine similarity can coexist with 30x selectivity because
binding-site-level divergence is invisible to global representations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kira.physics.core.parser import Residue, Structure


@dataclass
class BindingSite:
    """Extracted binding site from a structure.

    Attributes
    ----------
    residues : list[Residue]
        Residues within the binding pocket.
    residue_indices : list[int]
        Global residue indices in the parent structure.
    atom_indices : np.ndarray
        Global atom indices of all atoms in pocket residues.
    coords : np.ndarray
        (M, 3) coordinates of pocket atoms.
    centroid : np.ndarray
        (3,) center of mass of the pocket.
    res_names : list[str]
        Residue names in pocket order.
    """

    residues: list[Residue]
    residue_indices: list[int]
    atom_indices: np.ndarray
    coords: np.ndarray
    centroid: np.ndarray
    res_names: list[str]

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_atoms(self) -> int:
        return len(self.atom_indices)

    @property
    def sequence(self) -> str:
        """One-letter sequence of pocket residues."""
        three_to_one = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        return "".join(three_to_one.get(r, "X") for r in self.res_names)


def extract_binding_site(
    structure: Structure,
    ligand_coords: np.ndarray,
    cutoff: float = 5.0,
    protein_only: bool = True,
) -> BindingSite:
    """Extract binding pocket residues within cutoff distance of ligand coordinates.

    Parameters
    ----------
    structure : Structure
        Full protein structure.
    ligand_coords : np.ndarray
        (L, 3) coordinates of ligand atoms or a representative centroid.
    cutoff : float
        Distance cutoff in Angstroms. Residues with any atom within
        this distance of any ligand atom are included.
    protein_only : bool
        If True, only include protein residues in the pocket.

    Returns
    -------
    BindingSite
        Extracted binding pocket.
    """
    # Find all protein atoms within cutoff of any ligand atom
    pocket_residue_indices = set()

    for i, atom in enumerate(structure.atoms):
        if protein_only and not atom.is_protein:
            continue
        if atom.is_hydrogen:
            continue

        atom_coord = structure.coords[i]

        # Minimum distance to any ligand atom
        dists = np.linalg.norm(ligand_coords - atom_coord, axis=1)
        min_dist = np.min(dists)

        if min_dist <= cutoff:
            pocket_residue_indices.add(int(structure.res_indices[i]))

    # Build the BindingSite
    pocket_residue_indices = sorted(pocket_residue_indices)

    # Get residue objects and atom indices
    rid_list = list(structure.residues.keys())
    residues = []
    atom_indices = []
    res_names = []

    for ridx in pocket_residue_indices:
        if ridx < len(rid_list):
            rid = rid_list[ridx]
            residue = structure.residues[rid]
            residues.append(residue)
            res_names.append(residue.res_name)

            # Get atom indices for this residue
            for j in range(structure.n_atoms):
                if structure.res_indices[j] == ridx:
                    atom_indices.append(j)

    atom_indices = np.array(atom_indices, dtype=np.int32) if atom_indices else np.array([], dtype=np.int32)
    coords = structure.coords[atom_indices] if len(atom_indices) > 0 else np.empty((0, 3))
    centroid = np.mean(coords, axis=0) if len(coords) > 0 else np.zeros(3)

    return BindingSite(
        residues=residues,
        residue_indices=pocket_residue_indices,
        atom_indices=atom_indices,
        coords=coords,
        centroid=centroid,
        res_names=res_names,
    )


@dataclass
class PocketComparison:
    """Result of comparing two binding pockets.

    Summarizes which aligned pocket positions differ between two structures.
    """

    n_aligned_positions: int
    sequence_identity: float  # Fraction of identical residues at aligned positions
    identity_at_positions: list[bool]  # Per-position: same residue or not
    divergent_positions: list[int]  # Indices of non-identical positions
    pocket1_residues: list[str]
    pocket2_residues: list[str]

    @property
    def n_divergent(self) -> int:
        return len(self.divergent_positions)

    @property
    def divergence_fraction(self) -> float:
        if self.n_aligned_positions == 0:
            return 0.0
        return self.n_divergent / self.n_aligned_positions


def compare_binding_sites(
    site1: BindingSite,
    site2: BindingSite,
    alignment: list[tuple[int, int]] | None = None,
) -> PocketComparison:
    """Compare two binding pockets at aligned positions.

    For two homologous proteins (e.g., parasite target vs human ortholog),
    identify which binding-site residues differ and quantify the divergence.

    Parameters
    ----------
    site1 : BindingSite
        First binding pocket (e.g., parasite target).
    site2 : BindingSite
        Second binding pocket (e.g., human ortholog).
    alignment : list[tuple[int, int]] or None
        Pre-computed alignment mapping site1 position → site2 position.
        If None, uses a simple sequential alignment (assumes same-length
        pockets from homologous structures). For production use, this
        should come from a proper structural alignment.

    Returns
    -------
    PocketComparison
        Detailed comparison of the two pockets.
    """
    if alignment is None:
        # Simple sequential alignment — works for closely homologous structures
        n = min(site1.n_residues, site2.n_residues)
        alignment = [(i, i) for i in range(n)]

    identity_at_positions = []
    divergent_positions = []
    pocket1_res = []
    pocket2_res = []

    for pos, (i, j) in enumerate(alignment):
        r1 = site1.res_names[i] if i < len(site1.res_names) else "X"
        r2 = site2.res_names[j] if j < len(site2.res_names) else "X"
        pocket1_res.append(r1)
        pocket2_res.append(r2)

        is_same = r1 == r2
        identity_at_positions.append(is_same)
        if not is_same:
            divergent_positions.append(pos)

    n_aligned = len(alignment)
    n_identical = sum(identity_at_positions)
    seq_identity = n_identical / max(n_aligned, 1)

    return PocketComparison(
        n_aligned_positions=n_aligned,
        sequence_identity=seq_identity,
        identity_at_positions=identity_at_positions,
        divergent_positions=divergent_positions,
        pocket1_residues=pocket1_res,
        pocket2_residues=pocket2_res,
    )
