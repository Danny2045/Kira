"""Local vs global similarity metrics.

Computes divergence at binding-site level vs whole-protein level.
Demonstrates that global sequence/embedding similarity can mask
critical local differences (the SmDHODH-HsDHODH insight).

The central observation from Kira's ESM-2 analysis: SmDHODH and HsDHODH
have 99.97% cosine similarity in their global embeddings, yet compounds
can achieve 30x selectivity between them. This module quantifies that
gap by computing local (binding-site) vs global (whole-protein) metrics
and identifying which positions are responsible.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kira.causality.binding_site import (
    BindingSite,
    compare_binding_sites,
)

# Amino acid property groups for physicochemical divergence analysis
AA_PROPERTIES: dict[str, dict[str, float]] = {
    # (hydrophobicity, charge_at_pH7, volume_A3, H-bond_donors, H-bond_acceptors)
    "ALA": {"hydrophobicity": 1.8, "charge": 0.0, "volume": 88.6, "hbd": 0, "hba": 0},
    "ARG": {"hydrophobicity": -4.5, "charge": 1.0, "volume": 173.4, "hbd": 5, "hba": 1},
    "ASN": {"hydrophobicity": -3.5, "charge": 0.0, "volume": 114.1, "hbd": 2, "hba": 2},
    "ASP": {"hydrophobicity": -3.5, "charge": -1.0, "volume": 111.1, "hbd": 0, "hba": 3},
    "CYS": {"hydrophobicity": 2.5, "charge": 0.0, "volume": 108.5, "hbd": 1, "hba": 1},
    "GLN": {"hydrophobicity": -3.5, "charge": 0.0, "volume": 143.8, "hbd": 2, "hba": 2},
    "GLU": {"hydrophobicity": -3.5, "charge": -1.0, "volume": 138.4, "hbd": 0, "hba": 3},
    "GLY": {"hydrophobicity": -0.4, "charge": 0.0, "volume": 60.1, "hbd": 0, "hba": 0},
    "HIS": {"hydrophobicity": -3.2, "charge": 0.5, "volume": 153.2, "hbd": 1, "hba": 2},
    "ILE": {"hydrophobicity": 4.5, "charge": 0.0, "volume": 166.7, "hbd": 0, "hba": 0},
    "LEU": {"hydrophobicity": 3.8, "charge": 0.0, "volume": 166.7, "hbd": 0, "hba": 0},
    "LYS": {"hydrophobicity": -3.9, "charge": 1.0, "volume": 168.6, "hbd": 3, "hba": 1},
    "MET": {"hydrophobicity": 1.9, "charge": 0.0, "volume": 162.9, "hbd": 0, "hba": 1},
    "PHE": {"hydrophobicity": 2.8, "charge": 0.0, "volume": 189.9, "hbd": 0, "hba": 0},
    "PRO": {"hydrophobicity": -1.6, "charge": 0.0, "volume": 112.7, "hbd": 0, "hba": 0},
    "SER": {"hydrophobicity": -0.8, "charge": 0.0, "volume": 89.0, "hbd": 1, "hba": 2},
    "THR": {"hydrophobicity": -0.7, "charge": 0.0, "volume": 116.1, "hbd": 1, "hba": 2},
    "TRP": {"hydrophobicity": -0.9, "charge": 0.0, "volume": 227.8, "hbd": 1, "hba": 1},
    "TYR": {"hydrophobicity": -1.3, "charge": 0.0, "volume": 193.6, "hbd": 1, "hba": 1},
    "VAL": {"hydrophobicity": 4.2, "charge": 0.0, "volume": 140.0, "hbd": 0, "hba": 0},
}

DEFAULT_PROPERTIES = {"hydrophobicity": 0.0, "charge": 0.0, "volume": 100.0, "hbd": 0, "hba": 0}


@dataclass
class DivergenceProfile:
    """Local vs global divergence analysis for two homologous proteins.

    Attributes
    ----------
    global_sequence_identity : float
        Fraction of identical residues across the full protein sequences.
    local_sequence_identity : float
        Fraction of identical residues within the binding pocket only.
    identity_gap : float
        global_sequence_identity - local_sequence_identity.
        Positive = pocket is MORE divergent than the protein overall.
    global_physicochemical_distance : float
        Mean physicochemical property distance across all aligned positions.
    local_physicochemical_distance : float
        Mean physicochemical property distance at binding-site positions.
    physicochemical_gap : float
        local - global physicochemical distance.
        Positive = pocket is MORE physicochemically divergent.
    per_position_divergence : list[PositionDivergence]
        Per-position divergence detail for binding-site residues.
    n_charge_changes : int
        Number of pocket positions with charge sign changes.
    n_hydrophobicity_flips : int
        Number of positions where hydrophobicity crosses the polar/nonpolar
        boundary (Kyte-Doolittle sign change).
    n_volume_changes : int
        Number of positions with >30 A^3 volume difference.
    """

    global_sequence_identity: float
    local_sequence_identity: float
    identity_gap: float
    global_physicochemical_distance: float
    local_physicochemical_distance: float
    physicochemical_gap: float
    per_position_divergence: list[PositionDivergence]
    n_charge_changes: int
    n_hydrophobicity_flips: int
    n_volume_changes: int

    @property
    def is_locally_divergent(self) -> bool:
        """True if the binding site is more divergent than the global protein."""
        return self.identity_gap > 0.05  # >5% more divergent locally

    def summary(self) -> str:
        """One-paragraph summary of the divergence analysis."""
        if self.is_locally_divergent:
            qualifier = "MORE divergent"
            gap_pct = self.identity_gap * 100
        else:
            qualifier = "similarly conserved"
            gap_pct = abs(self.identity_gap) * 100

        return (
            f"Global sequence identity: {self.global_sequence_identity:.1%}, "
            f"pocket identity: {self.local_sequence_identity:.1%} "
            f"({qualifier}, gap={gap_pct:.1f}%). "
            f"Pocket has {self.n_charge_changes} charge changes, "
            f"{self.n_hydrophobicity_flips} hydrophobicity flips, "
            f"and {self.n_volume_changes} volume changes."
        )


@dataclass
class PositionDivergence:
    """Divergence detail for a single aligned position.

    Attributes
    ----------
    position : int
        Aligned position index.
    res1 : str
        Residue in structure 1.
    res2 : str
        Residue in structure 2.
    is_identical : bool
        Whether the residues are the same.
    hydrophobicity_delta : float
        Difference in Kyte-Doolittle hydrophobicity.
    charge_delta : float
        Difference in charge at pH 7.
    volume_delta : float
        Difference in sidechain volume (A^3).
    physicochemical_distance : float
        Euclidean distance in normalized property space.
    """

    position: int
    res1: str
    res2: str
    is_identical: bool
    hydrophobicity_delta: float
    charge_delta: float
    volume_delta: float
    physicochemical_distance: float


def _get_properties(res_name: str) -> dict[str, float]:
    """Get amino acid properties, falling back to defaults for unknowns."""
    return AA_PROPERTIES.get(res_name, DEFAULT_PROPERTIES)


def _physicochemical_distance(res1: str, res2: str) -> float:
    """Compute normalized physicochemical distance between two residues.

    Uses Kyte-Doolittle hydrophobicity, charge, and volume,
    each normalized to [0, 1] range before computing Euclidean distance.
    """
    p1 = _get_properties(res1)
    p2 = _get_properties(res2)

    # Normalization ranges (from AA_PROPERTIES extremes)
    hydro_range = 9.0  # -4.5 to 4.5
    charge_range = 2.0  # -1 to 1
    volume_range = 167.7  # 60.1 to 227.8

    dh = (p1["hydrophobicity"] - p2["hydrophobicity"]) / hydro_range
    dc = (p1["charge"] - p2["charge"]) / charge_range
    dv = (p1["volume"] - p2["volume"]) / volume_range

    return float(np.sqrt(dh**2 + dc**2 + dv**2))


def compute_divergence_profile(
    site1: BindingSite,
    site2: BindingSite,
    global_seq1: str | None = None,
    global_seq2: str | None = None,
    alignment: list[tuple[int, int]] | None = None,
) -> DivergenceProfile:
    """Compute local vs global divergence between two homologous proteins.

    Parameters
    ----------
    site1 : BindingSite
        Binding pocket from structure 1.
    site2 : BindingSite
        Binding pocket from structure 2.
    global_seq1 : str or None
        Full protein sequence for structure 1 (one-letter codes).
        If None, only local divergence is computed; global is set to 0.
    global_seq2 : str or None
        Full protein sequence for structure 2.
    alignment : list[tuple[int, int]] or None
        Residue-level alignment. If None, uses sequential alignment.

    Returns
    -------
    DivergenceProfile
        Complete local vs global divergence analysis.
    """
    # Local (pocket) comparison
    pocket_cmp = compare_binding_sites(site1, site2, alignment)
    local_identity = pocket_cmp.sequence_identity

    # Global sequence identity (if sequences provided)
    if global_seq1 is not None and global_seq2 is not None:
        n_global = min(len(global_seq1), len(global_seq2))
        if n_global > 0:
            n_identical = sum(
                1 for a, b in zip(global_seq1[:n_global], global_seq2[:n_global])
                if a == b
            )
            global_identity = n_identical / n_global
        else:
            global_identity = 0.0
    else:
        global_identity = 0.0

    identity_gap = global_identity - local_identity

    # Per-position physicochemical divergence
    per_pos: list[PositionDivergence] = []
    n_charge_changes = 0
    n_hydro_flips = 0
    n_volume_changes = 0
    local_distances: list[float] = []

    for pos in range(pocket_cmp.n_aligned_positions):
        r1 = pocket_cmp.pocket1_residues[pos]
        r2 = pocket_cmp.pocket2_residues[pos]
        p1 = _get_properties(r1)
        p2 = _get_properties(r2)

        h_delta = p1["hydrophobicity"] - p2["hydrophobicity"]
        c_delta = p1["charge"] - p2["charge"]
        v_delta = p1["volume"] - p2["volume"]
        pc_dist = _physicochemical_distance(r1, r2)

        local_distances.append(pc_dist)

        # Count significant changes
        if abs(c_delta) >= 1.0:
            n_charge_changes += 1
        if p1["hydrophobicity"] * p2["hydrophobicity"] < 0 and abs(h_delta) > 2.0:
            n_hydro_flips += 1
        if abs(v_delta) > 30.0:
            n_volume_changes += 1

        per_pos.append(PositionDivergence(
            position=pos,
            res1=r1,
            res2=r2,
            is_identical=pocket_cmp.identity_at_positions[pos],
            hydrophobicity_delta=h_delta,
            charge_delta=c_delta,
            volume_delta=v_delta,
            physicochemical_distance=pc_dist,
        ))

    local_pc_dist = float(np.mean(local_distances)) if local_distances else 0.0

    # Global physicochemical distance (sample from full sequences)
    global_distances: list[float] = []
    if global_seq1 is not None and global_seq2 is not None:
        three_letter = {
            "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
            "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
            "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
            "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        }
        n_global = min(len(global_seq1), len(global_seq2))
        for k in range(n_global):
            r1_3 = three_letter.get(global_seq1[k], "ALA")
            r2_3 = three_letter.get(global_seq2[k], "ALA")
            global_distances.append(_physicochemical_distance(r1_3, r2_3))

    global_pc_dist = float(np.mean(global_distances)) if global_distances else 0.0
    pc_gap = local_pc_dist - global_pc_dist

    return DivergenceProfile(
        global_sequence_identity=global_identity,
        local_sequence_identity=local_identity,
        identity_gap=identity_gap,
        global_physicochemical_distance=global_pc_dist,
        local_physicochemical_distance=local_pc_dist,
        physicochemical_gap=pc_gap,
        per_position_divergence=per_pos,
        n_charge_changes=n_charge_changes,
        n_hydrophobicity_flips=n_hydro_flips,
        n_volume_changes=n_volume_changes,
    )


def rank_divergent_positions(
    profile: DivergenceProfile,
    min_distance: float = 0.1,
) -> list[PositionDivergence]:
    """Rank pocket positions by physicochemical divergence.

    Returns positions that are most different between the two proteins,
    filtered by a minimum distance threshold.

    Parameters
    ----------
    profile : DivergenceProfile
        Computed divergence profile.
    min_distance : float
        Minimum physicochemical distance to include.

    Returns
    -------
    list[PositionDivergence]
        Positions sorted by physicochemical_distance descending,
        filtered to those above min_distance.
    """
    filtered = [
        p for p in profile.per_position_divergence
        if p.physicochemical_distance >= min_distance and not p.is_identical
    ]
    return sorted(filtered, key=lambda p: p.physicochemical_distance, reverse=True)
