"""Binding-site divergence features for selectivity prediction.

Replaces ESM-2 global embeddings with LOCAL pocket-level features.
The central hypothesis: per-residue physicochemical divergence at the
binding site predicts selectivity better than global protein similarity.

Kira Script 19 showed ESM-2 fails at LODO (0.38, 0.36, 0.55 AUROC).
These features should transfer across diseases because they capture
the MECHANISM of selectivity (pocket divergence), not just correlations
in compound fingerprints.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kira.causality.divergence import AA_PROPERTIES, _physicochemical_distance


@dataclass
class PocketFeatures:
    """Binding-site divergence features for a single target pair.

    These are the features that replace ESM-2 cosine similarity.
    Each feature captures a specific aspect of HOW the binding pocket
    differs between parasite and human — the mechanistic basis of selectivity.
    """

    pair_name: str
    disease: str
    pocket_size: int  # Number of aligned pocket positions
    pocket_identity: float  # Fraction identical residues in pocket
    pocket_divergence: float  # 1 - pocket_identity
    mean_physicochemical_distance: float  # Mean per-position distance
    max_physicochemical_distance: float  # Worst-case position distance
    std_physicochemical_distance: float  # Variability in divergence
    n_charge_changes: int  # Positions with charge sign change
    n_hydrophobicity_flips: int  # Positions crossing polar/nonpolar
    n_volume_changes: int  # Positions with >30 A^3 volume difference
    fraction_charge_changes: float  # n_charge_changes / pocket_size
    fraction_hydro_flips: float
    fraction_volume_changes: float
    total_hydrophobicity_delta: float  # Sum of absolute hydro differences
    total_charge_delta: float  # Sum of absolute charge differences
    total_volume_delta: float  # Sum of absolute volume differences

    def to_array(self) -> np.ndarray:
        """Convert to a feature vector for ML.

        Returns
        -------
        np.ndarray
            (15,) float32 feature vector.
        """
        return np.array([
            self.pocket_size,
            self.pocket_identity,
            self.pocket_divergence,
            self.mean_physicochemical_distance,
            self.max_physicochemical_distance,
            self.std_physicochemical_distance,
            self.n_charge_changes,
            self.n_hydrophobicity_flips,
            self.n_volume_changes,
            self.fraction_charge_changes,
            self.fraction_hydro_flips,
            self.fraction_volume_changes,
            self.total_hydrophobicity_delta,
            self.total_charge_delta,
            self.total_volume_delta,
        ], dtype=np.float32)

    @staticmethod
    def feature_names() -> list[str]:
        """Names for each feature dimension."""
        return [
            "pocket_size",
            "pocket_identity",
            "pocket_divergence",
            "mean_physicochemical_distance",
            "max_physicochemical_distance",
            "std_physicochemical_distance",
            "n_charge_changes",
            "n_hydrophobicity_flips",
            "n_volume_changes",
            "fraction_charge_changes",
            "fraction_hydro_flips",
            "fraction_volume_changes",
            "total_hydrophobicity_delta",
            "total_charge_delta",
            "total_volume_delta",
        ]


def _get_props(res_1letter: str) -> dict[str, float]:
    """Get amino acid properties from one-letter code."""
    one_to_three = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    three = one_to_three.get(res_1letter.upper(), "ALA")
    return AA_PROPERTIES.get(three, AA_PROPERTIES["ALA"])


def compute_pocket_features(
    parasite_pocket_seq: str,
    human_pocket_seq: str,
    pair_name: str = "",
    disease: str = "",
) -> PocketFeatures:
    """Compute binding-site divergence features from pocket sequences.

    This is the core computation that replaces ESM-2 global embeddings.
    Instead of computing similarity across the ENTIRE protein, we compute
    physicochemical property differences at each POCKET position.

    Parameters
    ----------
    parasite_pocket_seq : str
        One-letter amino acid codes for parasite binding pocket residues.
    human_pocket_seq : str
        Aligned one-letter codes for human pocket residues.
        Must be same length as parasite_pocket_seq.
    pair_name : str
        Target pair identifier.
    disease : str
        Disease context.

    Returns
    -------
    PocketFeatures
        Complete binding-site divergence feature set.
    """
    n = min(len(parasite_pocket_seq), len(human_pocket_seq))
    if n == 0:
        return PocketFeatures(
            pair_name=pair_name, disease=disease, pocket_size=0,
            pocket_identity=0.0, pocket_divergence=1.0,
            mean_physicochemical_distance=0.0,
            max_physicochemical_distance=0.0,
            std_physicochemical_distance=0.0,
            n_charge_changes=0, n_hydrophobicity_flips=0,
            n_volume_changes=0, fraction_charge_changes=0.0,
            fraction_hydro_flips=0.0, fraction_volume_changes=0.0,
            total_hydrophobicity_delta=0.0, total_charge_delta=0.0,
            total_volume_delta=0.0,
        )

    n_identical = 0
    n_charge = 0
    n_hydro = 0
    n_volume = 0
    pc_distances = []
    hydro_deltas = []
    charge_deltas = []
    volume_deltas = []

    for i in range(n):
        r1 = parasite_pocket_seq[i]
        r2 = human_pocket_seq[i]

        p1 = _get_props(r1)
        p2 = _get_props(r2)

        # Identity
        if r1 == r2:
            n_identical += 1

        # Property deltas
        h_delta = abs(p1["hydrophobicity"] - p2["hydrophobicity"])
        c_delta = abs(p1["charge"] - p2["charge"])
        v_delta = abs(p1["volume"] - p2["volume"])

        hydro_deltas.append(h_delta)
        charge_deltas.append(c_delta)
        volume_deltas.append(v_delta)

        # Physicochemical distance (normalized Euclidean)
        one_to_three = {
            "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
            "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
            "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
            "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        }
        r1_3 = one_to_three.get(r1.upper(), "ALA")
        r2_3 = one_to_three.get(r2.upper(), "ALA")
        pc_dist = _physicochemical_distance(r1_3, r2_3)
        pc_distances.append(pc_dist)

        # Count significant changes
        if c_delta >= 1.0:
            n_charge += 1
        if (p1["hydrophobicity"] * p2["hydrophobicity"] < 0
                and h_delta > 2.0):
            n_hydro += 1
        if v_delta > 30.0:
            n_volume += 1

    pocket_identity = n_identical / n
    pc_arr = np.array(pc_distances)

    return PocketFeatures(
        pair_name=pair_name,
        disease=disease,
        pocket_size=n,
        pocket_identity=pocket_identity,
        pocket_divergence=1.0 - pocket_identity,
        mean_physicochemical_distance=float(np.mean(pc_arr)),
        max_physicochemical_distance=float(np.max(pc_arr)),
        std_physicochemical_distance=float(np.std(pc_arr)),
        n_charge_changes=n_charge,
        n_hydrophobicity_flips=n_hydro,
        n_volume_changes=n_volume,
        fraction_charge_changes=n_charge / n,
        fraction_hydro_flips=n_hydro / n,
        fraction_volume_changes=n_volume / n,
        total_hydrophobicity_delta=float(sum(hydro_deltas)),
        total_charge_delta=float(sum(charge_deltas)),
        total_volume_delta=float(sum(volume_deltas)),
    )
