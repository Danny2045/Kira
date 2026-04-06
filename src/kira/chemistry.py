"""Chemistry utilities — ADMET filtering and molecular property computation.

Wraps RDKit functionality for computing drug-likeness properties,
Lipinski Rule of Five checks, and structural similarity via Morgan
fingerprints.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

# ---------------------------------------------------------------------------
# ADMET property computation
# ---------------------------------------------------------------------------

@dataclass
class ADMETProfile:
    """Drug-likeness profile for a single compound.

    Attributes
    ----------
    molecule_chembl_id : str
        ChEMBL compound identifier.
    smiles : str
        Canonical SMILES string.
    mw : float
        Molecular weight (Lipinski: <500).
    logp : float
        Partition coefficient (Lipinski: <5).
    hbd : int
        Hydrogen bond donors (Lipinski: <5).
    hba : int
        Hydrogen bond acceptors (Lipinski: <10).
    tpsa : float
        Topological polar surface area (<140 for oral absorption).
    rotatable_bonds : int
        Number of rotatable bonds (<10 preferred).
    qed : float or None
        Quantitative Estimate of Drug-likeness (0-1).
    lipinski_violations : int
        Count of Lipinski Rule of Five violations.
    drug_likeness : str
        Overall assessment: "excellent", "acceptable", "marginal", "poor".
    """

    molecule_chembl_id: str
    smiles: str
    mw: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    rotatable_bonds: int
    qed: float | None
    lipinski_violations: int
    drug_likeness: str


def count_lipinski_violations(
    mw: float, logp: float, hbd: int, hba: int
) -> int:
    """Count Lipinski Rule of Five violations.

    Parameters
    ----------
    mw : float
        Molecular weight.
    logp : float
        LogP partition coefficient.
    hbd : int
        Hydrogen bond donors.
    hba : int
        Hydrogen bond acceptors.

    Returns
    -------
    int
        Number of violations (0-4).
    """
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1
    return violations


def classify_drug_likeness(violations: int) -> str:
    """Classify drug-likeness based on Lipinski violation count.

    Parameters
    ----------
    violations : int
        Number of Lipinski violations.

    Returns
    -------
    str
        One of "excellent", "acceptable", "marginal", "poor".
    """
    if violations == 0:
        return "excellent"
    elif violations == 1:
        return "acceptable"
    elif violations == 2:
        return "marginal"
    else:
        return "poor"


def compute_admet_multiplier(
    lipinski_violations: float | None, tpsa: float | None
) -> float:
    """Compute ADMET penalty multiplier for final score adjustment.

    Parameters
    ----------
    lipinski_violations : float or None
        Number of Lipinski violations. None = unknown.
    tpsa : float or None
        Topological polar surface area.

    Returns
    -------
    float
        Multiplier between 0.0 and 1.0.
    """
    if pd.isna(lipinski_violations) or lipinski_violations is None:
        return 0.8  # Unknown = slight penalty

    multiplier = 1.0

    if lipinski_violations == 0:
        multiplier *= 1.0
    elif lipinski_violations == 1:
        multiplier *= 0.85
    elif lipinski_violations == 2:
        multiplier *= 0.60
    else:
        multiplier *= 0.35

    if tpsa is not None and not pd.isna(tpsa) and tpsa > 140:
        multiplier *= 0.75

    return multiplier
