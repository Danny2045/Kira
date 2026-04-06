"""Scoring functions for the Kira drug repurposing pipeline.

Converts raw measurements (IC50, phase, publication counts) into
normalized 0-1 scores that feed the composite ranking algorithm.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from kira.targets import DEFAULT_ESSENTIALITY, TARGET_ESSENTIALITY

# ---------------------------------------------------------------------------
# Signal weights — used by the composite scoring functions
# ---------------------------------------------------------------------------

# V1 weights (Script 04)
WEIGHTS_V1: dict[str, float] = {
    "potency": 0.40,
    "target": 0.25,
    "confidence": 0.10,
    "drug_stage": 0.15,
    "multitarget": 0.10,
}

# V2 weights (Script 05 — adds structural similarity)
WEIGHTS_V2: dict[str, float] = {
    "potency": 0.30,
    "target": 0.20,
    "confidence": 0.05,
    "drug_stage": 0.15,
    "multitarget": 0.05,
    "similarity": 0.25,
}

# V3 weights (Script 06 — adds whole-organism evidence)
WEIGHTS_V3: dict[str, float] = {
    "potency": 0.25,
    "target": 0.15,
    "confidence": 0.05,
    "drug_stage": 0.10,
    "multitarget": 0.05,
    "similarity": 0.15,
    "whole_org": 0.25,
}


# ---------------------------------------------------------------------------
# Signal 1: Potency Score
# ---------------------------------------------------------------------------

def compute_potency_score(ic50_nm: float) -> float:
    """Convert IC50 in nanomolar to a 0-1 potency score.

    Uses pIC50 (negative log of IC50 in molar), then scales to 0-1.

    pIC50 scale:
        IC50 = 1 nM       -> pIC50 = 9.0  -> score ~ 1.0
        IC50 = 10 nM      -> pIC50 = 8.0  -> score ~ 0.88
        IC50 = 100 nM     -> pIC50 = 7.0  -> score ~ 0.75
        IC50 = 1000 nM    -> pIC50 = 6.0  -> score ~ 0.50
        IC50 = 10000 nM   -> pIC50 = 5.0  -> score ~ 0.25
        IC50 = 100000 nM  -> pIC50 = 4.0  -> score ~ 0.0

    Parameters
    ----------
    ic50_nm : float
        IC50 in nanomolar. NaN or <=0 returns 0.0.

    Returns
    -------
    float
        Potency score between 0.0 and 1.0.
    """
    if pd.isna(ic50_nm) or ic50_nm <= 0:
        return 0.0

    ic50_molar = ic50_nm * 1e-9
    pic50 = -np.log10(ic50_molar)

    # Scale to 0-1 range (pIC50 of 4 = 0, pIC50 of 9 = 1)
    score = (pic50 - 4.0) / (9.0 - 4.0)
    return float(max(0.0, min(1.0, score)))


# ---------------------------------------------------------------------------
# Signal 2: Target Essentiality Score
# ---------------------------------------------------------------------------

def compute_target_score(target_name: str) -> float:
    """Look up the essentiality score for a target.

    Parameters
    ----------
    target_name : str
        Target name as it appears in ChEMBL.

    Returns
    -------
    float
        Essentiality score between 0 and 1.
    """
    return TARGET_ESSENTIALITY.get(target_name, DEFAULT_ESSENTIALITY)


# ---------------------------------------------------------------------------
# Signal 3: Data Confidence Score
# ---------------------------------------------------------------------------

def compute_confidence_score(n_measurements: int) -> float:
    """Score based on number of independent measurements.

    More independent measurements = more confidence.

    1 measurement  -> 0.2 (low confidence)
    2-3            -> 0.5
    4-9            -> 0.75
    10+            -> 1.0

    Parameters
    ----------
    n_measurements : int
        Number of independent measurements.

    Returns
    -------
    float
        Confidence score between 0.1 and 1.0.
    """
    if n_measurements <= 0:
        return 0.1

    score = min(1.0, 0.2 + 0.3 * np.log2(n_measurements))
    return float(max(0.1, score))


# ---------------------------------------------------------------------------
# Signal 4: Drug Stage Score
# ---------------------------------------------------------------------------

DRUG_STAGE_SCORES: dict[float, float] = {
    4.0: 1.0,   # Approved drug
    3.0: 0.7,   # Phase III
    2.0: 0.5,   # Phase II
    1.0: 0.3,   # Phase I
    0.0: 0.1,   # Preclinical
}

DEFAULT_DRUG_STAGE: float = 0.1


def compute_drug_stage_score(max_phase: float) -> float:
    """Score a compound based on its clinical development stage.

    Parameters
    ----------
    max_phase : float
        Maximum clinical phase (0-4). NaN returns DEFAULT_DRUG_STAGE.

    Returns
    -------
    float
        Drug stage score between 0.1 and 1.0.
    """
    if pd.isna(max_phase):
        return DEFAULT_DRUG_STAGE

    return DRUG_STAGE_SCORES.get(float(max_phase), DEFAULT_DRUG_STAGE)


# ---------------------------------------------------------------------------
# Selectivity classification
# ---------------------------------------------------------------------------

def classify_selectivity(ratio: float) -> str:
    """Classify a selectivity ratio into a human-readable category.

    Parameters
    ----------
    ratio : float
        human_IC50 / parasite_IC50. Higher = more selective for parasite.

    Returns
    -------
    str
        One of "SELECTIVE", "MODERATE", "POOR", "COUNTER-SELECTIVE".
    """
    if ratio >= 10:
        return "SELECTIVE"
    elif ratio >= 3:
        return "MODERATE"
    elif ratio >= 1:
        return "POOR"
    else:
        return "COUNTER-SELECTIVE"


# ---------------------------------------------------------------------------
# Novelty classification
# ---------------------------------------------------------------------------

def classify_novelty(pub_count: int) -> str:
    """Classify a compound based on its PubMed publication count.

    Parameters
    ----------
    pub_count : int
        Number of publications mentioning compound + disease.
        Negative values indicate a query error.

    Returns
    -------
    str
        One of "NOVEL", "EMERGING", "KNOWN", "WELL-KNOWN", "ERROR".
    """
    if pub_count < 0:
        return "ERROR"
    elif pub_count == 0:
        return "NOVEL"
    elif pub_count <= 3:
        return "EMERGING"
    elif pub_count <= 10:
        return "KNOWN"
    else:
        return "WELL-KNOWN"
