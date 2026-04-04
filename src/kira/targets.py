"""Target essentiality scores and orthologue mappings.

Domain knowledge about S. mansoni drug targets — how important each target
is to parasite survival and what the human orthologue is. These scores
encode expert judgment, not model output.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Target essentiality scores (0 to 1)
# ---------------------------------------------------------------------------
# These encode domain knowledge about how important each target is.
# Higher = more essential to the parasite = better drug target.

TARGET_ESSENTIALITY: dict[str, float] = {
    "Thioredoxin glutathione reductase": 1.0,
    # Single point of failure. No backup system in the worm.
    # Knockout is lethal. Best-validated target.

    "Histone deacetylase 8": 0.85,
    # Essential for gene regulation across life stages.
    # Well-validated. Selectivity window exists vs human HDAC8.

    "Dihydroorotate dehydrogenase (quinone), mitochondrial": 0.80,
    # Essential for pyrimidine synthesis (DNA/RNA building blocks).
    # Clear mechanism. Atovaquone already validates this target.

    "Cathepsin B1 isotype 1": 0.70,
    # Gut protease for hemoglobin digestion. Worm starves without it.
    # Less data available but strong biological rationale.

    "Voltage-activated calcium channel beta 1 subunit": 0.75,
    # Praziquantel's likely target family. Important but poorly characterized.

    "Voltage-activated calcium channel beta 2 subunit": 0.75,
    # Same as above.

    "Calcium channels (mechanism unclear)": 0.75,
    # Grouped calcium channel annotation.

    "Whole organism (phenotypic)": 0.70,
    # Phenotypic activity — target unknown but worm was killed.

    "Thioredoxin peroxidase": 0.65,
    # Downstream of SmTGR. Real target but no screening data.

    "Venus kinase receptor 2": 0.60,
    # No human equivalent (ideal selectivity). But only affects
    # reproduction, not survival. Worm lives but can't make eggs.

    "Sulfotransferase (prodrug activation)": 0.60,
    # Prodrug activation enzyme.

    "ATP-diphosphohydrolase 1": 0.55,
    # Immune evasion enzyme. Interesting biology but minimal data.

    "NAD-dependent protein deacetylase": 0.50,
    # Sirtuin. Stress response. Less validated. Weak data.

    "Unknown (whole-worm activity)": 0.50,
    # Whole-worm activity with unknown target.

    "None (negative control)": 0.0,
    # Negative control compounds.
}

DEFAULT_ESSENTIALITY: float = 0.3


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
# Human orthologues of S. mansoni drug targets in ChEMBL
# ---------------------------------------------------------------------------

ORTHOLOGUE_MAP: dict[str, dict[str, str]] = {
    "Histone deacetylase 8": {
        "parasite_id": "CHEMBL3797017",
        "human_name": "Histone deacetylase 8 (human)",
        "human_id": "CHEMBL3192",
        "notes": (
            "SmHDAC8 has structural differences in the active site "
            "loop region compared to human HDAC8. Selective inhibition "
            "is possible but requires careful compound design."
        ),
    },
    "Thioredoxin glutathione reductase": {
        "parasite_id": "CHEMBL6110",
        "human_name": "Thioredoxin reductase 1 (human)",
        "human_id": "CHEMBL3952",
        "notes": (
            "SmTGR is a fusion enzyme (TrxR + GR) unique to the parasite. "
            "Humans have separate TrxR and GR enzymes. Selectivity is "
            "theoretically favorable because the fusion creates a unique "
            "active site architecture."
        ),
    },
    "Dihydroorotate dehydrogenase (quinone), mitochondrial": {
        "parasite_id": "CHEMBL4523950",
        "human_name": "Dihydroorotate dehydrogenase (human)",
        "human_id": "CHEMBL1966",
        "notes": (
            "Both parasite and human DHODH are mitochondrial. "
            "Atovaquone shows some selectivity but this must be verified."
        ),
    },
}
