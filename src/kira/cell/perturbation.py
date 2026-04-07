"""Perturbation response modeling for cell-state transitions.

Models how chemical or genetic perturbations shift cell state,
enabling prediction of therapeutic and off-target effects at
the cellular level.

This is the computational analogue of what Arc's Virtual Cell
and CZ CELLxGENE are building: a model of how cells respond
to interventions. Here we implement a simplified version that
integrates with Kira's selectivity analysis.

The core abstraction:
    P(cell_state_t+1 | cell_state_t, perturbation)

For NTD drug repurposing, this answers:
  - If we inhibit SmDHODH with CHEMBL155771, what happens to
    parasite pyrimidine biosynthesis?
  - If the compound also partially inhibits HsDHODH (at 1/30th
    potency), what happens to human bone marrow cells?

References:
    - Dixit et al. (2016) Cell — Perturb-seq
    - Lotfollahi et al. (2023) Nature Methods — CPA
    - Roohani et al. (2024) Nature Biotechnology — GEARS
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class PerturbationType(Enum):
    """Types of cellular perturbation."""

    CHEMICAL = "chemical"  # Small molecule inhibitor/activator
    GENETIC_KO = "genetic_knockout"  # Gene knockout (CRISPR-Cas9)
    GENETIC_KD = "genetic_knockdown"  # Gene knockdown (siRNA/shRNA)
    GENETIC_OE = "genetic_overexpression"  # Gene overexpression
    EPIGENETIC = "epigenetic"  # Epigenetic modifier
    ENVIRONMENTAL = "environmental"  # pH, temperature, nutrient


class CellStateChange(Enum):
    """Qualitative cell state transitions."""

    UNCHANGED = "unchanged"
    GROWTH_ARREST = "growth_arrest"
    APOPTOSIS = "apoptosis"
    DIFFERENTIATION = "differentiation"
    SENESCENCE = "senescence"
    STRESS_RESPONSE = "stress_response"
    METABOLIC_SHIFT = "metabolic_shift"
    PROLIFERATION = "proliferation"


@dataclass
class Perturbation:
    """A defined cellular perturbation.

    Attributes
    ----------
    name : str
        Perturbation identifier.
    perturbation_type : PerturbationType
        Type of perturbation.
    target_gene : str
        Primary gene target.
    compound_id : str
        ChEMBL ID if chemical perturbation.
    concentration_nm : float
        Concentration in nM (for chemical perturbations).
    duration_hours : float
        Duration of perturbation.
    selectivity_ratio : float
        Selectivity for parasite vs human target (if applicable).
    """

    name: str
    perturbation_type: PerturbationType
    target_gene: str
    compound_id: str = ""
    concentration_nm: float = 0.0
    duration_hours: float = 24.0
    selectivity_ratio: float = 1.0


@dataclass
class PathwayEffect:
    """Effect of a perturbation on a specific biological pathway.

    Attributes
    ----------
    pathway_name : str
        Name of the affected pathway.
    direction : str
        'up', 'down', or 'mixed'.
    magnitude : float
        Effect magnitude (0-1, where 1 is complete shutdown/activation).
    confidence : float
        Confidence in the prediction (0-1).
    mechanism : str
        Brief mechanistic explanation.
    """

    pathway_name: str
    direction: str
    magnitude: float
    confidence: float
    mechanism: str

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "pathway_name": self.pathway_name,
            "direction": self.direction,
            "magnitude": self.magnitude,
            "confidence": self.confidence,
            "mechanism": self.mechanism,
        }


@dataclass
class PerturbationResponse:
    """Predicted cellular response to a perturbation.

    Attributes
    ----------
    perturbation : Perturbation
        The applied perturbation.
    cell_type : str
        Cell type being perturbed.
    organism : str
        Organism (parasite or human).
    predicted_state_change : CellStateChange
        Most likely cell state transition.
    state_change_probability : float
        Probability of the predicted state change (0-1).
    pathway_effects : list[PathwayEffect]
        Effects on individual pathways.
    viability : float
        Predicted cell viability (0-1, where 1 = fully viable).
    off_target_risk : float
        Risk of off-target effects (0-1).
    therapeutic_window : float
        Ratio of toxic concentration to effective concentration.
    summary : str
        Human-readable summary.
    """

    perturbation: Perturbation
    cell_type: str
    organism: str
    predicted_state_change: CellStateChange
    state_change_probability: float
    pathway_effects: list[PathwayEffect]
    viability: float
    off_target_risk: float
    therapeutic_window: float
    summary: str = ""

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "perturbation": {
                "name": self.perturbation.name,
                "type": self.perturbation.perturbation_type.value,
                "target_gene": self.perturbation.target_gene,
                "compound_id": self.perturbation.compound_id,
                "concentration_nm": self.perturbation.concentration_nm,
            },
            "cell_type": self.cell_type,
            "organism": self.organism,
            "predicted_state_change": self.predicted_state_change.value,
            "state_change_probability": self.state_change_probability,
            "viability": self.viability,
            "off_target_risk": self.off_target_risk,
            "therapeutic_window": self.therapeutic_window,
            "pathway_effects": [pe.to_dict() for pe in self.pathway_effects],
            "summary": self.summary,
        }


# ---------------------------------------------------------------------------
# Pathway knowledge base for Kira's NTD targets
# ---------------------------------------------------------------------------

# Maps target gene -> downstream pathways affected by inhibition
_TARGET_PATHWAYS: dict[str, list[dict]] = {
    "DHODH": [
        {
            "pathway": "de_novo_pyrimidine_biosynthesis",
            "direction": "down",
            "magnitude": 0.9,
            "confidence": 0.95,
            "mechanism": "DHODH catalyzes the fourth step in de novo pyrimidine synthesis; "
            "inhibition blocks UMP production",
        },
        {
            "pathway": "mitochondrial_electron_transport",
            "direction": "mixed",
            "magnitude": 0.3,
            "confidence": 0.6,
            "mechanism": "DHODH is coupled to Complex III via ubiquinone; "
            "inhibition may partially reduce electron flow",
        },
        {
            "pathway": "dna_replication",
            "direction": "down",
            "magnitude": 0.7,
            "confidence": 0.85,
            "mechanism": "Pyrimidine depletion limits dNTP pools, slowing DNA replication",
        },
        {
            "pathway": "rna_transcription",
            "direction": "down",
            "magnitude": 0.5,
            "confidence": 0.75,
            "mechanism": "Reduced UTP/CTP pools limit RNA synthesis",
        },
    ],
    "HDAC8": [
        {
            "pathway": "chromatin_remodeling",
            "direction": "mixed",
            "magnitude": 0.7,
            "confidence": 0.85,
            "mechanism": "HDAC8 inhibition increases histone acetylation, "
            "altering gene expression programs",
        },
        {
            "pathway": "smc3_deacetylation",
            "direction": "down",
            "magnitude": 0.8,
            "confidence": 0.9,
            "mechanism": "HDAC8 is the specific deacetylase for cohesin subunit SMC3; "
            "inhibition disrupts chromosome segregation",
        },
        {
            "pathway": "smooth_muscle_contraction",
            "direction": "down",
            "magnitude": 0.4,
            "confidence": 0.5,
            "mechanism": "HDAC8 regulates smooth muscle contractility via "
            "cortactin deacetylation",
        },
    ],
    "TGR": [
        {
            "pathway": "thioredoxin_glutathione_redox",
            "direction": "down",
            "magnitude": 0.95,
            "confidence": 0.95,
            "mechanism": "SmTGR is the sole enzyme maintaining both thioredoxin and "
            "glutathione systems in S. mansoni; inhibition is lethal",
        },
        {
            "pathway": "reactive_oxygen_species",
            "direction": "up",
            "magnitude": 0.8,
            "confidence": 0.85,
            "mechanism": "Loss of antioxidant defense leads to ROS accumulation",
        },
        {
            "pathway": "protein_folding",
            "direction": "down",
            "magnitude": 0.5,
            "confidence": 0.6,
            "mechanism": "Thioredoxin system assists disulfide bond formation; "
            "disruption impairs protein folding",
        },
    ],
    "TXNRD1": [
        {
            "pathway": "thioredoxin_redox",
            "direction": "down",
            "magnitude": 0.7,
            "confidence": 0.9,
            "mechanism": "TXNRD1 maintains the thioredoxin system in human cells; "
            "partial inhibition reduces antioxidant capacity",
        },
        {
            "pathway": "reactive_oxygen_species",
            "direction": "up",
            "magnitude": 0.5,
            "confidence": 0.8,
            "mechanism": "Reduced TXNRD1 activity leads to moderate ROS increase",
        },
        {
            "pathway": "selenoprotein_function",
            "direction": "down",
            "magnitude": 0.4,
            "confidence": 0.7,
            "mechanism": "TXNRD1 is a selenoprotein; systemic selenoprotein function "
            "may be partially affected",
        },
    ],
}


def predict_perturbation_response(
    perturbation: Perturbation,
    cell_type: str = "generic",
    organism: str = "Homo sapiens",
) -> PerturbationResponse:
    """Predict cellular response to a perturbation.

    Uses pathway knowledge and dose-response modeling to predict
    how a cell will respond to target inhibition.

    Parameters
    ----------
    perturbation : Perturbation
        The perturbation to model.
    cell_type : str
        Cell type being perturbed.
    organism : str
        Organism.

    Returns
    -------
    PerturbationResponse
        Predicted response including pathway effects and viability.
    """
    target = perturbation.target_gene.upper()

    # Strip common prefixes (Sm, Hs, etc.)
    clean_target = target
    for prefix in ["SM", "HS"]:
        if clean_target.startswith(prefix) and len(clean_target) > 2:
            clean_target = clean_target[len(prefix):]
            break

    # Look up pathway effects
    pathway_data = _TARGET_PATHWAYS.get(clean_target, [])

    # Compute dose-dependent magnitude scaling
    # Using a simple Hill equation: effect = conc^n / (EC50^n + conc^n)
    conc = perturbation.concentration_nm
    ec50 = 100.0  # Assume EC50 of 100 nM for generic targets
    hill_n = 1.5
    if conc > 0:
        dose_fraction = (conc ** hill_n) / (ec50 ** hill_n + conc ** hill_n)
    else:
        dose_fraction = 0.0

    # Apply selectivity correction for off-target prediction
    # If selectivity_ratio = 30, the effective concentration at the
    # human target is 1/30th of the applied concentration
    if perturbation.selectivity_ratio > 1.0 and "sapiens" in organism.lower():
        effective_conc = conc / perturbation.selectivity_ratio
        if effective_conc > 0:
            dose_fraction = (effective_conc ** hill_n) / (
                ec50 ** hill_n + effective_conc ** hill_n
            )
        else:
            dose_fraction = 0.0

    # Build pathway effects scaled by dose
    pathway_effects: list[PathwayEffect] = []
    for pd in pathway_data:
        scaled_magnitude = pd["magnitude"] * dose_fraction
        pathway_effects.append(PathwayEffect(
            pathway_name=pd["pathway"],
            direction=pd["direction"],
            magnitude=round(scaled_magnitude, 3),
            confidence=pd["confidence"],
            mechanism=pd["mechanism"],
        ))

    # Predict viability based on pathway effects
    # Strong downregulation of essential pathways reduces viability
    viability = 1.0
    for pe in pathway_effects:
        if pe.direction == "down" and pe.magnitude > 0.3:
            viability *= 1.0 - (pe.magnitude * 0.5)
        elif pe.direction == "up" and pe.pathway_name == "reactive_oxygen_species":
            viability *= 1.0 - (pe.magnitude * 0.3)
    viability = max(0.0, min(1.0, viability))

    # Determine state change
    if viability < 0.2:
        state = CellStateChange.APOPTOSIS
        state_prob = 0.8
    elif viability < 0.5:
        state = CellStateChange.GROWTH_ARREST
        state_prob = 0.7
    elif any(
        pe.pathway_name == "de_novo_pyrimidine_biosynthesis" and pe.magnitude > 0.5
        for pe in pathway_effects
    ):
        state = CellStateChange.GROWTH_ARREST
        state_prob = 0.75
    elif dose_fraction < 0.1:
        state = CellStateChange.UNCHANGED
        state_prob = 0.9
    else:
        state = CellStateChange.STRESS_RESPONSE
        state_prob = 0.5

    # Off-target risk
    off_target_risk = 0.0
    if perturbation.selectivity_ratio < 5.0:
        off_target_risk = 1.0 / max(perturbation.selectivity_ratio, 0.1)
    elif perturbation.selectivity_ratio < 10.0:
        off_target_risk = 0.3
    else:
        off_target_risk = 0.1

    # Therapeutic window
    therapeutic_window = perturbation.selectivity_ratio

    # Summary
    summary = (
        f"Perturbation of {perturbation.target_gene} at "
        f"{perturbation.concentration_nm:.0f} nM in {cell_type} ({organism}): "
        f"predicted {state.value} (p={state_prob:.2f}), "
        f"viability={viability:.2f}, "
        f"{len(pathway_effects)} pathways affected."
    )

    return PerturbationResponse(
        perturbation=perturbation,
        cell_type=cell_type,
        organism=organism,
        predicted_state_change=state,
        state_change_probability=state_prob,
        pathway_effects=pathway_effects,
        viability=viability,
        off_target_risk=off_target_risk,
        therapeutic_window=therapeutic_window,
        summary=summary,
    )


def compare_on_off_target(
    compound_name: str,
    compound_id: str,
    parasite_target: str,
    human_orthologue: str,
    parasite_ic50_nm: float,
    selectivity_ratio: float,
    parasite_cell_type: str = "schistosomulum",
    human_cell_type: str = "hepatocyte",
) -> dict:
    """Compare on-target (parasite) vs off-target (human) effects.

    This is the cell-level analogue of Kira's selectivity analysis:
    given a selective compound, predict what happens to parasite cells
    vs human cells at therapeutic concentration.

    Parameters
    ----------
    compound_name : str
        Compound name or identifier.
    compound_id : str
        ChEMBL ID.
    parasite_target : str
        Parasite target gene.
    human_orthologue : str
        Human orthologue gene.
    parasite_ic50_nm : float
        IC50 against parasite target.
    selectivity_ratio : float
        IC50(human) / IC50(parasite).
    parasite_cell_type : str
        Parasite cell type to model.
    human_cell_type : str
        Human cell type to model.

    Returns
    -------
    dict
        Comparison of on-target vs off-target cellular responses.
    """
    # Therapeutic concentration = 10x IC50 for efficacy
    therapeutic_conc = parasite_ic50_nm * 10

    # On-target perturbation (parasite)
    on_target = Perturbation(
        name=f"{compound_name}_on_target",
        perturbation_type=PerturbationType.CHEMICAL,
        target_gene=parasite_target,
        compound_id=compound_id,
        concentration_nm=therapeutic_conc,
        selectivity_ratio=1.0,  # Full potency at parasite target
    )

    # Off-target perturbation (human)
    off_target = Perturbation(
        name=f"{compound_name}_off_target",
        perturbation_type=PerturbationType.CHEMICAL,
        target_gene=human_orthologue,
        compound_id=compound_id,
        concentration_nm=therapeutic_conc,
        selectivity_ratio=selectivity_ratio,
    )

    on_response = predict_perturbation_response(
        on_target, parasite_cell_type, "Schistosoma mansoni"
    )
    off_response = predict_perturbation_response(
        off_target, human_cell_type, "Homo sapiens"
    )

    return {
        "compound": compound_name,
        "compound_id": compound_id,
        "therapeutic_concentration_nm": therapeutic_conc,
        "selectivity_ratio": selectivity_ratio,
        "on_target": on_response.to_dict(),
        "off_target": off_response.to_dict(),
        "therapeutic_window": selectivity_ratio,
        "on_target_viability": on_response.viability,
        "off_target_viability": off_response.viability,
        "selectivity_at_cell_level": (
            "GOOD" if off_response.viability > 0.8 and on_response.viability < 0.5
            else "MODERATE" if off_response.viability > 0.6
            else "POOR"
        ),
    }
