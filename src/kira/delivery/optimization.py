"""Multi-objective optimization for LNP formulation design.

Optimizes formulation parameters across competing objectives:
efficacy (on-target delivery), safety (off-target minimization),
manufacturability, and stability.

Uses Pareto-front exploration to identify non-dominated formulations,
since no single formulation maximizes all objectives simultaneously.

This is the delivery layer's equivalent of the selectivity problem:
you cannot optimize for one tissue without affecting others.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kira.delivery.biodistribution import (
    DeliveryProbabilityChain,
    TissueDistribution,
    estimate_delivery_chain,
    predict_biodistribution,
)
from kira.delivery.formulation import (
    AdministrationRoute,
    HelperLipidClass,
    IonizableLipidClass,
    LNPFormulation,
    PEGLipidClass,
)


@dataclass
class FormulationObjectives:
    """Evaluated objectives for a single formulation.

    Attributes
    ----------
    formulation : LNPFormulation
        The formulation being evaluated.
    efficacy : float
        On-target delivery probability (0-1, higher is better).
    safety : float
        Safety score based on off-target exposure (0-1, higher is better).
    manufacturability : float
        Ease of manufacturing score (0-1, higher is better).
    stability : float
        Predicted formulation stability (0-1, higher is better).
    overall : float
        Weighted overall score.
    biodistribution : TissueDistribution
        Full tissue distribution prediction.
    delivery_chain : DeliveryProbabilityChain
        Step-by-step delivery probabilities.
    """

    formulation: LNPFormulation
    efficacy: float
    safety: float
    manufacturability: float
    stability: float
    overall: float
    biodistribution: TissueDistribution
    delivery_chain: DeliveryProbabilityChain

    @property
    def is_viable(self) -> bool:
        """Whether the formulation meets minimum viability thresholds."""
        return (
            self.efficacy > 0.001
            and self.safety > 0.3
            and self.manufacturability > 0.3
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "formulation": self.formulation.to_dict(),
            "efficacy": self.efficacy,
            "safety": self.safety,
            "manufacturability": self.manufacturability,
            "stability": self.stability,
            "overall": self.overall,
            "is_viable": self.is_viable,
            "biodistribution": self.biodistribution.to_dict(),
            "delivery_chain": self.delivery_chain.to_dict(),
        }


def _score_manufacturability(formulation: LNPFormulation) -> float:
    """Score how easy a formulation is to manufacture.

    Factors:
    - Standard lipids are easier than custom
    - Moderate N:P ratios are easier to control
    - Standard sizes are easier to achieve
    - Fewer targeting peptides = simpler manufacturing
    """
    score = 0.8  # Base score for standard LNP

    # Lipid availability
    if formulation.ionizable_lipid == IonizableLipidClass.CUSTOM:
        score *= 0.5  # Custom lipids require synthesis
    elif formulation.ionizable_lipid in {
        IonizableLipidClass.SM102,
        IonizableLipidClass.ALC0315,
    }:
        score *= 1.1  # Well-characterized, commercially available

    # N:P ratio controllability
    if 4 <= formulation.np_ratio <= 8:
        score *= 1.0
    elif 2 <= formulation.np_ratio <= 12:
        score *= 0.9
    else:
        score *= 0.7

    # Size achievability
    if 60 <= formulation.target_size_nm <= 120:
        score *= 1.0
    elif 40 <= formulation.target_size_nm <= 150:
        score *= 0.85
    else:
        score *= 0.6

    # Targeting peptide complexity
    n_peptides = len(formulation.targeting_peptides)
    if n_peptides == 0:
        score *= 1.0
    elif n_peptides == 1:
        score *= 0.8
    elif n_peptides <= 3:
        score *= 0.6
    else:
        score *= 0.4

    # Payload complexity
    if formulation.payload_size_nt > 5000:
        score *= 0.85  # Large payloads harder to encapsulate

    return min(1.0, max(0.0, score))


def _score_stability(formulation: LNPFormulation) -> float:
    """Score predicted formulation stability.

    Factors:
    - PEG density affects colloidal stability
    - Cholesterol content affects membrane rigidity
    - Size affects Ostwald ripening tendency
    """
    score = 0.7

    # PEG stabilization
    if 1.0 <= formulation.peg_lipid_mol_pct <= 3.0:
        score *= 1.1
    elif formulation.peg_lipid_mol_pct < 0.5:
        score *= 0.7  # Insufficient steric stabilization

    # Cholesterol for membrane rigidity
    if 35 <= formulation.cholesterol_mol_pct <= 45:
        score *= 1.1
    elif formulation.cholesterol_mol_pct < 25:
        score *= 0.75

    # Size stability (smaller particles have higher surface energy)
    if formulation.target_size_nm < 50:
        score *= 0.85
    elif formulation.target_size_nm > 150:
        score *= 0.9

    return min(1.0, max(0.0, score))


def _score_safety(
    biodist: TissueDistribution,
    target_tissue: str,
) -> float:
    """Score safety based on off-target tissue exposure.

    Lower off-target exposure = higher safety score.
    Special penalties for brain and heart exposure when not targeted.
    """
    fractions = biodist.tissue_fractions
    target_frac = fractions.get(target_tissue, 0.0)

    # Base safety from on-target fraction
    score = 0.3 + 0.7 * target_frac

    # Penalty for critical organ exposure when not targeted
    sensitive_organs = {"brain", "heart", "kidney"}
    for organ in sensitive_organs:
        if organ != target_tissue:
            organ_frac = fractions.get(organ, 0.0)
            if organ_frac > 0.1:
                score *= 0.7
            elif organ_frac > 0.05:
                score *= 0.85

    # High liver exposure is expected for IV and less concerning
    # (liver toxicity is manageable at moderate doses)

    return min(1.0, max(0.0, score))


def evaluate_formulation(
    formulation: LNPFormulation,
    target_tissue: str = "liver",
    weights: dict[str, float] | None = None,
) -> FormulationObjectives:
    """Evaluate a formulation across all objectives.

    Parameters
    ----------
    formulation : LNPFormulation
        Formulation to evaluate.
    target_tissue : str
        Intended target tissue.
    weights : dict[str, float] or None
        Objective weights for overall score.
        Keys: 'efficacy', 'safety', 'manufacturability', 'stability'.
        Defaults to equal weights.

    Returns
    -------
    FormulationObjectives
        Complete multi-objective evaluation.
    """
    if weights is None:
        weights = {
            "efficacy": 0.4,
            "safety": 0.3,
            "manufacturability": 0.15,
            "stability": 0.15,
        }

    # Predict biodistribution
    biodist = predict_biodistribution(formulation, target_tissue)

    # Estimate delivery chain
    chain = estimate_delivery_chain(formulation, target_tissue)

    # Compute individual objectives
    efficacy = chain.p_effective
    safety = _score_safety(biodist, target_tissue)
    manuf = _score_manufacturability(formulation)
    stab = _score_stability(formulation)

    # Weighted overall score
    overall = (
        weights.get("efficacy", 0.25) * efficacy
        + weights.get("safety", 0.25) * safety
        + weights.get("manufacturability", 0.25) * manuf
        + weights.get("stability", 0.25) * stab
    )

    return FormulationObjectives(
        formulation=formulation,
        efficacy=efficacy,
        safety=safety,
        manufacturability=manuf,
        stability=stab,
        overall=overall,
        biodistribution=biodist,
        delivery_chain=chain,
    )


def _dominates(a: FormulationObjectives, b: FormulationObjectives) -> bool:
    """Check if formulation a Pareto-dominates formulation b.

    a dominates b if a is at least as good in all objectives
    and strictly better in at least one.
    """
    objectives_a = [a.efficacy, a.safety, a.manufacturability, a.stability]
    objectives_b = [b.efficacy, b.safety, b.manufacturability, b.stability]

    at_least_as_good = all(oa >= ob for oa, ob in zip(objectives_a, objectives_b))
    strictly_better = any(oa > ob for oa, ob in zip(objectives_a, objectives_b))

    return at_least_as_good and strictly_better


def find_pareto_front(
    evaluations: list[FormulationObjectives],
) -> list[FormulationObjectives]:
    """Find the Pareto-optimal formulations.

    A formulation is Pareto-optimal if no other formulation
    is better in all objectives simultaneously.

    Parameters
    ----------
    evaluations : list[FormulationObjectives]
        Evaluated formulations.

    Returns
    -------
    list[FormulationObjectives]
        Non-dominated (Pareto-optimal) formulations.
    """
    if not evaluations:
        return []

    pareto: list[FormulationObjectives] = []

    for candidate in evaluations:
        is_dominated = False
        for other in evaluations:
            if other is not candidate and _dominates(other, candidate):
                is_dominated = True
                break
        if not is_dominated:
            pareto.append(candidate)

    # Sort by overall score descending
    pareto.sort(key=lambda e: e.overall, reverse=True)
    return pareto


def grid_search_formulations(
    target_tissue: str = "liver",
    payload_type: str = "mRNA",
    payload_size_nt: int = 1000,
    routes: list[AdministrationRoute] | None = None,
    n_lipid_steps: int = 3,
) -> list[FormulationObjectives]:
    """Grid search over formulation parameter space.

    Generates a grid of formulations varying key parameters
    and evaluates each one.

    Parameters
    ----------
    target_tissue : str
        Target tissue for delivery.
    payload_type : str
        Type of payload.
    payload_size_nt : int
        Payload size in nucleotides.
    routes : list[AdministrationRoute] or None
        Routes to explore. Defaults to [IV, IM].
    n_lipid_steps : int
        Number of steps for lipid percentage grid.

    Returns
    -------
    list[FormulationObjectives]
        All evaluated formulations, sorted by overall score.
    """
    if routes is None:
        routes = [AdministrationRoute.IV, AdministrationRoute.IM]

    ionizable_lipids = [
        IonizableLipidClass.SM102,
        IonizableLipidClass.ALC0315,
        IonizableLipidClass.DLIN_MC3,
    ]
    helper_lipids = [HelperLipidClass.DSPC, HelperLipidClass.DOPE]
    peg_lipids = [PEGLipidClass.DMG_PEG2000, PEGLipidClass.DSPE_PEG2000]

    ionizable_pcts = np.linspace(35, 55, n_lipid_steps)

    evaluations: list[FormulationObjectives] = []
    idx = 0

    for route in routes:
        for il in ionizable_lipids:
            for hl in helper_lipids:
                for pl in peg_lipids:
                    for il_pct in ionizable_pcts:
                        # Balance the remaining lipid budget
                        helper_pct = 10.0
                        peg_pct = 1.5
                        chol_pct = 100.0 - il_pct - helper_pct - peg_pct

                        if chol_pct < 20 or chol_pct > 55:
                            continue

                        formulation = LNPFormulation(
                            name=f"grid_{idx:04d}",
                            ionizable_lipid=il,
                            ionizable_lipid_mol_pct=float(il_pct),
                            helper_lipid=hl,
                            helper_lipid_mol_pct=helper_pct,
                            cholesterol_mol_pct=chol_pct,
                            peg_lipid=pl,
                            peg_lipid_mol_pct=peg_pct,
                            np_ratio=6.0,
                            target_size_nm=80.0,
                            route=route,
                            payload_type=payload_type,
                            payload_size_nt=payload_size_nt,
                        )

                        evaluation = evaluate_formulation(
                            formulation, target_tissue
                        )
                        evaluations.append(evaluation)
                        idx += 1

    evaluations.sort(key=lambda e: e.overall, reverse=True)
    return evaluations
