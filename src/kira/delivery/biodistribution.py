"""Biodistribution prediction for LNP formulations.

Models tissue-level distribution of lipid nanoparticles based on
formulation parameters, administration route, and targeting ligands.

The core insight: standard LNPs are liver-tropic because serum ApoE
adsorbs onto the particle surface and mediates hepatocyte uptake via
LDLR. Breaking this default requires either:
  1. Surface modification (targeting peptides, antibodies)
  2. Lipid composition changes (SORT mechanism — Cheng et al. 2020)
  3. Route changes (IM/SC/IN instead of IV)
  4. PEG density tuning (affects protein corona composition)

This module implements a first-principles model of these effects.

References:
    - Akinc et al. (2019) Nature Nanotechnology — ApoE-mediated uptake
    - Cheng et al. (2020) Nature Nanotechnology — SORT lipids
    - Dilliard et al. (2021) Molecular Therapy — extrahepatic delivery
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from kira.delivery.formulation import (
    AdministrationRoute,
    HelperLipidClass,
    IonizableLipidClass,
    LNPFormulation,
    PEGLipidClass,
)

# ---------------------------------------------------------------------------
# Tissue distribution model
# ---------------------------------------------------------------------------

# Default tissue names tracked by the model
TRACKED_TISSUES = [
    "liver",
    "spleen",
    "lung",
    "heart",
    "kidney",
    "brain",
    "muscle",
    "lymph_node",
    "tumor",
    "bone_marrow",
]


@dataclass
class TissueDistribution:
    """Predicted tissue-level distribution of an LNP formulation.

    Attributes
    ----------
    formulation_name : str
        Name of the formulation being evaluated.
    tissue_fractions : dict[str, float]
        Tissue name -> predicted fraction of total dose (0-1).
        These are relative distributions, not absolute concentrations.
    dominant_tissue : str
        Tissue receiving the largest fraction.
    liver_fraction : float
        Fraction going to liver (the default sink for most IV LNPs).
    extrahepatic_fraction : float
        Total fraction reaching non-liver tissues.
    therapeutic_index_estimate : float
        Ratio of on-target to off-target exposure (higher is better).
    target_tissue : str
        Intended target tissue for therapeutic effect.
    warnings : list[str]
        Any warnings about the prediction.
    """

    formulation_name: str
    tissue_fractions: dict[str, float]
    dominant_tissue: str
    liver_fraction: float
    extrahepatic_fraction: float
    therapeutic_index_estimate: float
    target_tissue: str
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Human-readable summary of the distribution."""
        top3 = sorted(
            self.tissue_fractions.items(), key=lambda x: x[1], reverse=True
        )[:3]
        top3_str = ", ".join(f"{t}: {f:.1%}" for t, f in top3)
        return (
            f"Distribution for {self.formulation_name}: "
            f"dominant={self.dominant_tissue} ({self.tissue_fractions.get(self.dominant_tissue, 0):.1%}), "
            f"liver={self.liver_fraction:.1%}, "
            f"extrahepatic={self.extrahepatic_fraction:.1%}. "
            f"Top 3: {top3_str}. "
            f"TI estimate: {self.therapeutic_index_estimate:.2f}"
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "formulation_name": self.formulation_name,
            "tissue_fractions": self.tissue_fractions,
            "dominant_tissue": self.dominant_tissue,
            "liver_fraction": self.liver_fraction,
            "extrahepatic_fraction": self.extrahepatic_fraction,
            "therapeutic_index_estimate": self.therapeutic_index_estimate,
            "target_tissue": self.target_tissue,
            "warnings": self.warnings,
        }


# ---------------------------------------------------------------------------
# Route-dependent baseline distributions
# ---------------------------------------------------------------------------

# These baseline distributions represent the "default" tissue uptake
# for each route before formulation-specific adjustments.
# Values are approximate and derived from published biodistribution studies.

_ROUTE_BASELINES: dict[AdministrationRoute, dict[str, float]] = {
    AdministrationRoute.IV: {
        "liver": 0.70,
        "spleen": 0.15,
        "lung": 0.05,
        "heart": 0.02,
        "kidney": 0.03,
        "brain": 0.001,
        "muscle": 0.01,
        "lymph_node": 0.01,
        "tumor": 0.02,
        "bone_marrow": 0.019,
    },
    AdministrationRoute.IM: {
        "liver": 0.10,
        "spleen": 0.05,
        "lung": 0.02,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.001,
        "muscle": 0.50,
        "lymph_node": 0.25,
        "tumor": 0.01,
        "bone_marrow": 0.019,
    },
    AdministrationRoute.SC: {
        "liver": 0.08,
        "spleen": 0.04,
        "lung": 0.02,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.001,
        "muscle": 0.15,
        "lymph_node": 0.60,
        "tumor": 0.01,
        "bone_marrow": 0.019,
    },
    AdministrationRoute.IN: {
        "liver": 0.05,
        "spleen": 0.03,
        "lung": 0.65,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.05,
        "muscle": 0.02,
        "lymph_node": 0.10,
        "tumor": 0.01,
        "bone_marrow": 0.02,
    },
    AdministrationRoute.IT: {
        "liver": 0.05,
        "spleen": 0.02,
        "lung": 0.02,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.75,
        "muscle": 0.01,
        "lymph_node": 0.02,
        "tumor": 0.01,
        "bone_marrow": 0.09,
    },
    AdministrationRoute.INTRATUMORAL: {
        "liver": 0.10,
        "spleen": 0.05,
        "lung": 0.02,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.001,
        "muscle": 0.02,
        "lymph_node": 0.05,
        "tumor": 0.70,
        "bone_marrow": 0.019,
    },
    AdministrationRoute.NEBULIZED: {
        "liver": 0.05,
        "spleen": 0.03,
        "lung": 0.70,
        "heart": 0.01,
        "kidney": 0.02,
        "brain": 0.01,
        "muscle": 0.02,
        "lymph_node": 0.10,
        "tumor": 0.01,
        "bone_marrow": 0.05,
    },
}


# ---------------------------------------------------------------------------
# Lipid-specific modifiers
# ---------------------------------------------------------------------------

def _compute_liver_tropism_modifier(formulation: LNPFormulation) -> float:
    """Compute how much the lipid composition shifts liver tropism.

    Returns a multiplier (>1 increases liver fraction, <1 decreases it).

    Key factors:
    - MC3 and C12-200 are strongly hepatocyte-selective
    - Higher ionizable lipid % increases liver uptake
    - DOPE helper lipid slightly reduces liver specificity
    - Higher PEG density reduces ApoE adsorption -> less liver
    """
    modifier = 1.0

    # Ionizable lipid class effect
    lipid_liver_affinity = {
        IonizableLipidClass.DLIN_MC3: 1.3,  # Strong ApoE binder
        IonizableLipidClass.C12_200: 1.4,  # Very liver-tropic
        IonizableLipidClass.SM102: 1.0,  # Neutral baseline
        IonizableLipidClass.ALC0315: 0.95,  # Slightly less liver-tropic
        IonizableLipidClass.LIPID5: 0.9,  # Designed for lower liver uptake
        IonizableLipidClass.CUSTOM: 1.0,
    }
    modifier *= lipid_liver_affinity.get(formulation.ionizable_lipid, 1.0)

    # Ionizable lipid percentage effect (higher % -> more liver)
    # Baseline is 50%. Each 10% above/below shifts modifier by ~10%
    pct_effect = 1.0 + (formulation.ionizable_lipid_mol_pct - 50.0) * 0.01
    modifier *= max(0.5, min(1.5, pct_effect))

    # Helper lipid effect
    if formulation.helper_lipid == HelperLipidClass.DOPE:
        modifier *= 0.85  # Fusogenic lipid reduces liver specificity
    elif formulation.helper_lipid == HelperLipidClass.DSPC:
        modifier *= 1.05  # Standard, slightly liver-promoting

    # PEG density effect (more PEG -> less ApoE -> less liver)
    # Baseline is 1.5%
    peg_effect = 1.0 - (formulation.peg_lipid_mol_pct - 1.5) * 0.1
    modifier *= max(0.5, min(1.3, peg_effect))

    # PEG type effect
    if formulation.peg_lipid == PEGLipidClass.DSPE_PEG2000:
        modifier *= 0.85  # Longer acyl chain, more stable PEG coat, less ApoE

    return modifier


def _compute_targeting_modifier(
    formulation: LNPFormulation,
) -> dict[str, float]:
    """Compute tissue-level modifiers from targeting peptides.

    Returns a dict of tissue -> multiplier. Tissues with targeting
    peptide affinity get increased fractions.
    """
    modifiers: dict[str, float] = {t: 1.0 for t in TRACKED_TISSUES}

    if not formulation.targeting_peptides:
        return modifiers

    for peptide in formulation.targeting_peptides:
        for tissue, affinity in peptide.tissue_affinity.items():
            if tissue in modifiers:
                # Affinity of 1.0 doubles the tissue fraction
                # Affinity of 0.5 increases by 50%
                modifiers[tissue] *= 1.0 + affinity

    return modifiers


def _compute_size_modifier(formulation: LNPFormulation) -> dict[str, float]:
    """Compute tissue-level modifiers from particle size.

    Smaller particles (<100 nm) penetrate more tissues.
    Larger particles (>150 nm) are trapped more by spleen/liver.
    """
    modifiers: dict[str, float] = {t: 1.0 for t in TRACKED_TISSUES}
    size = formulation.target_size_nm

    if size < 60:
        # Very small — better tissue penetration, some kidney clearance
        modifiers["kidney"] *= 1.5
        modifiers["tumor"] *= 1.3
        modifiers["liver"] *= 0.8
    elif size < 100:
        # Optimal for most applications
        modifiers["tumor"] *= 1.1  # EPR effect
    elif size < 150:
        # Standard, slightly more spleen uptake
        modifiers["spleen"] *= 1.2
    else:
        # Large — trapped by liver sinusoids and spleen
        modifiers["liver"] *= 1.3
        modifiers["spleen"] *= 1.5
        modifiers["tumor"] *= 0.7

    return modifiers


# ---------------------------------------------------------------------------
# Main prediction function
# ---------------------------------------------------------------------------

def predict_biodistribution(
    formulation: LNPFormulation,
    target_tissue: str = "liver",
) -> TissueDistribution:
    """Predict tissue-level biodistribution of an LNP formulation.

    Combines route-dependent baselines with formulation-specific
    modifiers (lipid composition, PEG density, particle size,
    targeting peptides) to estimate tissue fractions.

    Parameters
    ----------
    formulation : LNPFormulation
        Complete formulation specification.
    target_tissue : str
        Intended target tissue (used for therapeutic index calculation).

    Returns
    -------
    TissueDistribution
        Predicted tissue distribution with therapeutic index estimate.
    """
    warnings: list[str] = []

    # Start with route-dependent baseline
    baseline = _ROUTE_BASELINES.get(
        formulation.route,
        _ROUTE_BASELINES[AdministrationRoute.IV],
    )
    fractions = dict(baseline)

    # Apply liver tropism modifier
    liver_mod = _compute_liver_tropism_modifier(formulation)
    fractions["liver"] *= liver_mod

    # Apply targeting peptide modifiers
    targeting_mods = _compute_targeting_modifier(formulation)
    for tissue, mod in targeting_mods.items():
        if tissue in fractions:
            fractions[tissue] *= mod

    # Apply size modifiers
    size_mods = _compute_size_modifier(formulation)
    for tissue, mod in size_mods.items():
        if tissue in fractions:
            fractions[tissue] *= mod

    # Renormalize to sum to 1.0
    total = sum(fractions.values())
    if total > 0:
        fractions = {t: f / total for t, f in fractions.items()}
    else:
        warnings.append("All tissue fractions are zero — check formulation parameters")
        fractions = {t: 1.0 / len(fractions) for t in fractions}

    # Find dominant tissue
    dominant = max(fractions, key=lambda t: fractions[t])

    # Compute therapeutic index
    liver_frac = fractions.get("liver", 0.0)
    extrahepatic = 1.0 - liver_frac
    target_frac = fractions.get(target_tissue, 0.0)

    # TI = on-target / average off-target
    off_target_sum = sum(
        f for t, f in fractions.items() if t != target_tissue
    )
    n_off = len(fractions) - 1
    avg_off_target = off_target_sum / max(n_off, 1)
    ti = target_frac / max(avg_off_target, 1e-6)

    # Generate warnings
    if target_tissue != "liver" and liver_frac > 0.5:
        warnings.append(
            f"High liver uptake ({liver_frac:.1%}) despite targeting {target_tissue}. "
            f"Consider targeting peptides or route change."
        )

    if target_frac < 0.1 and target_tissue != "brain":
        warnings.append(
            f"Low target tissue uptake ({target_frac:.1%} to {target_tissue}). "
            f"Consider formulation optimization."
        )

    if formulation.route == AdministrationRoute.IV and target_tissue == "brain":
        warnings.append(
            "IV route has very low brain penetration due to BBB. "
            "Consider intrathecal or intranasal route."
        )

    formulation_warnings = formulation.validate()
    warnings.extend(formulation_warnings)

    return TissueDistribution(
        formulation_name=formulation.name,
        tissue_fractions=fractions,
        dominant_tissue=dominant,
        liver_fraction=liver_frac,
        extrahepatic_fraction=extrahepatic,
        therapeutic_index_estimate=ti,
        target_tissue=target_tissue,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Delivery probability chain
# ---------------------------------------------------------------------------

@dataclass
class DeliveryProbabilityChain:
    """Probability chain for effective delivery.

    Models the sequential probability of a payload reaching its target:
    P(effective) = P(circulation) * P(tissue_access) * P(cell_binding)
                   * P(internalization) * P(endosomal_escape) * P(action)

    Each probability is estimated from formulation parameters.
    """

    p_circulation: float
    p_tissue_access: float
    p_cell_binding: float
    p_internalization: float
    p_endosomal_escape: float
    p_intracellular_action: float

    @property
    def p_effective(self) -> float:
        """Overall probability of effective delivery."""
        return (
            self.p_circulation
            * self.p_tissue_access
            * self.p_cell_binding
            * self.p_internalization
            * self.p_endosomal_escape
            * self.p_intracellular_action
        )

    @property
    def bottleneck(self) -> tuple[str, float]:
        """Identify the weakest link in the delivery chain."""
        probs = {
            "circulation": self.p_circulation,
            "tissue_access": self.p_tissue_access,
            "cell_binding": self.p_cell_binding,
            "internalization": self.p_internalization,
            "endosomal_escape": self.p_endosomal_escape,
            "intracellular_action": self.p_intracellular_action,
        }
        weakest = min(probs, key=lambda k: probs[k])
        return weakest, probs[weakest]

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        bottleneck_name, bottleneck_val = self.bottleneck
        return {
            "p_circulation": self.p_circulation,
            "p_tissue_access": self.p_tissue_access,
            "p_cell_binding": self.p_cell_binding,
            "p_internalization": self.p_internalization,
            "p_endosomal_escape": self.p_endosomal_escape,
            "p_intracellular_action": self.p_intracellular_action,
            "p_effective": self.p_effective,
            "bottleneck": bottleneck_name,
            "bottleneck_value": bottleneck_val,
        }


def estimate_delivery_chain(
    formulation: LNPFormulation,
    target_tissue: str = "liver",
) -> DeliveryProbabilityChain:
    """Estimate the delivery probability chain for a formulation.

    Each step probability is estimated from formulation parameters
    using simplified models derived from published structure-activity
    relationships.

    Parameters
    ----------
    formulation : LNPFormulation
        Complete formulation specification.
    target_tissue : str
        Target tissue for delivery.

    Returns
    -------
    DeliveryProbabilityChain
        Estimated probabilities for each delivery step.
    """
    # P(circulation): depends on PEG density, size, route
    peg_pct = formulation.peg_lipid_mol_pct
    # PEG density between 1-3% is optimal for circulation
    p_circ = 0.5 + 0.3 * math.exp(-((peg_pct - 2.0) ** 2) / 2.0)

    # Size effect: 60-100 nm optimal for circulation
    size = formulation.target_size_nm
    p_circ *= 0.5 + 0.5 * math.exp(-((size - 80) ** 2) / 2000.0)

    # Route effect: IV gives full systemic exposure
    route_circ = {
        AdministrationRoute.IV: 1.0,
        AdministrationRoute.IM: 0.3,
        AdministrationRoute.SC: 0.2,
        AdministrationRoute.IN: 0.15,
        AdministrationRoute.IT: 0.1,
        AdministrationRoute.INTRATUMORAL: 0.2,
        AdministrationRoute.NEBULIZED: 0.15,
    }
    p_circ *= route_circ.get(formulation.route, 0.5)
    p_circ = min(p_circ, 0.95)

    # P(tissue_access): depends on target tissue and biodistribution
    biodist = predict_biodistribution(formulation, target_tissue)
    p_tissue = biodist.tissue_fractions.get(target_tissue, 0.01)
    # Scale: tissue fraction of 0.7 -> p_access of 0.8
    p_tissue = min(0.95, p_tissue * 1.2)

    # P(cell_binding): depends on targeting and tissue type
    p_binding = 0.4  # Base probability for untargeted LNPs
    if formulation.has_targeting:
        # Targeting peptides improve binding
        max_affinity = max(
            p.tissue_affinity.get(target_tissue, 0.0)
            for p in formulation.targeting_peptides
        ) if formulation.targeting_peptides else 0.0
        p_binding = min(0.9, 0.4 + max_affinity * 0.4)

    # Liver is special — ApoE-mediated uptake is very efficient
    if target_tissue == "liver" and formulation.route == AdministrationRoute.IV:
        p_binding = max(p_binding, 0.7)

    # P(internalization): depends on lipid composition
    p_internal = 0.5
    if formulation.helper_lipid == HelperLipidClass.DOPE:
        p_internal *= 1.3  # Fusogenic
    if formulation.ionizable_lipid_mol_pct > 40:
        p_internal *= 1.1  # More ionizable lipid aids uptake

    p_internal = min(0.9, p_internal)

    # P(endosomal_escape): the hardest step, typically 1-5%
    p_escape = 0.02  # Base for standard LNPs
    if formulation.helper_lipid == HelperLipidClass.DOPE:
        p_escape *= 3.0  # DOPE significantly aids escape
    if formulation.ionizable_lipid in {
        IonizableLipidClass.SM102,
        IonizableLipidClass.ALC0315,
    }:
        p_escape *= 2.0  # Optimized lipids
    if formulation.np_ratio > 6:
        p_escape *= 1.2  # Higher N:P can improve escape

    p_escape = min(0.20, p_escape)  # Even best LNPs: ~15-20%

    # P(intracellular_action): depends on payload type
    payload_action = {
        "mRNA": 0.7,  # mRNA just needs ribosome access
        "siRNA": 0.6,  # siRNA needs RISC loading
        "CRISPR-Cas9 RNP": 0.3,  # Needs nuclear import
        "plasmid": 0.1,  # Needs nuclear import + transcription
        "base_editor": 0.25,  # Needs nuclear import
        "prime_editor": 0.2,  # Needs nuclear import + complex mechanism
    }
    p_action = payload_action.get(formulation.payload_type, 0.5)

    return DeliveryProbabilityChain(
        p_circulation=round(p_circ, 4),
        p_tissue_access=round(p_tissue, 4),
        p_cell_binding=round(p_binding, 4),
        p_internalization=round(p_internal, 4),
        p_endosomal_escape=round(p_escape, 4),
        p_intracellular_action=round(p_action, 4),
    )
