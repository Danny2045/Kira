"""Selectivity attribution maps.

For two homologous pockets and their per-residue LJ decompositions,
compute residue-level selectivity attribution: which specific residue
differences are associated with the residue-energy differential.

This module combines binding-site extraction, per-residue energy
decomposition, and pocket comparison to produce a mechanistic
hypothesis for selectivity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kira.causality.binding_site import (
    BindingSite,
    PocketComparison,
    compare_binding_sites,
)
from kira.causality.energy_decomp import (
    BindingSiteDecomposition,
    compare_residue_energies,
)


@dataclass
class ResidueAttribution:
    """Selectivity attribution for a single aligned position.

    Attributes
    ----------
    position : int
        Aligned position index.
    res1_name : str
        Residue name in structure 1 (e.g., parasite target).
    res2_name : str
        Residue name in structure 2 (e.g., human orthologue).
    is_divergent : bool
        True if the residue identity differs at this position.
    energy1 : float
        Per-residue LJ energy in structure 1 (kcal/mol).
    energy2 : float
        Per-residue LJ energy in structure 2 (kcal/mol).
    delta_energy : float
        energy1 - energy2. Negative means structure 1 is more favorable.
    attribution_score : float
        Normalized contribution to overall selectivity (0-1 scale).
    """

    position: int
    res1_name: str
    res2_name: str
    is_divergent: bool
    energy1: float
    energy2: float
    delta_energy: float
    attribution_score: float


@dataclass
class SelectivityMap:
    """Complete residue-level selectivity attribution between two homologous pockets.

    Attributes
    ----------
    structure1_name : str
        Name/identifier for structure 1 (e.g., "SmTGR").
    structure2_name : str
        Name/identifier for structure 2 (e.g., "HsTrxR1").
    total_delta_energy : float
        Sum of all per-residue delta energies (kcal/mol).
        Negative = structure 1 is more favorable overall.
    pocket_comparison : PocketComparison
        Residue-level pocket alignment and identity comparison.
    attributions : list[ResidueAttribution]
        Per-position residue-level attributions, sorted by |delta_energy|.
    n_positions : int
        Total number of aligned positions.
    n_divergent : int
        Number of positions where residue identity differs.
    selectivity_direction : str
        "target1" if compound prefers structure 1, "target2" otherwise,
        "neutral" if |total_delta_energy| < threshold.
    """

    structure1_name: str
    structure2_name: str
    total_delta_energy: float
    pocket_comparison: PocketComparison
    attributions: list[ResidueAttribution]
    n_positions: int
    n_divergent: int
    selectivity_direction: str

    @property
    def divergent_attributions(self) -> list[ResidueAttribution]:
        """Return only attributions at divergent (non-conserved) positions."""
        return [a for a in self.attributions if a.is_divergent]

    @property
    def conserved_attributions(self) -> list[ResidueAttribution]:
        """Return only attributions at conserved (identical) positions."""
        return [a for a in self.attributions if not a.is_divergent]

    @property
    def divergent_energy_fraction(self) -> float:
        """Fraction of total |delta_energy| attributable to divergent residues."""
        total_abs = sum(abs(a.delta_energy) for a in self.attributions)
        if total_abs == 0:
            return 0.0
        div_abs = sum(abs(a.delta_energy) for a in self.attributions if a.is_divergent)
        return div_abs / total_abs

    def top_drivers(self, n: int = 5) -> list[ResidueAttribution]:
        """Return the n positions contributing most to selectivity."""
        return self.attributions[:n]

    def summary(self) -> str:
        """One-paragraph heuristic summary of selectivity attribution."""
        if not self.attributions:
            return "No aligned positions available for selectivity analysis."

        direction = "favors" if self.total_delta_energy < 0 else "disfavors"
        top = self.attributions[:3]
        drivers = ", ".join(
            f"{a.res1_name}{a.position}->{a.res2_name} ({a.delta_energy:+.1f} kcal/mol)"
            for a in top
        )
        return (
            f"Net energy difference {direction} {self.structure1_name} by "
            f"{abs(self.total_delta_energy):.1f} kcal/mol. "
            f"{self.n_divergent}/{self.n_positions} pocket positions diverge, "
            f"accounting for {self.divergent_energy_fraction:.0%} of the "
            f"selectivity signal. Top drivers: {drivers}."
        )


def build_selectivity_map(
    site1: BindingSite,
    site2: BindingSite,
    decomp1: BindingSiteDecomposition,
    decomp2: BindingSiteDecomposition,
    structure1_name: str = "target",
    structure2_name: str = "ortholog",
    alignment: list[tuple[int, int]] | None = None,
    neutrality_threshold: float = 1.0,
) -> SelectivityMap:
    """Build a residue-level selectivity attribution map between two homologous pockets.

    Combines pocket comparison (which residues differ) with energy
    decomposition (how much each residue contributes) to highlight
    residue-level differences associated with the residue-energy differential.

    Parameters
    ----------
    site1 : BindingSite
        Binding pocket from structure 1 (e.g., parasite target).
    site2 : BindingSite
        Binding pocket from structure 2 (e.g., human orthologue).
    decomp1 : BindingSiteDecomposition
        Energy decomposition for structure 1.
    decomp2 : BindingSiteDecomposition
        Energy decomposition for structure 2.
    structure1_name : str
        Label for structure 1.
    structure2_name : str
        Label for structure 2.
    alignment : list[tuple[int, int]] or None
        Residue-level alignment between pockets. If None, sequential.
    neutrality_threshold : float
        |total_delta_energy| below this is classified as "neutral".

    Returns
    -------
    SelectivityMap
        Complete selectivity attribution.
    """
    # Compare pocket residue identities
    pocket_cmp = compare_binding_sites(site1, site2, alignment)

    # Compare per-residue energies
    energy_cmp = compare_residue_energies(decomp1, decomp2, alignment)

    # Build per-position attributions
    # We need to merge pocket comparison (identity) with energy comparison
    total_abs_delta = sum(abs(c["delta_energy"]) for c in energy_cmp)

    attributions: list[ResidueAttribution] = []
    for c in energy_cmp:
        pos = c["pos"]
        is_div = (
            pos in pocket_cmp.divergent_positions
            if pos < pocket_cmp.n_aligned_positions
            else True
        )

        attr_score = abs(c["delta_energy"]) / total_abs_delta if total_abs_delta > 0 else 0.0

        attributions.append(ResidueAttribution(
            position=pos,
            res1_name=c["res1_name"],
            res2_name=c["res2_name"],
            is_divergent=is_div,
            energy1=c["energy1"],
            energy2=c["energy2"],
            delta_energy=c["delta_energy"],
            attribution_score=attr_score,
        ))

    # Already sorted by |delta_energy| from compare_residue_energies
    total_delta = sum(a.delta_energy for a in attributions)
    n_divergent = sum(1 for a in attributions if a.is_divergent)

    if abs(total_delta) < neutrality_threshold:
        direction = "neutral"
    elif total_delta < 0:
        direction = "target1"
    else:
        direction = "target2"

    return SelectivityMap(
        structure1_name=structure1_name,
        structure2_name=structure2_name,
        total_delta_energy=total_delta,
        pocket_comparison=pocket_cmp,
        attributions=attributions,
        n_positions=len(attributions),
        n_divergent=n_divergent,
        selectivity_direction=direction,
    )


def compute_selectivity_score(sel_map: SelectivityMap) -> float:
    """Compute a scalar heuristic selectivity score from a SelectivityMap.

    Higher scores indicate stronger, more interpretable residue-level signal.
    The score considers both the magnitude of energy differences and
    the fraction attributable to divergent residues (which are more
    heuristically interpretable).

    Parameters
    ----------
    sel_map : SelectivityMap
        Completed selectivity map.

    Returns
    -------
    float
        Score in [0, 1]. Higher = stronger, more interpretable heuristic signal.
    """
    if sel_map.n_positions == 0:
        return 0.0

    # Component 1: Energy magnitude (sigmoid-like scaling)
    # 10 kcal/mol difference maps to ~0.9
    energy_mag = abs(sel_map.total_delta_energy)
    energy_score = 1.0 - np.exp(-energy_mag / 5.0)

    # Component 2: Fraction of selectivity from divergent residues
    # Higher = more mechanistically interpretable
    div_fraction = sel_map.divergent_energy_fraction

    # Component 3: Concentration — are a few residues responsible,
    # or is it spread across many? Concentrated is more actionable.
    if len(sel_map.attributions) >= 3:
        top3_frac = sum(a.attribution_score for a in sel_map.attributions[:3])
    else:
        top3_frac = 1.0

    # Weighted combination
    score = 0.5 * energy_score + 0.3 * div_fraction + 0.2 * top3_frac
    return float(np.clip(score, 0.0, 1.0))
