"""Per-residue energy decomposition with binding-site focus.

Decomposes LJ energy to per-residue contributions, focused on
binding site residues. Identifies which residues contribute
favorable vs unfavorable interactions and ranks them by magnitude.

This is the bridge between raw physics (energy.py) and mechanistic
explanation (selectivity_map.py): it answers "which residues are
responsible for the energy landscape at a binding site?"
"""

from __future__ import annotations

from dataclasses import dataclass

import jax.numpy as jnp
import numpy as np

from kira.causality.binding_site import BindingSite
from kira.physics.core.energy import (
    compute_lj_energy_matrix,
    compute_per_atom_lj_energy,
    get_lj_params_arrays,
)


@dataclass
class ResidueEnergyContribution:
    """Energy contribution of a single residue at a binding site.

    Attributes
    ----------
    residue_index : int
        Global residue index in the parent structure.
    res_name : str
        Three-letter residue name.
    total_energy : float
        Sum of LJ contributions for all atoms in this residue (kcal/mol).
    n_atoms : int
        Number of atoms in this residue contributing to the energy.
    favorable : bool
        True if total_energy < 0 (stabilizing interaction).
    per_atom_energies : np.ndarray
        Per-atom energy breakdown within this residue.
    """

    residue_index: int
    res_name: str
    total_energy: float
    n_atoms: int
    favorable: bool
    per_atom_energies: np.ndarray


@dataclass
class BindingSiteDecomposition:
    """Complete energy decomposition for a binding site.

    Attributes
    ----------
    total_energy : float
        Total LJ energy across all binding-site residues (kcal/mol).
    residue_contributions : list[ResidueEnergyContribution]
        Per-residue breakdown, sorted by absolute energy (largest first).
    n_favorable : int
        Number of residues with favorable (negative) energy.
    n_unfavorable : int
        Number of residues with unfavorable (positive) energy.
    hotspot_indices : list[int]
        Residue indices of the top energy contributors (|E| > threshold).
    pocket_atom_indices : np.ndarray
        Global atom indices of all pocket atoms used in the computation.
    """

    total_energy: float
    residue_contributions: list[ResidueEnergyContribution]
    n_favorable: int
    n_unfavorable: int
    hotspot_indices: list[int]
    pocket_atom_indices: np.ndarray

    @property
    def n_residues(self) -> int:
        return len(self.residue_contributions)

    @property
    def mean_residue_energy(self) -> float:
        if not self.residue_contributions:
            return 0.0
        return self.total_energy / len(self.residue_contributions)

    def top_favorable(self, n: int = 5) -> list[ResidueEnergyContribution]:
        """Return the n most favorable (most negative energy) residues."""
        favorable = [r for r in self.residue_contributions if r.favorable]
        return sorted(favorable, key=lambda r: r.total_energy)[:n]

    def top_unfavorable(self, n: int = 5) -> list[ResidueEnergyContribution]:
        """Return the n most unfavorable (most positive energy) residues."""
        unfavorable = [r for r in self.residue_contributions if not r.favorable]
        return sorted(unfavorable, key=lambda r: r.total_energy, reverse=True)[:n]


def decompose_binding_site_energy(
    coords: np.ndarray,
    elements: np.ndarray,
    res_indices: np.ndarray,
    nonbonded_mask: np.ndarray,
    site: BindingSite,
    energy_cap: float = 1000.0,
    hotspot_threshold: float = 5.0,
) -> BindingSiteDecomposition:
    """Decompose LJ energy at a binding site into per-residue contributions.

    Computes the full LJ energy matrix for the structure, then extracts
    and aggregates energies for residues within the binding pocket.

    Parameters
    ----------
    coords : np.ndarray
        (N, 3) coordinates of ALL atoms in the structure.
    elements : np.ndarray
        (N,) element symbols for all atoms.
    res_indices : np.ndarray
        (N,) residue index for each atom.
    nonbonded_mask : np.ndarray
        (N, N) boolean mask -- True for non-bonded pairs.
    site : BindingSite
        Extracted binding site specifying which residues to decompose.
    energy_cap : float
        Cap per-pair LJ energy (prevents infinity from overlapping atoms).
    hotspot_threshold : float
        Residues with |energy| > this threshold are flagged as hotspots.

    Returns
    -------
    BindingSiteDecomposition
        Complete per-residue energy breakdown for the binding site.
    """
    pocket_atom_idx = site.atom_indices
    if len(pocket_atom_idx) == 0:
        return BindingSiteDecomposition(
            total_energy=0.0,
            residue_contributions=[],
            n_favorable=0,
            n_unfavorable=0,
            hotspot_indices=[],
            pocket_atom_indices=pocket_atom_idx,
        )

    # Get LJ parameters for all atoms
    sigma, epsilon = get_lj_params_arrays(elements)
    sigma_jnp = jnp.array(sigma)
    epsilon_jnp = jnp.array(epsilon)

    # Compute full pairwise distance matrix
    coords_jnp = jnp.array(coords)
    diff = coords_jnp[:, None, :] - coords_jnp[None, :, :]
    dist_matrix = jnp.sqrt(jnp.sum(diff**2, axis=-1) + 1e-10)

    # Compute full LJ energy matrix
    mask_jnp = jnp.array(nonbonded_mask.astype(np.float32))
    energy_matrix = compute_lj_energy_matrix(
        dist_matrix, sigma_jnp, sigma_jnp,
        epsilon_jnp, epsilon_jnp, mask_jnp, energy_cap,
    )

    # Compute per-atom energy for ALL atoms
    per_atom_energy = np.asarray(compute_per_atom_lj_energy(energy_matrix))

    # Aggregate to per-residue for pocket residues only
    residue_contributions: list[ResidueEnergyContribution] = []

    for list_pos, ridx in enumerate(site.residue_indices):
        # Find atoms belonging to this residue within the pocket
        atom_mask = res_indices[pocket_atom_idx] == ridx
        residue_atom_global = pocket_atom_idx[atom_mask]

        if len(residue_atom_global) == 0:
            continue

        atom_energies = per_atom_energy[residue_atom_global]
        total_e = float(np.sum(atom_energies))

        # Get residue name from the binding site
        res_name = site.res_names[list_pos] if list_pos < len(site.res_names) else "UNK"

        residue_contributions.append(ResidueEnergyContribution(
            residue_index=ridx,
            res_name=res_name,
            total_energy=total_e,
            n_atoms=len(residue_atom_global),
            favorable=total_e < 0,
            per_atom_energies=atom_energies,
        ))

    # Sort by absolute energy (largest contributors first)
    residue_contributions.sort(key=lambda r: abs(r.total_energy), reverse=True)

    total_energy = sum(r.total_energy for r in residue_contributions)
    n_favorable = sum(1 for r in residue_contributions if r.favorable)
    n_unfavorable = sum(1 for r in residue_contributions if not r.favorable)

    # Identify hotspots
    hotspot_indices = [
        r.residue_index for r in residue_contributions
        if abs(r.total_energy) > hotspot_threshold
    ]

    return BindingSiteDecomposition(
        total_energy=total_energy,
        residue_contributions=residue_contributions,
        n_favorable=n_favorable,
        n_unfavorable=n_unfavorable,
        hotspot_indices=hotspot_indices,
        pocket_atom_indices=pocket_atom_idx,
    )


def compare_residue_energies(
    decomp1: BindingSiteDecomposition,
    decomp2: BindingSiteDecomposition,
    alignment: list[tuple[int, int]] | None = None,
) -> list[dict]:
    """Compare per-residue energies between two binding-site decompositions.

    For selectivity analysis: given a compound docked into a parasite target
    and its human orthologue, identify which residue positions show the
    largest energy differences.

    Parameters
    ----------
    decomp1 : BindingSiteDecomposition
        Energy decomposition for structure 1 (e.g., parasite target).
    decomp2 : BindingSiteDecomposition
        Energy decomposition for structure 2 (e.g., human orthologue).
    alignment : list[tuple[int, int]] or None
        Mapping of (index_in_contributions_1, index_in_contributions_2).
        If None, uses sequential alignment by contribution list order.

    Returns
    -------
    list[dict]
        Per-position comparison with keys: 'pos', 'res1_name', 'res2_name',
        'energy1', 'energy2', 'delta_energy', 'same_residue'.
        Sorted by |delta_energy| descending.
    """
    r1 = decomp1.residue_contributions
    r2 = decomp2.residue_contributions

    if alignment is None:
        n = min(len(r1), len(r2))
        alignment = [(i, i) for i in range(n)]

    comparisons: list[dict] = []
    for pos, (i, j) in enumerate(alignment):
        if i >= len(r1) or j >= len(r2):
            continue

        e1 = r1[i].total_energy
        e2 = r2[j].total_energy
        n1 = r1[i].res_name
        n2 = r2[j].res_name

        comparisons.append({
            "pos": pos,
            "res1_name": n1,
            "res2_name": n2,
            "res1_index": r1[i].residue_index,
            "res2_index": r2[j].residue_index,
            "energy1": e1,
            "energy2": e2,
            "delta_energy": e1 - e2,
            "same_residue": n1 == n2,
        })

    comparisons.sort(key=lambda c: abs(c["delta_energy"]), reverse=True)
    return comparisons
