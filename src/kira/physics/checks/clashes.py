"""Steric clash detection.

Identifies non-bonded atom pairs closer than the sum of their
van der Waals radii (minus a tolerance). This is the most common
physics hallucination in AI-generated structures.
"""

from __future__ import annotations

from dataclasses import dataclass

import jax.numpy as jnp
import numpy as np

from kira.physics.config import ClashConfig


@dataclass
class ClashResult:
    """Result of steric clash analysis."""

    n_clashes: int
    n_severe_clashes: int
    clashscore: float  # Clashes per 1000 atoms
    worst_overlap: float  # Angstroms, most severe overlap
    clashing_pairs: list[tuple[int, int, float]]  # (atom_i, atom_j, overlap)
    per_residue_clashes: np.ndarray  # (R,) count of clashes per residue
    subscore: float  # 0-1 normalized score (1 = no clashes)


def check_clashes(
    dist_matrix: jnp.ndarray,
    elements: np.ndarray,
    nonbonded_mask: jnp.ndarray,
    res_indices: np.ndarray,
    n_residues: int,
    config: ClashConfig | None = None,
) -> ClashResult:
    """Check for steric clashes in a structure.

    A clash occurs when the distance between two non-bonded atoms
    is less than the sum of their vdW radii minus a tolerance.

    Parameters
    ----------
    dist_matrix : jnp.ndarray
        (N, N) pairwise distance matrix.
    elements : np.ndarray
        (N,) element symbols.
    nonbonded_mask : jnp.ndarray
        (N, N) bool — True for non-bonded pairs.
    res_indices : np.ndarray
        (N,) residue index per atom.
    n_residues : int
        Total number of residues.
    config : ClashConfig or None
        Clash detection parameters. Uses defaults if None.

    Returns
    -------
    ClashResult
        Complete clash analysis results.
    """
    if config is None:
        config = ClashConfig()

    n_atoms = len(elements)
    radii = get_vdw_radii_array_from_elements(elements)

    # Compute vdW sum matrix
    radii_jnp = jnp.array(radii)
    vdw_sum = radii_jnp[:, None] + radii_jnp[None, :]

    # Clash threshold: vdW sum - tolerance
    threshold = vdw_sum - config.vdw_tolerance

    # Overlap: how much closer than threshold
    overlap = threshold - dist_matrix

    # Only consider non-bonded pairs in upper triangle
    upper_mask = jnp.triu(jnp.ones((n_atoms, n_atoms), dtype=bool), k=1)
    check_mask = nonbonded_mask & upper_mask

    # Clashing pairs: overlap > 0 and in the check mask
    is_clash = (overlap > 0) & check_mask
    is_severe = (overlap > config.severe_clash_threshold) & check_mask

    n_clashes = int(jnp.sum(is_clash))
    n_severe = int(jnp.sum(is_severe))

    # Extract clashing pairs for reporting
    clash_indices = np.argwhere(np.asarray(is_clash))
    overlaps = np.asarray(overlap)

    clashing_pairs = []
    for idx in clash_indices:
        i, j = int(idx[0]), int(idx[1])
        ov = float(overlaps[i, j])
        clashing_pairs.append((i, j, ov))

    # Sort by overlap severity (worst first)
    clashing_pairs.sort(key=lambda x: x[2], reverse=True)

    worst_overlap = clashing_pairs[0][2] if clashing_pairs else 0.0

    # Per-residue clash count
    per_res = np.zeros(n_residues, dtype=np.int32)
    for i, j, _ in clashing_pairs:
        per_res[res_indices[i]] += 1
        per_res[res_indices[j]] += 1

    # Clashscore: clashes per 1000 atoms (standard metric)
    clashscore = (n_clashes / max(n_atoms, 1)) * 1000.0

    # Subscore: 1.0 = perfect, decays with clashscore
    # A clashscore of 0 = 1.0, clashscore of 40+ = ~0.0
    # Using sigmoid-like decay
    subscore = max(0.0, 1.0 - clashscore / 40.0)
    subscore = min(1.0, subscore)

    return ClashResult(
        n_clashes=n_clashes,
        n_severe_clashes=n_severe,
        clashscore=clashscore,
        worst_overlap=worst_overlap,
        clashing_pairs=clashing_pairs[:100],  # Cap for reporting
        per_residue_clashes=per_res,
        subscore=subscore,
    )


def get_vdw_radii_array_from_elements(elements: np.ndarray) -> np.ndarray:
    """Get vdW radii array from element symbols.

    Parameters
    ----------
    elements : np.ndarray
        (N,) element symbols.

    Returns
    -------
    np.ndarray
        (N,) float32 vdW radii.
    """
    from kira.physics.core.topology import get_vdw_radius
    return np.array([get_vdw_radius(str(e).strip()) for e in elements], dtype=np.float32)
