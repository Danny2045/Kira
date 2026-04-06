"""JAX-accelerated geometry computations.

All functions are JIT-compiled for maximum throughput.
Distance matrix is computed once and reused across all checks.
Dihedral and angle computations are fully vectorized — no Python loops
over atom quads.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np


@jax.jit
def compute_distance_matrix(coords: jnp.ndarray) -> jnp.ndarray:
    """Compute pairwise Euclidean distance matrix.

    Parameters
    ----------
    coords : jnp.ndarray
        (N, 3) coordinate array.

    Returns
    -------
    jnp.ndarray
        (N, N) symmetric distance matrix. Diagonal is ~0 (epsilon for stability).
    """
    diff = coords[:, None, :] - coords[None, :, :]
    return jnp.sqrt(jnp.sum(diff**2, axis=-1) + 1e-10)


@jax.jit
def compute_dihedral_angles(p0: jnp.ndarray, p1: jnp.ndarray,
                            p2: jnp.ndarray, p3: jnp.ndarray) -> jnp.ndarray:
    """Compute dihedral angles for batches of atom quads.

    Uses the atan2 formulation for signed angles in [-π, π].

    Parameters
    ----------
    p0, p1, p2, p3 : jnp.ndarray
        Each (M, 3) — M sets of four atoms defining dihedrals.

    Returns
    -------
    jnp.ndarray
        (M,) dihedral angles in radians, range [-π, π].
    """
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2

    # Normal vectors to planes
    n1 = jnp.cross(b1, b2)
    n2 = jnp.cross(b2, b3)

    # Normalize b2 for the projection
    b2_norm = b2 / (jnp.linalg.norm(b2, axis=-1, keepdims=True) + 1e-10)

    # m1 is n1 × b2_hat (lies in the plane of n1 and perpendicular to b2)
    m1 = jnp.cross(n1, b2_norm)

    # x = n1 · n2, y = m1 · n2
    x = jnp.sum(n1 * n2, axis=-1)
    y = jnp.sum(m1 * n2, axis=-1)

    return jnp.arctan2(y, x)


@jax.jit
def compute_bond_angles(p0: jnp.ndarray, p1: jnp.ndarray,
                        p2: jnp.ndarray) -> jnp.ndarray:
    """Compute bond angles for batches of atom triples.

    Parameters
    ----------
    p0, p1, p2 : jnp.ndarray
        Each (M, 3) — M sets of three atoms. p1 is the central atom.

    Returns
    -------
    jnp.ndarray
        (M,) bond angles in radians, range [0, π].
    """
    v1 = p0 - p1
    v2 = p2 - p1

    cos_angle = jnp.sum(v1 * v2, axis=-1) / (
        jnp.linalg.norm(v1, axis=-1) * jnp.linalg.norm(v2, axis=-1) + 1e-10
    )
    # Clamp to [-1, 1] for numerical stability
    cos_angle = jnp.clip(cos_angle, -1.0, 1.0)

    return jnp.arccos(cos_angle)


@jax.jit
def compute_distances(p0: jnp.ndarray, p1: jnp.ndarray) -> jnp.ndarray:
    """Compute pairwise distances between corresponding atoms.

    Parameters
    ----------
    p0, p1 : jnp.ndarray
        Each (M, 3) — M pairs of atoms.

    Returns
    -------
    jnp.ndarray
        (M,) distances.
    """
    diff = p1 - p0
    return jnp.sqrt(jnp.sum(diff**2, axis=-1) + 1e-10)


def extract_backbone_dihedrals(coords: np.ndarray, atom_names: np.ndarray,
                               res_indices: np.ndarray, is_protein: np.ndarray,
                               chain_ids: np.ndarray) -> dict[str, dict]:
    """Extract phi, psi, and omega dihedral angles from backbone atoms.

    Parameters
    ----------
    coords : np.ndarray
        (N, 3) all-atom coordinates.
    atom_names : np.ndarray
        (N,) atom names.
    res_indices : np.ndarray
        (N,) residue indices.
    is_protein : np.ndarray
        (N,) protein mask.
    chain_ids : np.ndarray
        (N,) chain IDs.

    Returns
    -------
    dict
        Keys: 'phi', 'psi', 'omega'. Each maps to:
        - 'angles': jnp.ndarray of angles in degrees
        - 'res_indices': list of residue indices these angles correspond to
    """
    # Build per-residue backbone atom coordinate lookup
    # residue_idx → {atom_name: coord_index}
    res_backbone: dict[int, dict[str, int]] = {}

    protein_mask = is_protein.astype(bool)
    for i in range(len(coords)):
        if not protein_mask[i]:
            continue
        aname = str(atom_names[i]).strip()
        if aname not in ("N", "CA", "C"):
            continue
        ridx = int(res_indices[i])
        if ridx not in res_backbone:
            res_backbone[ridx] = {}
        res_backbone[ridx][aname] = i

    # Get ordered residue indices per chain
    chain_residues: dict[str, list[int]] = {}
    for i in range(len(coords)):
        if not protein_mask[i]:
            continue
        cid = str(chain_ids[i]).strip()
        ridx = int(res_indices[i])
        if cid not in chain_residues:
            chain_residues[cid] = []
        if ridx not in chain_residues[cid]:
            chain_residues[cid].append(ridx)

    # Sort residue indices within each chain
    for cid in chain_residues:
        chain_residues[cid] = sorted(set(chain_residues[cid]))

    # Collect dihedral quads
    phi_quads = []  # C(i-1), N(i), CA(i), C(i)
    phi_res = []
    psi_quads = []  # N(i), CA(i), C(i), N(i+1)
    psi_res = []
    omega_quads = []  # CA(i-1), C(i-1), N(i), CA(i)
    omega_res = []

    for cid, res_list in chain_residues.items():
        for k in range(len(res_list)):
            ridx = res_list[k]
            bb = res_backbone.get(ridx, {})

            if k > 0:
                prev_ridx = res_list[k - 1]
                prev_bb = res_backbone.get(prev_ridx, {})

                # Phi: C(i-1), N(i), CA(i), C(i)
                if "C" in prev_bb and "N" in bb and "CA" in bb and "C" in bb:
                    phi_quads.append((prev_bb["C"], bb["N"], bb["CA"], bb["C"]))
                    phi_res.append(ridx)

                # Omega: CA(i-1), C(i-1), N(i), CA(i)
                if "CA" in prev_bb and "C" in prev_bb and "N" in bb and "CA" in bb:
                    omega_quads.append((prev_bb["CA"], prev_bb["C"], bb["N"], bb["CA"]))
                    omega_res.append(ridx)

            if k < len(res_list) - 1:
                next_ridx = res_list[k + 1]
                next_bb = res_backbone.get(next_ridx, {})

                # Psi: N(i), CA(i), C(i), N(i+1)
                if "N" in bb and "CA" in bb and "C" in bb and "N" in next_bb:
                    psi_quads.append((bb["N"], bb["CA"], bb["C"], next_bb["N"]))
                    psi_res.append(ridx)

    results = {}

    for name, quads, res_list in [
        ("phi", phi_quads, phi_res),
        ("psi", psi_quads, psi_res),
        ("omega", omega_quads, omega_res),
    ]:
        if not quads:
            results[name] = {"angles": jnp.array([]), "res_indices": []}
            continue

        idx = np.array(quads, dtype=np.int32)
        p0 = jnp.array(coords[idx[:, 0]])
        p1 = jnp.array(coords[idx[:, 1]])
        p2 = jnp.array(coords[idx[:, 2]])
        p3 = jnp.array(coords[idx[:, 3]])

        angles_rad = compute_dihedral_angles(p0, p1, p2, p3)
        angles_deg = jnp.degrees(angles_rad)

        results[name] = {"angles": angles_deg, "res_indices": res_list}

    return results
