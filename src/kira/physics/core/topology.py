"""Topology builder — bond graph, atom typing, and mask generation.

Infers covalent bonds from atomic distances and residue templates,
assigns atom types for force field parameters, and generates the
bonded-pair mask used by LJ and clash computations.
"""

from __future__ import annotations

import numpy as np

from kira.physics.core.parser import Structure

# Maximum covalent bond distances by element pair (Angstroms)
# Generous thresholds to catch all real bonds without false positives
COVALENT_BOND_CUTOFFS: dict[tuple[str, str], float] = {
    ("C", "C"): 1.9,
    ("C", "N"): 1.8,
    ("C", "O"): 1.8,
    ("C", "S"): 2.1,
    ("C", "H"): 1.3,
    ("N", "H"): 1.3,
    ("N", "N"): 1.7,
    ("O", "H"): 1.2,
    ("S", "S"): 2.6,  # Disulfide bonds
    ("S", "H"): 1.5,
    ("P", "O"): 1.9,
    ("C", "F"): 1.6,
    ("C", "CL"): 2.0,
    ("C", "BR"): 2.1,
    ("C", "SE"): 2.2,
    ("N", "O"): 1.7,
    ("O", "P"): 1.9,
    ("S", "O"): 1.8,
}

# Minimum bond distance — anything closer is likely an error, not a bond
MIN_BOND_DISTANCE = 0.4  # Angstroms

# Standard backbone connectivity: N → CA → C → O, and C → next N
BACKBONE_BONDS = [("N", "CA"), ("CA", "C"), ("C", "O")]

# Simplified sidechain connectivity templates for standard residues
# Maps residue name → list of (atom1, atom2) bonds
# Only heavy atoms; hydrogens handled by distance-based detection
SIDECHAIN_BONDS: dict[str, list[tuple[str, str]]] = {
    "ALA": [("CA", "CB")],
    "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
            ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")],
    "ASN": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
    "ASP": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
    "CYS": [("CA", "CB"), ("CB", "SG")],
    "GLN": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")],
    "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
    "GLY": [],
    "HIS": [("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"),
            ("ND1", "CE1"), ("CD2", "NE2"), ("CE1", "NE2")],
    "ILE": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")],
    "LEU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")],
    "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "MET": [("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "PHE": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "PRO": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
    "SER": [("CA", "CB"), ("CB", "OG")],
    "THR": [("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")],
    "TRP": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "NE1"), ("CD2", "CE2"), ("CD2", "CE3"), ("NE1", "CE2"),
            ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ2", "CH2"), ("CZ3", "CH2")],
    "TYR": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH")],
    "VAL": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")],
    "MSE": [("CA", "CB"), ("CB", "CG"), ("CG", "SE"), ("SE", "CE")],
}


# Van der Waals radii by element (Angstroms) — from Bondi (1964) / AMBER
VDW_RADII: dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "SE": 1.90, "FE": 1.95, "ZN": 1.39, "MG": 1.73, "CA": 1.95,
    "MN": 1.95, "CO": 1.95, "CU": 1.40, "NI": 1.63,
}

DEFAULT_VDW_RADIUS = 1.70  # Fallback


def get_vdw_radius(element: str) -> float:
    """Get van der Waals radius for an element."""
    return VDW_RADII.get(element.upper(), DEFAULT_VDW_RADIUS)


@np.errstate(invalid="ignore")
def _compute_distance_matrix_np(coords: np.ndarray) -> np.ndarray:
    """Compute pairwise distance matrix using numpy (for topology building).

    Parameters
    ----------
    coords : np.ndarray
        (N, 3) coordinate array.

    Returns
    -------
    np.ndarray
        (N, N) distance matrix.
    """
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    return np.sqrt(np.sum(diff**2, axis=-1))


def infer_bonds_from_topology(structure: Structure) -> list[tuple[int, int]]:
    """Infer covalent bonds using residue templates + distance fallback.

    For protein residues: use backbone + sidechain templates.
    For non-protein atoms (ligands, waters): use distance-based detection.
    Peptide bonds (C → next N) inferred from sequential residues.

    Parameters
    ----------
    structure : Structure
        Parsed structure.

    Returns
    -------
    list[tuple[int, int]]
        List of bonded atom index pairs (i, j) where i < j.
    """
    bonds: set[tuple[int, int]] = set()

    # Build atom lookup: for each residue, map atom_name → global index
    residue_atom_indices: dict[tuple, dict[str, int]] = {}
    for i, atom in enumerate(structure.atoms):
        rid = atom.residue_id
        if rid not in residue_atom_indices:
            residue_atom_indices[rid] = {}
        residue_atom_indices[rid][atom.name] = i

    # 1. Template-based bonds for protein residues
    for rid, atom_map in residue_atom_indices.items():
        residue = structure.residues.get(rid)
        if residue is None or not residue.is_protein:
            continue

        res_name = residue.res_name

        # Backbone bonds
        for a1, a2 in BACKBONE_BONDS:
            if a1 in atom_map and a2 in atom_map:
                bond = tuple(sorted((atom_map[a1], atom_map[a2])))
                bonds.add(bond)

        # Sidechain bonds
        sc_bonds = SIDECHAIN_BONDS.get(res_name, [])
        for a1, a2 in sc_bonds:
            if a1 in atom_map and a2 in atom_map:
                bond = tuple(sorted((atom_map[a1], atom_map[a2])))
                bonds.add(bond)

    # 2. Peptide bonds between sequential residues in the same chain
    for chain in structure.chains.values():
        protein_residues = [r for r in chain.residues if r.is_protein]
        for k in range(len(protein_residues) - 1):
            r1 = protein_residues[k]
            r2 = protein_residues[k + 1]
            rid1 = r1.residue_id
            rid2 = r2.residue_id

            atoms1 = residue_atom_indices.get(rid1, {})
            atoms2 = residue_atom_indices.get(rid2, {})

            if "C" in atoms1 and "N" in atoms2:
                idx_c = atoms1["C"]
                idx_n = atoms2["N"]
                # Verify distance is reasonable for a peptide bond (< 1.8 Å)
                dist = np.linalg.norm(
                    structure.coords[idx_c] - structure.coords[idx_n]
                )
                if dist < 1.8:
                    bonds.add(tuple(sorted((idx_c, idx_n))))

    # 3. Distance-based bond detection for non-protein atoms (ligands, etc.)
    non_protein_indices = [
        i for i, atom in enumerate(structure.atoms)
        if not atom.is_protein
    ]

    if non_protein_indices:
        # Also check bonds between non-protein atoms and nearby protein atoms
        all_indices_to_check = set(non_protein_indices)
        # Include protein atoms near non-protein atoms for cross-bonds
        if len(non_protein_indices) < structure.n_atoms:
            for np_idx in non_protein_indices:
                dists = np.linalg.norm(
                    structure.coords - structure.coords[np_idx], axis=1
                )
                nearby = np.where(dists < 2.5)[0]
                all_indices_to_check.update(nearby.tolist())

        check_list = sorted(all_indices_to_check)
        if len(check_list) > 1:
            sub_coords = structure.coords[check_list]
            sub_dists = _compute_distance_matrix_np(sub_coords)

            for ii in range(len(check_list)):
                for jj in range(ii + 1, len(check_list)):
                    real_i = check_list[ii]
                    real_j = check_list[jj]
                    d = sub_dists[ii, jj]

                    if d < MIN_BOND_DISTANCE:
                        continue

                    e1 = structure.elements[real_i].upper()
                    e2 = structure.elements[real_j].upper()

                    pair = tuple(sorted((e1, e2)))
                    cutoff = COVALENT_BOND_CUTOFFS.get(pair, 1.9)

                    if d <= cutoff:
                        bonds.add(tuple(sorted((real_i, real_j))))

    return sorted(bonds)


def build_bonded_mask(n_atoms: int, bonds: list[tuple[int, int]],
                      exclude_neighbors: int = 3) -> np.ndarray:
    """Build a boolean mask for non-bonded atom pairs.

    For LJ and clash computations, we need to exclude bonded pairs
    and their close neighbors (1-2 and 1-3 interactions).

    Parameters
    ----------
    n_atoms : int
        Total number of atoms.
    bonds : list[tuple[int, int]]
        List of bonded atom index pairs.
    exclude_neighbors : int
        Exclude pairs up to this many bonds apart. 3 = exclude 1-2 and 1-3.

    Returns
    -------
    np.ndarray
        (N, N) bool array. True = non-bonded (include in LJ/clash).
        Diagonal is False.
    """
    # Build adjacency list
    adjacency: dict[int, set[int]] = {i: set() for i in range(n_atoms)}
    for i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)

    # BFS to find all pairs within `exclude_neighbors` bonds
    exclude = set()
    for start in range(n_atoms):
        # BFS from start
        visited = {start: 0}
        queue = [start]
        head = 0
        while head < len(queue):
            current = queue[head]
            head += 1
            current_dist = visited[current]
            if current_dist >= exclude_neighbors:
                continue
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited[neighbor] = current_dist + 1
                    queue.append(neighbor)

        for node, dist in visited.items():
            if node != start and dist <= exclude_neighbors:
                pair = (min(start, node), max(start, node))
                exclude.add(pair)

    # Build mask
    mask = np.ones((n_atoms, n_atoms), dtype=bool)
    np.fill_diagonal(mask, False)  # No self-interaction

    for i, j in exclude:
        mask[i, j] = False
        mask[j, i] = False

    return mask


def build_1_4_mask(n_atoms: int, bonds: list[tuple[int, int]]) -> np.ndarray:
    """Build mask for 1-4 interactions (exactly 3 bonds apart).

    1-4 pairs get scaled LJ/electrostatic interactions in most force fields.

    Parameters
    ----------
    n_atoms : int
        Total number of atoms.
    bonds : list[tuple[int, int]]
        List of bonded atom index pairs.

    Returns
    -------
    np.ndarray
        (N, N) bool array. True = 1-4 pair.
    """
    adjacency: dict[int, set[int]] = {i: set() for i in range(n_atoms)}
    for i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)

    mask = np.zeros((n_atoms, n_atoms), dtype=bool)

    for start in range(n_atoms):
        # Find atoms exactly 3 bonds away
        visited = {start: 0}
        queue = [start]
        head = 0
        while head < len(queue):
            current = queue[head]
            head += 1
            current_dist = visited[current]
            if current_dist >= 3:
                continue
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited[neighbor] = current_dist + 1
                    queue.append(neighbor)

        for node, dist in visited.items():
            if dist == 3 and node > start:
                mask[start, node] = True
                mask[node, start] = True

    return mask


def get_vdw_radii_array(structure: Structure) -> np.ndarray:
    """Get van der Waals radii for all atoms in a structure.

    Parameters
    ----------
    structure : Structure
        Parsed structure.

    Returns
    -------
    np.ndarray
        (N,) float32 array of vdW radii.
    """
    radii = np.array([
        get_vdw_radius(elem) for elem in structure.elements
    ], dtype=np.float32)
    return radii
