"""PDB file parser.

Reads PDB and minimal mmCIF files into a Structure dataclass
optimized for downstream JAX computation.

The parser extracts ATOM/HETATM records and builds a flat array
representation suitable for vectorized distance/energy calculations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Standard amino acid 3-letter codes
STANDARD_RESIDUES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    # Common modifications treated as standard
    "MSE",  # selenomethionine → methionine
})

# Map modified residues to their parent
MODIFIED_RESIDUE_MAP = {
    "MSE": "MET",
}

# Backbone atom names (for dihedral computation)
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})

# Element to atomic number (for common biologically relevant elements)
ELEMENT_ATOMIC_NUM = {
    "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
    "CL": 17, "SE": 34, "BR": 35, "I": 53, "FE": 26, "ZN": 30,
    "MG": 12, "CA": 20, "MN": 25, "CO": 27, "CU": 29, "NI": 28,
}


@dataclass
class Atom:
    """Single atom parsed from a PDB record."""

    serial: int
    name: str  # Atom name (e.g., "CA", "CB", "OG")
    alt_loc: str  # Alternate location indicator
    res_name: str  # Residue name (e.g., "ALA", "HOH")
    chain_id: str  # Chain identifier
    res_seq: int  # Residue sequence number
    insertion_code: str  # Insertion code
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    element: str  # Element symbol
    is_hetatm: bool  # True if from HETATM record

    @property
    def is_protein(self) -> bool:
        """Whether this atom belongs to a standard amino acid."""
        return self.res_name in STANDARD_RESIDUES

    @property
    def is_backbone(self) -> bool:
        """Whether this is a backbone atom (N, CA, C, O)."""
        return self.name in BACKBONE_ATOMS

    @property
    def is_hydrogen(self) -> bool:
        """Whether this is a hydrogen atom."""
        return self.element == "H"

    @property
    def residue_id(self) -> tuple[str, int, str]:
        """Unique residue identifier: (chain_id, res_seq, insertion_code)."""
        return (self.chain_id, self.res_seq, self.insertion_code)

    @property
    def canonical_res_name(self) -> str:
        """Residue name mapped through modifications (e.g., MSE → MET)."""
        return MODIFIED_RESIDUE_MAP.get(self.res_name, self.res_name)


@dataclass
class Residue:
    """A residue (amino acid, ligand, water, etc.) grouping its atoms."""

    chain_id: str
    res_seq: int
    insertion_code: str
    res_name: str
    atoms: dict[str, Atom] = field(default_factory=dict)  # atom_name → Atom

    @property
    def residue_id(self) -> tuple[str, int, str]:
        return (self.chain_id, self.res_seq, self.insertion_code)

    @property
    def is_protein(self) -> bool:
        return self.res_name in STANDARD_RESIDUES

    def get_atom(self, name: str) -> Atom | None:
        """Get atom by name, returning None if not found."""
        return self.atoms.get(name)

    def get_coord(self, name: str) -> np.ndarray | None:
        """Get (3,) coordinate array for an atom by name."""
        atom = self.atoms.get(name)
        if atom is None:
            return None
        return np.array([atom.x, atom.y, atom.z], dtype=np.float32)


@dataclass
class Chain:
    """A polymer or entity chain."""

    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    @property
    def is_protein(self) -> bool:
        """Whether this chain is primarily protein (>50% standard residues)."""
        if not self.residues:
            return False
        n_protein = sum(1 for r in self.residues if r.is_protein)
        return n_protein / len(self.residues) > 0.5

    @property
    def sequence(self) -> str:
        """One-letter amino acid sequence for protein residues."""
        three_to_one = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "MSE": "M",
        }
        return "".join(
            three_to_one.get(r.res_name, "X")
            for r in self.residues
            if r.is_protein
        )


@dataclass
class Structure:
    """Parsed protein structure — the central data object.

    Contains both the rich atom-level data (for queries and reporting)
    and flat numpy arrays (for JAX computation).

    Attributes
    ----------
    name : str
        Structure identifier (typically the filename or PDB ID).
    atoms : list[Atom]
        All atoms in the structure.
    residues : dict[tuple, Residue]
        Residues keyed by (chain_id, res_seq, insertion_code).
    chains : dict[str, Chain]
        Chains keyed by chain_id.
    coords : np.ndarray
        (N, 3) float32 coordinate array for all atoms.
    elements : np.ndarray
        (N,) array of element symbols.
    atom_names : np.ndarray
        (N,) array of atom names.
    res_indices : np.ndarray
        (N,) int array mapping each atom to its residue index.
    chain_ids : np.ndarray
        (N,) array of chain IDs per atom.
    is_protein_mask : np.ndarray
        (N,) bool array — True for protein atoms.
    is_backbone_mask : np.ndarray
        (N,) bool array — True for backbone atoms (N, CA, C, O).
    is_hydrogen_mask : np.ndarray
        (N,) bool array — True for hydrogen atoms.
    """

    name: str
    atoms: list[Atom]
    residues: dict[tuple, Residue]
    chains: dict[str, Chain]

    # Flat arrays for computation (built by _build_arrays)
    coords: np.ndarray = field(repr=False, default_factory=lambda: np.empty((0, 3)))
    elements: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    atom_names: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    res_indices: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0, dtype=int))
    res_names: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    chain_ids_array: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    is_protein_mask: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0, dtype=bool))
    is_backbone_mask: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0, dtype=bool))
    is_hydrogen_mask: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0, dtype=bool))

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_chains(self) -> int:
        return len(self.chains)

    @property
    def protein_chains(self) -> list[Chain]:
        """Return only protein chains."""
        return [c for c in self.chains.values() if c.is_protein]

    def heavy_atom_coords(self) -> np.ndarray:
        """Return (M, 3) coordinates for non-hydrogen atoms only."""
        mask = ~self.is_hydrogen_mask
        return self.coords[mask]

    def protein_coords(self) -> np.ndarray:
        """Return (M, 3) coordinates for protein atoms only."""
        return self.coords[self.is_protein_mask]

    def backbone_coords(self) -> np.ndarray:
        """Return (M, 3) coordinates for backbone atoms only."""
        return self.coords[self.is_backbone_mask]


def _parse_element(element_field: str, atom_name: str) -> str:
    """Extract element symbol from PDB columns 77-78 or infer from atom name.

    Parameters
    ----------
    element_field : str
        Columns 77-78 of PDB ATOM record (may be blank).
    atom_name : str
        Atom name field (columns 13-16).

    Returns
    -------
    str
        Uppercase element symbol.
    """
    elem = element_field.strip().upper()
    if elem and elem in ELEMENT_ATOMIC_NUM:
        return elem

    # Infer from atom name: first non-digit, non-space character
    clean = atom_name.strip()
    if not clean:
        return "X"

    # Standard PDB convention: element is right-justified in columns 13-14
    # For common cases: CA → C (carbon alpha), FE → FE (iron)
    if len(clean) >= 2 and clean[:2].upper() in ELEMENT_ATOMIC_NUM:
        return clean[:2].upper()
    if clean[0].upper() in ("C", "N", "O", "S", "H", "P", "F", "I"):
        return clean[0].upper()

    return "X"


def _parse_pdb_line(line: str, is_hetatm: bool) -> Atom | None:
    """Parse a single ATOM or HETATM record from a PDB file.

    Parameters
    ----------
    line : str
        Raw PDB line (must start with ATOM or HETATM).
    is_hetatm : bool
        Whether this is a HETATM record.

    Returns
    -------
    Atom or None
        Parsed atom, or None if parsing fails.
    """
    if len(line) < 54:
        return None

    try:
        serial = int(line[6:11].strip()) if line[6:11].strip() else 0
        name = line[12:16].strip()
        alt_loc = line[16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21].strip() or "A"
        res_seq = int(line[22:26].strip()) if line[22:26].strip() else 0
        insertion_code = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = float(line[54:60].strip()) if len(line) >= 60 and line[54:60].strip() else 1.0
        b_factor = float(line[60:66].strip()) if len(line) >= 66 and line[60:66].strip() else 0.0
        element = _parse_element(
            line[76:78] if len(line) >= 78 else "",
            line[12:16],
        )
    except (ValueError, IndexError):
        return None

    return Atom(
        serial=serial,
        name=name,
        alt_loc=alt_loc,
        res_name=res_name,
        chain_id=chain_id,
        res_seq=res_seq,
        insertion_code=insertion_code,
        x=x,
        y=y,
        z=z,
        occupancy=occupancy,
        b_factor=b_factor,
        element=element,
        is_hetatm=is_hetatm,
    )


def _build_structure(atoms: list[Atom], name: str) -> Structure:
    """Assemble atoms into a Structure with residues, chains, and flat arrays.

    Parameters
    ----------
    atoms : list[Atom]
        All parsed atoms.
    name : str
        Structure identifier.

    Returns
    -------
    Structure
        Complete structure with both rich objects and flat numpy arrays.
    """
    # Build residue and chain lookups
    residues: dict[tuple, Residue] = {}
    chains: dict[str, Chain] = {}
    residue_order: list[tuple] = []  # Ordered unique residue IDs

    for atom in atoms:
        rid = atom.residue_id

        if rid not in residues:
            residues[rid] = Residue(
                chain_id=atom.chain_id,
                res_seq=atom.res_seq,
                insertion_code=atom.insertion_code,
                res_name=atom.res_name,
            )
            residue_order.append(rid)

        residues[rid].atoms[atom.name] = atom

        if atom.chain_id not in chains:
            chains[atom.chain_id] = Chain(chain_id=atom.chain_id)

    # Add residues to chains in order
    chain_residues_added: dict[str, set] = {cid: set() for cid in chains}
    for rid in residue_order:
        cid = rid[0]
        if rid not in chain_residues_added[cid]:
            chains[cid].residues.append(residues[rid])
            chain_residues_added[cid].add(rid)

    # Build residue_id → index mapping
    rid_to_index = {rid: i for i, rid in enumerate(residue_order)}

    # Build flat arrays
    n = len(atoms)
    coords = np.zeros((n, 3), dtype=np.float32)
    elements = np.empty(n, dtype="U2")
    atom_names_arr = np.empty(n, dtype="U4")
    res_indices = np.zeros(n, dtype=np.int32)
    res_names_arr = np.empty(n, dtype="U4")
    chain_ids_arr = np.empty(n, dtype="U2")
    is_protein = np.zeros(n, dtype=bool)
    is_backbone = np.zeros(n, dtype=bool)
    is_hydrogen = np.zeros(n, dtype=bool)

    for i, atom in enumerate(atoms):
        coords[i] = [atom.x, atom.y, atom.z]
        elements[i] = atom.element
        atom_names_arr[i] = atom.name[:4]
        res_indices[i] = rid_to_index[atom.residue_id]
        res_names_arr[i] = atom.res_name[:4]
        chain_ids_arr[i] = atom.chain_id
        is_protein[i] = atom.is_protein
        is_backbone[i] = atom.is_backbone
        is_hydrogen[i] = atom.is_hydrogen

    return Structure(
        name=name,
        atoms=atoms,
        residues=residues,
        chains=chains,
        coords=coords,
        elements=elements,
        atom_names=atom_names_arr,
        res_indices=res_indices,
        res_names=res_names_arr,
        chain_ids_array=chain_ids_arr,
        is_protein_mask=is_protein,
        is_backbone_mask=is_backbone,
        is_hydrogen_mask=is_hydrogen,
    )


def parse_pdb(path: str | Path, keep_hydrogens: bool = False, keep_altloc: str = "A") -> Structure:
    """Parse a PDB file into a Structure.

    Parameters
    ----------
    path : str or Path
        Path to the PDB file.
    keep_hydrogens : bool
        Whether to retain hydrogen atoms (default: False for heavy-atom analysis).
    keep_altloc : str
        Which alternate location to keep. "A" keeps the first alt loc.
        Empty string keeps all conformations.

    Returns
    -------
    Structure
        Parsed structure ready for analysis.

    Raises
    ------
    FileNotFoundError
        If the PDB file doesn't exist.
    ValueError
        If no atoms could be parsed from the file.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")

    atoms: list[Atom] = []

    with open(path) as f:
        for line in f:
            record = line[:6].strip()

            if record == "ATOM":
                atom = _parse_pdb_line(line, is_hetatm=False)
            elif record == "HETATM":
                atom = _parse_pdb_line(line, is_hetatm=True)
            else:
                continue

            if atom is None:
                continue

            # Filter hydrogens
            if not keep_hydrogens and atom.is_hydrogen:
                continue

            # Handle alternate locations: keep only the specified one
            if keep_altloc and atom.alt_loc and atom.alt_loc != keep_altloc:
                continue

            atoms.append(atom)

    if not atoms:
        raise ValueError(f"No atoms parsed from {path}")

    name = path.stem
    return _build_structure(atoms, name)


def parse_pdb_string(pdb_string: str, name: str = "unknown",
                     keep_hydrogens: bool = False, keep_altloc: str = "A") -> Structure:
    """Parse a PDB-format string into a Structure.

    Parameters
    ----------
    pdb_string : str
        PDB file content as a string.
    name : str
        Structure identifier.
    keep_hydrogens : bool
        Whether to retain hydrogen atoms.
    keep_altloc : str
        Which alternate location to keep.

    Returns
    -------
    Structure
        Parsed structure.
    """
    atoms: list[Atom] = []

    for line in pdb_string.splitlines():
        record = line[:6].strip()

        if record == "ATOM":
            atom = _parse_pdb_line(line, is_hetatm=False)
        elif record == "HETATM":
            atom = _parse_pdb_line(line, is_hetatm=True)
        else:
            continue

        if atom is None:
            continue

        if not keep_hydrogens and atom.is_hydrogen:
            continue

        if keep_altloc and atom.alt_loc and atom.alt_loc != keep_altloc:
            continue

        atoms.append(atom)

    if not atoms:
        raise ValueError("No atoms parsed from PDB string")

    return _build_structure(atoms, name)
