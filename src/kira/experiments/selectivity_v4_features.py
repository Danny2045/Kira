from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors, rdmolops

from kira.causality.divergence import AA_PROPERTIES, DEFAULT_PROPERTIES
from kira.experiments import TARGET_PAIRS
from kira.experiments.selectivity_features import PocketFeatures, compute_pocket_features


ONE_TO_THREE = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

AROMATIC_ONE_LETTER = {"F", "W", "Y", "H"}
HYDROPHOBIC_ONE_LETTER = {"A", "V", "I", "L", "M", "F", "W", "Y", "C", "P"}

V4_BLOCK_ORDER = [
    "compound_fp",
    "compound_desc",
    "pair_delta",
    "parasite_side",
    "human_side",
    "compat_diff",
]


def _coerce_mol(smiles_or_mol: str | Chem.Mol) -> Chem.Mol:
    if isinstance(smiles_or_mol, Chem.Mol):
        return smiles_or_mol

    if not isinstance(smiles_or_mol, str) or not smiles_or_mol.strip():
        raise ValueError("canonical_smiles must be a non-empty string")

    mol = Chem.MolFromSmiles(smiles_or_mol)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles_or_mol}")

    return mol


def compound_descriptor_names() -> list[str]:
    return [
        "compound_desc__mol_wt",
        "compound_desc__logp",
        "compound_desc__tpsa",
        "compound_desc__hbd",
        "compound_desc__hba",
        "compound_desc__heavy_atom_count",
        "compound_desc__aromatic_ring_count",
        "compound_desc__rotatable_bonds",
        "compound_desc__formal_charge",
    ]


def compute_compound_descriptors(smiles_or_mol: str | Chem.Mol) -> np.ndarray:
    mol = _coerce_mol(smiles_or_mol)

    values = [
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        float(rdMolDescriptors.CalcTPSA(mol)),
        float(rdMolDescriptors.CalcNumHBD(mol)),
        float(rdMolDescriptors.CalcNumHBA(mol)),
        float(mol.GetNumHeavyAtoms()),
        float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        float(rdmolops.GetFormalCharge(mol)),
    ]
    return np.asarray(values, dtype=np.float32)


def compound_fingerprint_names(n_bits: int = 256) -> list[str]:
    return [f"compound_fp__bit_{i:03d}" for i in range(n_bits)]


def compute_compound_fingerprint(
    smiles_or_mol: str | Chem.Mol,
    *,
    radius: int = 2,
    n_bits: int = 256,
) -> np.ndarray:
    mol = _coerce_mol(smiles_or_mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)


def _residue_properties_from_one_letter(residue: str) -> dict[str, float]:
    three = ONE_TO_THREE.get(str(residue).upper(), None)
    if three is None:
        return DEFAULT_PROPERTIES
    return AA_PROPERTIES.get(three, DEFAULT_PROPERTIES)


def side_summary_feature_names(prefix: str | None = None) -> list[str]:
    base = [
        "mean_hydrophobicity",
        "mean_charge",
        "mean_volume",
        "mean_hbd",
        "mean_hba",
        "frac_positive",
        "frac_negative",
        "frac_aromatic",
        "frac_hydrophobic",
    ]
    if prefix is None:
        return base
    return [f"{prefix}__{name}" for name in base]


def compute_pocket_side_summary(pocket_sequence: str) -> np.ndarray:
    seq = str(pocket_sequence or "")
    if not seq:
        return np.zeros((len(side_summary_feature_names()),), dtype=np.float32)

    props = [_residue_properties_from_one_letter(res) for res in seq]
    hydro = np.array([p["hydrophobicity"] for p in props], dtype=np.float32)
    charge = np.array([p["charge"] for p in props], dtype=np.float32)
    volume = np.array([p["volume"] for p in props], dtype=np.float32)
    hbd = np.array([p["hbd"] for p in props], dtype=np.float32)
    hba = np.array([p["hba"] for p in props], dtype=np.float32)

    residues = [str(r).upper() for r in seq]

    out = np.array(
        [
            float(hydro.mean()),
            float(charge.mean()),
            float(volume.mean()),
            float(hbd.mean()),
            float(hba.mean()),
            float(np.mean(charge > 0.0)),
            float(np.mean(charge < 0.0)),
            float(np.mean([r in AROMATIC_ONE_LETTER for r in residues])),
            float(np.mean([r in HYDROPHOBIC_ONE_LETTER for r in residues])),
        ],
        dtype=np.float32,
    )
    return out


def compute_pair_delta(target_pair_id: str) -> np.ndarray:
    pair = TARGET_PAIRS[target_pair_id]
    pf = compute_pocket_features(
        pair.parasite_pocket_sequence,
        pair.human_pocket_sequence,
        pair_name=pair.name,
        disease=pair.disease,
    )
    return pf.to_array().astype(np.float32)


def compute_pair_side_summaries(target_pair_id: str) -> tuple[np.ndarray, np.ndarray]:
    pair = TARGET_PAIRS[target_pair_id]
    parasite_side = compute_pocket_side_summary(pair.parasite_pocket_sequence)
    human_side = compute_pocket_side_summary(pair.human_pocket_sequence)
    return parasite_side, human_side


def compatibility_feature_names() -> list[str]:
    desc = compound_descriptor_names()
    side = side_summary_feature_names()
    return [f"compat_diff__{d}__x__{s}" for d in desc for s in side]


def compute_compatibility_diff(
    compound_desc: np.ndarray,
    parasite_side: np.ndarray,
    human_side: np.ndarray,
) -> np.ndarray:
    side_diff = parasite_side - human_side
    compat = np.outer(compound_desc, side_diff).ravel()
    return compat.astype(np.float32)


def block_dimensions(n_bits: int = 256) -> dict[str, int]:
    desc_n = len(compound_descriptor_names())
    side_n = len(side_summary_feature_names())
    return {
        "compound_fp": n_bits,
        "compound_desc": desc_n,
        "pair_delta": len(PocketFeatures.feature_names()),
        "parasite_side": side_n,
        "human_side": side_n,
        "compat_diff": desc_n * side_n,
    }


def block_slices(n_bits: int = 256) -> dict[str, slice]:
    dims = block_dimensions(n_bits=n_bits)
    out: dict[str, slice] = {}
    start = 0
    for block in V4_BLOCK_ORDER:
        width = dims[block]
        out[block] = slice(start, start + width)
        start += width
    return out


def build_v4_blocks(
    canonical_smiles: str,
    target_pair_id: str,
    *,
    n_bits: int = 256,
) -> dict[str, np.ndarray]:
    if target_pair_id not in TARGET_PAIRS:
        raise KeyError(f"Unknown target_pair_id: {target_pair_id}")

    compound_fp = compute_compound_fingerprint(canonical_smiles, n_bits=n_bits)
    compound_desc = compute_compound_descriptors(canonical_smiles)
    pair_delta = compute_pair_delta(target_pair_id)
    parasite_side, human_side = compute_pair_side_summaries(target_pair_id)
    compat_diff = compute_compatibility_diff(compound_desc, parasite_side, human_side)

    return {
        "compound_fp": compound_fp,
        "compound_desc": compound_desc,
        "pair_delta": pair_delta,
        "parasite_side": parasite_side,
        "human_side": human_side,
        "compat_diff": compat_diff,
    }


def build_feature_names(n_bits: int = 256) -> list[str]:
    return (
        compound_fingerprint_names(n_bits=n_bits)
        + compound_descriptor_names()
        + PocketFeatures.feature_names()
        + side_summary_feature_names(prefix="parasite_side")
        + side_summary_feature_names(prefix="human_side")
        + compatibility_feature_names()
    )


def build_feature_vector(
    canonical_smiles: str,
    target_pair_id: str,
    *,
    n_bits: int = 256,
) -> np.ndarray:
    blocks = build_v4_blocks(canonical_smiles, target_pair_id, n_bits=n_bits)
    vec = np.concatenate([blocks[name] for name in V4_BLOCK_ORDER]).astype(np.float32)
    return vec


def featurize_dataframe(
    df: pd.DataFrame,
    *,
    n_bits: int = 256,
) -> tuple[np.ndarray, list[str]]:
    required = {"canonical_smiles", "target_pair_id"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"DataFrame missing required columns: {sorted(missing)}")

    names = build_feature_names(n_bits=n_bits)
    rows: list[np.ndarray] = []

    for row in df.itertuples(index=False):
        rows.append(
            build_feature_vector(
                canonical_smiles=getattr(row, "canonical_smiles"),
                target_pair_id=getattr(row, "target_pair_id"),
                n_bits=n_bits,
            )
        )

    if rows:
        X = np.vstack(rows).astype(np.float32)
    else:
        X = np.zeros((0, len(names)), dtype=np.float32)

    return X, names


def select_feature_blocks(
    X: np.ndarray,
    block_names: list[str],
    *,
    n_bits: int = 256,
) -> np.ndarray:
    if X.ndim != 2:
        raise ValueError("X must be a 2D feature matrix")

    slices = block_slices(n_bits=n_bits)
    parts = [X[:, slices[name]] for name in block_names]
    if not parts:
        return np.zeros((X.shape[0], 0), dtype=X.dtype)
    return np.concatenate(parts, axis=1)
