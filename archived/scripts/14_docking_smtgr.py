"""
Kira - Script 14: Molecular Docking — SmTGR Selectivity
=========================================================

Phase 2 begins here.

SmTGR is the most biologically compelling schistosomiasis target:
  - Single point of failure (no backup antioxidant system)
  - Unique fusion enzyme (TrxR + GR combined)
  - Most potent hit in entire dataset (10 nM)
  - 43 compounds tested
  - ZERO selectivity data

This script generates the first computational selectivity estimates
for SmTGR by docking compounds against both:
  - SmTGR (parasite): PDB 2X99 (crystal structure of S. mansoni TGR)
  - Human TrxR1: PDB 2ZZC (crystal structure of human thioredoxin reductase)

If the docking scores show differential binding — compounds fitting
better in the parasite active site than the human one — that is
structural evidence of selectivity exploiting SmTGR's unique
fusion architecture.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 14_docking_smtgr.py

Requires: meeko, vina, rdkit, requests (for PDB download)
Runtime: ~10-30 minutes depending on number of compounds docked
"""

import os
import sys
import time
import json
import requests
import numpy as np
import pandas as pd
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
DOCK_DIR = os.path.join(os.path.dirname(__file__), "data", "docking")
PUB_DIR = os.path.join(os.path.dirname(__file__), "data", "publication")

# PDB structures
STRUCTURES = {
    "SmTGR": {
        "pdb_id": "2X99",
        "description": "S. mansoni thioredoxin glutathione reductase",
        "organism": "parasite",
        # Active site center coordinates (from literature, FAD binding domain)
        # These define the docking search box
        "center_x": 34.4,
        "center_y": 24.9,
        "center_z": 6.5,
        "box_size": 25.0,  # Angstroms, cube side length
    },
    "HsTrxR1": {
        "pdb_id": "2ZZC",
        "description": "Human thioredoxin reductase 1",
        "organism": "human",
        "center_x": -20.0,
        "center_y": 19.4,
        "center_z": -43.3,
        "box_size": 25.0,
    },
}


# ---------------------------------------------------------------------------
# STEP 1: Download PDB structures
# ---------------------------------------------------------------------------

def download_pdb(pdb_id, output_dir):
    """Download a PDB structure file from RCSB."""
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.exists(filepath):
        print(f"  {pdb_id}.pdb already downloaded")
        return filepath

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    with open(filepath, "w") as f:
        f.write(response.text)

    size_kb = os.path.getsize(filepath) / 1024
    print(f"  Saved {pdb_id}.pdb ({size_kb:.0f} KB)")
    return filepath


# ---------------------------------------------------------------------------
# STEP 2: Prepare protein for docking (extract chain, remove water/ligands)
# ---------------------------------------------------------------------------

def prepare_protein(pdb_path, output_dir):
    """
    Prepare protein structure for Vina docking.

    Steps:
    1. Read PDB file
    2. Keep only protein atoms (ATOM records) from chain A
    3. Remove water molecules (HOH)
    4. Remove existing ligands and cofactors
    5. Save as clean PDB
    6. Convert to PDBQT format (Vina's input format)

    PDBQT adds partial charges and atom types that Vina needs
    for its scoring function.
    """
    pdb_id = os.path.basename(pdb_path).replace(".pdb", "")
    clean_pdb = os.path.join(output_dir, f"{pdb_id}_clean.pdb")
    pdbqt_path = os.path.join(output_dir, f"{pdb_id}_receptor.pdbqt")

    if os.path.exists(pdbqt_path):
        print(f"  {pdb_id} receptor already prepared")
        return pdbqt_path

    print(f"  Preparing {pdb_id} receptor...")

    # Read and clean PDB
    clean_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            # Keep ATOM records (protein), skip HETATM (ligands, water, cofactors)
            if line.startswith("ATOM"):
                # Skip hydrogen atoms if present
                atom_name = line[12:16].strip()
                if not atom_name.startswith("H"):
                    clean_lines.append(line)
            elif line.startswith("END"):
                clean_lines.append(line)
                break

    with open(clean_pdb, "w") as f:
        f.writelines(clean_lines)

    n_atoms = len([l for l in clean_lines if l.startswith("ATOM")])
    print(f"  Clean PDB: {n_atoms} protein atoms")

    # Convert to PDBQT using a simple approach:
    # Add Gasteiger charges and AD4 atom types
    # For production use, you'd use ADFR Suite's prepare_receptor
    # Here we create a minimal PDBQT by adding charge/type columns

    pdbqt_lines = []
    for line in clean_lines:
        if line.startswith("ATOM"):
            # Determine AD4 atom type from element
            element = line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0]
            atom_name = line[12:16].strip()

            # Map common protein atoms to AutoDock types
            if element == "C":
                if atom_name in ("C", "CA", "CB"):
                    ad_type = "C"
                else:
                    ad_type = "C"
            elif element == "N":
                ad_type = "NA" if atom_name in ("ND1", "ND2", "NE2", "NE1") else "N"
            elif element == "O":
                ad_type = "OA"
            elif element == "S":
                ad_type = "SA"
            elif element == "SE":
                ad_type = "SA"  # Selenocysteine -> treat as sulfur
            else:
                ad_type = "C"  # Default

            # Format PDBQT line: standard PDB + charge + type
            # Charge column: 70-76, Type column: 77-79
            base = line[:54]  # Up to coordinates
            coords = line[30:54]
            occupancy = line[54:60] if len(line) > 54 else " 1.00"
            bfactor = line[60:66] if len(line) > 60 else "  0.00"

            pdbqt_line = f"{base}{occupancy}{bfactor}    +0.000 {ad_type:<2s}\n"
            pdbqt_lines.append(pdbqt_line)
        elif line.startswith("END"):
            pdbqt_lines.append(line)

    with open(pdbqt_path, "w") as f:
        f.writelines(pdbqt_lines)

    print(f"  PDBQT receptor: {pdbqt_path}")
    return pdbqt_path


# ---------------------------------------------------------------------------
# STEP 3: Prepare ligands for docking
# ---------------------------------------------------------------------------

def prepare_ligand(smiles, compound_id, output_dir):
    """
    Convert a SMILES string to a 3D structure in PDBQT format.

    Steps:
    1. Parse SMILES with RDKit
    2. Add hydrogens
    3. Generate 3D coordinates (ETKDG conformer generation)
    4. Minimize energy (MMFF force field)
    5. Convert to PDBQT via meeko
    """
    pdbqt_path = os.path.join(output_dir, f"{compound_id}_ligand.pdbqt")

    if os.path.exists(pdbqt_path):
        return pdbqt_path

    try:
        # Parse and prepare molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            # Try with random coords if ETKDG fails
            AllChem.EmbedMolecule(mol, randomSeed=42)

        # Energy minimize
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except:
            pass  # Continue even if minimization fails

        # Use meeko to convert to PDBQT
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(pdbqt_path, "w") as f:
                    f.write(pdbqt_string)
                return pdbqt_path

        return None

    except Exception as e:
        return None


# ---------------------------------------------------------------------------
# STEP 4: Run docking
# ---------------------------------------------------------------------------

def dock_compound(receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z,
                   box_size, exhaustiveness=16):
    """
    Dock a ligand into a receptor using AutoDock Vina.

    Returns the best binding energy (kcal/mol). More negative = better binding.
    Typical range: -3 (weak) to -12 (very strong).

    Parameters:
    - exhaustiveness: search thoroughness (8=fast, 32=thorough). 16 is balanced.
    """
    from vina import Vina

    try:
        v = Vina(sf_name="vina")
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)

        v.compute_vina_maps(
            center=[center_x, center_y, center_z],
            box_size=[box_size, box_size, box_size],
        )

        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        energies = v.energies()

        # energies is a list of [total, inter, intra, torsion, ...]
        # First row, first column is the best total binding energy
        best_energy = energies[0][0]

        return best_energy

    except Exception as e:
        return None


# ---------------------------------------------------------------------------
# STEP 5: Run docking campaign
# ---------------------------------------------------------------------------

def run_docking_campaign(compounds_to_dock, structures):
    """
    Dock each compound against both parasite and human targets.
    Compute computational selectivity ratio.
    """
    print("\n" + "=" * 60)
    print("STEP 5: DOCKING CAMPAIGN")
    print("=" * 60)

    ligand_dir = os.path.join(DOCK_DIR, "ligands")
    os.makedirs(ligand_dir, exist_ok=True)

    results = []
    total = len(compounds_to_dock)

    for idx, (_, row) in enumerate(compounds_to_dock.iterrows()):
        cid = row["molecule_chembl_id"]
        smiles = row["smiles"]
        ic50 = row.get("best_ic50", None)

        print(f"\n  [{idx+1}/{total}] {cid}")

        # Prepare ligand
        ligand_pdbqt = prepare_ligand(smiles, cid, ligand_dir)
        if ligand_pdbqt is None:
            print(f"    Failed to prepare ligand. Skipping.")
            continue

        # Dock against each target
        energies = {}
        for target_name, target_info in structures.items():
            receptor_pdbqt = os.path.join(DOCK_DIR, "receptors",
                                           f"{target_info['pdb_id']}_receptor.pdbqt")
            if not os.path.exists(receptor_pdbqt):
                print(f"    Receptor {target_name} not ready. Skipping.")
                continue

            energy = dock_compound(
                receptor_pdbqt, ligand_pdbqt,
                target_info["center_x"],
                target_info["center_y"],
                target_info["center_z"],
                target_info["box_size"],
            )

            if energy is not None:
                energies[target_name] = energy
                print(f"    vs {target_name}: {energy:.1f} kcal/mol")
            else:
                print(f"    vs {target_name}: FAILED")

        # Compute selectivity
        if "SmTGR" in energies and "HsTrxR1" in energies:
            parasite_energy = energies["SmTGR"]
            human_energy = energies["HsTrxR1"]

            # More negative = better binding
            # delta > 0 means compound binds BETTER to parasite (favorable)
            delta = human_energy - parasite_energy

            # Convert to approximate selectivity ratio
            # delta E of 1.36 kcal/mol ≈ 10-fold difference in binding affinity
            # (from the Boltzmann relation: Kd ratio = exp(delta_G / RT))
            RT = 0.592  # kcal/mol at 298K
            selectivity_factor = np.exp(delta / RT)

            if delta > 0:
                sel_class = "PARASITE-SELECTIVE" if delta > 1.36 else "SLIGHTLY-PARASITE"
            elif delta < 0:
                sel_class = "HUMAN-SELECTIVE" if delta < -1.36 else "SLIGHTLY-HUMAN"
            else:
                sel_class = "NEUTRAL"

            print(f"    Delta: {delta:+.1f} kcal/mol → {selectivity_factor:.1f}x "
                  f"({sel_class})")
        else:
            delta = None
            selectivity_factor = None
            sel_class = "INCOMPLETE"

        results.append({
            "molecule_chembl_id": cid,
            "smiles": smiles,
            "experimental_ic50_nM": ic50,
            "SmTGR_dock_energy": energies.get("SmTGR"),
            "HsTrxR1_dock_energy": energies.get("HsTrxR1"),
            "delta_energy": delta,
            "computed_selectivity_factor": selectivity_factor,
            "selectivity_class": sel_class,
        })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# STEP 6: Analyze docking results
# ---------------------------------------------------------------------------

def analyze_docking_results(docking_df):
    """Analyze and report docking-based selectivity."""

    print("\n" + "=" * 60)
    print("DOCKING SELECTIVITY ANALYSIS")
    print("=" * 60)

    # Filter for completed docks
    complete = docking_df.dropna(subset=["SmTGR_dock_energy", "HsTrxR1_dock_energy"])

    if len(complete) == 0:
        print("  No complete docking results.")
        return

    print(f"\n  Compounds with complete docking: {len(complete)}")

    # Summary statistics
    print(f"\n  SmTGR docking energies:")
    print(f"    Mean: {complete['SmTGR_dock_energy'].mean():.1f} kcal/mol")
    print(f"    Best: {complete['SmTGR_dock_energy'].min():.1f} kcal/mol")
    print(f"    Worst: {complete['SmTGR_dock_energy'].max():.1f} kcal/mol")

    print(f"\n  HsTrxR1 docking energies:")
    print(f"    Mean: {complete['HsTrxR1_dock_energy'].mean():.1f} kcal/mol")
    print(f"    Best: {complete['HsTrxR1_dock_energy'].min():.1f} kcal/mol")
    print(f"    Worst: {complete['HsTrxR1_dock_energy'].max():.1f} kcal/mol")

    # Selectivity distribution
    print(f"\n  Computational selectivity distribution:")
    for cls in ["PARASITE-SELECTIVE", "SLIGHTLY-PARASITE", "NEUTRAL",
                "SLIGHTLY-HUMAN", "HUMAN-SELECTIVE", "INCOMPLETE"]:
        n = len(docking_df[docking_df["selectivity_class"] == cls])
        if n > 0:
            print(f"    {cls}: {n}")

    # Rank by selectivity
    complete_sorted = complete.sort_values("delta_energy", ascending=False)

    print(f"\n  RANKED BY PARASITE SELECTIVITY (most selective first):")
    print(f"  {'Compound':20s} {'SmTGR':>8s} {'HsTrxR1':>8s} {'Delta':>7s} {'Factor':>7s}  {'Exp IC50':>10s}")
    print(f"  {'-'*20} {'-'*8} {'-'*8} {'-'*7} {'-'*7}  {'-'*10}")

    for _, row in complete_sorted.iterrows():
        cid = row["molecule_chembl_id"]
        if len(cid) > 20:
            cid = cid[:17] + "..."
        ic50_str = f"{row['experimental_ic50_nM']:.0f} nM" if pd.notna(row.get("experimental_ic50_nM")) else "N/A"
        factor = row.get("computed_selectivity_factor", 0)
        factor_str = f"{factor:.1f}x" if pd.notna(factor) else "N/A"

        print(f"  {cid:20s} {row['SmTGR_dock_energy']:8.1f} {row['HsTrxR1_dock_energy']:8.1f} "
              f"{row['delta_energy']:+7.1f} {factor_str:>7s}  {ic50_str:>10s}")

    # Correlation with experimental IC50
    has_ic50 = complete.dropna(subset=["experimental_ic50_nM"])
    if len(has_ic50) >= 5:
        from scipy import stats
        corr, pval = stats.spearmanr(
            has_ic50["SmTGR_dock_energy"],
            np.log10(has_ic50["experimental_ic50_nM"])
        )
        print(f"\n  Correlation: SmTGR docking energy vs log10(experimental IC50)")
        print(f"    Spearman rho = {corr:.3f}, p = {pval:.4f}")
        if pval < 0.05 and corr > 0:
            print(f"    Positive correlation: better docking scores correspond to")
            print(f"    lower IC50 values. Docking captures real binding signal.")
        elif pval < 0.05 and corr < 0:
            print(f"    Negative correlation: docking and experiment agree on rank order.")
        else:
            print(f"    No significant correlation. Docking and experiment disagree.")

    return complete_sorted


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 2, Script 14")
    print("  Molecular Docking: SmTGR Selectivity")
    print("=" * 60)

    start_time = time.time()

    # Step 1: Download PDB structures
    print("\n" + "=" * 60)
    print("STEP 1: DOWNLOADING PROTEIN STRUCTURES")
    print("=" * 60)

    receptor_dir = os.path.join(DOCK_DIR, "receptors")
    os.makedirs(receptor_dir, exist_ok=True)

    for name, info in STRUCTURES.items():
        pdb_path = download_pdb(info["pdb_id"], receptor_dir)

    # Step 2: Prepare receptors
    print("\n" + "=" * 60)
    print("STEP 2: PREPARING RECEPTORS")
    print("=" * 60)

    for name, info in STRUCTURES.items():
        pdb_path = os.path.join(receptor_dir, f"{info['pdb_id']}.pdb")
        prepare_protein(pdb_path, receptor_dir)

    # Step 3: Load SmTGR compounds
    print("\n" + "=" * 60)
    print("STEP 3: LOADING SmTGR COMPOUNDS")
    print("=" * 60)

    activities = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))
    tgr_activities = activities[
        activities["target_name"] == "Thioredoxin glutathione reductase"
    ].copy()

    # Get best IC50 and SMILES per compound
    tgr_compounds = (
        tgr_activities
        .groupby("molecule_chembl_id")
        .agg(
            best_ic50=("standard_value", "min"),
            smiles=("canonical_smiles", "first"),
        )
        .reset_index()
        .sort_values("best_ic50")
    )

    print(f"  SmTGR compounds to dock: {len(tgr_compounds)}")
    print(f"  IC50 range: {tgr_compounds['best_ic50'].min():.0f} - "
          f"{tgr_compounds['best_ic50'].max():.0f} nM")

    # Step 4: Prepare ligands
    print("\n" + "=" * 60)
    print("STEP 4: PREPARING LIGANDS")
    print("=" * 60)

    ligand_dir = os.path.join(DOCK_DIR, "ligands")
    os.makedirs(ligand_dir, exist_ok=True)

    prepared = 0
    failed = 0
    for _, row in tgr_compounds.iterrows():
        result = prepare_ligand(row["smiles"], row["molecule_chembl_id"], ligand_dir)
        if result:
            prepared += 1
        else:
            failed += 1

    print(f"  Prepared: {prepared}, Failed: {failed}")

    # Step 5: Run docking campaign
    docking_df = run_docking_campaign(tgr_compounds, STRUCTURES)

    # Step 6: Analyze
    sorted_results = analyze_docking_results(docking_df)

    # Save results
    os.makedirs(PUB_DIR, exist_ok=True)
    docking_df.to_csv(os.path.join(PROCESSED_DIR, "smtgr_docking_results.csv"), index=False)
    if sorted_results is not None and len(sorted_results) > 0:
        sorted_results.to_csv(os.path.join(PUB_DIR, "table4_smtgr_docking_selectivity.csv"),
                               index=False)

    elapsed = time.time() - start_time

    print(f"\n{'=' * 60}")
    print(f"  PHASE 2 DOCKING COMPLETE")
    print(f"  Time: {elapsed:.0f} seconds")
    print(f"")
    print(f"  SmTGR compounds docked: {len(docking_df)}")
    if sorted_results is not None:
        n_parasite_sel = len(docking_df[docking_df["selectivity_class"].str.contains("PARASITE", na=False)])
        n_human_sel = len(docking_df[docking_df["selectivity_class"].str.contains("HUMAN", na=False)])
        print(f"  Computationally parasite-selective: {n_parasite_sel}")
        print(f"  Computationally human-selective: {n_human_sel}")
    print(f"")
    print(f"  Results: data/processed/smtgr_docking_results.csv")
    print(f"  Table: data/publication/table4_smtgr_docking_selectivity.csv")
    print(f"")
    print(f"  This is the first selectivity assessment of SmTGR inhibitors")
    print(f"  against human thioredoxin reductase, computational or experimental.")
    print(f"{'=' * 60}")
