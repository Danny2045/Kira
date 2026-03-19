#!/bin/bash
# Run this from ~/kira with: bash create_script05.sh
# It creates 05_structural_similarity.py for you.

cat > 05_structural_similarity.py << 'PYTHONSCRIPT'
"""
Kira - Script 05: Structural Similarity + Improved Ranking
============================================================

Script 04 revealed a critical flaw: praziquantel (the actual standard
of care) ranked 193/194 because it has no target-based IC50 in ChEMBL.

This script fixes that by adding a second evidence pathway: structural
similarity. Instead of asking only "does this compound have measured
activity?", we now also ask "does this compound LOOK LIKE compounds
that have measured activity?"

NEW SIGNAL: Structural similarity
  - Compute Morgan fingerprints for all compounds
  - For each compound, measure Tanimoto similarity to known actives
  - Structurally similar compounds score higher even without IC50 data

This also properly rescues curated drugs (praziquantel, oxamniquine,
mefloquine) by giving them molecular structure information.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 05_structural_similarity.py
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

# Suppress RDKit warnings (noisy but not dangerous)
RDLogger.logger().setLevel(RDLogger.ERROR)

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")

# --- SMILES for curated drugs (not in ChEMBL activity data) ---
# These are the standard SMILES representations from ChEMBL/PubChem.
CURATED_SMILES = {
    "PRAZIQUANTEL": "O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1",
    "OXAMNIQUINE": "CC(CO)Nc1ccc2c(c1)[C@@H](C)C[C@H](C)N2O",
    "MEFLOQUINE": "OC(c1cc(C(F)(F)F)nc2c(C(F)(F)F)cccc12)C1CCCCN1",
    "METFORMIN": "CN(C)C(=N)NC(=N)N",
    "LISINOPRIL": "NCCCC[C@@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCCC1C(=O)O",
    "SERTRALINE": "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21",
    "OMEPRAZOLE": "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    "ATORVASTATIN": "CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(-c2ccccc2)c(-c2ccc(F)cc2)c1C(=O)Nc1ccccc1",
    "AMLODIPINE": "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl",
    "LEVOTHYROXINE": "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O",
    "MONTELUKAST": "CC(C)(O)c1ccccc1CC[C@@H](SCC1(CC(=O)O)CC1)c1cccc(-c2cccc3ccc(Cl)cc23)c1",
    "CLOPIDOGREL": "COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1",
    "TAMSULOSIN": "CCOc1ccc(CC(C)NCC[C@@H](O)c2ccc(OC)c(S(N)(=O)=O)c2)cc1",
    "ESCITALOPRAM": "N#Cc1ccc2c(c1)C(CCCCN1CCC1)(OC2)c1ccc(F)cc1",
    "GABAPENTIN": "NCC1(CC(=O)O)CCCCC1",
    "PANTOPRAZOLE": "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC",
    "ROSUVASTATIN": "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O",
    "TRAMADOL": "COc1cccc(C2(O)CCCCC2CN(C)C)c1",
}

# Potency threshold for "reference active" set
ACTIVE_IC50_THRESHOLD = 1000  # nM

# Signal weights for v2 composite score
WEIGHT_POTENCY = 0.30        # Reduced from 0.40 (was dominating)
WEIGHT_TARGET = 0.20         # Reduced from 0.25
WEIGHT_CONFIDENCE = 0.05     # Reduced from 0.10
WEIGHT_DRUG_STAGE = 0.15     # Same
WEIGHT_MULTITARGET = 0.05    # Reduced from 0.10
WEIGHT_SIMILARITY = 0.25     # NEW: structural similarity to known actives

# Target essentiality (same as Script 04)
TARGET_ESSENTIALITY = {
    "Thioredoxin glutathione reductase": 1.0,
    "Histone deacetylase 8": 0.85,
    "Dihydroorotate dehydrogenase (quinone), mitochondrial": 0.80,
    "Cathepsin B1 isotype 1": 0.70,
    "Venus kinase receptor 2": 0.60,
    "NAD-dependent protein deacetylase": 0.50,
    "ATP-diphosphohydrolase 1": 0.55,
    "Thioredoxin peroxidase": 0.65,
    "Voltage-activated calcium channel beta 1 subunit": 0.75,
    "Voltage-activated calcium channel beta 2 subunit": 0.75,
    "Calcium channels (mechanism unclear)": 0.75,
    "Sulfotransferase (prodrug activation)": 0.60,
    "Unknown (whole-worm activity)": 0.50,
    "None (negative control)": 0.0,
}
DEFAULT_ESSENTIALITY = 0.3

DRUG_STAGE_SCORES = {
    4.0: 1.0, 3.0: 0.7, 2.0: 0.5, 1.0: 0.3, 0.0: 0.1,
}


# ---------------------------------------------------------------------------
# STEP 1: Collect SMILES for all compounds
# ---------------------------------------------------------------------------

def collect_smiles(eval_df, activities_df):
    """
    Gather SMILES strings for every compound in the evaluation set.
    Sources: ChEMBL activity data (has canonical_smiles) + curated lookup.
    """
    print("\n" + "=" * 60)
    print("COLLECTING MOLECULAR STRUCTURES (SMILES)")
    print("=" * 60)

    smiles_map = {}

    # From ChEMBL activity data
    if "canonical_smiles" in activities_df.columns:
        for _, row in activities_df.iterrows():
            cid = row["molecule_chembl_id"]
            smi = row.get("canonical_smiles")
            if pd.notna(smi) and cid not in smiles_map:
                smiles_map[cid] = smi

    print(f"  SMILES from ChEMBL activity data: {len(smiles_map)}")

    # From curated lookup
    curated_count = 0
    for _, row in eval_df.iterrows():
        cid = row["molecule_chembl_id"]
        name = row.get("pref_name", "")
        if cid not in smiles_map and pd.notna(name) and name in CURATED_SMILES:
            smiles_map[cid] = CURATED_SMILES[name]
            curated_count += 1

    print(f"  SMILES from curated lookup: {curated_count}")
    print(f"  Total compounds with SMILES: {len(smiles_map)}")

    # Check which eval compounds are missing SMILES
    missing = []
    for _, row in eval_df.iterrows():
        if row["molecule_chembl_id"] not in smiles_map:
            missing.append(row.get("pref_name", row["molecule_chembl_id"]))
    if missing:
        print(f"  Missing SMILES: {len(missing)} compounds")
        for m in missing[:5]:
            print(f"    {m}")
    else:
        print(f"  All evaluation set compounds have SMILES.")

    return smiles_map


# ---------------------------------------------------------------------------
# STEP 2: Compute Morgan fingerprints
# ---------------------------------------------------------------------------

def compute_fingerprints(smiles_map):
    """
    Convert SMILES strings to Morgan fingerprints (radius=2, 2048 bits).

    Morgan fingerprints encode circular substructures around each atom.
    Radius=2 means it captures features up to 2 bonds away from each atom.
    2048 bits is the standard fingerprint length.

    This is the standard representation used in medicinal chemistry for
    similarity searching.
    """
    print("\n" + "=" * 60)
    print("COMPUTING MOLECULAR FINGERPRINTS")
    print("=" * 60)

    fingerprints = {}
    failures = []

    for cid, smiles in smiles_map.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            fingerprints[cid] = fp
        else:
            failures.append(cid)

    print(f"  Successfully computed: {len(fingerprints)} fingerprints")
    if failures:
        print(f"  Failed (invalid SMILES): {len(failures)}")

    return fingerprints


# ---------------------------------------------------------------------------
# STEP 3: Compute similarity to known actives
# ---------------------------------------------------------------------------

def compute_similarity_scores(eval_df, activities_df, fingerprints):
    """
    For each compound, compute its maximum Tanimoto similarity to any
    known potent active (IC50 < 1000 nM).

    Tanimoto similarity = |A intersection B| / |A union B|
    where A and B are the bit sets of two fingerprints.

    Range: 0 (completely different) to 1 (identical).
    In drug discovery, Tanimoto > 0.4 is often considered "similar".
    """
    print("\n" + "=" * 60)
    print("COMPUTING STRUCTURAL SIMILARITY TO KNOWN ACTIVES")
    print("=" * 60)

    # Identify reference actives (potent compounds with fingerprints)
    potent_actives = set(
        activities_df[
            activities_df["standard_value"] < ACTIVE_IC50_THRESHOLD
        ]["molecule_chembl_id"].unique()
    )

    reference_fps = {
        cid: fp for cid, fp in fingerprints.items()
        if cid in potent_actives
    }

    print(f"  Reference active compounds: {len(reference_fps)}")
    print(f"  Total compounds to score: {len(fingerprints)}")

    if not reference_fps:
        print("  WARNING: No reference actives with fingerprints found!")
        return {}

    # Compute max similarity for each compound
    similarity_scores = {}
    ref_fps_list = list(reference_fps.values())
    ref_ids_list = list(reference_fps.keys())

    for cid, fp in fingerprints.items():
        # Compute Tanimoto to all reference actives
        similarities = DataStructs.BulkTanimotoSimilarity(fp, ref_fps_list)
        max_sim = max(similarities)
        max_idx = similarities.index(max_sim)
        nearest_active = ref_ids_list[max_idx]

        similarity_scores[cid] = {
            "max_similarity": max_sim,
            "nearest_active": nearest_active,
        }

    # Show some examples
    print(f"\n  Similarity examples:")

    # Known drugs
    for drug_name in ["PRAZIQUANTEL", "OXAMNIQUINE", "MEFLOQUINE", "ATOVAQUONE"]:
        for _, row in eval_df.iterrows():
            if row.get("pref_name") == drug_name:
                cid = row["molecule_chembl_id"]
                if cid in similarity_scores:
                    sim = similarity_scores[cid]["max_similarity"]
                    nearest = similarity_scores[cid]["nearest_active"]
                    print(f"    {drug_name:20s}  max similarity = {sim:.3f}  (nearest: {nearest})")
                break

    # Negative controls
    print()
    neg_sims = []
    for _, row in eval_df.iterrows():
        if row.get("source") == "curated_negative":
            cid = row["molecule_chembl_id"]
            if cid in similarity_scores:
                sim = similarity_scores[cid]["max_similarity"]
                neg_sims.append(sim)
                name = row.get("pref_name", cid)

    if neg_sims:
        print(f"    Negative controls: mean similarity = {np.mean(neg_sims):.3f}, "
              f"max = {max(neg_sims):.3f}")

    return similarity_scores


# ---------------------------------------------------------------------------
# STEP 4: Rebuild composite score with similarity signal
# ---------------------------------------------------------------------------

def compute_v2_scores(eval_df, activities_df, fingerprints, similarity_scores):
    """
    Compute v2 composite score with structural similarity added.
    """
    print("\n" + "=" * 60)
    print("COMPUTING V2 COMPOSITE SCORES")
    print("=" * 60)

    # Multi-target scores
    target_counts = (
        activities_df
        .groupby("molecule_chembl_id")["target_name"]
        .nunique()
    )
    multitarget_scores = {}
    for cid, n in target_counts.items():
        if n >= 3:
            multitarget_scores[cid] = 1.0
        elif n == 2:
            multitarget_scores[cid] = 0.7
        else:
            multitarget_scores[cid] = 0.3

    results = []

    for _, row in eval_df.iterrows():
        cid = row["molecule_chembl_id"]

        # Signal 1: Potency
        ic50 = row.get("best_value")
        if pd.notna(ic50) and ic50 > 0:
            pic50 = -np.log10(ic50 * 1e-9)
            potency = max(0.0, min(1.0, (pic50 - 4.0) / 5.0))
        else:
            potency = 0.0

        # Signal 2: Target essentiality
        target_name = row.get("best_target", "")
        target_score = TARGET_ESSENTIALITY.get(target_name, DEFAULT_ESSENTIALITY)

        # Signal 3: Confidence
        n_meas = row.get("n_measurements", 0)
        if pd.isna(n_meas):
            n_meas = 0
        confidence = min(1.0, 0.2 + 0.3 * np.log2(max(1, int(n_meas))))

        # Signal 4: Drug stage
        phase = row.get("max_phase")
        if pd.notna(phase):
            drug_stage = DRUG_STAGE_SCORES.get(float(phase), 0.1)
        else:
            drug_stage = 0.1

        # Signal 5: Multi-target
        multitarget = multitarget_scores.get(cid, 0.1)

        # Signal 6: Structural similarity (NEW)
        if cid in similarity_scores:
            similarity = similarity_scores[cid]["max_similarity"]
            nearest = similarity_scores[cid]["nearest_active"]
        else:
            similarity = 0.0
            nearest = "N/A"

        # V2 Composite
        composite = (
            WEIGHT_POTENCY * potency
            + WEIGHT_TARGET * target_score
            + WEIGHT_CONFIDENCE * confidence
            + WEIGHT_DRUG_STAGE * drug_stage
            + WEIGHT_MULTITARGET * multitarget
            + WEIGHT_SIMILARITY * similarity
        )

        results.append({
            "molecule_chembl_id": cid,
            "pref_name": row.get("pref_name"),
            "label": row.get("label"),
            "source": row.get("source"),
            "best_value": row.get("best_value"),
            "best_target": row.get("best_target"),
            "score_potency": round(potency, 3),
            "score_target": round(target_score, 3),
            "score_confidence": round(confidence, 3),
            "score_drug_stage": round(drug_stage, 3),
            "score_multitarget": round(multitarget, 3),
            "score_similarity": round(similarity, 3),
            "nearest_active": nearest,
            "composite_score": round(composite, 4),
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("composite_score", ascending=False)

    # Print top 25
    print(f"\n  TOP 25 COMPOUNDS BY V2 COMPOSITE SCORE:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'Potency':>7s}  {'Target':>6s}  {'Simil':>5s}  {'Stage':>5s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*5}  {'-'*5}  {'-'*8}  {'-'*30}")

    for rank, (_, row) in enumerate(results_df.head(25).iterrows(), 1):
        name = row["pref_name"] if pd.notna(row["pref_name"]) else row["molecule_chembl_id"]
        if len(str(name)) > 30:
            name = str(name)[:27] + "..."
        print(f"  {rank:4d}  {row['composite_score']:6.4f}  {row['score_potency']:7.3f}  "
              f"{row['score_target']:6.3f}  {row['score_similarity']:5.3f}  "
              f"{row['score_drug_stage']:5.3f}  {row['label']:>8s}  {name}")

    # Bottom 10
    print(f"\n  BOTTOM 10 COMPOUNDS:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'Simil':>5s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*5}  {'-'*8}  {'-'*30}")
    total = len(results_df)
    for rank, (_, row) in enumerate(results_df.tail(10).iterrows()):
        name = row["pref_name"] if pd.notna(row["pref_name"]) else row["molecule_chembl_id"]
        actual_rank = total - 9 + rank
        print(f"  {actual_rank:4d}  {row['composite_score']:6.4f}  {row['score_similarity']:5.3f}  "
              f"{row['label']:>8s}  {name}")

    return results_df


# ---------------------------------------------------------------------------
# EVALUATION
# ---------------------------------------------------------------------------

def evaluate_v2(results_df):
    """Evaluate v2 ranking against ground truth."""
    from sklearn.metrics import roc_auc_score, average_precision_score

    print("\n" + "=" * 60)
    print("V2 EVALUATION AGAINST GROUND TRUTH")
    print("=" * 60)

    # Binary: active vs inactive
    print("\n  TEST 1: Active vs Inactive separation")
    print("  " + "-" * 50)

    binary_df = results_df[results_df["label"].isin(["active", "inactive"])].copy()
    binary_labels = (binary_df["label"] == "active").astype(int).values
    binary_scores = binary_df["composite_score"].values

    n_active = sum(binary_labels == 1)
    n_inactive = sum(binary_labels == 0)
    print(f"  Actives: {n_active}, Inactives: {n_inactive}")

    auroc = None
    if n_active > 0 and n_inactive > 0:
        auroc = roc_auc_score(binary_labels, binary_scores)
        auprc = average_precision_score(binary_labels, binary_scores)
        print(f"  AUROC:  {auroc:.4f}  (v1 was 0.9855)")
        print(f"  AUPRC:  {auprc:.4f}")

    # Three-class ordering
    print(f"\n  TEST 2: Three-class score distributions")
    print("  " + "-" * 50)
    for label in ["active", "weak", "inactive"]:
        subset = results_df[results_df["label"] == label]["composite_score"]
        if len(subset) > 0:
            print(f"  {label:10s}  n={len(subset):3d}  "
                  f"mean={subset.mean():.4f}  median={subset.median():.4f}  "
                  f"min={subset.min():.4f}  max={subset.max():.4f}")

    active_mean = results_df[results_df["label"] == "active"]["composite_score"].mean()
    weak_mean = results_df[results_df["label"] == "weak"]["composite_score"].mean()
    inactive_mean = results_df[results_df["label"] == "inactive"]["composite_score"].mean()

    if active_mean > weak_mean > inactive_mean:
        print(f"\n  --> CORRECT ordering: active ({active_mean:.4f}) > weak ({weak_mean:.4f}) > inactive ({inactive_mean:.4f})")

    # Known drug rankings
    print(f"\n  TEST 3: Known drug rankings (CRITICAL TEST)")
    print("  " + "-" * 50)

    ranked = results_df.reset_index(drop=True)
    ranked["rank"] = range(1, len(ranked) + 1)
    total = len(ranked)

    known_drugs = ["PRAZIQUANTEL", "OXAMNIQUINE", "MEFLOQUINE", "ATOVAQUONE",
                   "PANOBINOSTAT", "VORINOSTAT", "IDEBENONE"]
    for drug_name in known_drugs:
        match = ranked[ranked["pref_name"] == drug_name]
        if len(match) > 0:
            row = match.iloc[0]
            percentile = (1 - row["rank"] / total) * 100
            sim = row.get("score_similarity", 0)
            print(f"  {drug_name:20s}  rank {int(row['rank']):3d}/{total}  "
                  f"score={row['composite_score']:.4f}  "
                  f"similarity={sim:.3f}  "
                  f"(top {percentile:.0f}%)")

    # Negative control check
    print(f"\n  TEST 4: Negative control analysis")
    print("  " + "-" * 50)
    negatives = results_df[results_df["source"] == "curated_negative"]
    if len(negatives) > 0:
        active_median = results_df[results_df["label"] == "active"]["composite_score"].median()
        neg_above = negatives[negatives["composite_score"] > active_median]
        print(f"  Negative controls: {len(negatives)}")
        print(f"  Active median score: {active_median:.4f}")
        print(f"  Negatives above active median: {len(neg_above)}")

    return auroc


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 05")
    print("  Structural Similarity + V2 Ranking")
    print("=" * 60)

    # Load data
    eval_df = pd.read_csv(os.path.join(EVAL_DIR, "evaluation_set_v1.csv"))
    eval_df["max_phase"] = pd.to_numeric(eval_df.get("max_phase", pd.Series()), errors="coerce")
    print(f"\n  Evaluation set: {len(eval_df)} compounds")

    activities_df = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))
    print(f"  Activity data: {len(activities_df)} measurements")

    # Step 1: Collect SMILES
    smiles_map = collect_smiles(eval_df, activities_df)

    # Step 2: Compute fingerprints
    fingerprints = compute_fingerprints(smiles_map)

    # Step 3: Compute similarity to known actives
    similarity_scores = compute_similarity_scores(eval_df, activities_df, fingerprints)

    # Step 4: Compute v2 composite scores
    results_df = compute_v2_scores(eval_df, activities_df, fingerprints, similarity_scores)

    # Step 5: Evaluate
    auroc = evaluate_v2(results_df)

    # Save
    results_df.to_csv(
        os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v2.csv"), index=False
    )

    report_path = os.path.join(PROCESSED_DIR, "kira_evaluation_report_v2.txt")
    with open(report_path, "w") as f:
        f.write(f"KIRA V2 EVALUATION REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"AUROC: {auroc:.4f}\n" if auroc else "AUROC: N/A\n")
        f.write(f"\nKey improvement: Added structural similarity signal.\n")
        f.write(f"Praziquantel now scored based on molecular structure\n")
        f.write(f"similarity to known active compounds.\n")

    print(f"\n{'=' * 60}")
    print(f"  DONE")
    print(f"  AUROC: {auroc:.4f}" if auroc else "  AUROC: N/A")
    print(f"")
    print(f"  Key question: Where does praziquantel rank now?")
    print(f"")
    print(f"  Saved: kira_ranked_candidates_v2.csv")
    print(f"         kira_evaluation_report_v2.txt")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 05_structural_similarity.py successfully."
echo "Now run: python 05_structural_similarity.py"
