"""
Kira - Script 06: Whole-Organism Activity Data + V3 Ranking
=============================================================

Scripts 04 and 05 revealed that praziquantel ranks low because its
mechanism is not captured in target-based IC50 assays.

This script adds a third evidence pathway: whole-organism activity data.
These are experiments where researchers exposed live Schistosoma mansoni
worms to compounds and measured whether the worms died.

This directly captures drugs like praziquantel, oxamniquine, and
mefloquine that work but whose molecular targets are poorly characterized.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 06_whole_organism.py

Requires internet connection. May take a few minutes.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

# Import shared modules (canonical source for scoring data and drug SMILES)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from kira.scoring import WEIGHTS_V3, DRUG_STAGE_SCORES  # noqa: E402
from kira.targets import TARGET_ESSENTIALITY, DEFAULT_ESSENTIALITY  # noqa: E402
from kira.drugs import CURATED_SMILES  # noqa: E402

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")

ORGANISM_TARGET_ID = "CHEMBL612893"  # Schistosoma mansoni whole organism

# We'll accept multiple activity types for whole-organism assays
WHOLE_ORG_ACTIVITY_TYPES = [
    "IC50", "EC50", "AC50", "GI50", "LC50",
    "ED50", "CC50", "Potency",
]

ACTIVE_IC50_THRESHOLD = 1000  # nM for target-based actives

# V3 signal weights (from shared module)
WEIGHT_POTENCY = WEIGHTS_V3["potency"]
WEIGHT_TARGET = WEIGHTS_V3["target"]
WEIGHT_CONFIDENCE = WEIGHTS_V3["confidence"]
WEIGHT_DRUG_STAGE = WEIGHTS_V3["drug_stage"]
WEIGHT_MULTITARGET = WEIGHTS_V3["multitarget"]
WEIGHT_SIMILARITY = WEIGHTS_V3["similarity"]
WEIGHT_WHOLE_ORG = WEIGHTS_V3["whole_org"]


# ---------------------------------------------------------------------------
# STEP 1: Query whole-organism assay data from ChEMBL
# ---------------------------------------------------------------------------


def fetch_whole_organism_data():
    from chembl_webresource_client.new_client import new_client

    print()
    print('=' * 60)
    print('FETCHING WHOLE-ORGANISM ACTIVITY DATA')
    print('=' * 60)

    cache_path = os.path.join(PROCESSED_DIR, 'schisto_whole_organism_raw.csv')
    if os.path.exists(cache_path):
        import pandas as pd
        print(f'  Loading cached data from {cache_path}')
        df = pd.read_csv(cache_path)
        print(f'  Loaded {len(df)} cached records')
        return df

    print('  Strategy: query per-molecule to avoid API timeout')
    print('  Looking up known drugs + eval set compounds...')

    activity_api = new_client.activity
    import pandas as pd

    eval_df = pd.read_csv(os.path.join(EVAL_DIR, 'evaluation_set_v1.csv'))
    compound_ids = eval_df['molecule_chembl_id'].unique().tolist()

    # Also add key known drugs explicitly
    key_drugs = ['CHEMBL976', 'CHEMBL1135', 'CHEMBL1594', 'CHEMBL1450',
                 'CHEMBL483254', 'CHEMBL98', 'CHEMBL252556']
    for d in key_drugs:
        if d not in compound_ids:
            compound_ids.append(d)

    all_records = []
    checked = 0

    for cid in compound_ids:
        checked += 1
        try:
            acts = activity_api.filter(
                molecule_chembl_id=cid,
                target_chembl_id=ORGANISM_TARGET_ID,
            ).only([
                'molecule_chembl_id', 'canonical_smiles',
                'standard_type', 'standard_value', 'standard_units',
                'standard_relation', 'pchembl_value',
                'assay_chembl_id', 'assay_type',
            ])
            records = list(acts)
            if records:
                all_records.extend(records)
                print(f'    {cid}: {len(records)} records')
        except Exception as e:
            pass

        if checked % 50 == 0:
            print(f'  Checked {checked}/{len(compound_ids)} compounds...')

    print(f'  Total records found: {len(all_records)}')

    if not all_records:
        print('  No whole-organism data found for eval set compounds.')
        return pd.DataFrame()

    df = pd.DataFrame(all_records)
    df.to_csv(cache_path, index=False)
    print(f'  Cached to {cache_path}')
    return df


def process_whole_organism_data(raw_df):
    """
    Filter whole-organism data for quality and extract usable signals.
    """
    print("\n" + "=" * 60)
    print("PROCESSING WHOLE-ORGANISM DATA")
    print("=" * 60)

    df = raw_df.copy()
    print(f"\n  Raw records: {len(df)}")

    # Show what activity types exist
    print(f"\n  Activity types present:")
    for atype, count in df["standard_type"].value_counts().head(15).items():
        print(f"    {atype}: {count}")

    # Show units
    print(f"\n  Units present:")
    for unit, count in df["standard_units"].value_counts().head(10).items():
        print(f"    {unit}: {count}")

    # Filter for dose-response measurements
    type_mask = df["standard_type"].isin(WHOLE_ORG_ACTIVITY_TYPES)
    df_filtered = df[type_mask].copy()
    print(f"\n  After filtering for dose-response types: {len(df_filtered)}")

    # Convert values to numeric
    df_filtered["standard_value"] = pd.to_numeric(
        df_filtered["standard_value"], errors="coerce"
    )
    df_filtered = df_filtered.dropna(subset=["standard_value"])
    print(f"  After removing non-numeric values: {len(df_filtered)}")

    # Normalize units to nM
    # Keep nM directly, convert uM to nM
    nm_mask = df_filtered["standard_units"] == "nM"
    um_mask = df_filtered["standard_units"] == "uM"
    ug_ml_mask = df_filtered["standard_units"] == "ug.mL-1"

    df_nm = df_filtered[nm_mask].copy()
    df_um = df_filtered[um_mask].copy()
    df_um["standard_value"] = df_um["standard_value"] * 1000  # uM -> nM
    df_um["standard_units"] = "nM"

    df_combined = pd.concat([df_nm, df_um], ignore_index=True)
    print(f"  After normalizing to nM (nM + uM converted): {len(df_combined)}")

    if len(df_combined) == 0:
        print("  WARNING: No usable dose-response data after filtering.")
        # Fall back to any numeric data
        print("  Attempting to use all numeric measurements...")
        df_combined = df_filtered.copy()

    # Get best value per compound
    best = (
        df_combined
        .groupby("molecule_chembl_id")
        .agg(
            best_whole_org_value=("standard_value", "min"),
            n_whole_org_measurements=("standard_value", "count"),
            whole_org_smiles=("canonical_smiles", "first"),
        )
        .reset_index()
    )

    print(f"\n  Unique compounds with whole-organism data: {len(best)}")

    # Show compounds with best activity
    top_compounds = best.nsmallest(15, "best_whole_org_value")
    print(f"\n  Most potent compounds (whole-organism):")
    for _, row in top_compounds.iterrows():
        print(f"    {row['molecule_chembl_id']:20s}  "
              f"best = {row['best_whole_org_value']:.0f} nM  "
              f"({int(row['n_whole_org_measurements'])} measurements)")

    return best


# ---------------------------------------------------------------------------
# STEP 3: Check if known drugs appear in whole-organism data
# ---------------------------------------------------------------------------

def check_known_drugs_in_whole_org(whole_org_df, eval_df):
    """
    Critical check: do praziquantel, oxamniquine, and mefloquine appear?
    """
    print("\n" + "=" * 60)
    print("CHECKING KNOWN DRUGS IN WHOLE-ORGANISM DATA")
    print("=" * 60)

    # Get ChEMBL IDs for known drugs from eval set
    known_drugs = {
        "PRAZIQUANTEL": "CHEMBL976",
        "OXAMNIQUINE": "CHEMBL1135",
        "MEFLOQUINE": "CHEMBL1594",
        "ATOVAQUONE": "CHEMBL1450",
    }

    for name, cid in known_drugs.items():
        match = whole_org_df[whole_org_df["molecule_chembl_id"] == cid]
        if len(match) > 0:
            row = match.iloc[0]
            print(f"  {name:20s}  FOUND  "
                  f"best = {row['best_whole_org_value']:.0f} nM  "
                  f"({int(row['n_whole_org_measurements'])} measurements)")
        else:
            print(f"  {name:20s}  NOT FOUND in whole-organism data")


# ---------------------------------------------------------------------------
# STEP 4: Merge all data and compute V3 scores
# ---------------------------------------------------------------------------

def compute_v3_scores(eval_df, activities_df, whole_org_df, fingerprints, similarity_scores):
    """
    V3 composite score: target-based + structural similarity + whole-organism.
    """
    print("\n" + "=" * 60)
    print("COMPUTING V3 COMPOSITE SCORES")
    print("=" * 60)

    # Multi-target scores
    target_counts = activities_df.groupby("molecule_chembl_id")["target_name"].nunique()
    multitarget_scores = {}
    for cid, n in target_counts.items():
        if n >= 3:
            multitarget_scores[cid] = 1.0
        elif n == 2:
            multitarget_scores[cid] = 0.7
        else:
            multitarget_scores[cid] = 0.3

    # Build whole-organism lookup
    whole_org_lookup = {}
    for _, row in whole_org_df.iterrows():
        whole_org_lookup[row["molecule_chembl_id"]] = {
            "value": row["best_whole_org_value"],
            "n_measurements": row["n_whole_org_measurements"],
        }

    results = []

    for _, row in eval_df.iterrows():
        cid = row["molecule_chembl_id"]

        # Signal 1: Target-based potency
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

        # Signal 6: Structural similarity
        if cid in similarity_scores:
            similarity = similarity_scores[cid]["max_similarity"]
        else:
            similarity = 0.0

        # Signal 7: Whole-organism activity (NEW)
        if cid in whole_org_lookup:
            wo_value = whole_org_lookup[cid]["value"]
            if wo_value > 0:
                wo_pic50 = -np.log10(wo_value * 1e-9)
                whole_org_score = max(0.0, min(1.0, (wo_pic50 - 4.0) / 5.0))
            else:
                whole_org_score = 0.0
        else:
            whole_org_score = 0.0

        # If we have whole-organism evidence, also boost target score
        if whole_org_score > 0 and target_score < 0.7:
            target_score = max(target_score, 0.70)

        # V3 Composite
        composite = (
            WEIGHT_POTENCY * potency
            + WEIGHT_TARGET * target_score
            + WEIGHT_CONFIDENCE * confidence
            + WEIGHT_DRUG_STAGE * drug_stage
            + WEIGHT_MULTITARGET * multitarget
            + WEIGHT_SIMILARITY * similarity
            + WEIGHT_WHOLE_ORG * whole_org_score
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
            "score_whole_org": round(whole_org_score, 3),
            "composite_score": round(composite, 4),
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("composite_score", ascending=False)

    # Top 25
    print(f"\n  TOP 25 BY V3 COMPOSITE SCORE:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'Potency':>7s}  {'Target':>6s}  {'Simil':>5s}  {'WOrg':>5s}  {'Stage':>5s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*5}  {'-'*5}  {'-'*5}  {'-'*8}  {'-'*30}")

    for rank, (_, r) in enumerate(results_df.head(25).iterrows(), 1):
        name = r["pref_name"] if pd.notna(r["pref_name"]) else r["molecule_chembl_id"]
        if len(str(name)) > 30:
            name = str(name)[:27] + "..."
        print(f"  {rank:4d}  {r['composite_score']:6.4f}  {r['score_potency']:7.3f}  "
              f"{r['score_target']:6.3f}  {r['score_similarity']:5.3f}  "
              f"{r['score_whole_org']:5.3f}  {r['score_drug_stage']:5.3f}  "
              f"{r['label']:>8s}  {name}")

    # Bottom 10
    print(f"\n  BOTTOM 10:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'WOrg':>5s}  {'Simil':>5s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*5}  {'-'*5}  {'-'*8}  {'-'*30}")
    total = len(results_df)
    for rank, (_, r) in enumerate(results_df.tail(10).iterrows()):
        name = r["pref_name"] if pd.notna(r["pref_name"]) else r["molecule_chembl_id"]
        actual_rank = total - 9 + rank
        print(f"  {actual_rank:4d}  {r['composite_score']:6.4f}  {r['score_whole_org']:5.3f}  "
              f"{r['score_similarity']:5.3f}  {r['label']:>8s}  {name}")

    return results_df


# ---------------------------------------------------------------------------
# EVALUATION
# ---------------------------------------------------------------------------

def evaluate_v3(results_df):
    """Evaluate v3 ranking."""
    from sklearn.metrics import roc_auc_score, average_precision_score

    print("\n" + "=" * 60)
    print("V3 EVALUATION")
    print("=" * 60)

    binary_df = results_df[results_df["label"].isin(["active", "inactive"])].copy()
    binary_labels = (binary_df["label"] == "active").astype(int).values
    binary_scores = binary_df["composite_score"].values

    n_active = sum(binary_labels == 1)
    n_inactive = sum(binary_labels == 0)

    auroc = None
    if n_active > 0 and n_inactive > 0:
        auroc = roc_auc_score(binary_labels, binary_scores)
        auprc = average_precision_score(binary_labels, binary_scores)
        print(f"\n  AUROC:  {auroc:.4f}  (v1: 0.9855, v2: 0.9989)")
        print(f"  AUPRC:  {auprc:.4f}")

    # Three-class
    print(f"\n  Three-class distributions:")
    for label in ["active", "weak", "inactive"]:
        subset = results_df[results_df["label"] == label]["composite_score"]
        if len(subset) > 0:
            print(f"    {label:10s}  n={len(subset):3d}  "
                  f"mean={subset.mean():.4f}  median={subset.median():.4f}")

    # Known drugs - THE CRITICAL TEST
    print(f"\n  KNOWN DRUG RANKINGS (v1 -> v2 -> v3):")
    print("  " + "-" * 60)

    ranked = results_df.reset_index(drop=True)
    ranked["rank"] = range(1, len(ranked) + 1)
    total = len(ranked)

    v1_ranks = {"PRAZIQUANTEL": 193, "OXAMNIQUINE": 192, "MEFLOQUINE": 191,
                "ATOVAQUONE": 4, "PANOBINOSTAT": 5, "VORINOSTAT": 2}
    v2_ranks = {"PRAZIQUANTEL": 132, "OXAMNIQUINE": 169, "MEFLOQUINE": 170,
                "ATOVAQUONE": 3, "PANOBINOSTAT": 4, "VORINOSTAT": 1}

    for drug_name in ["PRAZIQUANTEL", "OXAMNIQUINE", "MEFLOQUINE",
                      "ATOVAQUONE", "PANOBINOSTAT", "VORINOSTAT"]:
        match = ranked[ranked["pref_name"] == drug_name]
        if len(match) > 0:
            row = match.iloc[0]
            v3_rank = int(row["rank"])
            v1_r = v1_ranks.get(drug_name, "?")
            v2_r = v2_ranks.get(drug_name, "?")
            wo = row.get("score_whole_org", 0)
            print(f"  {drug_name:20s}  v1:{v1_r:>3}  v2:{v2_r:>3}  v3:{v3_rank:>3}  "
                  f"score={row['composite_score']:.4f}  whole_org={wo:.3f}")

    # Negative controls
    negatives = results_df[results_df["source"] == "curated_negative"]
    active_median = results_df[results_df["label"] == "active"]["composite_score"].median()
    neg_above = negatives[negatives["composite_score"] > active_median]
    print(f"\n  Negative controls above active median: {len(neg_above)}/{len(negatives)}")

    return auroc


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 06")
    print("  Whole-Organism Data + V3 Ranking")
    print("=" * 60)

    start_time = time.time()

    # Load existing data
    eval_df = pd.read_csv(os.path.join(EVAL_DIR, "evaluation_set_v1.csv"))
    eval_df["max_phase"] = pd.to_numeric(eval_df.get("max_phase", pd.Series()), errors="coerce")
    activities_df = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))

    print(f"\n  Evaluation set: {len(eval_df)} compounds")
    print(f"  Target-based activities: {len(activities_df)} measurements")

    # Step 1: Fetch whole-organism data
    raw_whole_org = fetch_whole_organism_data()

    if len(raw_whole_org) == 0:
        print("No whole-organism data available. Exiting.")
        sys.exit(1)

    # Step 2: Process it
    whole_org_df = process_whole_organism_data(raw_whole_org)

    # Step 3: Check known drugs
    check_known_drugs_in_whole_org(whole_org_df, eval_df)

    # Step 4: Recompute fingerprints and similarity (from Script 05)
    print("\n  Recomputing fingerprints and similarity...")
    smiles_map = {}
    if "canonical_smiles" in activities_df.columns:
        for _, row in activities_df.iterrows():
            cid = row["molecule_chembl_id"]
            smi = row.get("canonical_smiles")
            if pd.notna(smi) and cid not in smiles_map:
                smiles_map[cid] = smi

    # Add whole-organism SMILES
    if "whole_org_smiles" in whole_org_df.columns:
        for _, row in whole_org_df.iterrows():
            cid = row["molecule_chembl_id"]
            smi = row.get("whole_org_smiles")
            if pd.notna(smi) and cid not in smiles_map:
                smiles_map[cid] = smi

    # Add curated SMILES
    for _, row in eval_df.iterrows():
        cid = row["molecule_chembl_id"]
        name = row.get("pref_name", "")
        if cid not in smiles_map and pd.notna(name) and name in CURATED_SMILES:
            smiles_map[cid] = CURATED_SMILES[name]

    fingerprints = {}
    for cid, smi in smiles_map.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fingerprints[cid] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

    # Similarity to target-based actives
    potent_ids = set(
        activities_df[activities_df["standard_value"] < ACTIVE_IC50_THRESHOLD]["molecule_chembl_id"].unique()
    )
    ref_fps = {cid: fp for cid, fp in fingerprints.items() if cid in potent_ids}
    ref_fps_list = list(ref_fps.values())
    ref_ids_list = list(ref_fps.keys())

    similarity_scores = {}
    if ref_fps_list:
        for cid, fp in fingerprints.items():
            sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps_list)
            max_sim = max(sims)
            similarity_scores[cid] = {"max_similarity": max_sim}

    print(f"  Fingerprints: {len(fingerprints)}, Similarity scores: {len(similarity_scores)}")

    # Step 5: Compute v3 scores
    results_df = compute_v3_scores(eval_df, activities_df, whole_org_df,
                                    fingerprints, similarity_scores)

    # Step 6: Evaluate
    auroc = evaluate_v3(results_df)

    # Save
    results_df.to_csv(
        os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv"), index=False
    )

    elapsed = time.time() - start_time
    print(f"\n{'=' * 60}")
    print(f"  DONE in {elapsed:.0f} seconds")
    print(f"  AUROC: {auroc:.4f}" if auroc else "  AUROC: N/A")
    print(f"")
    print(f"  Three evidence pathways now active:")
    print(f"    1. Target-based IC50 (ChEMBL)")
    print(f"    2. Structural similarity (RDKit)")
    print(f"    3. Whole-organism activity (ChEMBL)")
    print(f"")
    print(f"  Saved: kira_ranked_candidates_v3.csv")
    print(f"{'=' * 60}")
