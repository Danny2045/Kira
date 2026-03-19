#!/bin/bash
# Run this from ~/kira with: bash create_script03.sh
# It creates 03_build_eval_set.py for you.

cat > 03_build_eval_set.py << 'PYTHONSCRIPT'
"""
Kira - Script 03: Build Ground-Truth Evaluation Set
=====================================================

This script creates the dataset that lets us measure whether Kira works.

It builds a labeled set of compounds:
  - "active"   = known to potently inhibit a parasite target (IC50 < 1000 nM)
  - "weak"     = measurable but weak activity (1000-10000 nM)
  - "inactive" = tested but not active, or unrelated approved drugs

Later, when we build a ranking algorithm, we test it against this set.
If it ranks actives above inactives, the algorithm has learned something real.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 03_build_eval_set.py

Requires: Script 02 must have been run first (needs data/processed/ files).
"""

import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")

# Activity thresholds (in nanomolar)
ACTIVE_THRESHOLD = 1000      # IC50 < 1000 nM = "active"
WEAK_THRESHOLD = 10000       # 1000 <= IC50 < 10000 = "weak"
                             # IC50 >= 10000 = "inactive"

# ---------------------------------------------------------------------------
# STEP 1: Load ChEMBL activity data from Script 02
# ---------------------------------------------------------------------------

def load_chembl_data():
    """Load the filtered activity data and repurposing candidates from Script 02."""

    print("\n" + "=" * 60)
    print("LOADING DATA FROM SCRIPT 02")
    print("=" * 60)

    # Load filtered activities (has activity values)
    filtered_path = os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv")
    if not os.path.exists(filtered_path):
        print(f"ERROR: {filtered_path} not found.")
        print("Run Script 02 first: python 02_query_chembl.py")
        sys.exit(1)

    filtered_df = pd.read_csv(filtered_path)
    print(f"\n  Loaded {len(filtered_df)} filtered activity measurements")

    # Load repurposing candidates (has drug approval status)
    candidates_path = os.path.join(PROCESSED_DIR, "schisto_repurposing_candidates.csv")
    if os.path.exists(candidates_path):
        candidates_df = pd.read_csv(candidates_path)
        # Fix max_phase type
        if "max_phase" in candidates_df.columns:
            candidates_df["max_phase"] = pd.to_numeric(
                candidates_df["max_phase"], errors="coerce"
            )
        print(f"  Loaded {len(candidates_df)} repurposing candidate records")
    else:
        candidates_df = filtered_df.copy()
        print("  No repurposing candidates file found, using filtered data only")

    return filtered_df, candidates_df


# ---------------------------------------------------------------------------
# STEP 2: Label compounds by activity level
# ---------------------------------------------------------------------------

def label_compounds(filtered_df, candidates_df):
    """
    Assign each compound a label based on its measured activity.

    For compounds tested against multiple targets or in multiple assays,
    we use the BEST (lowest) IC50 value. This is standard practice because
    a compound only needs to hit one target effectively to be interesting.
    """

    print("\n" + "=" * 60)
    print("LABELING COMPOUNDS BY ACTIVITY")
    print("=" * 60)

    # Get best (lowest) activity value per compound
    best_activity = (
        filtered_df
        .groupby("molecule_chembl_id")
        .agg(
            best_value=("standard_value", "min"),
            best_target=("target_name", "first"),
            n_measurements=("standard_value", "count"),
        )
        .reset_index()
    )

    # Merge in drug names and approval status
    if "pref_name" in candidates_df.columns and "max_phase" in candidates_df.columns:
        drug_info = (
            candidates_df[["molecule_chembl_id", "pref_name", "max_phase"]]
            .drop_duplicates(subset=["molecule_chembl_id"])
        )
        best_activity = best_activity.merge(drug_info, on="molecule_chembl_id", how="left")

    # Assign labels
    def assign_label(value):
        if value < ACTIVE_THRESHOLD:
            return "active"
        elif value < WEAK_THRESHOLD:
            return "weak"
        else:
            return "inactive"

    best_activity["label"] = best_activity["best_value"].apply(assign_label)

    # Summary
    print(f"\n  Label distribution:")
    for label, count in best_activity["label"].value_counts().items():
        print(f"    {label:10s}: {count} compounds")

    print(f"\n  Active compounds (IC50 < {ACTIVE_THRESHOLD} nM):")
    actives = best_activity[best_activity["label"] == "active"].sort_values("best_value")
    for _, row in actives.iterrows():
        name = row.get("pref_name", None)
        if pd.isna(name) or name is None:
            name = row["molecule_chembl_id"]
        phase = row.get("max_phase", None)
        phase_str = ""
        if pd.notna(phase) and phase > 0:
            phase_labels = {4: "APPROVED", 3: "Phase III", 2: "Phase II", 1: "Phase I"}
            phase_str = f" [{phase_labels.get(phase, '')}]"
        print(f"    {name:40s} IC50 = {row['best_value']:8.1f} nM  vs {row['best_target']}{phase_str}")

    return best_activity


# ---------------------------------------------------------------------------
# STEP 3: Add known drugs as curated entries
# ---------------------------------------------------------------------------

def add_curated_entries(eval_df):
    """
    Add manually curated entries for known schistosomiasis drugs.

    These are drugs where we KNOW the clinical outcome but the mechanism
    may not be fully captured in ChEMBL's target-based assays.

    This is domain knowledge — the kind of information that makes the
    evaluation set trustworthy beyond what any database query can provide.
    """

    print("\n" + "=" * 60)
    print("ADDING CURATED ENTRIES")
    print("=" * 60)

    curated = [
        {
            "molecule_chembl_id": "CHEMBL976",
            "pref_name": "PRAZIQUANTEL",
            "best_value": None,  # No target-based IC50 available
            "best_target": "Calcium channels (mechanism unclear)",
            "n_measurements": 0,
            "max_phase": 4.0,
            "label": "active",
            "source": "curated",
            "notes": "Current standard of care. Mechanism involves calcium "
                     "influx causing worm paralysis. No target-based IC50 in "
                     "ChEMBL but clinically proven effective.",
        },
        {
            "molecule_chembl_id": "CHEMBL1135",
            "pref_name": "OXAMNIQUINE",
            "best_value": None,
            "best_target": "Sulfotransferase (prodrug activation)",
            "n_measurements": 0,
            "max_phase": 4.0,
            "label": "active",
            "source": "curated",
            "notes": "Older drug effective against S. mansoni only. Prodrug "
                     "activated by parasite sulfotransferase. Largely replaced "
                     "by praziquantel but still a known active.",
        },
        {
            "molecule_chembl_id": "CHEMBL1594",
            "pref_name": "MEFLOQUINE",
            "best_value": None,
            "best_target": "Unknown (whole-worm activity)",
            "n_measurements": 0,
            "max_phase": 4.0,
            "label": "active",
            "source": "curated",
            "notes": "Antimalarial with demonstrated in vitro and in vivo "
                     "activity against juvenile and adult schistosomes. "
                     "Explored as combination partner with praziquantel.",
        },
        {
            "molecule_chembl_id": "CHEMBL1450",
            "pref_name": "ATOVAQUONE",
            "best_value": 430.0,
            "best_target": "Dihydroorotate dehydrogenase (quinone), mitochondrial",
            "n_measurements": 2,
            "max_phase": 4.0,
            "label": "active",
            "source": "chembl+curated",
            "notes": "Approved antiparasitic (malaria). IC50 430 nM against "
                     "SmDHODH. Available in sub-Saharan Africa. Strong "
                     "repurposing candidate.",
        },
    ]

    curated_df = pd.DataFrame(curated)

    # Mark existing entries as "chembl" source
    if "source" not in eval_df.columns:
        eval_df["source"] = "chembl"
    if "notes" not in eval_df.columns:
        eval_df["notes"] = ""

    # Remove any duplicates (if atovaquone is already in eval_df from ChEMBL)
    existing_ids = set(eval_df["molecule_chembl_id"].values)
    new_curated = curated_df[~curated_df["molecule_chembl_id"].isin(existing_ids)]

    # Update existing entries that overlap (like atovaquone)
    for _, curated_row in curated_df.iterrows():
        if curated_row["molecule_chembl_id"] in existing_ids:
            mask = eval_df["molecule_chembl_id"] == curated_row["molecule_chembl_id"]
            eval_df.loc[mask, "source"] = "chembl+curated"
            eval_df.loc[mask, "notes"] = curated_row["notes"]
            print(f"  Updated existing entry: {curated_row['pref_name']}")

    if len(new_curated) > 0:
        eval_df = pd.concat([eval_df, new_curated], ignore_index=True)
        for _, row in new_curated.iterrows():
            print(f"  Added curated entry: {row['pref_name']} (label={row['label']})")

    return eval_df


# ---------------------------------------------------------------------------
# STEP 4: Add negative controls (assumed-inactive approved drugs)
# ---------------------------------------------------------------------------

def add_negative_controls(eval_df):
    """
    Add approved drugs with NO known connection to schistosomiasis as
    negative controls.

    These are drugs used for completely unrelated conditions (hypertension,
    depression, diabetes, etc.) that have no biological reason to be active
    against schistosome targets. We assume they are inactive.

    Why this matters: an evaluation set with only positives cannot measure
    false positive rate. If your ranking algorithm puts everything high,
    you need negatives to expose that failure.
    """
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("ADDING NEGATIVE CONTROLS")
    print("=" * 60)

    # Well-known drugs for unrelated conditions that we can reasonably
    # assume are inactive against schistosome targets
    negative_drugs = [
        {"name": "METFORMIN", "id": "CHEMBL1431", "indication": "Type 2 diabetes"},
        {"name": "LISINOPRIL", "id": "CHEMBL1237", "indication": "Hypertension"},
        {"name": "SERTRALINE", "id": "CHEMBL809", "indication": "Depression"},
        {"name": "OMEPRAZOLE", "id": "CHEMBL1503", "indication": "Acid reflux"},
        {"name": "ATORVASTATIN", "id": "CHEMBL1487", "indication": "High cholesterol"},
        {"name": "AMLODIPINE", "id": "CHEMBL1491", "indication": "Hypertension"},
        {"name": "LEVOTHYROXINE", "id": "CHEMBL1365", "indication": "Hypothyroidism"},
        {"name": "MONTELUKAST", "id": "CHEMBL787", "indication": "Asthma"},
        {"name": "CLOPIDOGREL", "id": "CHEMBL1771", "indication": "Blood clotting"},
        {"name": "TAMSULOSIN", "id": "CHEMBL953", "indication": "Benign prostatic hyperplasia"},
        {"name": "ESCITALOPRAM", "id": "CHEMBL1508", "indication": "Depression/Anxiety"},
        {"name": "GABAPENTIN", "id": "CHEMBL940", "indication": "Neuropathic pain"},
        {"name": "PANTOPRAZOLE", "id": "CHEMBL1502", "indication": "Acid reflux"},
        {"name": "ROSUVASTATIN", "id": "CHEMBL1496", "indication": "High cholesterol"},
        {"name": "TRAMADOL", "id": "CHEMBL1232048", "indication": "Pain"},
    ]

    existing_ids = set(eval_df["molecule_chembl_id"].values)
    negatives = []

    for drug in negative_drugs:
        if drug["id"] not in existing_ids:
            negatives.append({
                "molecule_chembl_id": drug["id"],
                "pref_name": drug["name"],
                "best_value": None,
                "best_target": "None (negative control)",
                "n_measurements": 0,
                "max_phase": 4.0,
                "label": "inactive",
                "source": "curated_negative",
                "notes": f"Approved for {drug['indication']}. No known "
                         f"schistosomiasis activity. Assumed negative control.",
            })

    neg_df = pd.DataFrame(negatives)
    eval_df = pd.concat([eval_df, neg_df], ignore_index=True)

    print(f"  Added {len(negatives)} negative control drugs")
    print(f"  These are approved drugs for unrelated conditions")
    print(f"  (diabetes, hypertension, depression, etc.)")

    return eval_df


# ---------------------------------------------------------------------------
# STEP 5: Final summary and save
# ---------------------------------------------------------------------------

def save_evaluation_set(eval_df):
    """Save the evaluation set and print a comprehensive summary."""

    os.makedirs(EVAL_DIR, exist_ok=True)

    print("\n" + "=" * 60)
    print("EVALUATION SET SUMMARY")
    print("=" * 60)

    # Overall stats
    print(f"\n  Total compounds: {len(eval_df)}")
    print(f"\n  By label:")
    for label in ["active", "weak", "inactive"]:
        subset = eval_df[eval_df["label"] == label]
        print(f"    {label:10s}: {len(subset):3d} compounds")

    print(f"\n  By source:")
    for source, count in eval_df["source"].value_counts().items():
        print(f"    {source:20s}: {count}")

    # Show all actives
    print(f"\n  {'=' * 50}")
    print(f"  ALL ACTIVE COMPOUNDS")
    print(f"  {'=' * 50}")
    actives = eval_df[eval_df["label"] == "active"].sort_values(
        "best_value", na_position="last"
    )
    for _, row in actives.iterrows():
        name = row.get("pref_name", row["molecule_chembl_id"])
        if pd.isna(name):
            name = row["molecule_chembl_id"]
        value_str = f"IC50 = {row['best_value']:.0f} nM" if pd.notna(row["best_value"]) else "clinical evidence"
        target = row.get("best_target", "unknown")
        print(f"    {name}")
        print(f"      {value_str} vs {target}")
        if pd.notna(row.get("notes", None)) and row["notes"]:
            # Print first 80 chars of notes
            note_preview = row["notes"][:100]
            print(f"      Note: {note_preview}")
        print()

    # Save
    eval_path = os.path.join(EVAL_DIR, "evaluation_set_v1.csv")
    eval_df.to_csv(eval_path, index=False)
    print(f"  Saved to: {eval_path}")

    # Also save a clean version with just the key columns
    clean_cols = ["molecule_chembl_id", "pref_name", "label", "best_value",
                  "best_target", "source", "notes"]
    clean_cols = [c for c in clean_cols if c in eval_df.columns]
    clean_df = eval_df[clean_cols].copy()
    clean_path = os.path.join(EVAL_DIR, "evaluation_set_v1_clean.csv")
    clean_df.to_csv(clean_path, index=False)
    print(f"  Clean version: {clean_path}")

    # Save a text report
    report_lines = [
        "=" * 60,
        "KIRA EVALUATION SET v1",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 60,
        "",
        f"Total compounds: {len(eval_df)}",
        f"  Active:   {len(eval_df[eval_df['label'] == 'active'])}",
        f"  Weak:     {len(eval_df[eval_df['label'] == 'weak'])}",
        f"  Inactive: {len(eval_df[eval_df['label'] == 'inactive'])}",
        "",
        "PURPOSE:",
        "  This set defines ground truth for evaluating Kira's ranking algorithm.",
        "  A correct ranking places actives above weak, and weak above inactive.",
        "",
        "METHODOLOGY:",
        f"  - Active: IC50 < {ACTIVE_THRESHOLD} nM against S. mansoni targets,",
        "    plus curated known effective drugs",
        f"  - Weak: {ACTIVE_THRESHOLD} <= IC50 < {WEAK_THRESHOLD} nM",
        "  - Inactive: curated approved drugs for unrelated conditions",
        "",
        "LIMITATIONS:",
        "  - Negative controls are ASSUMED inactive (not experimentally verified)",
        "  - Some actives have clinical evidence but no target-based IC50",
        "  - Evaluation set is small; results should be interpreted with caution",
        "  - Does not include whole-organism assay data (future improvement)",
    ]
    report_path = os.path.join(EVAL_DIR, "evaluation_set_v1_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    print(f"  Report: {report_path}")

    return eval_df


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 03")
    print("  Building Ground-Truth Evaluation Set")
    print("=" * 60)

    # Step 1: Load data from Script 02
    filtered_df, candidates_df = load_chembl_data()

    # Step 2: Label compounds by activity
    eval_df = label_compounds(filtered_df, candidates_df)

    # Step 3: Add curated known drugs
    eval_df = add_curated_entries(eval_df)

    # Step 4: Add negative controls
    eval_df = add_negative_controls(eval_df)

    # Step 5: Save and summarize
    eval_df = save_evaluation_set(eval_df)

    print(f"\n{'=' * 60}")
    print(f"  DONE")
    print(f"")
    print(f"  Your evaluation set is ready at:")
    print(f"    data/eval/evaluation_set_v1.csv")
    print(f"")
    print(f"  This is Kira's ground truth. Every future ranking")
    print(f"  algorithm will be tested against this set.")
    print(f"")
    print(f"  Next: 04_build_ranking.py -- build the composite")
    print(f"  scoring function that ranks repurposing candidates")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 03_build_eval_set.py successfully."
echo "Now run: python 03_build_eval_set.py"
