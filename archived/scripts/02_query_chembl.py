"""
Kira - Script 02: Query ChEMBL for Schistosoma mansoni Targets
================================================================

This script queries ChEMBL to find:
1. All protein targets from Schistosoma mansoni (the parasite)
2. Compounds tested against those targets
3. How potently those compounds inhibit each target
4. Which compounds are already approved drugs

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 02_query_chembl.py
"""

import os
import sys
import time
import pandas as pd
from datetime import datetime

ORGANISM = "Schistosoma mansoni"
ACTIVITY_TYPES = ["IC50", "EC50", "Ki", "Kd", "AC50"]
ACTIVITY_THRESHOLD_NM = 10000
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")


def find_parasite_targets(organism):
    from chembl_webresource_client.new_client import new_client
    print(f"\nSearching ChEMBL for targets from: {organism}")
    print("-" * 50)
    target_api = new_client.target
    targets = target_api.filter(organism=organism).only([
        "target_chembl_id", "pref_name", "target_type", "organism",
    ])
    target_list = list(targets)
    # Filter out ORGANISM and NUCLEIC-ACID targets (too broad, causes API timeout)
    target_list = [t for t in target_list if t.get('target_type') == 'SINGLE PROTEIN']
    if not target_list:
        print(f"WARNING: No targets found for '{organism}'")
        return pd.DataFrame()
    df = pd.DataFrame(target_list)
    print(f"\nFound {len(df)} targets from {organism}")
    print(f"\nTarget types:")
    for ttype, count in df["target_type"].value_counts().items():
        print(f"  {ttype}: {count}")
    print(f"\nAll targets:")
    for _, row in df.iterrows():
        name = row["pref_name"] if pd.notna(row["pref_name"]) else "(unnamed)"
        print(f"  {row['target_chembl_id']:15s}  {row['target_type']:25s}  {name}")
    return df


def get_bioactivity_for_target(target_chembl_id, target_name):
    from chembl_webresource_client.new_client import new_client
    activity_api = new_client.activity
    activities = activity_api.filter(
        target_chembl_id=target_chembl_id,
        standard_type__in=ACTIVITY_TYPES,
    ).only([
        "molecule_chembl_id", "canonical_smiles", "standard_type",
        "standard_value", "standard_units", "standard_relation",
        "pchembl_value", "assay_chembl_id", "assay_type", "target_chembl_id",
    ])
    return list(activities)


def get_all_bioactivities(targets_df):
    print(f"\n{'=' * 60}")
    print("RETRIEVING BIOACTIVITY DATA")
    print(f"{'=' * 60}")
    all_activities = []
    for idx, row in targets_df.iterrows():
        target_id = row["target_chembl_id"]
        target_name = row["pref_name"] if pd.notna(row["pref_name"]) else "(unnamed)"
        print(f"\n  [{idx + 1}/{len(targets_df)}] {target_name} ({target_id})...")
        try:
            activities = get_bioactivity_for_target(target_id, target_name)
            print(f"    -> {len(activities)} activity measurements found")
            for act in activities:
                act["target_name"] = target_name
            all_activities.extend(activities)
        except Exception as e:
            print(f"    -> ERROR: {e}")
            continue
    if not all_activities:
        print("\nNo bioactivity data found for any target.")
        return pd.DataFrame()
    df = pd.DataFrame(all_activities)
    print(f"\nTotal raw activity measurements: {len(df):,}")
    return df


def filter_activities(activities_df):
    print(f"\n{'=' * 60}")
    print("FILTERING BIOACTIVITY DATA")
    print(f"{'=' * 60}")
    df = activities_df.copy()
    initial_count = len(df)
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    print(f"\n  After removing non-numeric values: {len(df):,} (dropped {initial_count - len(df):,})")
    df = df[df["standard_units"] == "nM"]
    print(f"  After keeping only nM units: {len(df):,}")
    if "standard_relation" in df.columns:
        exact_mask = (
            df["standard_relation"].isna() |
            (df["standard_relation"] == "=") |
            (df["standard_relation"] == "'='")
        )
        df = df[exact_mask]
        print(f"  After keeping only exact (=) measurements: {len(df):,}")
    df = df[df["standard_value"] <= ACTIVITY_THRESHOLD_NM]
    print(f"  After activity threshold (<= {ACTIVITY_THRESHOLD_NM} nM): {len(df):,}")
    if len(df) == 0:
        print("\n  WARNING: No measurements survived filtering.")
        return df
    df = df.sort_values("standard_value", ascending=True)
    print(f"\n  Measurements by target:")
    for target, group in df.groupby("target_name"):
        n_compounds = group["molecule_chembl_id"].nunique()
        best_value = group["standard_value"].min()
        print(f"    {target}")
        print(f"      {len(group)} measurements, {n_compounds} unique compounds")
        print(f"      Best (most potent) value: {best_value:.1f} nM")
    return df


def check_approved_drugs(activities_df):
    from chembl_webresource_client.new_client import new_client
    print(f"\n{'=' * 60}")
    print("CHECKING WHICH COMPOUNDS ARE APPROVED DRUGS")
    print(f"{'=' * 60}")
    molecule_api = new_client.molecule
    compound_ids = activities_df["molecule_chembl_id"].unique().tolist()
    print(f"\n  Checking {len(compound_ids)} unique compounds...")
    BATCH_SIZE = 50
    all_molecules = []
    for i in range(0, len(compound_ids), BATCH_SIZE):
        batch = compound_ids[i:i + BATCH_SIZE]
        try:
            mols = molecule_api.filter(
                molecule_chembl_id__in=batch
            ).only([
                "molecule_chembl_id", "pref_name", "max_phase",
                "molecule_type", "first_approval",
            ])
            all_molecules.extend(list(mols))
        except Exception as e:
            print(f"    Batch {i//BATCH_SIZE + 1} error: {e}")
            continue
        if (i // BATCH_SIZE + 1) % 5 == 0:
            print(f"    Processed {min(i + BATCH_SIZE, len(compound_ids))}/{len(compound_ids)} compounds...")
    if not all_molecules:
        print("  Could not retrieve molecule information.")
        return activities_df
    mol_df = pd.DataFrame(all_molecules)
    merged = activities_df.merge(
        mol_df[["molecule_chembl_id", "pref_name", "max_phase", "first_approval"]],
        on="molecule_chembl_id",
        how="left",
    )
    # Fix: ChEMBL sometimes returns max_phase as string
    merged["max_phase"] = pd.to_numeric(merged["max_phase"], errors="coerce")
    phase_labels = {
        4: "Approved drug", 3: "Phase III", 2: "Phase II",
        1: "Phase I", 0.5: "Early Phase I", 0: "Preclinical / Research",
    }
    if "max_phase" in merged.columns:
        print(f"\n  Compound development phases:")
        for phase in sorted(merged["max_phase"].dropna().unique(), reverse=True):
            n = merged[merged["max_phase"] == phase]["molecule_chembl_id"].nunique()
            label = phase_labels.get(phase, f"Phase {phase}")
            print(f"    {label} (phase {phase}): {n} compounds")
        approved = merged[merged["max_phase"] == 4].copy()
        if len(approved) > 0:
            print(f"\n  {'=' * 50}")
            print(f"  APPROVED DRUGS WITH ACTIVITY AGAINST PARASITE TARGETS")
            print(f"  {'=' * 50}")
            for compound_id in approved["molecule_chembl_id"].unique():
                compound_data = approved[approved["molecule_chembl_id"] == compound_id]
                name = compound_data["pref_name"].iloc[0]
                if pd.isna(name):
                    name = compound_id
                print(f"\n    {name} ({compound_id})")
                for _, row in compound_data.iterrows():
                    print(f"      vs {row['target_name']}: "
                          f"{row['standard_type']} = {row['standard_value']:.0f} nM")
        else:
            print(f"\n  No approved drugs found with activity against parasite targets")
            print(f"  at the current threshold (<= {ACTIVITY_THRESHOLD_NM} nM).")
            clinical = merged[merged["max_phase"] >= 1].copy()
            if len(clinical) > 0:
                print(f"\n  CLINICAL-STAGE COMPOUNDS:")
                for compound_id in clinical["molecule_chembl_id"].unique():
                    compound_data = clinical[clinical["molecule_chembl_id"] == compound_id]
                    name = compound_data["pref_name"].iloc[0]
                    phase = compound_data["max_phase"].iloc[0]
                    if pd.isna(name):
                        name = compound_id
                    label = phase_labels.get(phase, f"Phase {phase}")
                    print(f"\n    {name} ({compound_id}) -- {label}")
                    for _, row in compound_data.iterrows():
                        print(f"      vs {row['target_name']}: "
                              f"{row['standard_type']} = {row['standard_value']:.0f} nM")
    return merged


def generate_report(targets_df, activities_df, filtered_df, final_df):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if len(targets_df) > 0:
        targets_df.to_csv(os.path.join(OUTPUT_DIR, "schisto_parasite_targets.csv"), index=False)
    if len(filtered_df) > 0:
        filtered_df.to_csv(os.path.join(OUTPUT_DIR, "schisto_filtered_activities.csv"), index=False)
    if len(final_df) > 0:
        final_df.to_csv(os.path.join(OUTPUT_DIR, "schisto_repurposing_candidates.csv"), index=False)
    report_lines = [
        "=" * 60,
        "KIRA -- ChEMBL Schistosomiasis Target Report",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 60, "",
        f"Organism: {ORGANISM}",
        f"Targets found: {len(targets_df)}", "",
    ]
    if len(targets_df) > 0:
        report_lines.append("PARASITE TARGETS:")
        for _, row in targets_df.iterrows():
            name = row["pref_name"] if pd.notna(row["pref_name"]) else "(unnamed)"
            report_lines.append(f"  {row['target_chembl_id']:15s}  {name}")
    if len(filtered_df) > 0:
        report_lines.append("")
        report_lines.append(f"FILTERED ACTIVITY MEASUREMENTS: {len(filtered_df)}")
        report_lines.append(f"  Unique compounds: {filtered_df['molecule_chembl_id'].nunique()}")
    if len(final_df) > 0 and "max_phase" in final_df.columns:
        approved = final_df[final_df["max_phase"] == 4]
        if len(approved) > 0:
            report_lines.append("")
            report_lines.append("APPROVED DRUGS WITH PARASITE ACTIVITY:")
            for cid in approved["molecule_chembl_id"].unique():
                cdata = approved[approved["molecule_chembl_id"] == cid]
                name = cdata["pref_name"].iloc[0]
                if pd.isna(name):
                    name = cid
                report_lines.append(f"  {name}")
    report_text = "\n".join(report_lines)
    report_path = os.path.join(OUTPUT_DIR, "schisto_chembl_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"\n  Report saved to: {report_path}")
    print(f"  Data files saved to: {OUTPUT_DIR}/")


if __name__ == "__main__":
    print("=" * 60)
    print("  KIRA -- Phase 1, Script 02")
    print("  ChEMBL Query: Schistosoma mansoni Targets & Bioactivity")
    print("=" * 60)
    start_time = time.time()
    targets_df = find_parasite_targets(ORGANISM)
    if len(targets_df) == 0:
        print("\nNo targets found. Cannot continue.")
        sys.exit(1)
    activities_df = get_all_bioactivities(targets_df)
    if len(activities_df) == 0:
        print("\nNo bioactivity data found. Cannot continue.")
        generate_report(targets_df, activities_df, pd.DataFrame(), pd.DataFrame())
        sys.exit(1)
    filtered_df = filter_activities(activities_df)
    if len(filtered_df) == 0:
        print("\nNo activities passed filtering. Saving raw data for inspection.")
        activities_df.to_csv(
            os.path.join(OUTPUT_DIR, "schisto_raw_activities.csv"), index=False,
        )
        print(f"  Raw activities saved to: {OUTPUT_DIR}/schisto_raw_activities.csv")
        generate_report(targets_df, activities_df, filtered_df, pd.DataFrame())
        sys.exit(0)
    final_df = check_approved_drugs(filtered_df)
    generate_report(targets_df, activities_df, filtered_df, final_df)
    elapsed = time.time() - start_time
    print(f"\n{'=' * 60}")
    print(f"  DONE in {elapsed:.0f} seconds")
    print(f"")
    print(f"  Files saved to data/processed/:")
    print(f"    schisto_parasite_targets.csv      -- parasite protein targets")
    print(f"    schisto_filtered_activities.csv    -- quality-filtered activity data")
    print(f"    schisto_repurposing_candidates.csv -- activities + drug approval status")
    print(f"    schisto_chembl_report.txt          -- human-readable summary")
    print(f"")
    print(f"  Next: 03_build_evaluation_set.py -- create ground-truth test set")
    print(f"{'=' * 60}")
