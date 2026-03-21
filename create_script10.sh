#!/bin/bash
cat > 10_selectivity_analysis.py << 'PYTHONSCRIPT'
"""
Kira - Script 10: Selectivity Analysis Against Human Orthologues
==================================================================

The most critical translational question: does this compound hit the
PARASITE target without also hitting the HUMAN version?

For each parasite target with data (SmHDAC8, SmTGR, SmDHODH), we:
1. Identify the human orthologue in ChEMBL
2. Query activity data for our top compounds against the human target
3. Compute selectivity ratio: human_IC50 / parasite_IC50
   - Ratio > 10: SELECTIVE (compound prefers parasite target)
   - Ratio 3-10: MODERATE selectivity
   - Ratio 1-3: POOR selectivity (hits both similarly)
   - Ratio < 1: COUNTER-SELECTIVE (prefers human target — toxic)

A compound that is potent AND selective is a real drug candidate.
A compound that is potent but NOT selective is a research tool at best.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 10_selectivity_analysis.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
REPORT_DIR = os.path.join(os.path.dirname(__file__), "data", "reports")

# ---------------------------------------------------------------------------
# Human orthologues of S. mansoni drug targets in ChEMBL
# ---------------------------------------------------------------------------
# These were looked up manually from ChEMBL.
# Each parasite target has a human equivalent (orthologue).

ORTHOLOGUE_MAP = {
    "Histone deacetylase 8": {
        "parasite_id": "CHEMBL3797017",
        "human_name": "Histone deacetylase 8 (human)",
        "human_id": "CHEMBL3192",  # Human HDAC8
        "notes": "SmHDAC8 has structural differences in the active site "
                 "loop region compared to human HDAC8. Selective inhibition "
                 "is possible but requires careful compound design.",
    },
    "Thioredoxin glutathione reductase": {
        "parasite_id": "CHEMBL6110",
        "human_name": "Thioredoxin reductase 1 (human)",
        "human_id": "CHEMBL3952",  # Human TrxR1 (closest orthologue)
        "notes": "SmTGR is a fusion enzyme (TrxR + GR) unique to the parasite. "
                 "Humans have separate TrxR and GR enzymes. Selectivity is "
                 "theoretically favorable because the fusion creates a unique "
                 "active site architecture.",
    },
    "Dihydroorotate dehydrogenase (quinone), mitochondrial": {
        "parasite_id": "CHEMBL4523950",
        "human_name": "Dihydroorotate dehydrogenase (human)",
        "human_id": "CHEMBL1966",  # Human DHODH
        "notes": "Both parasite and human DHODH are mitochondrial. "
                 "Atovaquone shows some selectivity but this must be verified.",
    },
}


# ---------------------------------------------------------------------------
# STEP 1: Query human target activity for our compounds
# ---------------------------------------------------------------------------

def query_human_orthologues(activities_df):
    """
    For each compound in our parasite activity data, check if ChEMBL
    has activity data against the corresponding human orthologue.
    """
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("STEP 1: QUERYING HUMAN ORTHOLOGUE DATA")
    print("=" * 60)

    activity_api = new_client.activity
    all_human_data = []

    for parasite_target, info in ORTHOLOGUE_MAP.items():
        human_id = info["human_id"]
        human_name = info["human_name"]

        print(f"\n  {parasite_target}")
        print(f"  Human orthologue: {human_name} ({human_id})")

        # Get compounds tested against this parasite target
        parasite_compounds = set(
            activities_df[activities_df["target_name"] == parasite_target]["molecule_chembl_id"].unique()
        )
        print(f"  Compounds tested against parasite target: {len(parasite_compounds)}")

        if not parasite_compounds:
            continue

        # Query each compound against the human target
        found = 0
        for cid in parasite_compounds:
            try:
                acts = activity_api.filter(
                    molecule_chembl_id=cid,
                    target_chembl_id=human_id,
                    standard_type__in=["IC50", "EC50", "Ki", "Kd"],
                ).only([
                    "molecule_chembl_id", "standard_type",
                    "standard_value", "standard_units",
                    "standard_relation",
                ])
                records = list(acts)
                if records:
                    found += 1
                    for r in records:
                        r["parasite_target"] = parasite_target
                        r["human_target"] = human_name
                        r["human_target_id"] = human_id
                    all_human_data.extend(records)
            except Exception:
                pass

            time.sleep(0.35)  # Rate limit

        print(f"  Compounds with human orthologue data: {found}")

    if not all_human_data:
        print("\n  No human orthologue data found for any compounds.")
        return pd.DataFrame()

    human_df = pd.DataFrame(all_human_data)
    print(f"\n  Total human orthologue activity records: {len(human_df)}")

    # Cache
    cache_path = os.path.join(PROCESSED_DIR, "human_orthologue_activities.csv")
    human_df.to_csv(cache_path, index=False)
    print(f"  Cached to {cache_path}")

    return human_df


# ---------------------------------------------------------------------------
# STEP 2: Compute selectivity ratios
# ---------------------------------------------------------------------------

def compute_selectivity(activities_df, human_df):
    """
    For each compound with data against both parasite and human targets,
    compute the selectivity ratio: human_IC50 / parasite_IC50.

    Higher ratio = more selective for parasite = better drug candidate.
    """
    print("\n" + "=" * 60)
    print("STEP 2: COMPUTING SELECTIVITY RATIOS")
    print("=" * 60)

    if len(human_df) == 0:
        print("  No human data available.")
        return pd.DataFrame()

    # Clean human data
    human_clean = human_df.copy()
    human_clean["standard_value"] = pd.to_numeric(
        human_clean["standard_value"], errors="coerce"
    )
    human_clean = human_clean.dropna(subset=["standard_value"])

    # Normalize to nM
    um_mask = human_clean["standard_units"] == "uM"
    human_clean.loc[um_mask, "standard_value"] = human_clean.loc[um_mask, "standard_value"] * 1000

    # Keep only nM
    human_clean = human_clean[human_clean["standard_units"].isin(["nM", "uM"])]

    # Best human IC50 per compound per target
    human_best = (
        human_clean
        .groupby(["molecule_chembl_id", "parasite_target"])
        .agg(human_best_ic50=("standard_value", "min"),
             human_n_measurements=("standard_value", "count"))
        .reset_index()
    )

    # Best parasite IC50 per compound per target
    parasite_best = (
        activities_df
        .groupby(["molecule_chembl_id", "target_name"])
        .agg(parasite_best_ic50=("standard_value", "min"),
             parasite_n_measurements=("standard_value", "count"))
        .reset_index()
        .rename(columns={"target_name": "parasite_target"})
    )

    # Merge
    merged = parasite_best.merge(
        human_best,
        on=["molecule_chembl_id", "parasite_target"],
        how="inner",
    )

    if len(merged) == 0:
        print("  No compounds with both parasite and human data.")
        return pd.DataFrame()

    # Compute selectivity ratio
    merged["selectivity_ratio"] = merged["human_best_ic50"] / merged["parasite_best_ic50"]

    # Classify
    def classify_selectivity(ratio):
        if ratio >= 10:
            return "SELECTIVE"
        elif ratio >= 3:
            return "MODERATE"
        elif ratio >= 1:
            return "POOR"
        else:
            return "COUNTER-SELECTIVE"

    merged["selectivity_class"] = merged["selectivity_ratio"].apply(classify_selectivity)

    # Sort by selectivity
    merged = merged.sort_values("selectivity_ratio", ascending=False)

    # Get compound names
    v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv")
    if os.path.exists(v3_path):
        v3 = pd.read_csv(v3_path)
        name_map = dict(zip(v3["molecule_chembl_id"], v3.get("pref_name", pd.Series())))
        merged["pref_name"] = merged["molecule_chembl_id"].map(name_map)

    # Print results
    print(f"\n  Compounds with selectivity data: {len(merged)}")
    print(f"\n  Selectivity distribution:")
    for cls in ["SELECTIVE", "MODERATE", "POOR", "COUNTER-SELECTIVE"]:
        n = len(merged[merged["selectivity_class"] == cls])
        if n > 0:
            print(f"    {cls:20s}: {n}")

    print(f"\n  {'Compound':35s} {'Target':20s} {'Para IC50':>10s} {'Human IC50':>11s} {'Ratio':>7s} {'Class'}")
    print(f"  {'-'*35} {'-'*20} {'-'*10} {'-'*11} {'-'*7} {'-'*20}")

    for _, row in merged.iterrows():
        name = row.get("pref_name")
        if pd.isna(name):
            name = row["molecule_chembl_id"]
        if len(str(name)) > 35:
            name = str(name)[:32] + "..."

        target_short = row["parasite_target"][:20]
        print(f"  {name:35s} {target_short:20s} "
              f"{row['parasite_best_ic50']:10.0f} {row['human_best_ic50']:11.0f} "
              f"{row['selectivity_ratio']:7.1f} {row['selectivity_class']}")

    return merged


# ---------------------------------------------------------------------------
# STEP 3: Integrate selectivity into candidate assessment
# ---------------------------------------------------------------------------

def assess_candidates_with_selectivity(selectivity_df):
    """
    Provide translational assessment based on selectivity data.
    """
    print("\n" + "=" * 60)
    print("STEP 3: TRANSLATIONAL ASSESSMENT")
    print("=" * 60)

    if len(selectivity_df) == 0:
        print("  No selectivity data available for assessment.")
        return

    # By target
    for target in selectivity_df["parasite_target"].unique():
        target_data = selectivity_df[selectivity_df["parasite_target"] == target]
        info = ORTHOLOGUE_MAP.get(target, {})

        print(f"\n  TARGET: {target}")
        print(f"  Human orthologue: {info.get('human_name', 'unknown')}")
        print(f"  Biology: {info.get('notes', '')}")
        print(f"  Compounds with selectivity data: {len(target_data)}")

        selective = target_data[target_data["selectivity_class"] == "SELECTIVE"]
        moderate = target_data[target_data["selectivity_class"] == "MODERATE"]
        poor = target_data[target_data["selectivity_class"] == "POOR"]
        counter = target_data[target_data["selectivity_class"] == "COUNTER-SELECTIVE"]

        if len(selective) > 0:
            print(f"\n  SELECTIVE compounds (ratio >= 10x):")
            for _, row in selective.iterrows():
                name = row.get("pref_name", row["molecule_chembl_id"])
                print(f"    {name}: {row['selectivity_ratio']:.1f}x selective for parasite")
                print(f"      Parasite IC50: {row['parasite_best_ic50']:.0f} nM")
                print(f"      Human IC50: {row['human_best_ic50']:.0f} nM")

        if len(counter) > 0:
            print(f"\n  WARNING — COUNTER-SELECTIVE compounds (prefer human target):")
            for _, row in counter.iterrows():
                name = row.get("pref_name", row["molecule_chembl_id"])
                print(f"    {name}: {row['selectivity_ratio']:.2f}x (PREFERS HUMAN TARGET)")


# ---------------------------------------------------------------------------
# STEP 4: Save results
# ---------------------------------------------------------------------------

def save_selectivity_results(selectivity_df):
    """Save selectivity analysis."""

    os.makedirs(REPORT_DIR, exist_ok=True)

    if len(selectivity_df) > 0:
        selectivity_df.to_csv(
            os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv"),
            index=False,
        )

    report = []
    report.append("=" * 70)
    report.append("KIRA v1 — SELECTIVITY ANALYSIS REPORT")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("=" * 70)
    report.append("")
    report.append("QUESTION: Do our top compounds hit the parasite target")
    report.append("selectively, or do they also hit the human orthologue?")
    report.append("")

    if len(selectivity_df) > 0:
        for cls in ["SELECTIVE", "MODERATE", "POOR", "COUNTER-SELECTIVE"]:
            n = len(selectivity_df[selectivity_df["selectivity_class"] == cls])
            report.append(f"  {cls:20s}: {n} compounds")

        report.append("")
        report.append("DETAIL:")
        for _, row in selectivity_df.iterrows():
            name = row.get("pref_name", row["molecule_chembl_id"])
            report.append(f"  {name}")
            report.append(f"    vs {row['parasite_target']}: {row['parasite_best_ic50']:.0f} nM")
            report.append(f"    vs human orthologue: {row['human_best_ic50']:.0f} nM")
            report.append(f"    Selectivity: {row['selectivity_ratio']:.1f}x ({row['selectivity_class']})")
            report.append("")
    else:
        report.append("  No selectivity data was obtained.")
        report.append("  This means either:")
        report.append("    - Our compounds haven't been tested against human targets")
        report.append("    - The ChEMBL query found no matching data")
        report.append("  This is itself a finding: selectivity is UNKNOWN, not confirmed.")

    report_path = os.path.join(REPORT_DIR, "kira_v1_selectivity_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report))
    print(f"\n  Report saved to: {report_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 10")
    print("  Selectivity Analysis Against Human Orthologues")
    print("=" * 60)

    # Load parasite activity data
    activities_df = pd.read_csv(
        os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv")
    )
    print(f"\n  Parasite activity data: {len(activities_df)} measurements")
    print(f"  Unique compounds: {activities_df['molecule_chembl_id'].nunique()}")

    # Check for cached human data
    cache_path = os.path.join(PROCESSED_DIR, "human_orthologue_activities.csv")
    if os.path.exists(cache_path):
        print(f"\n  Loading cached human orthologue data...")
        human_df = pd.read_csv(cache_path)
        print(f"  Loaded {len(human_df)} records")
    else:
        # Step 1: Query human orthologues
        human_df = query_human_orthologues(activities_df)

    # Step 2: Compute selectivity
    selectivity_df = compute_selectivity(activities_df, human_df)

    # Step 3: Translational assessment
    assess_candidates_with_selectivity(selectivity_df)

    # Step 4: Save
    save_selectivity_results(selectivity_df)

    n_compounds = len(selectivity_df) if len(selectivity_df) > 0 else 0
    print(f"\n{'=' * 60}")
    print(f"  SELECTIVITY ANALYSIS COMPLETE")
    print(f"  Compounds with selectivity data: {n_compounds}")
    if n_compounds > 0:
        n_sel = len(selectivity_df[selectivity_df["selectivity_class"] == "SELECTIVE"])
        n_counter = len(selectivity_df[selectivity_df["selectivity_class"] == "COUNTER-SELECTIVE"])
        print(f"  Selective (>10x): {n_sel}")
        print(f"  Counter-selective (<1x): {n_counter}")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 10_selectivity_analysis.py successfully."
echo "Now run: python 10_selectivity_analysis.py"
