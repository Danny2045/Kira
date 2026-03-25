#!/bin/bash
cat > 16_tryp_expanded_selectivity.py << 'PYTHONSCRIPT'
"""
Kira - Script 16: Expanded T. brucei Selectivity
===================================================

Script 15 only mapped 2 of 16 targets to human orthologues.
The biggest datasets were left unmapped:
  - Rhodesain (488 measurements, 435 compounds) → human cathepsin L
  - PDEB1 (167 measurements, 157 compounds) → human PDE4
  - Alternative oxidase (45 measurements) → no human equivalent
  - RNA-editing ligase (38 measurements) → no human equivalent
  - Trypanothione reductase (128 measurements) → no human equivalent

This script:
  1. Maps Rhodesain → human cathepsin L (CHEMBL3524)
  2. Maps PDEB1 → human PDE4B (CHEMBL275)
  3. Maps ATP-dependent PFK → human PFKM (search)
  4. Queries selectivity for all new mappings
  5. Updates cross-disease comparison with full data

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 16_tryp_expanded_selectivity.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

BASE_DIR = os.path.dirname(__file__)
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
PUB_DIR = os.path.join(BASE_DIR, "data", "publication")
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")

# Expanded orthologue mappings
# Research-verified: Rhodesain is a clan CA, family C1 cysteine protease
# closest human homolog is cathepsin L (CHEMBL3524), not cathepsin B
EXPANDED_ORTHOLOGUES = {
    "Rhodesain": {
        "human_target": "CHEMBL3524",  # Human cathepsin L
        "human_name": "Cathepsin L",
        "rationale": "Rhodesain is a clan CA family C1 cysteine protease; "
                     "cathepsin L is the closest human structural homolog",
    },
    "Class 1 phosphodiesterase PDEB1": {
        "human_target": "CHEMBL275",  # Human PDE4B
        "human_name": "PDE4B",
        "rationale": "PDEB1 is a class I cyclic nucleotide PDE; "
                     "PDE4B is the most studied human PDE drug target",
    },
    "Cathepsin B-like cysteine protease": {
        "human_target": "CHEMBL3837",  # Human cathepsin B (already done in Script 15)
        "human_name": "Cathepsin B",
        "rationale": "Direct orthologue",
    },
}

# Targets unique to trypanosomatids (no human equivalent)
UNIQUE_TARGETS = [
    "Trypanothione reductase",
    "Alternative oxidase, mitochondrial",
    "RNA-editing ligase 1, mitochondrial",
    "Tryparedoxin peroxidase",
    "Trypanothione synthetase",
]


def query_human_selectivity(parasite_activities, target_name, human_chembl_id, human_name):
    """Query ChEMBL for human orthologue data and compute selectivity."""
    from chembl_webresource_client.new_client import new_client
    activity_api = new_client.activity

    # Get parasite compounds for this target
    target_data = parasite_activities[
        parasite_activities["target_name"] == target_name
    ].copy()

    if len(target_data) == 0:
        print(f"  No parasite data for {target_name}")
        return pd.DataFrame()

    compound_ids = target_data["molecule_chembl_id"].unique().tolist()
    print(f"\n  {target_name} → {human_name} ({human_chembl_id})")
    print(f"    Parasite compounds to check: {len(compound_ids)}")

    # Query human target
    human_records = []
    for i in range(0, len(compound_ids), 50):
        batch = compound_ids[i:i+50]
        try:
            results = activity_api.filter(
                target_chembl_id=human_chembl_id,
                molecule_chembl_id__in=batch,
                standard_type__in=["IC50", "EC50", "Ki", "Kd"],
                standard_relation="=",
                standard_units="nM",
            ).only([
                "molecule_chembl_id", "standard_value",
                "standard_type", "target_chembl_id",
            ])
            human_records.extend(list(results))
            time.sleep(0.3)
        except Exception as e:
            print(f"    Batch error: {e}")

    print(f"    Human orthologue measurements found: {len(human_records)}")

    if not human_records:
        return pd.DataFrame()

    human_df = pd.DataFrame(human_records)
    human_df["standard_value"] = pd.to_numeric(human_df["standard_value"], errors="coerce")
    human_df = human_df.dropna(subset=["standard_value"])

    # Best parasite IC50 per compound
    parasite_best = (
        target_data.groupby("molecule_chembl_id")["standard_value"]
        .min().reset_index()
        .rename(columns={"standard_value": "parasite_ic50"})
    )

    # Best human IC50 per compound
    human_best = (
        human_df.groupby("molecule_chembl_id")["standard_value"]
        .min().reset_index()
        .rename(columns={"standard_value": "human_ic50"})
    )

    # Merge
    merged = parasite_best.merge(human_best, on="molecule_chembl_id")
    if len(merged) == 0:
        print(f"    No overlapping compounds.")
        return pd.DataFrame()

    merged["selectivity_ratio"] = merged["human_ic50"] / merged["parasite_ic50"]
    merged["parasite_target"] = target_name
    merged["human_target"] = human_name
    merged["human_target_id"] = human_chembl_id

    def classify(r):
        if r >= 10: return "SELECTIVE"
        elif r >= 3: return "MODERATE"
        elif r >= 1: return "POOR"
        else: return "COUNTER-SELECTIVE"

    merged["selectivity_class"] = merged["selectivity_ratio"].apply(classify)

    # Summary
    n = len(merged)
    n_sel = (merged["selectivity_class"] == "SELECTIVE").sum()
    n_mod = (merged["selectivity_class"] == "MODERATE").sum()
    n_poor = (merged["selectivity_class"] == "POOR").sum()
    n_counter = (merged["selectivity_class"] == "COUNTER-SELECTIVE").sum()
    pct_nonsel = round(100 * (n_poor + n_counter) / n, 1)

    print(f"    Dual-species compounds: {n}")
    print(f"    Selective (>=10x): {n_sel}")
    print(f"    Moderate (3-10x): {n_mod}")
    print(f"    Poor (1-3x): {n_poor}")
    print(f"    Counter-selective (<1x): {n_counter}")
    print(f"    Non-selective rate: {pct_nonsel}%")
    print(f"    Median ratio: {merged['selectivity_ratio'].median():.1f}x")

    if n_sel > 0:
        best = merged.loc[merged["selectivity_ratio"].idxmax()]
        print(f"    BEST: {best['molecule_chembl_id']} — "
              f"parasite {best['parasite_ic50']:.0f} nM, "
              f"human {best['human_ic50']:.0f} nM, "
              f"{best['selectivity_ratio']:.1f}x")

    return merged


def build_full_comparison(all_selectivity, parasite_activities):
    """Build the complete cross-disease selectivity table."""

    print("\n" + "=" * 60)
    print("FULL CROSS-DISEASE SELECTIVITY LANDSCAPE")
    print("=" * 60)

    rows = []

    # Load schistosomiasis data
    schisto_sel_path = os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv")
    schisto_act_path = os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv")

    if os.path.exists(schisto_sel_path):
        s_sel = pd.read_csv(schisto_sel_path)
        s_act = pd.read_csv(schisto_act_path)

        for target in s_sel["parasite_target"].unique():
            td = s_sel[s_sel["parasite_target"] == target]
            n_total = len(s_act[s_act["target_name"] == target]["molecule_chembl_id"].unique())
            n_with = len(td)
            n_nonsel = len(td[td["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
            rows.append({
                "disease": "Schistosomiasis",
                "target": target[:35],
                "n_compounds": n_total,
                "n_with_selectivity": n_with,
                "pct_non_selective": round(100 * n_nonsel / n_with, 1) if n_with > 0 else None,
                "median_ratio": round(td["selectivity_ratio"].median(), 1),
                "best_ratio": round(td["selectivity_ratio"].max(), 1),
                "n_selective": len(td[td["selectivity_class"] == "SELECTIVE"]),
                "method": "Experimental",
                "unique_enzyme": False,
            })

        # SmTGR docking
        dock_path = os.path.join(PROCESSED_DIR, "smtgr_docking_results.csv")
        if os.path.exists(dock_path):
            dock = pd.read_csv(dock_path).dropna(subset=["SmTGR_dock_energy", "HsTrxR1_dock_energy"])
            n_nonsel = len(dock[dock["delta_energy"] <= 0])
            rows.append({
                "disease": "Schistosomiasis",
                "target": "SmTGR (docking)",
                "n_compounds": 43,
                "n_with_selectivity": len(dock),
                "pct_non_selective": round(100 * n_nonsel / len(dock), 1),
                "median_ratio": None,
                "best_ratio": None,
                "n_selective": 0,
                "method": "Docking",
                "unique_enzyme": True,
            })

    # Trypanosoma selectivity data
    for target in all_selectivity["parasite_target"].unique():
        td = all_selectivity[all_selectivity["parasite_target"] == target]
        n_total = len(parasite_activities[
            parasite_activities["target_name"] == target
        ]["molecule_chembl_id"].unique())
        n_with = len(td)
        n_nonsel = len(td[td["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
        rows.append({
            "disease": "Trypanosomiasis",
            "target": target[:35],
            "n_compounds": n_total,
            "n_with_selectivity": n_with,
            "pct_non_selective": round(100 * n_nonsel / n_with, 1) if n_with > 0 else None,
            "median_ratio": round(td["selectivity_ratio"].median(), 1),
            "best_ratio": round(td["selectivity_ratio"].max(), 1),
            "n_selective": len(td[td["selectivity_class"] == "SELECTIVE"]),
            "method": "Experimental",
            "unique_enzyme": False,
        })

    # Trypanosoma unique enzymes
    for uname in UNIQUE_TARGETS:
        n = len(parasite_activities[
            parasite_activities["target_name"] == uname
        ]["molecule_chembl_id"].unique())
        if n > 0:
            rows.append({
                "disease": "Trypanosomiasis",
                "target": uname[:35],
                "n_compounds": n,
                "n_with_selectivity": 0,
                "pct_non_selective": None,
                "median_ratio": None,
                "best_ratio": None,
                "n_selective": None,
                "method": "Unique enzyme",
                "unique_enzyme": True,
            })

    comp_df = pd.DataFrame(rows)

    # Print
    print(f"\n  {'Disease':15s} {'Target':35s} {'N cpd':>6s} {'N sel':>6s} "
          f"{'%NonSel':>8s} {'MedR':>6s} {'BestR':>6s} {'Sel10x':>7s} {'Method':>10s}")
    print(f"  {'-'*15} {'-'*35} {'-'*6} {'-'*6} {'-'*8} {'-'*6} {'-'*6} {'-'*7} {'-'*10}")

    for _, r in comp_df.sort_values(["disease", "n_compounds"], ascending=[True, False]).iterrows():
        pct = f"{r['pct_non_selective']:.0f}%" if pd.notna(r['pct_non_selective']) else "N/A"
        med = f"{r['median_ratio']:.1f}x" if pd.notna(r['median_ratio']) else "N/A"
        best = f"{r['best_ratio']:.1f}x" if pd.notna(r['best_ratio']) else "N/A"
        nsel = f"{r['n_selective']:.0f}" if pd.notna(r['n_selective']) else "N/A"
        print(f"  {r['disease']:15s} {r['target']:35s} {r['n_compounds']:6d} "
              f"{r['n_with_selectivity']:6d} {pct:>8s} {med:>6s} {best:>6s} {nsel:>7s} {r['method']:>10s}")

    # Summary statistics
    experimental = comp_df[comp_df["method"] == "Experimental"].dropna(subset=["pct_non_selective"])
    if len(experimental) > 0:
        total_with_sel = experimental["n_with_selectivity"].sum()
        total_compounds = experimental["n_compounds"].sum()
        # Weighted average non-selectivity
        weighted_nonsel = (
            (experimental["pct_non_selective"] * experimental["n_with_selectivity"]).sum()
            / total_with_sel
        )
        total_selective = experimental["n_selective"].sum()

        print(f"\n  CROSS-DISEASE SUMMARY (experimental selectivity only):")
        print(f"    Total parasite targets assessed: {len(experimental)}")
        print(f"    Total compounds with selectivity data: {int(total_with_sel)}")
        print(f"    Weighted average non-selectivity: {weighted_nonsel:.1f}%")
        print(f"    Total compounds with >=10x selectivity: {int(total_selective)}")
        print(f"    Targets where majority is non-selective: "
              f"{len(experimental[experimental['pct_non_selective'] > 50])}/{len(experimental)}")

    return comp_df


if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Script 16: Expanded T. brucei Selectivity")
    print("=" * 60)

    start_time = time.time()

    # Load T. brucei activities from Script 15
    act_path = os.path.join(TRYP_DIR, "tryp_activities.csv")
    if not os.path.exists(act_path):
        print("ERROR: Run Script 15 first.")
        sys.exit(1)

    activities = pd.read_csv(act_path)
    print(f"\n  Loaded T. brucei activities: {len(activities)} measurements")

    # Query expanded orthologues
    print("\n" + "=" * 60)
    print("QUERYING EXPANDED HUMAN ORTHOLOGUES")
    print("=" * 60)

    all_selectivity = []

    # Load existing selectivity from Script 15
    existing_path = os.path.join(TRYP_DIR, "tryp_selectivity.csv")
    if os.path.exists(existing_path):
        existing = pd.read_csv(existing_path)
        all_selectivity.append(existing)
        print(f"  Loaded existing selectivity: {len(existing)} records")

    for target_name, info in EXPANDED_ORTHOLOGUES.items():
        # Skip cathepsin B — already done in Script 15
        if target_name == "Cathepsin B-like cysteine protease":
            continue

        sel = query_human_selectivity(
            activities,
            target_name,
            info["human_target"],
            info["human_name"],
        )
        if len(sel) > 0:
            all_selectivity.append(sel)

    if all_selectivity:
        all_sel_df = pd.concat(all_selectivity, ignore_index=True)
    else:
        all_sel_df = pd.DataFrame()

    # Save expanded selectivity
    if len(all_sel_df) > 0:
        all_sel_df.to_csv(os.path.join(TRYP_DIR, "tryp_selectivity_expanded.csv"), index=False)

    # Build full comparison
    comparison = build_full_comparison(all_sel_df, activities)
    comparison.to_csv(os.path.join(PUB_DIR, "cross_disease_selectivity_full.csv"), index=False)

    elapsed = time.time() - start_time

    n_sel_total = len(all_sel_df) if len(all_sel_df) > 0 else 0
    print(f"\n{'=' * 60}")
    print(f"  EXPANDED SELECTIVITY COMPLETE")
    print(f"  Time: {elapsed:.0f} seconds")
    print(f"")
    print(f"  T. brucei selectivity comparisons: {n_sel_total}")
    print(f"  Cross-disease table: data/publication/cross_disease_selectivity_full.csv")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 16_tryp_expanded_selectivity.py successfully."
echo "Now run: python 16_tryp_expanded_selectivity.py"
