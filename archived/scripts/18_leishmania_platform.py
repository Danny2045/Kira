"""
Kira - Script 18: Third Disease — Leishmania
================================================

Leishmania causes leishmaniasis (~12 million infected, ~1 billion at risk).
Three clinical forms: visceral (fatal if untreated), cutaneous, mucocutaneous.
Endemic in East Africa, South Asia, Middle East, Latin America.

Why Leishmania matters for the platform:
  - Leishmania and Trypanosoma are BOTH trypanosomatids
  - They share trypanothione-based thiol metabolism
  - If selectivity patterns match across all three diseases,
    the platform finding extends to the entire kinetoplastid family

Key targets:
  - Trypanothione reductase (shared with T. brucei — unique to trypanosomatids)
  - Pteridine reductase (PTR1) — reduces folates, partially overlaps DHFR
  - Dihydrofolate reductase (DHFR) — human orthologue: CHEMBL202
  - Cysteine proteases — human orthologues: cathepsins
  - Topoisomerases — human orthologues exist
  - Sterol 14-alpha demethylase (CYP51) — human CYP51A1

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 18_leishmania_platform.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

BASE_DIR = os.path.dirname(__file__)
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
LEISH_DIR = os.path.join(BASE_DIR, "data", "leishmania")
PUB_DIR = os.path.join(BASE_DIR, "data", "publication")

# Human orthologue mappings for Leishmania targets
LEISH_ORTHOLOGUES = {
    "trypanothione reductase": {"human_id": None, "status": "UNIQUE", "human_name": "None (unique to trypanosomatids)"},
    "pteridine reductase": {"human_id": "CHEMBL202", "status": "QUERY", "human_name": "DHFR (closest functional equivalent)"},
    "dihydrofolate reductase": {"human_id": "CHEMBL202", "status": "QUERY", "human_name": "DHFR"},
    "sterol 14": {"human_id": "CHEMBL1978", "status": "QUERY", "human_name": "CYP51A1"},
    "topoisomerase": {"human_id": "CHEMBL1806", "status": "QUERY", "human_name": "TOP2A"},
    "cathepsin": {"human_id": "CHEMBL3524", "status": "QUERY", "human_name": "Cathepsin L"},
    "cysteine protease": {"human_id": "CHEMBL3524", "status": "QUERY", "human_name": "Cathepsin L"},
    "phosphodiesterase": {"human_id": "CHEMBL275", "status": "QUERY", "human_name": "PDE4B"},
    "n-myristoyltransferase": {"human_id": "CHEMBL2146302", "status": "QUERY", "human_name": "NMT1"},
    "miltefosine": {"human_id": None, "status": "SKIP", "human_name": "N/A"},
}


def query_leishmania_targets():
    """Find all Leishmania protein targets in ChEMBL."""
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("STEP 1: QUERYING LEISHMANIA TARGETS")
    print("=" * 60)

    target_api = new_client.target
    all_targets = []

    species = [
        "Leishmania donovani",      # Visceral leishmaniasis (kala-azar)
        "Leishmania major",         # Cutaneous leishmaniasis
        "Leishmania infantum",      # Visceral (Mediterranean)
        "Leishmania mexicana",      # Cutaneous (Americas)
        "Leishmania amazonensis",   # Cutaneous (Americas)
        "Leishmania tropica",       # Cutaneous (Old World)
        "Leishmania braziliensis",  # Mucocutaneous
        "Leishmania",               # Generic
    ]

    for organism in species:
        try:
            results = target_api.filter(
                organism__istartswith=organism,
                target_type="SINGLE PROTEIN",
            ).only(["target_chembl_id", "pref_name", "organism", "target_type"])
            batch = list(results)
            if batch:
                all_targets.extend(batch)
                print(f"  {organism}: {len(batch)} targets")
            time.sleep(0.5)
        except Exception as e:
            print(f"  {organism}: error — {e}")

    # Deduplicate
    seen = set()
    unique = []
    for t in all_targets:
        tid = t["target_chembl_id"]
        if tid not in seen:
            seen.add(tid)
            unique.append(t)

    targets_df = pd.DataFrame(unique)
    print(f"\n  Total unique Leishmania targets: {len(targets_df)}")

    for _, row in targets_df.iterrows():
        print(f"    {row['target_chembl_id']:15s}  {str(row['pref_name'])[:45]:45s}  {row['organism']}")

    return targets_df


def query_activities(targets_df):
    """Get activity data for each target."""
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("STEP 2: RETRIEVING ACTIVITY DATA")
    print("=" * 60)

    activity_api = new_client.activity
    all_activities = []

    for _, row in targets_df.iterrows():
        tid = row["target_chembl_id"]
        tname = row["pref_name"]

        try:
            results = activity_api.filter(
                target_chembl_id=tid,
                standard_type__in=["IC50", "EC50", "Ki", "Kd"],
                standard_relation="=",
                standard_units="nM",
            ).only([
                "molecule_chembl_id", "canonical_smiles",
                "standard_type", "standard_value", "standard_units",
                "target_chembl_id", "pref_name",
            ])

            batch = list(results)
            if batch:
                for r in batch:
                    r["target_name"] = tname
                    r["target_organism"] = row["organism"]
                all_activities.extend(batch)
                print(f"  {str(tname)[:40]:40s}  {len(batch):5d} measurements")

            time.sleep(0.3)

        except Exception as e:
            print(f"  {str(tname)[:40]:40s}  ERROR: {e}")

    if not all_activities:
        print("  No activities found!")
        return pd.DataFrame()

    df = pd.DataFrame(all_activities)
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value", "canonical_smiles"])

    print(f"\n  Total quality-filtered measurements: {len(df)}")
    print(f"  Unique compounds: {df['molecule_chembl_id'].nunique()}")

    # Per-target summary
    summary = (
        df.groupby("target_name")
        .agg(n_meas=("standard_value", "count"),
             n_cpd=("molecule_chembl_id", "nunique"),
             best=("standard_value", "min"),
             median=("standard_value", "median"))
        .sort_values("n_cpd", ascending=False)
    )

    print(f"\n  {'Target':40s} {'N meas':>7s} {'N cpd':>6s} {'Best nM':>8s}")
    print(f"  {'-'*40} {'-'*7} {'-'*6} {'-'*8}")
    for t, r in summary.iterrows():
        print(f"  {str(t)[:40]:40s} {r['n_meas']:7.0f} {r['n_cpd']:6.0f} {r['best']:8.0f}")

    return df


def query_selectivity(activities_df):
    """Query human orthologues and compute selectivity."""
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("STEP 3: SELECTIVITY ANALYSIS")
    print("=" * 60)

    activity_api = new_client.activity
    all_selectivity = []
    target_status = {}

    for target_name in activities_df["target_name"].unique():
        target_lower = str(target_name).lower()

        # Find matching orthologue
        matched = None
        for pattern, info in LEISH_ORTHOLOGUES.items():
            if pattern in target_lower:
                matched = info
                break

        if matched is None:
            target_status[target_name] = "UNMAPPED"
            continue

        if matched["status"] == "UNIQUE":
            print(f"  {str(target_name)[:40]:40s} → UNIQUE (no human equivalent)")
            target_status[target_name] = "UNIQUE"
            continue

        if matched["status"] == "SKIP":
            target_status[target_name] = "SKIP"
            continue

        human_id = matched["human_id"]
        human_name = matched["human_name"]

        # Get parasite compounds
        target_data = activities_df[activities_df["target_name"] == target_name]
        compound_ids = target_data["molecule_chembl_id"].unique().tolist()

        print(f"  {str(target_name)[:40]:40s} → {human_name} ({human_id})")
        print(f"    Compounds to check: {len(compound_ids)}")

        # Query human
        human_records = []
        for i in range(0, len(compound_ids), 50):
            batch = compound_ids[i:i+50]
            try:
                results = activity_api.filter(
                    target_chembl_id=human_id,
                    molecule_chembl_id__in=batch,
                    standard_type__in=["IC50", "EC50", "Ki", "Kd"],
                    standard_relation="=",
                    standard_units="nM",
                ).only(["molecule_chembl_id", "standard_value", "standard_type"])
                human_records.extend(list(results))
                time.sleep(0.3)
            except:
                pass

        print(f"    Human measurements: {len(human_records)}")

        if not human_records:
            target_status[target_name] = "NO_OVERLAP"
            continue

        human_df = pd.DataFrame(human_records)
        human_df["standard_value"] = pd.to_numeric(human_df["standard_value"], errors="coerce")
        human_df = human_df.dropna(subset=["standard_value"])

        # Best IC50 per compound
        para_best = (target_data.groupby("molecule_chembl_id")["standard_value"]
                     .min().reset_index().rename(columns={"standard_value": "parasite_ic50"}))
        human_best = (human_df.groupby("molecule_chembl_id")["standard_value"]
                      .min().reset_index().rename(columns={"standard_value": "human_ic50"}))

        merged = para_best.merge(human_best, on="molecule_chembl_id")
        if len(merged) == 0:
            target_status[target_name] = "NO_OVERLAP"
            continue

        merged["selectivity_ratio"] = merged["human_ic50"] / merged["parasite_ic50"]
        merged["parasite_target"] = target_name
        merged["human_target"] = human_name
        merged["disease"] = "Leishmaniasis"

        def classify(r):
            if r >= 10: return "SELECTIVE"
            elif r >= 3: return "MODERATE"
            elif r >= 1: return "POOR"
            else: return "COUNTER-SELECTIVE"

        merged["selectivity_class"] = merged["selectivity_ratio"].apply(classify)

        n = len(merged)
        n_sel = (merged["selectivity_class"] == "SELECTIVE").sum()
        n_counter = (merged["selectivity_class"] == "COUNTER-SELECTIVE").sum()
        n_nonsel = len(merged[merged["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
        pct = round(100 * n_nonsel / n, 1)

        print(f"    Dual-species: {n}, Selective: {n_sel}, Counter: {n_counter}, Non-sel: {pct}%")
        if n > 0:
            print(f"    Median ratio: {merged['selectivity_ratio'].median():.1f}x")

        all_selectivity.append(merged)
        target_status[target_name] = "ASSESSED"

    if all_selectivity:
        sel_df = pd.concat(all_selectivity, ignore_index=True)
    else:
        sel_df = pd.DataFrame()

    return sel_df, target_status


def three_disease_comparison(leish_activities, leish_selectivity, leish_status):
    """Build the definitive three-disease selectivity table."""

    print("\n" + "=" * 60)
    print("THREE-DISEASE SELECTIVITY LANDSCAPE")
    print("=" * 60)

    rows = []

    # Schistosomiasis
    s_sel_path = os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv")
    s_act_path = os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv")
    if os.path.exists(s_sel_path):
        s_sel = pd.read_csv(s_sel_path)
        s_act = pd.read_csv(s_act_path)
        for target in s_sel["parasite_target"].unique():
            td = s_sel[s_sel["parasite_target"] == target]
            n_total = len(s_act[s_act["target_name"] == target]["molecule_chembl_id"].unique())
            n_with = len(td)
            n_nonsel = len(td[td["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
            rows.append({
                "disease": "Schistosomiasis", "target": target[:30],
                "n_compounds": n_total, "n_with_sel": n_with,
                "pct_nonsel": round(100*n_nonsel/n_with, 1) if n_with > 0 else None,
                "median_ratio": round(td["selectivity_ratio"].median(), 1),
                "n_selective": len(td[td["selectivity_class"] == "SELECTIVE"]),
                "method": "Experimental",
            })

    # SmTGR docking
    dock_path = os.path.join(PROCESSED_DIR, "smtgr_docking_results.csv")
    if os.path.exists(dock_path):
        dock = pd.read_csv(dock_path).dropna(subset=["SmTGR_dock_energy", "HsTrxR1_dock_energy"])
        n_nonsel = len(dock[dock["delta_energy"] <= 0])
        rows.append({
            "disease": "Schistosomiasis", "target": "SmTGR (docking)",
            "n_compounds": 43, "n_with_sel": len(dock),
            "pct_nonsel": round(100*n_nonsel/len(dock), 1),
            "median_ratio": None, "n_selective": 0, "method": "Docking",
        })

    # Trypanosomiasis
    t_sel_path = os.path.join(TRYP_DIR, "tryp_selectivity_expanded.csv")
    t_act_path = os.path.join(TRYP_DIR, "tryp_activities.csv")
    if os.path.exists(t_sel_path):
        t_sel = pd.read_csv(t_sel_path)
        t_act = pd.read_csv(t_act_path)
        for target in t_sel["parasite_target"].unique():
            td = t_sel[t_sel["parasite_target"] == target]
            n_total = len(t_act[t_act["target_name"] == target]["molecule_chembl_id"].unique())
            n_with = len(td)
            n_nonsel = len(td[td["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
            rows.append({
                "disease": "Trypanosomiasis", "target": target[:30],
                "n_compounds": n_total, "n_with_sel": n_with,
                "pct_nonsel": round(100*n_nonsel/n_with, 1) if n_with > 0 else None,
                "median_ratio": round(td["selectivity_ratio"].median(), 1),
                "n_selective": len(td[td["selectivity_class"] == "SELECTIVE"]),
                "method": "Experimental",
            })

    # Leishmaniasis
    if len(leish_selectivity) > 0:
        for target in leish_selectivity["parasite_target"].unique():
            td = leish_selectivity[leish_selectivity["parasite_target"] == target]
            n_total = len(leish_activities[
                leish_activities["target_name"] == target
            ]["molecule_chembl_id"].unique())
            n_with = len(td)
            n_nonsel = len(td[td["selectivity_class"].isin(["POOR", "COUNTER-SELECTIVE"])])
            rows.append({
                "disease": "Leishmaniasis", "target": target[:30],
                "n_compounds": n_total, "n_with_sel": n_with,
                "pct_nonsel": round(100*n_nonsel/n_with, 1) if n_with > 0 else None,
                "median_ratio": round(td["selectivity_ratio"].median(), 1),
                "n_selective": len(td[td["selectivity_class"] == "SELECTIVE"]),
                "method": "Experimental",
            })

    # Add unique enzyme entries for Leishmania
    for tname, status in leish_status.items():
        if status == "UNIQUE":
            n = len(leish_activities[
                leish_activities["target_name"] == tname
            ]["molecule_chembl_id"].unique())
            if n > 0:
                rows.append({
                    "disease": "Leishmaniasis", "target": str(tname)[:30],
                    "n_compounds": n, "n_with_sel": 0,
                    "pct_nonsel": None, "median_ratio": None,
                    "n_selective": None, "method": "Unique enzyme",
                })

    comp = pd.DataFrame(rows)

    # Print
    print(f"\n  {'Disease':15s} {'Target':30s} {'N cpd':>6s} {'N sel':>6s} "
          f"{'%NonSel':>8s} {'MedR':>6s} {'Sel10x':>7s} {'Method':>10s}")
    print(f"  {'-'*15} {'-'*30} {'-'*6} {'-'*6} {'-'*8} {'-'*6} {'-'*7} {'-'*10}")

    for _, r in comp.sort_values(["disease", "n_compounds"], ascending=[True, False]).iterrows():
        pct = f"{r['pct_nonsel']:.0f}%" if pd.notna(r['pct_nonsel']) else "N/A"
        med = f"{r['median_ratio']:.1f}x" if pd.notna(r['median_ratio']) else "N/A"
        nsel = f"{r['n_selective']:.0f}" if pd.notna(r['n_selective']) else "N/A"
        print(f"  {r['disease']:15s} {r['target']:30s} {r['n_compounds']:6d} "
              f"{r['n_with_sel']:6d} {pct:>8s} {med:>6s} {nsel:>7s} {r['method']:>10s}")

    # Cross-disease summary
    exp = comp[comp["method"] == "Experimental"].dropna(subset=["pct_nonsel"])
    if len(exp) > 0:
        total_sel_data = exp["n_with_sel"].sum()
        weighted_nonsel = (exp["pct_nonsel"] * exp["n_with_sel"]).sum() / total_sel_data
        total_selective = exp["n_selective"].sum()
        n_diseases = comp["disease"].nunique()

        print(f"\n  ╔══════════════════════════════════════════════════════════╗")
        print(f"  ║  THREE-DISEASE PLATFORM SUMMARY                         ║")
        print(f"  ╠══════════════════════════════════════════════════════════╣")
        print(f"  ║  Diseases: {n_diseases}                                          ║")
        print(f"  ║  Experimental selectivity comparisons: {int(total_sel_data):>4d}             ║")
        print(f"  ║  Weighted non-selectivity: {weighted_nonsel:>5.1f}%                    ║")
        print(f"  ║  Compounds with >=10x selectivity: {int(total_selective):>3d}                 ║")
        print(f"  ║  Majority-non-selective targets: "
              f"{len(exp[exp['pct_nonsel'] > 50])}/{len(exp)}                    ║")
        print(f"  ╚══════════════════════════════════════════════════════════╝")

    return comp


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Script 18: Third Disease — Leishmania")
    print("  Three-Disease Platform Validation")
    print("=" * 60)

    start_time = time.time()
    os.makedirs(LEISH_DIR, exist_ok=True)

    # Step 1: Find targets
    targets_df = query_leishmania_targets()
    if len(targets_df) == 0:
        print("ERROR: No Leishmania targets found.")
        sys.exit(1)
    targets_df.to_csv(os.path.join(LEISH_DIR, "leish_targets.csv"), index=False)

    # Step 2: Activities
    activities_df = query_activities(targets_df)
    if len(activities_df) == 0:
        print("ERROR: No activity data.")
        sys.exit(1)
    activities_df.to_csv(os.path.join(LEISH_DIR, "leish_activities.csv"), index=False)

    # Step 3: Selectivity
    selectivity_df, target_status = query_selectivity(activities_df)
    if len(selectivity_df) > 0:
        selectivity_df.to_csv(os.path.join(LEISH_DIR, "leish_selectivity.csv"), index=False)

    # Step 4: Three-disease comparison
    comparison = three_disease_comparison(activities_df, selectivity_df, target_status)
    comparison.to_csv(os.path.join(PUB_DIR, "three_disease_selectivity.csv"), index=False)

    elapsed = time.time() - start_time

    n_targets = activities_df["target_name"].nunique()
    n_compounds = activities_df["molecule_chembl_id"].nunique()
    n_sel = len(selectivity_df) if len(selectivity_df) > 0 else 0

    print(f"\n{'=' * 60}")
    print(f"  LEISHMANIA COMPLETE — THREE-DISEASE PLATFORM")
    print(f"  Time: {elapsed:.0f} seconds")
    print(f"")
    print(f"  Leishmania:")
    print(f"    Targets: {n_targets}")
    print(f"    Compounds: {n_compounds}")
    print(f"    Selectivity comparisons: {n_sel}")
    print(f"")
    print(f"  KIRA NOW COVERS THREE DISEASES:")
    print(f"    1. Schistosomiasis (S. mansoni)")
    print(f"    2. Trypanosomiasis (T. brucei)")
    print(f"    3. Leishmaniasis (Leishmania spp.)")
    print(f"")
    print(f"  18 scripts. 3 diseases. 1 platform.")
    print(f"{'=' * 60}")
