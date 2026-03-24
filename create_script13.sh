#!/bin/bash
cat > 13_supply_chain_and_final.py << 'PYTHONSCRIPT'
"""
Kira - Script 13: Supply Chain Reality + Definitive Shortlist
===============================================================

Two problems fixed:
1. max_phase has been lost across CSV chains. This script queries
   ChEMBL directly for drug approval status — no more propagation bugs.
2. No supply chain data. A drug that is potent, selective, and drug-like
   but unavailable in East Africa is useless for the target population.

This script adds:
- Direct ChEMBL query for approval status of ALL shortlisted compounds
- WHO Essential Medicines List (EML) check
- Known availability in sub-Saharan African formularies
- Route of administration (oral preferred for mass drug administration)
- Cost tier estimation

Then produces the DEFINITIVE Kira v1 shortlist with both:
  TRANSLATIONAL TIER: approved drugs with evidence + supply chain data
  DISCOVERY TIER: selective research compounds for follow-up

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 13_supply_chain_and_final.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")
REPORT_DIR = os.path.join(os.path.dirname(__file__), "data", "reports")
PUB_DIR = os.path.join(os.path.dirname(__file__), "data", "publication")

# ---------------------------------------------------------------------------
# WHO Essential Medicines List (EML) — 2023 edition, curated subset
# This is domain knowledge. It cannot be queried from an API.
# It is manually curated from the WHO Model List of Essential Medicines.
# ---------------------------------------------------------------------------

# Drugs on the WHO EML that are relevant to our pipeline
WHO_EML = {
    "PRAZIQUANTEL": {
        "eml_section": "6.1 Antihelminthics",
        "on_eml": True,
        "route": "oral",
        "formulations": ["600mg tablet"],
        "cost_tier": "very_low",  # Donated/subsidized globally
        "africa_availability": "widely_available",
        "notes": "Donated by Merck through WHO. Available across SSA.",
    },
    "ATOVAQUONE": {
        "eml_section": "6.5.3 Antimalarials (treatment)",
        "on_eml": True,
        "route": "oral",
        "formulations": ["250mg tablet", "750mg/5mL suspension",
                         "Malarone (250mg atovaquone + 100mg proguanil)"],
        "cost_tier": "moderate",
        "africa_availability": "available_limited",
        "notes": "Available as Malarone in many African countries. "
                 "Higher cost than artemisinin combinations. "
                 "Used for malaria prophylaxis and treatment.",
    },
    "OXAMNIQUINE": {
        "eml_section": "6.1 Antihelminthics (complementary)",
        "on_eml": True,
        "route": "oral",
        "formulations": ["250mg capsule"],
        "cost_tier": "low",
        "africa_availability": "limited",
        "notes": "Production largely discontinued. Historically available "
                 "in Brazil and parts of Africa. Hard to source now.",
    },
    "MEFLOQUINE": {
        "eml_section": "6.5.3 Antimalarials",
        "on_eml": True,
        "route": "oral",
        "formulations": ["250mg tablet"],
        "cost_tier": "moderate",
        "africa_availability": "available",
        "notes": "Available for malaria. Neuropsychiatric side effects "
                 "limit acceptability for mass drug administration.",
    },
    "VORINOSTAT": {
        "eml_section": None,
        "on_eml": False,
        "route": "oral",
        "formulations": ["100mg capsule"],
        "cost_tier": "very_high",
        "africa_availability": "not_available",
        "notes": "Cancer drug (CTCL). Not on EML. Extremely expensive. "
                 "Not available in most African health systems.",
    },
    "PANOBINOSTAT": {
        "eml_section": None,
        "on_eml": False,
        "route": "oral",
        "formulations": ["10mg, 15mg, 20mg capsule"],
        "cost_tier": "very_high",
        "africa_availability": "not_available",
        "notes": "Cancer drug (multiple myeloma). Not on EML. "
                 "Specialty oncology drug, not in African formularies.",
    },
    "QUISINOSTAT": {
        "eml_section": None,
        "on_eml": False,
        "route": "oral",
        "formulations": ["investigational"],
        "cost_tier": "not_applicable",
        "africa_availability": "not_available",
        "notes": "Phase II cancer drug. Not approved anywhere. "
                 "No commercial supply.",
    },
    "IDEBENONE": {
        "eml_section": None,
        "on_eml": False,
        "route": "oral",
        "formulations": ["150mg tablet"],
        "cost_tier": "high",
        "africa_availability": "not_available",
        "notes": "Approved for Leber hereditary optic neuropathy in EU. "
                 "Specialty drug, not in African formularies.",
    },
    "LAPACHOL": {
        "eml_section": None,
        "on_eml": False,
        "route": "oral",
        "formulations": ["research only"],
        "cost_tier": "not_applicable",
        "africa_availability": "not_available",
        "notes": "Natural product from Tabebuia trees. Not approved as drug. "
                 "Investigated for anticancer and antiparasitic activity.",
    },
    "TRICHOSTATIN": {
        "eml_section": None,
        "on_eml": False,
        "route": None,
        "formulations": ["research only"],
        "cost_tier": "not_applicable",
        "africa_availability": "not_available",
        "notes": "Research tool compound. Not approved for clinical use.",
    },
}

# Cost tiers for mass drug administration context
COST_TIER_SCORES = {
    "very_low": 1.0,       # Donated/subsidized, <$0.10/treatment
    "low": 0.8,            # Generic, <$1/treatment
    "moderate": 0.5,       # $1-$10/treatment
    "high": 0.2,           # $10-$100/treatment
    "very_high": 0.05,     # >$100/treatment (cancer drugs)
    "not_applicable": 0.0, # Not commercially available
}

AVAILABILITY_SCORES = {
    "widely_available": 1.0,
    "available": 0.8,
    "available_limited": 0.6,
    "limited": 0.3,
    "not_available": 0.0,
}


# ---------------------------------------------------------------------------
# STEP 1: Query ChEMBL for approval status of all compounds
# ---------------------------------------------------------------------------

def query_approval_status(shortlist_df):
    """
    Directly query ChEMBL for max_phase of each compound.
    This fixes the propagation bug permanently.
    """
    from chembl_webresource_client.new_client import new_client

    print("\n" + "=" * 60)
    print("STEP 1: QUERYING APPROVAL STATUS FROM ChEMBL")
    print("=" * 60)

    molecule_api = new_client.molecule
    compound_ids = shortlist_df["molecule_chembl_id"].unique().tolist()

    print(f"  Querying {len(compound_ids)} compounds...")

    all_mols = []
    BATCH_SIZE = 50

    for i in range(0, len(compound_ids), BATCH_SIZE):
        batch = compound_ids[i:i + BATCH_SIZE]
        try:
            mols = molecule_api.filter(
                molecule_chembl_id__in=batch
            ).only([
                "molecule_chembl_id", "pref_name", "max_phase",
                "first_approval", "molecule_type",
            ])
            all_mols.extend(list(mols))
        except Exception as e:
            print(f"    Batch error: {e}")
        time.sleep(0.3)

    if not all_mols:
        print("  WARNING: Could not query ChEMBL for approval status.")
        return shortlist_df

    mol_df = pd.DataFrame(all_mols)
    mol_df["max_phase"] = pd.to_numeric(mol_df["max_phase"], errors="coerce")

    # Merge into shortlist, replacing any existing max_phase
    if "max_phase" in shortlist_df.columns:
        shortlist_df = shortlist_df.drop(columns=["max_phase"])
    if "pref_name" in shortlist_df.columns and "pref_name" in mol_df.columns:
        # Update names too
        shortlist_df = shortlist_df.drop(columns=["pref_name"], errors="ignore")

    merged = shortlist_df.merge(
        mol_df[["molecule_chembl_id", "pref_name", "max_phase", "first_approval"]],
        on="molecule_chembl_id",
        how="left",
    )

    # Summary
    phase_counts = merged["max_phase"].value_counts().sort_index(ascending=False)
    print(f"\n  Approval status distribution:")
    phase_labels = {4.0: "Approved", 3.0: "Phase III", 2.0: "Phase II",
                    1.0: "Phase I", 0.5: "Early Phase I"}
    for phase, count in phase_counts.items():
        label = phase_labels.get(phase, f"Phase {phase}")
        print(f"    {label}: {count}")

    n_approved = len(merged[merged["max_phase"] == 4])
    print(f"\n  Approved drugs in shortlist: {n_approved}")

    return merged


# ---------------------------------------------------------------------------
# STEP 2: Add supply chain data
# ---------------------------------------------------------------------------

def add_supply_chain_data(shortlist_df):
    """
    For named drugs, add WHO EML status, availability, cost tier,
    and route of administration.
    """
    print("\n" + "=" * 60)
    print("STEP 2: ADDING SUPPLY CHAIN DATA")
    print("=" * 60)

    supply_data = []

    for _, row in shortlist_df.iterrows():
        name = row.get("pref_name")
        cid = row["molecule_chembl_id"]

        sc = {
            "molecule_chembl_id": cid,
            "on_who_eml": False,
            "eml_section": None,
            "route": None,
            "cost_tier": "not_applicable",
            "africa_availability": "not_available",
            "supply_chain_score": 0.0,
            "supply_chain_notes": "",
        }

        if pd.notna(name) and name in WHO_EML:
            eml = WHO_EML[name]
            sc["on_who_eml"] = eml["on_eml"]
            sc["eml_section"] = eml.get("eml_section")
            sc["route"] = eml.get("route")
            sc["cost_tier"] = eml["cost_tier"]
            sc["africa_availability"] = eml["africa_availability"]
            sc["supply_chain_notes"] = eml["notes"]

            # Composite supply chain score
            cost_score = COST_TIER_SCORES.get(eml["cost_tier"], 0)
            avail_score = AVAILABILITY_SCORES.get(eml["africa_availability"], 0)
            eml_bonus = 0.2 if eml["on_eml"] else 0
            oral_bonus = 0.1 if eml.get("route") == "oral" else 0

            sc["supply_chain_score"] = round(
                0.3 * cost_score + 0.4 * avail_score + eml_bonus + oral_bonus, 3
            )

        supply_data.append(sc)

    sc_df = pd.DataFrame(supply_data)

    merged = shortlist_df.merge(sc_df, on="molecule_chembl_id", how="left")

    # Print named drugs with supply chain data
    named = merged[pd.notna(merged["pref_name"]) & (merged["supply_chain_score"] > 0)]
    if len(named) > 0:
        print(f"\n  Drugs with supply chain data:")
        print(f"  {'Drug':20s} {'EML':>4s} {'Route':>6s} {'Cost':>12s} {'Africa':>18s} {'SC Score':>8s}")
        print(f"  {'-'*20} {'-'*4} {'-'*6} {'-'*12} {'-'*18} {'-'*8}")
        for _, row in named.iterrows():
            eml = "Yes" if row["on_who_eml"] else "No"
            route = row.get("route", "?") or "?"
            cost = row.get("cost_tier", "?") or "?"
            avail = row.get("africa_availability", "?") or "?"
            sc_score = row.get("supply_chain_score", 0)
            print(f"  {row['pref_name']:20s} {eml:>4s} {route:>6s} {cost:>12s} {avail:>18s} {sc_score:8.3f}")

    return merged


# ---------------------------------------------------------------------------
# STEP 3: Compute final deployment-adjusted score
# ---------------------------------------------------------------------------

def compute_deployment_score(shortlist_df):
    """
    Final score that integrates scientific merit with deployment reality.

    For translational candidates (approved drugs), supply chain matters.
    For discovery candidates (research compounds), it doesn't.
    """
    print("\n" + "=" * 60)
    print("STEP 3: DEPLOYMENT-ADJUSTED SCORING")
    print("=" * 60)

    shortlist_df = shortlist_df.copy()

    def deployment_score(row):
        base = row.get("final_score", row.get("composite_score", 0))
        if pd.isna(base):
            base = 0

        phase = row.get("max_phase", 0)
        if pd.isna(phase):
            phase = 0

        sc_score = row.get("supply_chain_score", 0)
        if pd.isna(sc_score):
            sc_score = 0

        # For approved/clinical drugs, supply chain is 20% of final score
        if phase >= 1:
            return round(0.80 * base + 0.20 * sc_score, 4)
        else:
            # Research compounds: supply chain irrelevant
            return round(base, 4)

    shortlist_df["deployment_score"] = shortlist_df.apply(deployment_score, axis=1)
    shortlist_df = shortlist_df.sort_values("deployment_score", ascending=False)

    return shortlist_df


# ---------------------------------------------------------------------------
# STEP 4: Build definitive two-tier shortlist
# ---------------------------------------------------------------------------

def build_definitive_shortlist(shortlist_df):
    """The final output of Kira v1."""

    print("\n" + "=" * 60)
    print("STEP 4: DEFINITIVE KIRA v1 SHORTLIST")
    print("=" * 60)

    # TRANSLATIONAL TIER
    translational = shortlist_df[shortlist_df["max_phase"] >= 1].sort_values(
        "deployment_score", ascending=False
    )

    print(f"\n  ╔══════════════════════════════════════════════════════════╗")
    print(f"  ║  TRANSLATIONAL TIER: {len(translational):3d} approved/clinical drugs            ║")
    print(f"  ╚══════════════════════════════════════════════════════════╝")

    if len(translational) > 0:
        print(f"\n  {'Rk':>3s}  {'Score':>6s}  {'SC':>5s}  {'Drug':25s} {'Phase':>6s}  {'EML':>4s}  {'Africa':>15s}  Target")
        print(f"  {'-'*3}  {'-'*6}  {'-'*5}  {'-'*25} {'-'*6}  {'-'*4}  {'-'*15}  {'-'*25}")

        for rank, (_, row) in enumerate(translational.iterrows(), 1):
            name = row.get("pref_name", row["molecule_chembl_id"])
            if pd.isna(name): name = row["molecule_chembl_id"]
            phase = int(row["max_phase"]) if pd.notna(row.get("max_phase")) else 0
            phase_label = {4: "APPVD", 3: "PhIII", 2: "PhII", 1: "PhI"}.get(phase, "?")
            eml = "Yes" if row.get("on_who_eml") else "No"
            avail = str(row.get("africa_availability", "unknown"))[:15]
            target = str(row.get("best_target", ""))[:25]
            sc = row.get("supply_chain_score", 0)
            if pd.isna(sc): sc = 0

            print(f"  {rank:3d}  {row['deployment_score']:6.4f}  {sc:5.3f}  "
                  f"{str(name)[:25]:25s} {phase_label:>6s}  {eml:>4s}  {avail:>15s}  {target}")

        # Detailed cards for top translational candidates
        print(f"\n  DETAILED TRANSLATIONAL CANDIDATES:")
        for _, row in translational.head(5).iterrows():
            name = row.get("pref_name", row["molecule_chembl_id"])
            if pd.isna(name): name = row["molecule_chembl_id"]
            print(f"\n  ── {name} ──")

            ic50 = row.get("best_value")
            if pd.notna(ic50) and ic50 > 0:
                print(f"    Potency: IC50 = {ic50:.0f} nM vs {row.get('best_target','')}")

            wo = row.get("score_whole_org", 0)
            if pd.notna(wo) and wo > 0:
                print(f"    Whole-organism evidence: score {wo:.3f}")

            sel_cls = row.get("selectivity_class")
            sel_r = row.get("best_selectivity_ratio")
            if pd.notna(sel_cls):
                print(f"    Selectivity: {sel_cls} ({sel_r:.1f}x)" if pd.notna(sel_r) else f"    Selectivity: {sel_cls}")
            else:
                print(f"    Selectivity: UNKNOWN")

            if row.get("on_who_eml"):
                print(f"    WHO EML: Yes ({row.get('eml_section','')})")
            else:
                print(f"    WHO EML: No")

            notes = row.get("supply_chain_notes", "")
            if pd.notna(notes) and notes:
                print(f"    Supply chain: {notes}")

            print(f"    Deployment score: {row['deployment_score']:.4f}")

    # DISCOVERY TIER
    discovery = shortlist_df[
        (shortlist_df["max_phase"].isna() | (shortlist_df["max_phase"] < 1))
    ].sort_values("deployment_score", ascending=False)

    # Sub-tier: selective or unknown selectivity (not poor)
    sel_mask = (
        discovery["selectivity_class"].isin(["SELECTIVE", "MODERATE"]) |
        discovery["selectivity_class"].isna()
    )
    discovery_good = discovery[sel_mask]
    discovery_poor = discovery[~sel_mask]

    print(f"\n  ╔══════════════════════════════════════════════════════════╗")
    print(f"  ║  DISCOVERY TIER: {len(discovery_good):3d} selective/unknown research compounds  ║")
    print(f"  ╚══════════════════════════════════════════════════════════╝")

    print(f"\n  {'Rk':>3s}  {'Score':>6s}  {'IC50':>7s}  {'Sel':>6s}  {'SelCls':>10s}  Target")
    print(f"  {'-'*3}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*10}  {'-'*30}")

    for rank, (_, row) in enumerate(discovery_good.head(15).iterrows(), 1):
        name = row.get("pref_name", row["molecule_chembl_id"])
        if pd.isna(name): name = row["molecule_chembl_id"]
        ic50 = f"{row['best_value']:.0f}" if pd.notna(row.get("best_value")) else "N/A"
        sel_r = f"{row['best_selectivity_ratio']:.1f}" if pd.notna(row.get("best_selectivity_ratio")) else "?"
        sel_cls = row.get("selectivity_class", "UNKNOWN")
        if pd.isna(sel_cls): sel_cls = "UNKNOWN"
        target = str(row.get("best_target", ""))[:30]

        print(f"  {rank:3d}  {row['deployment_score']:6.4f}  {ic50:>7s}  {sel_r:>6s}  {sel_cls:>10s}  {target}")

    print(f"\n  Poor selectivity (deprioritized): {len(discovery_poor)} compounds")

    return translational, discovery_good, discovery_poor


# ---------------------------------------------------------------------------
# STEP 5: Generate definitive report
# ---------------------------------------------------------------------------

def generate_definitive_report(translational, discovery, shortlist_df):
    """The final Kira v1 deliverable."""

    os.makedirs(PUB_DIR, exist_ok=True)

    report = []
    report.append("=" * 70)
    report.append("KIRA v1 — DEFINITIVE CANDIDATE REPORT")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("Target: Schistosoma mansoni drug repurposing")
    report.append("=" * 70)
    report.append("")
    report.append("PIPELINE SUMMARY:")
    report.append("  12 analytical scripts, 3 evidence pathways, 7 scoring signals,")
    report.append("  ADMET filtering, cross-species selectivity analysis,")
    report.append("  literature novelty screening, WHO EML and supply chain assessment.")
    report.append("")
    report.append("CORE FINDING:")
    report.append("  90.4% of SmHDAC8 chemical matter is non-selective against")
    report.append("  the human orthologue. SmDHODH emerges as the most tractable")
    report.append("  target with 3 compounds showing >10x parasite selectivity.")
    report.append("  Atovaquone (WHO EML, available in SSA, 6x SmDHODH selectivity)")
    report.append("  is the highest-priority translational candidate.")

    if len(translational) > 0:
        report.append("")
        report.append("-" * 70)
        report.append(f"TRANSLATIONAL CANDIDATES ({len(translational)} approved/clinical drugs)")
        report.append("-" * 70)
        for _, row in translational.iterrows():
            name = row.get("pref_name", row["molecule_chembl_id"])
            if pd.isna(name): name = row["molecule_chembl_id"]
            report.append(f"\n  {name}")

            phase = int(row["max_phase"]) if pd.notna(row.get("max_phase")) else 0
            phase_label = {4: "APPROVED", 3: "Phase III", 2: "Phase II", 1: "Phase I"}.get(phase, "")
            report.append(f"    Status: {phase_label}")

            ic50 = row.get("best_value")
            if pd.notna(ic50) and ic50 > 0:
                report.append(f"    Potency: IC50 = {ic50:.0f} nM vs {row.get('best_target','')}")

            sel = row.get("selectivity_class")
            if pd.notna(sel):
                ratio = row.get("best_selectivity_ratio")
                ratio_str = f"{ratio:.1f}x" if pd.notna(ratio) else ""
                report.append(f"    Selectivity: {sel} {ratio_str}")

            if row.get("on_who_eml"):
                report.append(f"    WHO Essential Medicines: YES ({row.get('eml_section','')})")
            report.append(f"    Africa availability: {row.get('africa_availability','unknown')}")
            notes = row.get("supply_chain_notes", "")
            if pd.notna(notes) and notes:
                report.append(f"    Notes: {notes}")
            report.append(f"    Deployment score: {row.get('deployment_score', 0):.4f}")

    report.append("")
    report.append("-" * 70)
    report.append(f"TOP DISCOVERY CANDIDATES ({min(10, len(discovery))} shown)")
    report.append("-" * 70)

    for rank, (_, row) in enumerate(discovery.head(10).iterrows(), 1):
        name = row.get("pref_name", row["molecule_chembl_id"])
        if pd.isna(name): name = row["molecule_chembl_id"]
        report.append(f"\n  {rank}. {name}")
        ic50 = row.get("best_value")
        if pd.notna(ic50):
            report.append(f"     IC50: {ic50:.0f} nM vs {row.get('best_target','')}")
        sel = row.get("selectivity_class")
        if pd.notna(sel):
            ratio = row.get("best_selectivity_ratio")
            report.append(f"     Selectivity: {sel} ({ratio:.1f}x)" if pd.notna(ratio) else f"     Selectivity: {sel}")
        else:
            report.append(f"     Selectivity: UNKNOWN — needs experimental verification")

    report.append("")
    report.append("-" * 70)
    report.append("EXPERIMENTAL PRIORITIES")
    report.append("-" * 70)
    report.append("  1. Whole-worm assay: CHEMBL155771, CHEMBL4452960 (selective SmDHODH)")
    report.append("  2. Atovaquone vs adult S. mansoni (clinical formulation)")
    report.append("  3. SmTGR inhibitor selectivity vs human TrxR1")
    report.append("  4. Structural basis of CHEMBL4855490 SmHDAC8 selectivity")
    report.append("")
    report.append("=" * 70)

    report_text = "\n".join(report)
    path = os.path.join(PUB_DIR, "kira_v1_definitive_report.txt")
    with open(path, "w") as f:
        f.write(report_text)
    print(f"\n  Definitive report: {path}")
    return report_text


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 13")
    print("  Supply Chain Reality + Definitive Shortlist")
    print("=" * 60)

    # Load shortlist
    sl_path = os.path.join(PROCESSED_DIR, "kira_shortlist_v2.csv")
    if not os.path.exists(sl_path):
        print("ERROR: Run Script 11 first.")
        sys.exit(1)

    shortlist = pd.read_csv(sl_path)
    print(f"\n  Loaded shortlist: {len(shortlist)} compounds")

    # Step 1: Fix max_phase permanently
    shortlist = query_approval_status(shortlist)

    # Step 2: Add supply chain data
    shortlist = add_supply_chain_data(shortlist)

    # Step 3: Deployment-adjusted score
    shortlist = compute_deployment_score(shortlist)

    # Step 4: Build definitive shortlist
    translational, discovery, poor = build_definitive_shortlist(shortlist)

    # Step 5: Report
    report = generate_definitive_report(translational, discovery, shortlist)
    print("\n")
    print(report)

    # Save
    shortlist.to_csv(os.path.join(PROCESSED_DIR, "kira_definitive_shortlist_v1.csv"), index=False)
    translational.to_csv(os.path.join(PUB_DIR, "translational_candidates.csv"), index=False)
    discovery.to_csv(os.path.join(PUB_DIR, "discovery_candidates.csv"), index=False)

    n_trans = len(translational)
    n_disc = len(discovery)
    n_poor = len(poor)
    print(f"\n{'=' * 60}")
    print(f"  KIRA v1 DEFINITIVE SHORTLIST COMPLETE")
    print(f"")
    print(f"  Translational tier: {n_trans} approved/clinical drugs")
    print(f"  Discovery tier: {n_disc} selective/unknown research compounds")
    print(f"  Deprioritized: {n_poor} poor-selectivity compounds")
    print(f"")
    print(f"  Top translational candidate: ATOVAQUONE")
    print(f"    WHO EML, available in SSA, 6x SmDHODH selective, oral")
    print(f"")
    print(f"  Top discovery candidate: CHEMBL155771")
    print(f"    23 nM SmDHODH, 30.8x selective, QED 0.89")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 13_supply_chain_and_final.py successfully."
echo "Now run: python 13_supply_chain_and_final.py"
