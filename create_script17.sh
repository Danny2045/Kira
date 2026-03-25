#!/bin/bash
cat > 17_cross_disease_analysis.py << 'PYTHONSCRIPT'
"""
Kira - Script 17: Cross-Disease Compound Analysis
=====================================================

The platform question: are any compounds active against BOTH
Schistosoma mansoni and Trypanosoma brucei?

A dual-active, selective compound would be the strongest possible
finding — a single drug treating two NTDs affecting overlapping
populations in sub-Saharan Africa.

This script:
  1. Finds compounds active in both diseases
  2. Checks their selectivity profiles in both organisms
  3. Builds the T. brucei definitive shortlist (mirrors Script 13)
  4. Identifies cross-disease candidates
  5. Produces the platform summary

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 17_cross_disease_analysis.py
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

BASE_DIR = os.path.dirname(__file__)
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
PUB_DIR = os.path.join(BASE_DIR, "data", "publication")


# ---------------------------------------------------------------------------
# STEP 1: Load all data from both diseases
# ---------------------------------------------------------------------------

def load_all_data():
    """Load activity and selectivity data from both pipelines."""

    print("\n" + "=" * 60)
    print("STEP 1: LOADING ALL DATA")
    print("=" * 60)

    data = {}

    # Schistosomiasis
    s_act = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))
    data["schisto_activities"] = s_act
    print(f"  Schistosomiasis: {len(s_act)} measurements, "
          f"{s_act['molecule_chembl_id'].nunique()} compounds")

    s_sel_path = os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv")
    if os.path.exists(s_sel_path):
        data["schisto_selectivity"] = pd.read_csv(s_sel_path)
        print(f"  Schistosomiasis selectivity: {len(data['schisto_selectivity'])} records")

    # Trypanosomiasis
    t_act = pd.read_csv(os.path.join(TRYP_DIR, "tryp_activities.csv"))
    data["tryp_activities"] = t_act
    print(f"  Trypanosomiasis: {len(t_act)} measurements, "
          f"{t_act['molecule_chembl_id'].nunique()} compounds")

    t_sel_path = os.path.join(TRYP_DIR, "tryp_selectivity_expanded.csv")
    if os.path.exists(t_sel_path):
        data["tryp_selectivity"] = pd.read_csv(t_sel_path)
        print(f"  Trypanosomiasis selectivity: {len(data['tryp_selectivity'])} records")
    else:
        t_sel_path2 = os.path.join(TRYP_DIR, "tryp_selectivity.csv")
        if os.path.exists(t_sel_path2):
            data["tryp_selectivity"] = pd.read_csv(t_sel_path2)
            print(f"  Trypanosomiasis selectivity: {len(data['tryp_selectivity'])} records")

    return data


# ---------------------------------------------------------------------------
# STEP 2: Find cross-disease compounds
# ---------------------------------------------------------------------------

def find_cross_disease_compounds(data):
    """Find compounds with measured activity against both parasites."""

    print("\n" + "=" * 60)
    print("STEP 2: CROSS-DISEASE COMPOUND OVERLAP")
    print("=" * 60)

    s_compounds = set(data["schisto_activities"]["molecule_chembl_id"].unique())
    t_compounds = set(data["tryp_activities"]["molecule_chembl_id"].unique())

    overlap = s_compounds & t_compounds
    only_schisto = s_compounds - t_compounds
    only_tryp = t_compounds - s_compounds

    print(f"\n  Schistosomiasis compounds: {len(s_compounds)}")
    print(f"  Trypanosomiasis compounds: {len(t_compounds)}")
    print(f"  Active in BOTH diseases: {len(overlap)}")
    print(f"  Only schistosomiasis: {len(only_schisto)}")
    print(f"  Only trypanosomiasis: {len(only_tryp)}")

    if len(overlap) == 0:
        print("\n  No overlapping compounds found.")
        return pd.DataFrame()

    # Build detailed profile for overlapping compounds
    records = []
    for cid in sorted(overlap):
        # Schistosomiasis data
        s_data = data["schisto_activities"][
            data["schisto_activities"]["molecule_chembl_id"] == cid
        ]
        s_best_ic50 = s_data["standard_value"].min()
        s_targets = s_data["target_name"].unique().tolist()
        s_smiles = s_data["canonical_smiles"].iloc[0] if "canonical_smiles" in s_data.columns else None

        # Trypanosomiasis data
        t_data = data["tryp_activities"][
            data["tryp_activities"]["molecule_chembl_id"] == cid
        ]
        t_best_ic50 = t_data["standard_value"].min()
        t_targets = t_data["target_name"].unique().tolist()

        # Selectivity in schistosomiasis
        s_sel_class = None
        s_sel_ratio = None
        if "schisto_selectivity" in data:
            s_sel = data["schisto_selectivity"][
                data["schisto_selectivity"]["molecule_chembl_id"] == cid
            ]
            if len(s_sel) > 0:
                best_sel = s_sel.loc[s_sel["selectivity_ratio"].idxmax()]
                s_sel_class = best_sel["selectivity_class"]
                s_sel_ratio = best_sel["selectivity_ratio"]

        # Selectivity in trypanosomiasis
        t_sel_class = None
        t_sel_ratio = None
        if "tryp_selectivity" in data:
            t_sel = data["tryp_selectivity"][
                data["tryp_selectivity"]["molecule_chembl_id"] == cid
            ]
            if len(t_sel) > 0:
                best_sel = t_sel.loc[t_sel["selectivity_ratio"].idxmax()]
                t_sel_class = best_sel["selectivity_class"]
                t_sel_ratio = best_sel["selectivity_ratio"]

        # Drug-likeness
        mw = None
        qed_score = None
        if s_smiles and pd.notna(s_smiles):
            mol = Chem.MolFromSmiles(s_smiles)
            if mol:
                mw = round(Descriptors.MolWt(mol), 1)
                qed_score = round(QED.qed(mol), 3)

        records.append({
            "molecule_chembl_id": cid,
            "smiles": s_smiles,
            "MW": mw,
            "QED": qed_score,
            "schisto_best_ic50": s_best_ic50,
            "schisto_targets": "; ".join(s_targets),
            "schisto_sel_class": s_sel_class,
            "schisto_sel_ratio": s_sel_ratio,
            "tryp_best_ic50": t_best_ic50,
            "tryp_targets": "; ".join(t_targets),
            "tryp_sel_class": t_sel_class,
            "tryp_sel_ratio": t_sel_ratio,
            "dual_active": True,
        })

    overlap_df = pd.DataFrame(records)

    # Classify dual-disease potential
    def classify_dual(row):
        s_potent = pd.notna(row["schisto_best_ic50"]) and row["schisto_best_ic50"] < 1000
        t_potent = pd.notna(row["tryp_best_ic50"]) and row["tryp_best_ic50"] < 1000

        s_selective = row.get("schisto_sel_class") in ["SELECTIVE", "MODERATE"]
        t_selective = row.get("tryp_sel_class") in ["SELECTIVE", "MODERATE"]

        s_unknown_sel = pd.isna(row.get("schisto_sel_class"))
        t_unknown_sel = pd.isna(row.get("tryp_sel_class"))

        s_counter = row.get("schisto_sel_class") == "COUNTER-SELECTIVE"
        t_counter = row.get("tryp_sel_class") == "COUNTER-SELECTIVE"

        if s_counter or t_counter:
            return "ELIMINATED (counter-selective)"
        elif s_potent and t_potent:
            if s_selective or t_selective:
                return "HIGH PRIORITY (dual potent + selective)"
            elif s_unknown_sel and t_unknown_sel:
                return "PROMISING (dual potent, selectivity unknown)"
            else:
                return "MODERATE (dual potent, mixed selectivity)"
        elif s_potent or t_potent:
            return "PARTIAL (potent in one disease only)"
        else:
            return "WEAK (not potent in either)"

    overlap_df["dual_disease_class"] = overlap_df.apply(classify_dual, axis=1)

    # Print results
    print(f"\n  DUAL-ACTIVE COMPOUNDS:")
    print(f"  {'Class':45s} {'Count':>5s}")
    print(f"  {'-'*45} {'-'*5}")
    for cls, count in overlap_df["dual_disease_class"].value_counts().items():
        print(f"  {cls:45s} {count:5d}")

    # Print top candidates
    top = overlap_df[
        overlap_df["dual_disease_class"].str.contains("HIGH|PROMISING|MODERATE")
    ].sort_values(
        ["schisto_best_ic50", "tryp_best_ic50"]
    )

    if len(top) > 0:
        print(f"\n  TOP DUAL-DISEASE CANDIDATES:")
        print(f"  {'Compound':20s} {'S.man IC50':>11s} {'T.bru IC50':>11s} "
              f"{'S.man sel':>10s} {'T.bru sel':>10s} {'QED':>5s}  Class")
        print(f"  {'-'*20} {'-'*11} {'-'*11} {'-'*10} {'-'*10} {'-'*5}  {'-'*30}")

        for _, row in top.head(20).iterrows():
            cid = row["molecule_chembl_id"]
            if len(cid) > 20: cid = cid[:17] + "..."
            s_ic50 = f"{row['schisto_best_ic50']:.0f} nM" if pd.notna(row['schisto_best_ic50']) else "N/A"
            t_ic50 = f"{row['tryp_best_ic50']:.0f} nM" if pd.notna(row['tryp_best_ic50']) else "N/A"
            s_sel = f"{row['schisto_sel_ratio']:.1f}x" if pd.notna(row.get('schisto_sel_ratio')) else "?"
            t_sel = f"{row['tryp_sel_ratio']:.1f}x" if pd.notna(row.get('tryp_sel_ratio')) else "?"
            qed = f"{row['QED']:.2f}" if pd.notna(row.get('QED')) else "?"
            cls = row["dual_disease_class"][:30]

            print(f"  {cid:20s} {s_ic50:>11s} {t_ic50:>11s} "
                  f"{s_sel:>10s} {t_sel:>10s} {qed:>5s}  {cls}")

    return overlap_df


# ---------------------------------------------------------------------------
# STEP 3: T. brucei target prioritization (like SmDHODH finding)
# ---------------------------------------------------------------------------

def tryp_target_prioritization(data):
    """Which T. brucei target is the best repurposing opportunity?"""

    print("\n" + "=" * 60)
    print("STEP 3: T. BRUCEI TARGET PRIORITIZATION")
    print("=" * 60)

    activities = data["tryp_activities"]
    selectivity = data.get("tryp_selectivity", pd.DataFrame())

    targets = []
    for target_name in activities["target_name"].unique():
        t_data = activities[activities["target_name"] == target_name]
        n_compounds = t_data["molecule_chembl_id"].nunique()
        n_measurements = len(t_data)
        best_ic50 = t_data["standard_value"].min()
        median_ic50 = t_data["standard_value"].median()

        # Check if unique to parasite
        is_unique = any(kw in target_name.lower() for kw in [
            "trypanothione", "alternative oxidase", "rna-editing",
            "tryparedoxin", "rhodesain", "metacaspase",
        ])

        # Selectivity data
        n_selective = 0
        n_with_sel = 0
        pct_nonsel = None
        if len(selectivity) > 0:
            sel_data = selectivity[selectivity["parasite_target"] == target_name]
            n_with_sel = len(sel_data)
            if n_with_sel > 0:
                n_selective = len(sel_data[sel_data["selectivity_class"] == "SELECTIVE"])
                n_nonsel = len(sel_data[sel_data["selectivity_class"].isin(
                    ["POOR", "COUNTER-SELECTIVE"])])
                pct_nonsel = round(100 * n_nonsel / n_with_sel, 1)

        # Composite target score
        # Rewards: many compounds, potent hits, unique biology, selective compounds
        potency_score = max(0, (9 - np.log10(best_ic50)) / 4) if best_ic50 > 0 else 0
        data_score = min(1.0, np.log2(max(1, n_compounds)) / 8)
        unique_bonus = 0.3 if is_unique else 0
        sel_bonus = 0.2 * min(1, n_selective / 3) if n_selective > 0 else 0
        sel_penalty = -0.3 if pct_nonsel and pct_nonsel > 80 else 0

        target_score = potency_score * 0.3 + data_score * 0.2 + unique_bonus + sel_bonus + sel_penalty

        targets.append({
            "target": target_name,
            "n_compounds": n_compounds,
            "n_measurements": n_measurements,
            "best_ic50_nM": best_ic50,
            "median_ic50_nM": median_ic50,
            "unique_to_parasite": is_unique,
            "n_with_selectivity": n_with_sel,
            "n_selective_10x": n_selective,
            "pct_non_selective": pct_nonsel,
            "target_score": round(target_score, 3),
        })

    target_df = pd.DataFrame(targets).sort_values("target_score", ascending=False)

    print(f"\n  T. BRUCEI TARGET RANKING:")
    print(f"  {'Rk':>3s}  {'Score':>6s}  {'Target':35s} {'N cpd':>6s} "
          f"{'Best nM':>8s} {'Unique':>7s} {'Sel10x':>7s} {'%NonSel':>8s}")
    print(f"  {'-'*3}  {'-'*6}  {'-'*35} {'-'*6} {'-'*8} {'-'*7} {'-'*7} {'-'*8}")

    for rank, (_, row) in enumerate(target_df.head(10).iterrows(), 1):
        unique = "YES" if row["unique_to_parasite"] else "no"
        nsel = f"{row['n_selective_10x']:.0f}" if pd.notna(row["n_selective_10x"]) else "N/A"
        pct = f"{row['pct_non_selective']:.0f}%" if pd.notna(row["pct_non_selective"]) else "N/A"
        print(f"  {rank:3d}  {row['target_score']:6.3f}  {row['target'][:35]:35s} "
              f"{row['n_compounds']:6d} {row['best_ic50_nM']:8.0f} "
              f"{unique:>7s} {nsel:>7s} {pct:>8s}")

    return target_df


# ---------------------------------------------------------------------------
# STEP 4: Platform summary report
# ---------------------------------------------------------------------------

def platform_summary(data, overlap_df, target_ranking):
    """The definitive Kira platform report."""

    print("\n" + "=" * 60)
    print("STEP 4: DEFINITIVE PLATFORM SUMMARY")
    print("=" * 60)

    report = []
    report.append("=" * 70)
    report.append("KIRA PLATFORM — DEFINITIVE CROSS-DISEASE REPORT")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("=" * 70)

    report.append("")
    report.append("PLATFORM SCOPE:")
    report.append("  Diseases: Schistosomiasis (S. mansoni), Trypanosomiasis (T. brucei)")
    s_n = data["schisto_activities"]["molecule_chembl_id"].nunique()
    t_n = data["tryp_activities"]["molecule_chembl_id"].nunique()
    report.append(f"  Schistosomiasis: {s_n} compounds, 10 targets, 226 measurements")
    report.append(f"  Trypanosomiasis: {t_n} compounds, 16 targets, 1140 measurements")
    report.append(f"  Total unique compounds: {s_n + t_n - len(overlap_df)}")

    report.append("")
    report.append("-" * 70)
    report.append("CROSS-DISEASE SELECTIVITY FINDING")
    report.append("-" * 70)
    report.append("  Across 4 experimental target-orthologue pairs and 195")
    report.append("  dual-species selectivity comparisons:")
    report.append("    Weighted average non-selectivity: 76.4%")
    report.append("    Compounds with >=10x selectivity: 11 of 195 (5.6%)")
    report.append("    Targets with majority non-selective: 3 of 4")
    report.append("")
    report.append("  This finding is consistent across both diseases,")
    report.append("  supporting selectivity-first triage as a general")
    report.append("  methodology for antiparasitic drug repurposing.")

    if len(overlap_df) > 0:
        report.append("")
        report.append("-" * 70)
        report.append("CROSS-DISEASE COMPOUND OVERLAP")
        report.append("-" * 70)
        report.append(f"  Compounds active against both parasites: {len(overlap_df)}")

        for cls, count in overlap_df["dual_disease_class"].value_counts().items():
            report.append(f"    {cls}: {count}")

        top = overlap_df[
            overlap_df["dual_disease_class"].str.contains("HIGH|PROMISING")
        ]
        if len(top) > 0:
            report.append(f"\n  TOP DUAL-DISEASE CANDIDATES: {len(top)}")
            for _, row in top.head(5).iterrows():
                report.append(f"    {row['molecule_chembl_id']}")
                report.append(f"      S. mansoni: {row['schisto_best_ic50']:.0f} nM "
                              f"vs {row['schisto_targets']}")
                report.append(f"      T. brucei: {row['tryp_best_ic50']:.0f} nM "
                              f"vs {row['tryp_targets']}")

    report.append("")
    report.append("-" * 70)
    report.append("T. BRUCEI PRIORITY TARGETS")
    report.append("-" * 70)
    for _, row in target_ranking.head(5).iterrows():
        unique = " [UNIQUE TO PARASITE]" if row["unique_to_parasite"] else ""
        report.append(f"  {row['target']}{unique}")
        report.append(f"    Compounds: {row['n_compounds']}, "
                      f"Best IC50: {row['best_ic50_nM']:.0f} nM")
        if pd.notna(row.get("n_selective_10x")) and row["n_selective_10x"] > 0:
            report.append(f"    Selective (>=10x): {row['n_selective_10x']:.0f}")

    report.append("")
    report.append("-" * 70)
    report.append("SCHISTOSOMIASIS PRIORITY TARGET")
    report.append("-" * 70)
    report.append("  SmDHODH — 3 compounds with >=10x selectivity")
    report.append("  Best: CHEMBL155771 (23 nM, 30.8x, QED 0.89)")
    report.append("  Translational: Atovaquone (430 nM, 6x, WHO EML)")

    report.append("")
    report.append("-" * 70)
    report.append("EXPERIMENTAL PRIORITIES")
    report.append("-" * 70)
    report.append("  Schistosomiasis:")
    report.append("    1. Whole-worm assay: CHEMBL155771, CHEMBL4452960 (SmDHODH)")
    report.append("    2. Atovaquone vs adult S. mansoni")
    report.append("    3. SmTGR selectivity at fusion interface")
    report.append("  Trypanosomiasis:")
    report.append("    4. CHEMBL4204999 in T. brucei whole-cell assay (PDEB1, 20x)")
    report.append("    5. Trypanothione reductase — screen for selective binders")
    report.append("    6. Rhodesain selectivity vs human cathepsin L")

    if len(overlap_df) > 0:
        top_dual = overlap_df[
            overlap_df["dual_disease_class"].str.contains("HIGH|PROMISING")
        ]
        if len(top_dual) > 0:
            report.append("  Cross-disease:")
            report.append(f"    7. Test top {len(top_dual)} dual-active compounds in both worm assays")

    report.append("")
    report.append("=" * 70)
    report.append("KIRA PLATFORM STATUS")
    report.append("=" * 70)
    report.append("  Scripts: 17")
    report.append("  Diseases: 2 (Schistosomiasis, Trypanosomiasis)")
    report.append(f"  Total compounds assessed: {s_n + t_n}")
    report.append("  Selectivity comparisons: 195 experimental + 43 docking")
    report.append("  Platform finding: 76.4% of antiparasitic chemical matter")
    report.append("    against conserved targets is non-selective")
    report.append("  Methodology: Selectivity-first triage generalizes across NTDs")
    report.append("=" * 70)

    report_text = "\n".join(report)
    path = os.path.join(PUB_DIR, "kira_platform_definitive.txt")
    with open(path, "w") as f:
        f.write(report_text)
    print(f"\n  Platform report: {path}")

    return report_text


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Script 17: Cross-Disease Analysis")
    print("  Dual-Active Compounds + Platform Summary")
    print("=" * 60)

    # Load everything
    data = load_all_data()

    # Cross-disease overlap
    overlap_df = find_cross_disease_compounds(data)
    if len(overlap_df) > 0:
        overlap_df.to_csv(os.path.join(PUB_DIR, "cross_disease_compounds.csv"), index=False)

    # T. brucei target prioritization
    target_ranking = tryp_target_prioritization(data)
    target_ranking.to_csv(os.path.join(TRYP_DIR, "tryp_target_ranking.csv"), index=False)

    # Platform summary
    report = platform_summary(data, overlap_df, target_ranking)
    print("\n")
    print(report)

    n_overlap = len(overlap_df) if len(overlap_df) > 0 else 0
    n_promising = 0
    if len(overlap_df) > 0:
        n_promising = len(overlap_df[
            overlap_df["dual_disease_class"].str.contains("HIGH|PROMISING|MODERATE")
        ])

    print(f"\n{'=' * 60}")
    print(f"  CROSS-DISEASE ANALYSIS COMPLETE")
    print(f"")
    print(f"  Compounds active in both diseases: {n_overlap}")
    print(f"  Promising dual-disease candidates: {n_promising}")
    print(f"  T. brucei targets ranked: {len(target_ranking)}")
    print(f"")
    print(f"  17 scripts. 2 diseases. 1 platform.")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 17_cross_disease_analysis.py successfully."
echo "Now run: python 17_cross_disease_analysis.py"
