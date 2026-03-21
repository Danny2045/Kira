"""
Kira - Script 12: Publication-Ready Selectivity Analysis
==========================================================

This script packages the selectivity finding into publication form.

Core claim: "Systematic cross-species selectivity analysis of 91
anti-schistosomal compounds reveals that SmHDAC8 chemical matter is
overwhelmingly non-selective (79%), while SmDHODH harbors compounds
with >10x parasite selectivity."

What this script produces:
1. Resolved compound identities (no more "nan" names)
2. Per-target selectivity breakdown with statistics
3. Hard negatives: structurally similar but inactive/non-selective compounds
4. Sanity anchor analysis (praziquantel, atovaquone as reference points)
5. Discovery tier vs translational tier separation
6. Publication-ready summary tables
7. Formal limitations and scope statement

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 12_publication_analysis.py
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs, QED, Draw
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")
REPORT_DIR = os.path.join(os.path.dirname(__file__), "data", "reports")
PUB_DIR = os.path.join(os.path.dirname(__file__), "data", "publication")


# ---------------------------------------------------------------------------
# STEP 1: Resolve all compound identities
# ---------------------------------------------------------------------------

def resolve_compound_identities():
    """
    Fix the 'nan' name problem. For every compound in the pipeline,
    establish: ChEMBL ID, preferred name (if any), SMILES, molecular
    formula, and which targets it hits.
    """
    print("\n" + "=" * 60)
    print("STEP 1: RESOLVING COMPOUND IDENTITIES")
    print("=" * 60)

    activities = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))
    sel_path = os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv")
    selectivity = pd.read_csv(sel_path) if os.path.exists(sel_path) else pd.DataFrame()

    # Build master compound registry
    compounds = {}

    for _, row in activities.iterrows():
        cid = row["molecule_chembl_id"]
        if cid not in compounds:
            compounds[cid] = {
                "molecule_chembl_id": cid,
                "pref_name": None,
                "smiles": None,
                "targets": set(),
                "best_ic50": {},
                "selectivity": {},
            }

        target = row["target_name"]
        compounds[cid]["targets"].add(target)

        smi = row.get("canonical_smiles")
        if pd.notna(smi):
            compounds[cid]["smiles"] = smi

        # Track best IC50 per target
        val = row.get("standard_value")
        if pd.notna(val):
            current = compounds[cid]["best_ic50"].get(target, float("inf"))
            if val < current:
                compounds[cid]["best_ic50"][target] = val

    # Add names from various sources
    v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv")
    if os.path.exists(v3_path):
        v3 = pd.read_csv(v3_path)
        for _, row in v3.iterrows():
            cid = row["molecule_chembl_id"]
            name = row.get("pref_name")
            if cid in compounds and pd.notna(name):
                compounds[cid]["pref_name"] = name

    # Add selectivity data
    if len(selectivity) > 0:
        for _, row in selectivity.iterrows():
            cid = row["molecule_chembl_id"]
            if cid in compounds:
                target = row["parasite_target"]
                compounds[cid]["selectivity"][target] = {
                    "ratio": row["selectivity_ratio"],
                    "parasite_ic50": row["parasite_best_ic50"],
                    "human_ic50": row["human_best_ic50"],
                    "class": row["selectivity_class"],
                }

    # Compute molecular properties for unnamed compounds
    for cid, data in compounds.items():
        smi = data["smiles"]
        if smi and not data["pref_name"]:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mw = Descriptors.MolWt(mol)
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                data["molecular_formula"] = formula
                data["MW"] = round(mw, 1)

    # Create display name: prefer pref_name, else ChEMBL ID + formula
    for cid, data in compounds.items():
        if data["pref_name"]:
            data["display_name"] = data["pref_name"]
        elif data.get("molecular_formula"):
            data["display_name"] = f"{cid} ({data['molecular_formula']})"
        else:
            data["display_name"] = cid

    print(f"  Resolved {len(compounds)} compounds")
    named = sum(1 for c in compounds.values() if c["pref_name"])
    print(f"  With preferred names: {named}")
    print(f"  With SMILES: {sum(1 for c in compounds.values() if c['smiles'])}")
    print(f"  With selectivity data: {sum(1 for c in compounds.values() if c['selectivity'])}")

    return compounds


# ---------------------------------------------------------------------------
# STEP 2: Per-target selectivity statistics
# ---------------------------------------------------------------------------

def per_target_selectivity_analysis(compounds):
    """
    For each parasite target, compute:
    - Total compounds tested
    - Compounds with selectivity data
    - Distribution: selective / moderate / poor / counter-selective
    - Median selectivity ratio
    - Key statistics for publication
    """
    print("\n" + "=" * 60)
    print("STEP 2: PER-TARGET SELECTIVITY BREAKDOWN")
    print("=" * 60)

    target_stats = {}

    # Collect selectivity data by target
    target_sel_data = defaultdict(list)
    target_compounds = defaultdict(set)

    for cid, data in compounds.items():
        for target in data["targets"]:
            target_compounds[target].add(cid)
        for target, sel in data["selectivity"].items():
            target_sel_data[target].append({
                "compound": cid,
                "display_name": data["display_name"],
                "ratio": sel["ratio"],
                "class": sel["class"],
                "parasite_ic50": sel["parasite_ic50"],
                "human_ic50": sel["human_ic50"],
            })

    for target in sorted(target_compounds.keys()):
        n_total = len(target_compounds[target])
        sel_data = target_sel_data.get(target, [])
        n_with_sel = len(sel_data)

        stats = {
            "target": target,
            "n_total_compounds": n_total,
            "n_with_selectivity": n_with_sel,
            "n_selective": 0,
            "n_moderate": 0,
            "n_poor": 0,
            "n_counter": 0,
            "median_ratio": None,
            "mean_ratio": None,
            "best_ratio": None,
            "best_compound": None,
            "pct_non_selective": None,
        }

        if n_with_sel > 0:
            ratios = [s["ratio"] for s in sel_data]
            stats["median_ratio"] = round(np.median(ratios), 2)
            stats["mean_ratio"] = round(np.mean(ratios), 2)
            stats["best_ratio"] = round(max(ratios), 2)

            best_entry = max(sel_data, key=lambda x: x["ratio"])
            stats["best_compound"] = best_entry["display_name"]

            for s in sel_data:
                if s["class"] == "SELECTIVE":
                    stats["n_selective"] += 1
                elif s["class"] == "MODERATE":
                    stats["n_moderate"] += 1
                elif s["class"] == "POOR":
                    stats["n_poor"] += 1
                elif s["class"] == "COUNTER-SELECTIVE":
                    stats["n_counter"] += 1

            n_problematic = stats["n_poor"] + stats["n_counter"]
            stats["pct_non_selective"] = round(100 * n_problematic / n_with_sel, 1)

        target_stats[target] = stats

        # Print
        print(f"\n  {target}")
        print(f"    Total compounds: {n_total}")
        print(f"    With selectivity data: {n_with_sel}")
        if n_with_sel > 0:
            print(f"    Selective (>=10x): {stats['n_selective']}")
            print(f"    Moderate (3-10x): {stats['n_moderate']}")
            print(f"    Poor (1-3x): {stats['n_poor']}")
            print(f"    Counter-selective (<1x): {stats['n_counter']}")
            print(f"    % non-selective (poor+counter): {stats['pct_non_selective']}%")
            print(f"    Median ratio: {stats['median_ratio']}x")
            print(f"    Best: {stats['best_compound']} at {stats['best_ratio']}x")
        else:
            print(f"    NO SELECTIVITY DATA AVAILABLE")

    return target_stats


# ---------------------------------------------------------------------------
# STEP 3: Build hard negatives from structurally similar non-selective compounds
# ---------------------------------------------------------------------------

def build_hard_negatives(compounds):
    """
    Hard negatives are compounds that are structurally SIMILAR to actives
    but are non-selective or counter-selective. These are the hardest
    cases for a ranking algorithm to handle correctly.

    They test: "Can the pipeline tell the difference between a compound
    that looks right and a compound that actually IS right?"
    """
    print("\n" + "=" * 60)
    print("STEP 3: IDENTIFYING HARD NEGATIVES")
    print("=" * 60)

    # Get fingerprints for selective and counter-selective compounds
    selective_fps = {}
    counter_fps = {}

    for cid, data in compounds.items():
        smi = data.get("smiles")
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        has_selective = any(s["class"] == "SELECTIVE" for s in data["selectivity"].values())
        has_counter = any(s["class"] == "COUNTER-SELECTIVE" for s in data["selectivity"].values())

        if has_selective:
            selective_fps[cid] = fp
        elif has_counter:
            counter_fps[cid] = fp

    print(f"  Selective compounds with fingerprints: {len(selective_fps)}")
    print(f"  Counter-selective compounds with fingerprints: {len(counter_fps)}")

    if not selective_fps or not counter_fps:
        print("  Cannot compute hard negatives: need both selective and counter-selective compounds")
        return []

    # For each counter-selective compound, find its similarity to the nearest selective
    hard_negatives = []
    sel_fps_list = list(selective_fps.values())
    sel_ids_list = list(selective_fps.keys())

    for cid, fp in counter_fps.items():
        sims = DataStructs.BulkTanimotoSimilarity(fp, sel_fps_list)
        max_sim = max(sims)
        nearest_sel = sel_ids_list[sims.index(max_sim)]

        data = compounds[cid]
        nearest_data = compounds[nearest_sel]

        hard_negatives.append({
            "compound": cid,
            "display_name": data["display_name"],
            "similarity_to_nearest_selective": round(max_sim, 3),
            "nearest_selective": nearest_sel,
            "nearest_selective_name": nearest_data["display_name"],
            "selectivity_ratio": min(s["ratio"] for s in data["selectivity"].values()),
        })

    hard_negatives.sort(key=lambda x: x["similarity_to_nearest_selective"], reverse=True)

    print(f"\n  Hard negatives (counter-selective compounds similar to selectives):")
    print(f"  {'Compound':35s} {'Sim':>5s} {'Ratio':>6s}  Nearest selective")
    print(f"  {'-'*35} {'-'*5} {'-'*6}  {'-'*35}")

    for hn in hard_negatives[:15]:
        print(f"  {hn['display_name']:35s} {hn['similarity_to_nearest_selective']:5.3f} "
              f"{hn['selectivity_ratio']:6.2f}  {hn['nearest_selective_name']}")

    # Highlight the hardest cases: high similarity but terrible selectivity
    hardest = [h for h in hard_negatives if h["similarity_to_nearest_selective"] > 0.3]
    print(f"\n  Compounds with Tanimoto > 0.3 to a selective compound but counter-selective: {len(hardest)}")
    print(f"  These are the hardest cases for any ranking algorithm.")

    return hard_negatives


# ---------------------------------------------------------------------------
# STEP 4: Sanity anchor analysis
# ---------------------------------------------------------------------------

def sanity_anchor_analysis(compounds):
    """
    Praziquantel and atovaquone are sanity anchors: drugs we KNOW work.
    Any credible pipeline must handle them correctly.

    This section documents exactly how the pipeline treats them and why.
    """
    print("\n" + "=" * 60)
    print("STEP 4: SANITY ANCHOR ANALYSIS")
    print("=" * 60)

    anchors = {
        "PRAZIQUANTEL": {
            "expected": "Should rank in top 25% despite no target-based IC50",
            "mechanism": "Calcium channel disruption → worm paralysis",
            "evidence_type": "Whole-organism (77 nM), structural similarity (0.524)",
            "selectivity_status": "Unknown — no human orthologue data for target",
            "clinical_status": "Standard of care, 4030 publications, available globally",
        },
        "ATOVAQUONE": {
            "expected": "Should rank in top 10 among approved drugs",
            "mechanism": "SmDHODH inhibition → pyrimidine synthesis block",
            "evidence_type": "Target-based IC50 (430 nM), 6x selective over human DHODH",
            "selectivity_status": "MODERATE (6.0x)",
            "clinical_status": "Approved antimalarial, available in Africa as Malarone",
        },
    }

    # Load shortlist to find actual ranks
    shortlist_path = os.path.join(PROCESSED_DIR, "kira_shortlist_v2.csv")
    if os.path.exists(shortlist_path):
        shortlist = pd.read_csv(shortlist_path)
        shortlist = shortlist.sort_values("final_score", ascending=False).reset_index(drop=True)
        shortlist["rank"] = range(1, len(shortlist) + 1)
        total = len(shortlist)

        for name, info in anchors.items():
            match = shortlist[shortlist["pref_name"] == name]
            if len(match) > 0:
                row = match.iloc[0]
                rank = int(row["rank"])
                pct = round((1 - rank / total) * 100)
                info["actual_rank"] = f"{rank}/{total} (top {pct}%)"
                info["actual_score"] = round(row["final_score"], 4)
            else:
                info["actual_rank"] = "Not in shortlist"
                info["actual_score"] = None

    for name, info in anchors.items():
        print(f"\n  {name}")
        for key, val in info.items():
            print(f"    {key}: {val}")


# ---------------------------------------------------------------------------
# STEP 5: Discovery tier vs translational tier
# ---------------------------------------------------------------------------

def split_tiers(compounds):
    """
    Separate compounds into two tiers:

    DISCOVERY TIER: Research compounds with potent and selective activity
    against parasite targets. These need synthesis, optimization, and
    years of development. Audience: medicinal chemists, academic labs.

    TRANSLATIONAL TIER: Approved or clinical-stage drugs with evidence
    of anti-schistosomal activity. These can potentially enter clinical
    evaluation much faster. Audience: clinical researchers, WHO, NTD
    drug access organizations.
    """
    print("\n" + "=" * 60)
    print("STEP 5: DISCOVERY vs TRANSLATIONAL TIERS")
    print("=" * 60)

    # Load shortlist
    shortlist_path = os.path.join(PROCESSED_DIR, "kira_shortlist_v2.csv")
    shortlist = pd.read_csv(shortlist_path) if os.path.exists(shortlist_path) else pd.DataFrame()

    if len(shortlist) == 0:
        print("  No shortlist available.")
        return [], []

    shortlist["max_phase"] = pd.to_numeric(shortlist.get("max_phase", pd.Series()), errors="coerce")

    # Translational: approved or clinical-stage (phase >= 1)
    translational = shortlist[shortlist["max_phase"] >= 1].sort_values("final_score", ascending=False)

    # Discovery: research compounds with good selectivity or unknown (SmTGR)
    discovery = shortlist[
        (shortlist["max_phase"].isna() | (shortlist["max_phase"] < 1)) &
        (shortlist["selectivity_class"].isin(["SELECTIVE", "MODERATE"]) |
         shortlist["selectivity_class"].isna())
    ].sort_values("final_score", ascending=False)

    print(f"\n  TRANSLATIONAL TIER: {len(translational)} compounds (approved/clinical)")
    for rank, (_, row) in enumerate(translational.head(10).iterrows(), 1):
        name = row.get("pref_name", row["molecule_chembl_id"])
        if pd.isna(name): name = row["molecule_chembl_id"]
        phase = int(row["max_phase"]) if pd.notna(row.get("max_phase")) else 0
        phase_label = {4: "APPROVED", 3: "PhIII", 2: "PhII", 1: "PhI"}.get(phase, "?")
        sel = row.get("selectivity_class", "UNKNOWN")
        if pd.isna(sel): sel = "UNKNOWN"
        ic50 = f"{row['best_value']:.0f}" if pd.notna(row.get("best_value")) else "N/A"
        print(f"    {rank:2d}. {name:25s} {phase_label:8s} IC50={ic50:>7s}  sel={sel}")

    print(f"\n  DISCOVERY TIER: {len(discovery)} compounds (research, selective or unknown)")
    for rank, (_, row) in enumerate(discovery.head(10).iterrows(), 1):
        name = row.get("pref_name", row["molecule_chembl_id"])
        if pd.isna(name): name = row["molecule_chembl_id"]
        sel = row.get("selectivity_class", "UNKNOWN")
        if pd.isna(sel): sel = "UNKNOWN"
        ic50 = f"{row['best_value']:.0f}" if pd.notna(row.get("best_value")) else "N/A"
        target = str(row.get("best_target", ""))[:25]
        print(f"    {rank:2d}. {name:25s} IC50={ic50:>7s}  sel={sel:10s}  {target}")

    return translational, discovery


# ---------------------------------------------------------------------------
# STEP 6: Generate publication summary tables
# ---------------------------------------------------------------------------

def generate_publication_tables(target_stats, compounds, hard_negatives):
    """Generate clean tables suitable for a manuscript."""

    os.makedirs(PUB_DIR, exist_ok=True)

    # Table 1: Target selectivity summary
    rows = []
    for target, stats in target_stats.items():
        rows.append({
            "Target": target,
            "N_compounds": stats["n_total_compounds"],
            "N_with_selectivity": stats["n_with_selectivity"],
            "Selective_ge10x": stats["n_selective"],
            "Moderate_3_10x": stats["n_moderate"],
            "Poor_1_3x": stats["n_poor"],
            "Counter_lt1x": stats["n_counter"],
            "Pct_non_selective": stats["pct_non_selective"],
            "Median_ratio": stats["median_ratio"],
            "Best_ratio": stats["best_ratio"],
        })

    table1 = pd.DataFrame(rows)
    table1.to_csv(os.path.join(PUB_DIR, "table1_target_selectivity.csv"), index=False)

    # Table 2: Top selective compounds
    top_selective = []
    for cid, data in compounds.items():
        for target, sel in data["selectivity"].items():
            if sel["class"] in ["SELECTIVE", "MODERATE"]:
                mol = Chem.MolFromSmiles(data["smiles"]) if data.get("smiles") else None
                top_selective.append({
                    "ChEMBL_ID": cid,
                    "Name": data["display_name"],
                    "Target": target,
                    "Parasite_IC50_nM": sel["parasite_ic50"],
                    "Human_IC50_nM": sel["human_ic50"],
                    "Selectivity_ratio": round(sel["ratio"], 1),
                    "Class": sel["class"],
                    "MW": round(Descriptors.MolWt(mol), 1) if mol else None,
                    "QED": round(QED.qed(mol), 3) if mol else None,
                    "SMILES": data.get("smiles", ""),
                })

    table2 = pd.DataFrame(top_selective)
    table2 = table2.sort_values("Selectivity_ratio", ascending=False)
    table2.to_csv(os.path.join(PUB_DIR, "table2_selective_compounds.csv"), index=False)

    # Table 3: Hard negatives
    if hard_negatives:
        table3 = pd.DataFrame(hard_negatives)
        table3.to_csv(os.path.join(PUB_DIR, "table3_hard_negatives.csv"), index=False)

    print(f"\n  Publication tables saved to {PUB_DIR}/")
    print(f"    table1_target_selectivity.csv")
    print(f"    table2_selective_compounds.csv")
    print(f"    table3_hard_negatives.csv")

    return table1, table2


# ---------------------------------------------------------------------------
# STEP 7: Generate publication-ready report
# ---------------------------------------------------------------------------

def generate_publication_report(target_stats, compounds, hard_negatives,
                                 translational, discovery):
    """The report a supervisor or reviewer would read."""

    report = []
    report.append("=" * 70)
    report.append("KIRA — PUBLICATION-READY ANALYSIS")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("=" * 70)

    report.append("")
    report.append("TITLE (working):")
    report.append("  Systematic cross-species selectivity analysis reveals SmDHODH")
    report.append("  as the most tractable target for schistosomiasis drug repurposing")

    report.append("")
    report.append("ABSTRACT (working):")
    report.append("  Drug repurposing offers a fast path to new schistosomiasis")
    report.append("  treatments, but candidate prioritization requires assessment of")
    report.append("  selectivity against human orthologues — a step frequently omitted")
    report.append("  in computational repurposing studies. We built an integrated pipeline")
    report.append("  combining target-based bioactivity (ChEMBL), whole-organism evidence,")
    report.append("  structural similarity, ADMET properties, and cross-species selectivity")
    report.append("  to systematically evaluate the schistosomiasis repurposing landscape.")

    # Core finding
    hdac_stats = target_stats.get("Histone deacetylase 8", {})
    dhodh_stats = target_stats.get("Dihydroorotate dehydrogenase (quinone), mitochondrial", {})

    if hdac_stats.get("pct_non_selective"):
        report.append(f"  Of {hdac_stats['n_with_selectivity']} SmHDAC8 compounds with human orthologue")
        report.append(f"  data, {hdac_stats['pct_non_selective']}% showed poor or counter-selectivity,")
        report.append(f"  including all clinically advanced HDAC inhibitors (vorinostat,")
        report.append(f"  panobinostat, quisinostat).")

    if dhodh_stats.get("n_selective"):
        report.append(f"  In contrast, SmDHODH yielded {dhodh_stats['n_selective']} compounds with")
        report.append(f"  >10-fold parasite selectivity, with the best achieving")
        report.append(f"  {dhodh_stats['best_ratio']}x (IC50 = 23 nM, QED = 0.89).")

    report.append("  Atovaquone, an approved antimalarial available in sub-Saharan Africa,")
    report.append("  showed 6-fold SmDHODH selectivity, supporting its prioritization")
    report.append("  for anti-schistosomal evaluation.")

    report.append("")
    report.append("-" * 70)
    report.append("KEY RESULTS")
    report.append("-" * 70)

    report.append("")
    report.append("1. TARGET SELECTIVITY LANDSCAPE")
    for target, stats in sorted(target_stats.items()):
        report.append(f"")
        report.append(f"   {target}")
        report.append(f"   Compounds tested: {stats['n_total_compounds']}")
        report.append(f"   With selectivity data: {stats['n_with_selectivity']}")
        if stats["n_with_selectivity"] > 0:
            report.append(f"   Selective (>=10x): {stats['n_selective']}, "
                          f"Moderate (3-10x): {stats['n_moderate']}, "
                          f"Poor (1-3x): {stats['n_poor']}, "
                          f"Counter (<1x): {stats['n_counter']}")
            report.append(f"   Non-selective rate: {stats['pct_non_selective']}%")
            report.append(f"   Median ratio: {stats['median_ratio']}x, "
                          f"Best: {stats['best_ratio']}x ({stats['best_compound']})")
        else:
            report.append(f"   SELECTIVITY UNKNOWN — experimental gap identified")

    report.append("")
    report.append("2. TOP SELECTIVE COMPOUNDS (>=10x)")
    for cid, data in compounds.items():
        for target, sel in data["selectivity"].items():
            if sel["class"] == "SELECTIVE":
                report.append(f"   {data['display_name']}")
                report.append(f"     vs {target}: {sel['parasite_ic50']:.0f} nM parasite, "
                              f"{sel['human_ic50']:.0f} nM human, {sel['ratio']:.1f}x selective")

    if hard_negatives:
        n_hard = len([h for h in hard_negatives if h["similarity_to_nearest_selective"] > 0.3])
        report.append("")
        report.append(f"3. HARD NEGATIVES")
        report.append(f"   {n_hard} counter-selective compounds with Tanimoto > 0.3")
        report.append(f"   to a selective compound. These are structurally similar")
        report.append(f"   but non-selective — the hardest discrimination cases.")

    report.append("")
    report.append("-" * 70)
    report.append("SCOPE AND LIMITATIONS")
    report.append("-" * 70)
    report.append("  1. Selectivity assessed for SmHDAC8 and SmDHODH only; SmTGR")
    report.append("     (the most biologically compelling target) lacks human orthologue data")
    report.append("  2. Selectivity ratios computed from ChEMBL data, not matched assay conditions")
    report.append("     (different labs, assay formats, and conditions may introduce systematic bias)")
    report.append("  3. Only direct orthologue selectivity assessed; off-target pharmacology")
    report.append("     against other human proteins not evaluated")
    report.append("  4. Whole-organism validation data available for 6 compounds only")
    report.append("  5. No structural modeling or docking performed")
    report.append("  6. AUROC metrics reflect easy discrimination (known actives vs unrelated drugs)")
    report.append("  7. Pipeline is a prioritization tool, not a clinical recommendation engine")

    report.append("")
    report.append("-" * 70)
    report.append("RECOMMENDED EXPERIMENTAL FOLLOW-UP")
    report.append("-" * 70)
    report.append("  Priority 1: Test CHEMBL155771 and CHEMBL4452960 in S. mansoni")
    report.append("    whole-worm killing assay to confirm phenotypic activity")
    report.append("  Priority 2: Verify atovaquone activity against adult S. mansoni")
    report.append("    worms (clinical formulation, relevant concentrations)")
    report.append("  Priority 3: Test top SmTGR inhibitors against human TrxR1")
    report.append("    to establish selectivity for this target class")
    report.append("  Priority 4: Structural analysis of SmHDAC8 selective compound")
    report.append("    (CHEMBL4855490, 11x) to understand basis of selectivity")

    report.append("")
    report.append("=" * 70)

    report_text = "\n".join(report)
    report_path = os.path.join(PUB_DIR, "publication_analysis.txt")
    with open(report_path, "w") as f:
        f.write(report_text)

    print(f"\n  Publication report: {report_path}")
    return report_text


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 12")
    print("  Publication-Ready Analysis")
    print("=" * 60)

    # Step 1: Resolve identities
    compounds = resolve_compound_identities()

    # Step 2: Per-target selectivity
    target_stats = per_target_selectivity_analysis(compounds)

    # Step 3: Hard negatives
    hard_negatives = build_hard_negatives(compounds)

    # Step 4: Sanity anchors
    sanity_anchor_analysis(compounds)

    # Step 5: Tier split
    translational, discovery = split_tiers(compounds)

    # Step 6: Publication tables
    table1, table2 = generate_publication_tables(target_stats, compounds, hard_negatives)

    # Step 7: Publication report
    report = generate_publication_report(target_stats, compounds, hard_negatives,
                                          translational, discovery)

    # Print the report
    print("\n")
    print(report)

    print(f"\n{'=' * 60}")
    print(f"  PUBLICATION ANALYSIS COMPLETE")
    print(f"")
    print(f"  All outputs in data/publication/:")
    print(f"    table1_target_selectivity.csv")
    print(f"    table2_selective_compounds.csv")
    print(f"    table3_hard_negatives.csv")
    print(f"    publication_analysis.txt")
    print(f"")
    print(f"  Working title:")
    print(f"  'Systematic cross-species selectivity analysis reveals SmDHODH")
    print(f"   as the most tractable target for schistosomiasis drug repurposing'")
    print(f"{'=' * 60}")
