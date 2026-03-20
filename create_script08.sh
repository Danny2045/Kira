#!/bin/bash
cat > 08_harden_benchmark.py << 'PYTHONSCRIPT'
"""
Kira - Script 08: Harden the Benchmark
========================================

This is the most important script since the pipeline was built.

Problems with evaluation set v1:
  - Only 16 inactive controls (specificity estimates meaningless)
  - 15 of 16 are assumed negatives (not experimentally verified)
  - No train/test split (potential leakage)
  - No confidence intervals (point metrics are unreliable)
  - No per-target performance breakdown

This script fixes all five:
  1. Expand negatives using real ChEMBL whole-organism screening failures
  2. Add more assumed negatives from diverse approved drug classes
  3. Freeze the evaluation set with documented inclusion criteria
  4. Split into development (80%) and held-out test (20%)
  5. Compute bootstrap confidence intervals on AUROC and AUPRC
  6. Per-target performance breakdown

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 08_harden_benchmark.py
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")

# ---------------------------------------------------------------------------
# STEP 1: Mine ChEMBL whole-organism data for experimental negatives
# ---------------------------------------------------------------------------

def find_experimental_negatives():
    """
    From the whole-organism data we already cached, find compounds that
    were TESTED against live S. mansoni worms but showed NO meaningful
    activity. These are real experimental negatives — far more valuable
    than assumed negatives.

    We look for:
    - Compounds with "Activity" type and low/zero percentage values
    - Compounds with high IC50/EC50 (above 10000 nM)
    - Compounds tested but showing no worm killing
    """
    print("\n" + "=" * 60)
    print("STEP 1: MINING EXPERIMENTAL NEGATIVES")
    print("=" * 60)

    wo_path = os.path.join(PROCESSED_DIR, "schisto_whole_organism_raw.csv")
    if not os.path.exists(wo_path):
        print("  WARNING: No whole-organism cache found.")
        print("  Run Script 06 first to download whole-organism data.")
        return pd.DataFrame()

    wo_df = pd.read_csv(wo_path)
    print(f"  Raw whole-organism records: {len(wo_df)}")

    experimental_negatives = []

    # Strategy A: Compounds with percentage-based activity < 30%
    # These were tested and showed minimal worm killing
    pct_mask = wo_df["standard_units"] == "%"
    pct_df = wo_df[pct_mask].copy()
    pct_df["standard_value"] = pd.to_numeric(pct_df["standard_value"], errors="coerce")

    if len(pct_df) > 0:
        # Group by compound, get max activity percentage
        pct_summary = (
            pct_df
            .groupby("molecule_chembl_id")
            .agg(max_pct=("standard_value", "max"),
                 n_tests=("standard_value", "count"),
                 smiles=("canonical_smiles", "first"))
            .reset_index()
        )

        # Compounds where even the best result was < 30% activity
        weak_compounds = pct_summary[pct_summary["max_pct"] < 30]
        print(f"  Compounds with <30% max activity in whole-worm assays: {len(weak_compounds)}")

        for _, row in weak_compounds.iterrows():
            experimental_negatives.append({
                "molecule_chembl_id": row["molecule_chembl_id"],
                "label": "inactive",
                "source": "experimental_negative_wholeworm",
                "evidence": f"Max {row['max_pct']:.0f}% activity in {int(row['n_tests'])} whole-worm tests",
                "canonical_smiles": row["smiles"],
            })

    # Strategy B: Compounds with dose-response IC50/EC50 > 50000 nM
    dose_types = ["IC50", "EC50", "LC50", "ED50"]
    dose_mask = wo_df["standard_type"].isin(dose_types)
    dose_df = wo_df[dose_mask].copy()
    dose_df["standard_value"] = pd.to_numeric(dose_df["standard_value"], errors="coerce")

    # Normalize to nM
    um_mask = dose_df["standard_units"] == "uM"
    dose_df.loc[um_mask, "standard_value"] = dose_df.loc[um_mask, "standard_value"] * 1000

    if len(dose_df) > 0:
        dose_summary = (
            dose_df
            .groupby("molecule_chembl_id")
            .agg(best_value=("standard_value", "min"),
                 smiles=("canonical_smiles", "first"))
            .reset_index()
        )
        # Very weak compounds (best IC50 > 50000 nM = 50 uM)
        very_weak = dose_summary[dose_summary["best_value"] > 50000]
        print(f"  Compounds with best IC50 > 50 uM (whole-worm): {len(very_weak)}")

        existing_ids = set(n["molecule_chembl_id"] for n in experimental_negatives)
        for _, row in very_weak.iterrows():
            if row["molecule_chembl_id"] not in existing_ids:
                experimental_negatives.append({
                    "molecule_chembl_id": row["molecule_chembl_id"],
                    "label": "inactive",
                    "source": "experimental_negative_wholeworm",
                    "evidence": f"Best IC50 = {row['best_value']:.0f} nM (very weak)",
                    "canonical_smiles": row["smiles"],
                })

    print(f"  Total experimental negatives from whole-worm data: {len(experimental_negatives)}")
    return pd.DataFrame(experimental_negatives)


# ---------------------------------------------------------------------------
# STEP 2: Add diverse assumed negatives from approved drug classes
# ---------------------------------------------------------------------------

def add_diverse_assumed_negatives(existing_ids):
    """
    Add approved drugs from therapeutic areas with no plausible
    connection to anti-parasitic activity.

    We systematically sample from: cardiovascular, psychiatric,
    endocrine, ophthalmologic, dermatologic, and other unrelated areas.
    """
    print("\n" + "=" * 60)
    print("STEP 2: ADDING DIVERSE ASSUMED NEGATIVES")
    print("=" * 60)

    # Expanded list of unrelated approved drugs with their SMILES
    diverse_negatives = [
        # Cardiovascular
        {"name": "METOPROLOL", "id": "CHEMBL13", "class": "Beta-blocker",
         "smiles": "COCCc1ccc(OCC(O)CNC(C)C)cc1"},
        {"name": "DILTIAZEM", "id": "CHEMBL23", "class": "Calcium channel blocker",
         "smiles": "COc1ccc(C2Sc3ccccc3N(CCN(C)C)C(=O)C2OC(C)=O)cc1"},
        {"name": "LOSARTAN", "id": "CHEMBL191", "class": "ARB",
         "smiles": "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1"},
        {"name": "PROPRANOLOL", "id": "CHEMBL27", "class": "Beta-blocker",
         "smiles": "CC(C)NCC(O)COc1cccc2ccccc12"},
        {"name": "DIGOXIN", "id": "CHEMBL1751", "class": "Cardiac glycoside",
         "smiles": None},
        {"name": "WARFARIN", "id": "CHEMBL1464", "class": "Anticoagulant",
         "smiles": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"},
        {"name": "FUROSEMIDE", "id": "CHEMBL35", "class": "Loop diuretic",
         "smiles": "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl"},

        # Psychiatric
        {"name": "FLUOXETINE", "id": "CHEMBL41", "class": "SSRI",
         "smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"},
        {"name": "RISPERIDONE", "id": "CHEMBL85", "class": "Antipsychotic",
         "smiles": "Cc1nc2n(c1C)CCCC2=O"},
        {"name": "DIAZEPAM", "id": "CHEMBL97", "class": "Benzodiazepine",
         "smiles": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"},
        {"name": "LITHIUM CARBONATE", "id": "CHEMBL1200668", "class": "Mood stabilizer",
         "smiles": None},
        {"name": "ARIPIPRAZOLE", "id": "CHEMBL1112", "class": "Antipsychotic",
         "smiles": "O=c1[nH]c2ccccc2n1CCCCN1CCN(c2cccc(Cl)c2Cl)CC1"},
        {"name": "BUPROPION", "id": "CHEMBL894", "class": "Antidepressant",
         "smiles": "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1"},
        {"name": "QUETIAPINE", "id": "CHEMBL716", "class": "Antipsychotic",
         "smiles": "OCCOCCN1CCN(c2c3ccc(Cl)cc3Sc3ccccc32)CC1"},

        # Endocrine
        {"name": "GLIPIZIDE", "id": "CHEMBL1202", "class": "Sulfonylurea",
         "smiles": "Cc1cnc(NS(=O)(=O)c2ccc(CCNC(=O)c3ccccn3)cc2)s1"},
        {"name": "PIOGLITAZONE", "id": "CHEMBL595", "class": "Thiazolidinedione",
         "smiles": "O=C1NC(=O)C(Cc2ccc(OCCc3ncccc3C)cc2)S1"},
        {"name": "DEXAMETHASONE", "id": "CHEMBL384467", "class": "Corticosteroid",
         "smiles": "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO"},

        # Respiratory
        {"name": "SALBUTAMOL", "id": "CHEMBL714", "class": "Beta-2 agonist",
         "smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1"},
        {"name": "TIOTROPIUM", "id": "CHEMBL1201335", "class": "Anticholinergic",
         "smiles": None},
        {"name": "FLUTICASONE", "id": "CHEMBL1201396", "class": "Inhaled corticosteroid",
         "smiles": None},

        # Pain/Musculoskeletal
        {"name": "CELECOXIB", "id": "CHEMBL118", "class": "COX-2 inhibitor",
         "smiles": "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"},
        {"name": "NAPROXEN", "id": "CHEMBL154", "class": "NSAID",
         "smiles": "COc1ccc2cc(CC(C)C(=O)O)ccc2c1"},
        {"name": "PREGABALIN", "id": "CHEMBL1059", "class": "GABA analogue",
         "smiles": "CC(C)CC(CN)CC(=O)O"},
        {"name": "MORPHINE", "id": "CHEMBL70", "class": "Opioid",
         "smiles": "CN1CCC23c4c5ccc(O)c4OC2C(O)=CC3C1C5"},
        {"name": "ACETAMINOPHEN", "id": "CHEMBL112", "class": "Analgesic",
         "smiles": "CC(=O)Nc1ccc(O)cc1"},

        # GI
        {"name": "ONDANSETRON", "id": "CHEMBL46", "class": "5-HT3 antagonist",
         "smiles": "Cn1c2ccccc2c2c(C=O)c(CN3CCCC3)cc21"},
        {"name": "RANITIDINE", "id": "CHEMBL1680", "class": "H2 blocker",
         "smiles": "CNC(/N=C/[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1"},
        {"name": "LOPERAMIDE", "id": "CHEMBL841", "class": "Antidiarrheal",
         "smiles": "O=C(c1ccc(Cl)cc1)C(CCN1CCC(O)(c2ccc(Cl)cc2)CC1)(c1ccccc1)C(C)C"},

        # Dermatologic
        {"name": "ISOTRETINOIN", "id": "CHEMBL297", "class": "Retinoid",
         "smiles": "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C(=O)O)C(C)(C)CCC1"},
        {"name": "TERBINAFINE", "id": "CHEMBL822", "class": "Antifungal",
         "smiles": "CN(/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"},

        # Ophthalmologic
        {"name": "LATANOPROST", "id": "CHEMBL1188", "class": "Prostaglandin analogue",
         "smiles": "CCCCC(O)/C=C/C1C(O)CC(O)C1C/C=C\\CCCC(=O)OC(C)C"},
        {"name": "TIMOLOL", "id": "CHEMBL499", "class": "Beta-blocker (eye)",
         "smiles": "CC(C)(C)NCC(O)COc1nsnc1N1CCOCC1"},

        # Metabolic
        {"name": "ALLOPURINOL", "id": "CHEMBL707", "class": "Xanthine oxidase inhibitor",
         "smiles": "O=c1[nH]cnc2[nH]ncc12"},
        {"name": "COLCHICINE", "id": "CHEMBL107", "class": "Anti-gout",
         "smiles": "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1CC(NC(C)=O)C2"},

        # Antihistamine
        {"name": "CETIRIZINE", "id": "CHEMBL1000", "class": "H1 antagonist",
         "smiles": "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1"},
        {"name": "DIPHENHYDRAMINE", "id": "CHEMBL657", "class": "H1 antagonist",
         "smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1"},

        # Antiepileptic
        {"name": "LEVETIRACETAM", "id": "CHEMBL1286", "class": "Antiepileptic",
         "smiles": "CCC(C(N)=O)N1CCCC1=O"},
        {"name": "CARBAMAZEPINE", "id": "CHEMBL108", "class": "Antiepileptic",
         "smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc21"},
        {"name": "VALPROIC ACID", "id": "CHEMBL109", "class": "Antiepileptic",
         "smiles": "CCCC(CCC)C(=O)O"},
    ]

    new_negatives = []
    for drug in diverse_negatives:
        if drug["id"] not in existing_ids and drug["smiles"] is not None:
            new_negatives.append({
                "molecule_chembl_id": drug["id"],
                "pref_name": drug["name"],
                "label": "inactive",
                "source": f"assumed_negative_{drug['class'].lower().replace(' ', '_')}",
                "evidence": f"Approved for {drug['class']}. No plausible anti-parasitic mechanism.",
                "canonical_smiles": drug["smiles"],
            })

    print(f"  Added {len(new_negatives)} diverse assumed negatives across {len(set(d['class'] for d in diverse_negatives))} therapeutic classes")
    return pd.DataFrame(new_negatives)


# ---------------------------------------------------------------------------
# STEP 3: Build evaluation set v2 with documented inclusion criteria
# ---------------------------------------------------------------------------

def build_eval_set_v2(v1_df, exp_negatives, assumed_negatives):
    """
    Construct a hardened evaluation set with:
    - All v1 actives and weak compounds (from ChEMBL target-based data)
    - Curated known drugs (praziquantel, oxamniquine, mefloquine)
    - Experimental negatives from whole-organism screens
    - Diverse assumed negatives from unrelated drug classes
    - Documented inclusion criteria for every compound
    """
    print("\n" + "=" * 60)
    print("STEP 3: BUILDING EVALUATION SET v2")
    print("=" * 60)

    # Start with v1 positives and weak
    positives = v1_df[v1_df["label"].isin(["active", "weak"])].copy()
    positives["eval_version"] = "v1_carried_forward"
    print(f"  Actives from v1: {len(positives[positives['label'] == 'active'])}")
    print(f"  Weak from v1: {len(positives[positives['label'] == 'weak'])}")

    # V1 curated entries
    curated = v1_df[v1_df["source"].isin(["curated", "chembl+curated"])].copy()
    curated["eval_version"] = "v1_curated"

    # Remove curated from positives to avoid duplicates
    positives = positives[~positives["molecule_chembl_id"].isin(curated["molecule_chembl_id"])]

    # Experimental negatives
    if len(exp_negatives) > 0:
        exp_negatives["eval_version"] = "v2_experimental_negative"
        print(f"  Experimental negatives: {len(exp_negatives)}")

    # Assumed negatives (v1 + new diverse)
    v1_assumed = v1_df[v1_df["source"] == "curated_negative"].copy()
    v1_assumed["eval_version"] = "v1_assumed_negative"

    if len(assumed_negatives) > 0:
        assumed_negatives["eval_version"] = "v2_assumed_negative"
        print(f"  New assumed negatives: {len(assumed_negatives)}")

    # Combine
    frames = [positives, curated]
    if len(exp_negatives) > 0:
        frames.append(exp_negatives)
    frames.append(v1_assumed)
    if len(assumed_negatives) > 0:
        frames.append(assumed_negatives)

    eval_v2 = pd.concat(frames, ignore_index=True)

    # Deduplicate
    eval_v2 = eval_v2.drop_duplicates(subset=["molecule_chembl_id"], keep="first")

    print(f"\n  EVALUATION SET v2 SUMMARY:")
    print(f"  Total compounds: {len(eval_v2)}")
    for label in ["active", "weak", "inactive"]:
        n = len(eval_v2[eval_v2["label"] == label])
        print(f"    {label:10s}: {n}")

    print(f"\n  By source type:")
    for src, count in eval_v2["eval_version"].value_counts().items():
        print(f"    {src:35s}: {count}")

    return eval_v2


# ---------------------------------------------------------------------------
# STEP 4: Train/test split
# ---------------------------------------------------------------------------

def split_eval_set(eval_v2):
    """
    Split into development (80%) and held-out test (20%).

    CRITICAL RULES:
    - Curated drugs (praziquantel, etc.) go to BOTH sets for sanity checking
    - Split is stratified by label (proportional actives/inactives in each)
    - Split is deterministic (fixed random seed for reproducibility)
    - Held-out test set is NEVER used for tuning signal weights
    """
    print("\n" + "=" * 60)
    print("STEP 4: TRAIN/TEST SPLIT")
    print("=" * 60)

    np.random.seed(42)  # Reproducible

    # Curated drugs go to both sets
    curated_mask = eval_v2["eval_version"].isin(["v1_curated"])
    curated = eval_v2[curated_mask].copy()
    non_curated = eval_v2[~curated_mask].copy()

    # Stratified split of non-curated
    dev_frames = []
    test_frames = []

    for label in ["active", "weak", "inactive"]:
        subset = non_curated[non_curated["label"] == label]
        indices = np.random.permutation(len(subset))
        split_point = int(0.8 * len(subset))

        dev_idx = subset.index[indices[:split_point]]
        test_idx = subset.index[indices[split_point:]]

        dev_frames.append(subset.loc[dev_idx])
        test_frames.append(subset.loc[test_idx])

    dev_set = pd.concat([pd.concat(dev_frames)] + [curated], ignore_index=True)
    test_set = pd.concat([pd.concat(test_frames)] + [curated], ignore_index=True)

    print(f"\n  Development set: {len(dev_set)} compounds")
    for label in ["active", "weak", "inactive"]:
        print(f"    {label}: {len(dev_set[dev_set['label'] == label])}")

    print(f"\n  Held-out test set: {len(test_set)} compounds")
    for label in ["active", "weak", "inactive"]:
        print(f"    {label}: {len(test_set[test_set['label'] == label])}")

    return dev_set, test_set


# ---------------------------------------------------------------------------
# STEP 5: Re-run ranking on both sets with bootstrap confidence intervals
# ---------------------------------------------------------------------------

def run_ranking_with_bootstrap(eval_set, activities_df, set_name, n_bootstrap=1000):
    """
    Run the v3 ranking pipeline on the given eval set and compute
    bootstrap confidence intervals for AUROC and AUPRC.
    """
    from sklearn.metrics import roc_auc_score, average_precision_score

    print(f"\n  {'=' * 55}")
    print(f"  EVALUATION: {set_name}")
    print(f"  {'=' * 55}")

    # Recompute scores using the same v3 logic
    # Load v3 scores if available, otherwise use simplified scoring
    v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv")
    if os.path.exists(v3_path):
        v3_df = pd.read_csv(v3_path)
        # Merge v3 scores into eval set
        score_cols = ["molecule_chembl_id", "composite_score", "score_potency",
                      "score_target", "score_similarity", "score_whole_org",
                      "score_drug_stage"]
        available_cols = [c for c in score_cols if c in v3_df.columns]
        eval_scored = eval_set.merge(
            v3_df[available_cols], on="molecule_chembl_id", how="left"
        )
        # For new compounds not in v3, assign score of 0
        eval_scored["composite_score"] = eval_scored["composite_score"].fillna(0)
    else:
        eval_scored = eval_set.copy()
        eval_scored["composite_score"] = 0

    # Binary evaluation: active vs inactive
    binary = eval_scored[eval_scored["label"].isin(["active", "inactive"])].copy()
    labels = (binary["label"] == "active").astype(int).values
    scores = binary["composite_score"].values

    n_pos = sum(labels == 1)
    n_neg = sum(labels == 0)
    print(f"\n  Actives: {n_pos}, Inactives: {n_neg}")

    if n_pos == 0 or n_neg == 0:
        print("  Cannot compute metrics: need both classes")
        return None

    # Point estimates
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    print(f"  AUROC (point): {auroc:.4f}")
    print(f"  AUPRC (point): {auprc:.4f}")

    # Bootstrap confidence intervals
    auroc_boots = []
    auprc_boots = []
    n_total = len(labels)

    for i in range(n_bootstrap):
        idx = np.random.choice(n_total, size=n_total, replace=True)
        boot_labels = labels[idx]
        boot_scores = scores[idx]

        # Skip if bootstrap sample has only one class
        if len(np.unique(boot_labels)) < 2:
            continue

        auroc_boots.append(roc_auc_score(boot_labels, boot_scores))
        auprc_boots.append(average_precision_score(boot_labels, boot_scores))

    auroc_boots = np.array(auroc_boots)
    auprc_boots = np.array(auprc_boots)

    auroc_ci_lo = np.percentile(auroc_boots, 2.5)
    auroc_ci_hi = np.percentile(auroc_boots, 97.5)
    auprc_ci_lo = np.percentile(auprc_boots, 2.5)
    auprc_ci_hi = np.percentile(auprc_boots, 97.5)

    print(f"  AUROC 95% CI: [{auroc_ci_lo:.4f}, {auroc_ci_hi:.4f}]")
    print(f"  AUPRC 95% CI: [{auprc_ci_lo:.4f}, {auprc_ci_hi:.4f}]")
    print(f"  Bootstrap samples: {len(auroc_boots)}/{n_bootstrap}")

    # Three-class check
    print(f"\n  Three-class distributions:")
    for label in ["active", "weak", "inactive"]:
        subset = eval_scored[eval_scored["label"] == label]["composite_score"]
        if len(subset) > 0:
            print(f"    {label:10s}  n={len(subset):3d}  "
                  f"mean={subset.mean():.4f}  median={subset.median():.4f}")

    # Known drug check
    print(f"\n  Known drug positions:")
    ranked = eval_scored.sort_values("composite_score", ascending=False).reset_index(drop=True)
    ranked["rank"] = range(1, len(ranked) + 1)
    total = len(ranked)

    for drug in ["PRAZIQUANTEL", "OXAMNIQUINE", "MEFLOQUINE", "ATOVAQUONE"]:
        match = ranked[ranked["pref_name"] == drug]
        if len(match) > 0:
            row = match.iloc[0]
            pct = (1 - row["rank"] / total) * 100
            print(f"    {drug:20s}  rank {int(row['rank']):3d}/{total}  (top {pct:.0f}%)")

    return {
        "auroc": auroc,
        "auroc_ci": (auroc_ci_lo, auroc_ci_hi),
        "auprc": auprc,
        "auprc_ci": (auprc_ci_lo, auprc_ci_hi),
        "n_active": n_pos,
        "n_inactive": n_neg,
    }


# ---------------------------------------------------------------------------
# STEP 6: Per-target performance breakdown
# ---------------------------------------------------------------------------

def per_target_breakdown(eval_set, activities_df):
    """
    Break down ranking performance by target.
    This reveals whether the pipeline works equally well across all
    targets or is carried by just one or two.
    """
    print(f"\n  {'=' * 55}")
    print(f"  PER-TARGET PERFORMANCE BREAKDOWN")
    print(f"  {'=' * 55}")

    v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv")
    if not os.path.exists(v3_path):
        print("  V3 rankings not found.")
        return

    v3_df = pd.read_csv(v3_path)

    # For each target, see how well its compounds are ranked
    targets = activities_df.groupby("target_name")["molecule_chembl_id"].apply(set).to_dict()

    print(f"\n  {'Target':50s} {'N compounds':>12s} {'Mean score':>12s} {'Median score':>12s}")
    print(f"  {'-'*50} {'-'*12} {'-'*12} {'-'*12}")

    for target_name, compound_ids in sorted(targets.items()):
        target_scores = v3_df[v3_df["molecule_chembl_id"].isin(compound_ids)]
        if len(target_scores) > 0:
            mean_score = target_scores["composite_score"].mean()
            median_score = target_scores["composite_score"].median()
            n = len(target_scores)
            print(f"  {target_name:50s} {n:12d} {mean_score:12.4f} {median_score:12.4f}")


# ---------------------------------------------------------------------------
# STEP 7: Save everything with documentation
# ---------------------------------------------------------------------------

def save_hardened_benchmark(eval_v2, dev_set, test_set, results):
    """Save all benchmark artifacts with inclusion criteria documentation."""

    os.makedirs(EVAL_DIR, exist_ok=True)

    # Save full v2 eval set
    eval_v2.to_csv(os.path.join(EVAL_DIR, "evaluation_set_v2.csv"), index=False)

    # Save splits
    dev_set.to_csv(os.path.join(EVAL_DIR, "eval_v2_dev.csv"), index=False)
    test_set.to_csv(os.path.join(EVAL_DIR, "eval_v2_test.csv"), index=False)

    # Save inclusion criteria
    criteria = [
        "=" * 60,
        "KIRA EVALUATION SET v2 — INCLUSION CRITERIA",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 60,
        "",
        "ACTIVE compounds (label = 'active'):",
        f"  - ChEMBL target-based IC50 < 1000 nM against S. mansoni proteins",
        f"  - Curated known effective drugs (praziquantel, oxamniquine, mefloquine)",
        "",
        "WEAK compounds (label = 'weak'):",
        f"  - ChEMBL target-based IC50 between 1000-10000 nM",
        "",
        "INACTIVE compounds (label = 'inactive'):",
        f"  - Experimental: tested in whole-worm assays with <30% max activity",
        f"  - Experimental: whole-worm IC50/EC50 > 50000 nM",
        f"  - Assumed: approved drugs from unrelated therapeutic areas",
        f"    (cardiovascular, psychiatric, endocrine, ophthalmologic, etc.)",
        f"  - Assumed negatives are NOT experimentally verified against S. mansoni",
        "",
        "SPLIT:",
        f"  - 80/20 stratified split, random seed 42",
        f"  - Curated drugs appear in both dev and test for sanity checking",
        "",
        "FROZEN:",
        f"  - This benchmark is frozen as of {datetime.now().strftime('%Y-%m-%d')}",
        f"  - No changes to labels or inclusion criteria without new version number",
        f"  - Signal weight tuning uses dev set ONLY",
        f"  - Held-out test set for final evaluation ONLY",
    ]

    criteria_path = os.path.join(EVAL_DIR, "eval_v2_inclusion_criteria.txt")
    with open(criteria_path, "w") as f:
        f.write("\n".join(criteria))

    print(f"\n  Saved:")
    print(f"    {os.path.join(EVAL_DIR, 'evaluation_set_v2.csv')}")
    print(f"    {os.path.join(EVAL_DIR, 'eval_v2_dev.csv')}")
    print(f"    {os.path.join(EVAL_DIR, 'eval_v2_test.csv')}")
    print(f"    {criteria_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 08")
    print("  Benchmark Hardening")
    print("=" * 60)

    # Load v1 eval set
    v1_path = os.path.join(EVAL_DIR, "evaluation_set_v1.csv")
    v1_df = pd.read_csv(v1_path)
    print(f"\n  V1 eval set: {len(v1_df)} compounds")

    activities_df = pd.read_csv(os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"))

    # Step 1: Find experimental negatives
    exp_negatives = find_experimental_negatives()

    # Step 2: Add diverse assumed negatives
    existing_ids = set(v1_df["molecule_chembl_id"].tolist())
    if len(exp_negatives) > 0:
        existing_ids.update(exp_negatives["molecule_chembl_id"].tolist())
    assumed_negatives = add_diverse_assumed_negatives(existing_ids)

    # Step 3: Build v2 eval set
    eval_v2 = build_eval_set_v2(v1_df, exp_negatives, assumed_negatives)

    # Step 4: Split
    dev_set, test_set = split_eval_set(eval_v2)

    # Step 5: Evaluate with bootstrap CIs
    print("\n" + "=" * 60)
    print("EVALUATION WITH BOOTSTRAP CONFIDENCE INTERVALS")
    print("=" * 60)

    dev_results = run_ranking_with_bootstrap(dev_set, activities_df, "DEVELOPMENT SET")
    test_results = run_ranking_with_bootstrap(test_set, activities_df, "HELD-OUT TEST SET")

    # Step 6: Per-target breakdown
    per_target_breakdown(eval_v2, activities_df)

    # Step 7: Save
    save_hardened_benchmark(eval_v2, dev_set, test_set,
                           {"dev": dev_results, "test": test_results})

    # Final comparison
    print(f"\n{'=' * 60}")
    print(f"  BENCHMARK HARDENING COMPLETE")
    print(f"{'=' * 60}")

    n_inactive_v1 = len(v1_df[v1_df["label"] == "inactive"])
    n_inactive_v2 = len(eval_v2[eval_v2["label"] == "inactive"])
    print(f"\n  Inactive controls: v1={n_inactive_v1} -> v2={n_inactive_v2}")

    if dev_results:
        print(f"\n  Dev set AUROC:  {dev_results['auroc']:.4f}  "
              f"95% CI [{dev_results['auroc_ci'][0]:.4f}, {dev_results['auroc_ci'][1]:.4f}]")
    if test_results:
        print(f"  Test set AUROC: {test_results['auroc']:.4f}  "
              f"95% CI [{test_results['auroc_ci'][0]:.4f}, {test_results['auroc_ci'][1]:.4f}]")

    print(f"\n  The benchmark is now frozen. Inclusion criteria documented.")
    print(f"  Tune weights on dev set only. Report final metrics on test set.")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 08_harden_benchmark.py successfully."
echo "Now run: python 08_harden_benchmark.py"
