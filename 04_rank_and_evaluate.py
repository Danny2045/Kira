"""
Kira - Script 04: Composite Ranking Algorithm + Evaluation
============================================================

This is where Kira becomes an intelligence system.

It takes every compound in our dataset, computes a composite repurposing
score from multiple signals, and then tests whether that score actually
separates known actives from known inactives.

THE FIVE SIGNALS:
  1. Binding potency (pIC50) - how strongly it hits the target
  2. Target essentiality - how important that target is to the worm
  3. Data confidence - how many independent measurements support the claim
  4. Drug stage - approved > clinical > preclinical
  5. Multi-target bonus - hitting multiple parasite targets

EVALUATION:
  We test the ranking against our ground-truth evaluation set using:
  - AUROC (area under ROC curve) - does the algorithm rank actives above inactives?
  - AUPRC (area under precision-recall curve) - more sensitive to false positives
  - Rank inspection - where do known drugs actually land?

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 04_rank_and_evaluate.py
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

# Import shared modules (canonical source for scoring functions and data)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from kira.scoring import (  # noqa: E402
    WEIGHTS_V1,
    DRUG_STAGE_SCORES,
    DEFAULT_DRUG_STAGE,
    compute_potency_score,
    compute_target_score,
    compute_confidence_score,
    compute_drug_stage_score,
)
from kira.targets import TARGET_ESSENTIALITY, DEFAULT_ESSENTIALITY  # noqa: E402

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")

# --- Signal weights (from shared module) ---
WEIGHT_POTENCY = WEIGHTS_V1["potency"]
WEIGHT_TARGET = WEIGHTS_V1["target"]
WEIGHT_CONFIDENCE = WEIGHTS_V1["confidence"]
WEIGHT_DRUG_STAGE = WEIGHTS_V1["drug_stage"]
WEIGHT_MULTITARGET = WEIGHTS_V1["multitarget"]


# ---------------------------------------------------------------------------
# SIGNAL 1: Potency Score
# ---------------------------------------------------------------------------

# Note: compute_potency_score, compute_target_score,
# compute_confidence_score, compute_drug_stage_score are now
# imported from kira.scoring (shared module).


# ---------------------------------------------------------------------------
# SIGNAL 5: Multi-Target Score
# ---------------------------------------------------------------------------

def compute_multitarget_scores(activities_df):
    """
    Count how many distinct parasite targets each compound hits
    (with IC50 < 10000 nM). Compounds hitting multiple targets get a bonus.

    Returns a dict: molecule_chembl_id -> score (0 to 1)
    """
    target_counts = (
        activities_df
        .groupby("molecule_chembl_id")["target_name"]
        .nunique()
    )

    scores = {}
    for compound_id, n_targets in target_counts.items():
        if n_targets >= 3:
            scores[compound_id] = 1.0
        elif n_targets == 2:
            scores[compound_id] = 0.7
        else:
            scores[compound_id] = 0.3

    return scores


# ---------------------------------------------------------------------------
# COMPOSITE SCORE
# ---------------------------------------------------------------------------

def compute_composite_scores(eval_df, activities_df):
    """
    Compute the composite repurposing score for every compound.

    composite = w1 * potency + w2 * target + w3 * confidence
              + w4 * drug_stage + w5 * multitarget

    All signals are 0-1, weights sum to 1, so composite is 0-1.
    """

    print("\n" + "=" * 60)
    print("COMPUTING COMPOSITE SCORES")
    print("=" * 60)

    # Pre-compute multi-target scores
    multitarget_scores = compute_multitarget_scores(activities_df)

    results = []

    for _, row in eval_df.iterrows():
        compound_id = row["molecule_chembl_id"]

        # Signal 1: Potency
        potency = compute_potency_score(row.get("best_value"))

        # Signal 2: Target essentiality
        target = compute_target_score(row.get("best_target", ""))

        # Signal 3: Confidence
        n_meas = row.get("n_measurements", 0)
        if pd.isna(n_meas):
            n_meas = 0
        confidence = compute_confidence_score(int(n_meas))

        # Signal 4: Drug stage
        drug_stage = compute_drug_stage_score(row.get("max_phase"))

        # Signal 5: Multi-target
        multitarget = multitarget_scores.get(compound_id, 0.1)

        # Composite
        composite = (
            WEIGHT_POTENCY * potency
            + WEIGHT_TARGET * target
            + WEIGHT_CONFIDENCE * confidence
            + WEIGHT_DRUG_STAGE * drug_stage
            + WEIGHT_MULTITARGET * multitarget
        )

        results.append({
            "molecule_chembl_id": compound_id,
            "pref_name": row.get("pref_name"),
            "label": row.get("label"),
            "source": row.get("source"),
            "best_value": row.get("best_value"),
            "best_target": row.get("best_target"),
            "score_potency": round(potency, 3),
            "score_target": round(target, 3),
            "score_confidence": round(confidence, 3),
            "score_drug_stage": round(drug_stage, 3),
            "score_multitarget": round(multitarget, 3),
            "composite_score": round(composite, 4),
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("composite_score", ascending=False)

    # Print top 20
    print(f"\n  TOP 20 COMPOUNDS BY COMPOSITE SCORE:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'Potency':>7s}  {'Target':>6s}  {'Conf':>5s}  {'Stage':>5s}  {'Multi':>5s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*5}  {'-'*5}  {'-'*5}  {'-'*8}  {'-'*30}")

    for rank, (_, row) in enumerate(results_df.head(20).iterrows(), 1):
        name = row["pref_name"] if pd.notna(row["pref_name"]) else row["molecule_chembl_id"]
        if len(str(name)) > 30:
            name = str(name)[:27] + "..."
        print(f"  {rank:4d}  {row['composite_score']:6.4f}  {row['score_potency']:7.3f}  {row['score_target']:6.3f}  "
              f"{row['score_confidence']:5.3f}  {row['score_drug_stage']:5.3f}  {row['score_multitarget']:5.3f}  "
              f"{row['label']:>8s}  {name}")

    # Print bottom 10
    print(f"\n  BOTTOM 10 COMPOUNDS:")
    print(f"  {'Rank':>4s}  {'Score':>6s}  {'Label':>8s}  Name")
    print(f"  {'-'*4}  {'-'*6}  {'-'*8}  {'-'*30}")
    total = len(results_df)
    for rank, (_, row) in enumerate(results_df.tail(10).iterrows()):
        name = row["pref_name"] if pd.notna(row["pref_name"]) else row["molecule_chembl_id"]
        actual_rank = total - 9 + rank
        print(f"  {actual_rank:4d}  {row['composite_score']:6.4f}  {row['label']:>8s}  {name}")

    return results_df


# ---------------------------------------------------------------------------
# EVALUATION
# ---------------------------------------------------------------------------

def evaluate_ranking(results_df):
    """
    Test whether the composite score actually separates actives from inactives.

    Metrics:
    - AUROC: Does the ranking put actives above inactives?
      0.5 = random, 1.0 = perfect separation
    - AUPRC: Of the things scored highly, how many are truly active?
      More sensitive to false positives than AUROC.
    - Score distributions: What do actives vs inactives look like?
    """
    from sklearn.metrics import roc_auc_score, average_precision_score

    print("\n" + "=" * 60)
    print("EVALUATION AGAINST GROUND TRUTH")
    print("=" * 60)

    scores = results_df["composite_score"].values

    # --- Test 1: Active vs Inactive (binary) ---
    print("\n  TEST 1: Active vs Inactive separation")
    print("  " + "-" * 50)

    binary_df = results_df[results_df["label"].isin(["active", "inactive"])].copy()
    binary_labels = (binary_df["label"] == "active").astype(int).values
    binary_scores = binary_df["composite_score"].values

    n_active = sum(binary_labels == 1)
    n_inactive = sum(binary_labels == 0)
    print(f"  Actives: {n_active}, Inactives: {n_inactive}")

    if n_active > 0 and n_inactive > 0:
        auroc = roc_auc_score(binary_labels, binary_scores)
        auprc = average_precision_score(binary_labels, binary_scores)
        print(f"  AUROC:  {auroc:.4f}  (0.5 = random, 1.0 = perfect)")
        print(f"  AUPRC:  {auprc:.4f}  (baseline = {n_active/(n_active+n_inactive):.4f})")

        if auroc > 0.9:
            print(f"  --> STRONG separation. Algorithm clearly distinguishes actives from inactives.")
        elif auroc > 0.75:
            print(f"  --> GOOD separation. Algorithm has learned meaningful signal.")
        elif auroc > 0.6:
            print(f"  --> WEAK separation. Some signal but significant overlap.")
        else:
            print(f"  --> POOR separation. Algorithm is not much better than random.")
    else:
        print(f"  Cannot compute AUROC: need both actives and inactives.")

    # --- Test 2: Three-class ordering (active > weak > inactive) ---
    print(f"\n  TEST 2: Three-class score distributions")
    print("  " + "-" * 50)

    for label in ["active", "weak", "inactive"]:
        subset = results_df[results_df["label"] == label]["composite_score"]
        if len(subset) > 0:
            print(f"  {label:10s}  n={len(subset):3d}  "
                  f"mean={subset.mean():.4f}  "
                  f"median={subset.median():.4f}  "
                  f"min={subset.min():.4f}  "
                  f"max={subset.max():.4f}")

    # Check if ordering is correct: mean(active) > mean(weak) > mean(inactive)
    active_mean = results_df[results_df["label"] == "active"]["composite_score"].mean()
    weak_mean = results_df[results_df["label"] == "weak"]["composite_score"].mean()
    inactive_mean = results_df[results_df["label"] == "inactive"]["composite_score"].mean()

    if active_mean > weak_mean > inactive_mean:
        print(f"\n  --> CORRECT ordering: active ({active_mean:.4f}) > weak ({weak_mean:.4f}) > inactive ({inactive_mean:.4f})")
    elif active_mean > inactive_mean:
        print(f"\n  --> PARTIAL ordering: active > inactive but weak is misplaced")
    else:
        print(f"\n  --> INCORRECT ordering: algorithm needs revision")

    # --- Test 3: Where do known drugs rank? ---
    print(f"\n  TEST 3: Known drug rankings")
    print("  " + "-" * 50)

    known_drugs = ["PRAZIQUANTEL", "OXAMNIQUINE", "MEFLOQUINE", "ATOVAQUONE",
                   "PANOBINOSTAT", "VORINOSTAT", "IDEBENONE"]

    ranked = results_df.reset_index(drop=True)
    ranked["rank"] = range(1, len(ranked) + 1)
    total = len(ranked)

    for drug_name in known_drugs:
        match = ranked[ranked["pref_name"] == drug_name]
        if len(match) > 0:
            row = match.iloc[0]
            percentile = (1 - row["rank"] / total) * 100
            print(f"  {drug_name:20s}  rank {int(row['rank']):3d}/{total}  "
                  f"score={row['composite_score']:.4f}  "
                  f"(top {percentile:.0f}%)")
        else:
            print(f"  {drug_name:20s}  not found in dataset")

    # --- Test 4: False positive check ---
    print(f"\n  TEST 4: Negative control analysis")
    print("  " + "-" * 50)

    negatives = results_df[results_df["source"] == "curated_negative"]
    if len(negatives) > 0:
        active_median = results_df[results_df["label"] == "active"]["composite_score"].median()
        neg_above_active_median = negatives[negatives["composite_score"] > active_median]
        print(f"  Negative controls: {len(negatives)}")
        print(f"  Active median score: {active_median:.4f}")
        print(f"  Negatives scoring above active median: {len(neg_above_active_median)}")

        if len(neg_above_active_median) == 0:
            print(f"  --> GOOD: No negative controls are scored in the active range.")
        else:
            print(f"  --> WARNING: Some negative controls scored suspiciously high:")
            for _, row in neg_above_active_median.iterrows():
                print(f"      {row['pref_name']}: {row['composite_score']:.4f}")

    return auroc if (n_active > 0 and n_inactive > 0) else None


# ---------------------------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------------------------

def save_results(results_df, auroc):
    """Save ranked results and evaluation report."""

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save ranked compounds
    ranked_path = os.path.join(OUTPUT_DIR, "kira_ranked_candidates_v1.csv")
    results_df.to_csv(ranked_path, index=False)
    print(f"\n  Ranked candidates saved to: {ranked_path}")

    # Save evaluation report
    report_lines = [
        "=" * 60,
        "KIRA RANKING EVALUATION REPORT v1",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 60,
        "",
        "ALGORITHM: Weighted composite of 5 signals",
        f"  Potency weight:     {WEIGHT_POTENCY}",
        f"  Target weight:      {WEIGHT_TARGET}",
        f"  Confidence weight:  {WEIGHT_CONFIDENCE}",
        f"  Drug stage weight:  {WEIGHT_DRUG_STAGE}",
        f"  Multi-target weight:{WEIGHT_MULTITARGET}",
        "",
        f"AUROC: {auroc:.4f}" if auroc else "AUROC: could not compute",
        "",
        "INTERPRETATION:",
    ]

    if auroc and auroc > 0.9:
        report_lines.append("  Strong separation between actives and inactives.")
    elif auroc and auroc > 0.75:
        report_lines.append("  Good separation. Algorithm captures meaningful signal.")
    elif auroc and auroc > 0.6:
        report_lines.append("  Weak separation. Improvement needed.")
    else:
        report_lines.append("  Poor separation. Algorithm needs fundamental revision.")

    report_lines.extend([
        "",
        "LIMITATIONS:",
        "  - Only 16 inactive compounds (specificity estimates are noisy)",
        "  - Target essentiality weights are human estimates, not data-derived",
        "  - Signal weights are hand-tuned, not optimized",
        "  - No molecular structure features (future: docking scores, fingerprints)",
        "  - No whole-organism assay data incorporated",
        "",
        "NEXT STEPS:",
        "  - Expand negative control set for more reliable evaluation",
        "  - Add molecular docking scores as a new signal",
        "  - Learn signal weights from data instead of hand-tuning",
        "  - Incorporate ADMET predictions",
        "  - Add supply chain availability data",
    ])

    report_path = os.path.join(OUTPUT_DIR, "kira_evaluation_report_v1.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    print(f"  Evaluation report saved to: {report_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 04")
    print("  Composite Ranking Algorithm + Evaluation")
    print("=" * 60)

    # Load evaluation set
    eval_path = os.path.join(EVAL_DIR, "evaluation_set_v1.csv")
    if not os.path.exists(eval_path):
        print(f"ERROR: {eval_path} not found. Run Script 03 first.")
        sys.exit(1)

    eval_df = pd.read_csv(eval_path)
    eval_df["max_phase"] = pd.to_numeric(eval_df.get("max_phase", pd.Series()), errors="coerce")
    print(f"\n  Loaded evaluation set: {len(eval_df)} compounds")

    # Load activity data for multi-target computation
    activities_path = os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv")
    if os.path.exists(activities_path):
        activities_df = pd.read_csv(activities_path)
    else:
        activities_df = pd.DataFrame()

    # Compute composite scores
    results_df = compute_composite_scores(eval_df, activities_df)

    # Check if sklearn is available for evaluation
    try:
        from sklearn.metrics import roc_auc_score
        has_sklearn = True
    except ImportError:
        has_sklearn = False
        print("\n  WARNING: scikit-learn not installed. Installing now...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install",
                               "scikit-learn", "--quiet"])
        has_sklearn = True

    # Evaluate
    auroc = evaluate_ranking(results_df)

    # Save
    save_results(results_df, auroc)

    print(f"\n{'=' * 60}")
    print(f"  DONE")
    print(f"")
    print(f"  Kira v1 ranking complete.")
    print(f"  AUROC: {auroc:.4f}" if auroc else "  AUROC: could not compute")
    print(f"")
    print(f"  This is your first end-to-end pipeline:")
    print(f"    PrimeKG -> ChEMBL -> Evaluation Set -> Ranking -> Evaluation")
    print(f"")
    print(f"  Files in data/processed/:")
    print(f"    kira_ranked_candidates_v1.csv    -- all compounds ranked")
    print(f"    kira_evaluation_report_v1.txt    -- evaluation summary")
    print(f"{'=' * 60}")
