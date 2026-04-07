"""Selectivity Prediction v3: Binding-Site Features vs ESM-2 Baseline.

THE EXPERIMENT:
    Kira Script 19 showed that ESM-2 global embeddings fail at cross-disease
    selectivity prediction (LODO AUROC: 0.38, 0.36, 0.55). The features
    that work in 5-fold CV (AUROC 0.895) don't transfer because they learn
    compound-specific patterns, not the biological mechanism of selectivity.

    This experiment replaces ESM-2 global similarity with BINDING-SITE-LEVEL
    divergence features — physicochemical property differences at the pocket
    positions that actually contact the drug. The hypothesis: these features
    capture WHY selectivity exists (pocket divergence), not just THAT it exists,
    and therefore should transfer across diseases.

HOW TO RUN:
    cd ~/Kira
    conda activate bio-builder
    python -m kira.experiments.run_selectivity_v3

WHAT IT OUTPUTS:
    - Per-target-pair pocket feature summary
    - 5-fold CV AUROC for binding-site features
    - Leave-one-disease-out AUROC for binding-site features
    - Comparison against Script 19 ESM-2 baseline
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from kira.experiments import TARGET_PAIRS, TargetPair
from kira.experiments.selectivity_features import compute_pocket_features, PocketFeatures


# ---------------------------------------------------------------------------
# Script 19 baseline numbers (from Kira's published results)
# ---------------------------------------------------------------------------
BASELINE_CV_AUROC = 0.895
BASELINE_LODO = {
    "Schistosomiasis": 0.382,
    "Trypanosomiasis": 0.359,
    "Leishmaniasis": 0.547,
}


# ---------------------------------------------------------------------------
# Map from CSV target names to our target pair IDs
# ---------------------------------------------------------------------------
TARGET_NAME_MAP: dict[str, str] = {
    # Trypanosomiasis
    "Cathepsin B-like cysteine protease": "TbCathB",
    "Class 1 phosphodiesterase PDEB1": "TbPDEB1",
    # Leishmaniasis
    "Pteridine reductase 1": "LmPTR1",
    "Bifunctional dihydrofolate reductase-thymidylate synthase": "LmDHFR",
    "Class I phosphodiesterase PDEB1": "TbPDEB1",  # Same enzyme, shared
    # Schistosomiasis (from publication tables)
    "Histone deacetylase 8": "SmHDAC8",
    "Dihydroorotate dehydrogenase (quinone), mitochondrial": "SmDHODH",
    "Dihydroorotate dehydrogenase (quino": "SmDHODH",  # Truncated in CSV
}


def find_data_dir() -> Path:
    """Find the Kira data directory."""
    # Try relative to repo root
    candidates = [
        Path(__file__).resolve().parents[3] / "data",  # src/kira/experiments/ → repo/data
        Path.home() / "Kira" / "data",
        Path.home() / "kira" / "data",
        Path("data"),
    ]
    for p in candidates:
        if p.is_dir():
            return p
    raise FileNotFoundError(
        "Cannot find Kira data directory. Run from ~/Kira or set KIRA_DATA env var."
    )


def load_selectivity_data(data_dir: Path) -> pd.DataFrame:
    """Load compound-level selectivity data from all disease CSVs.

    Returns
    -------
    pd.DataFrame
        Columns: molecule_chembl_id, parasite_ic50, human_ic50,
                 selectivity_ratio, parasite_target, disease,
                 selectivity_class, target_pair_id, is_selective
    """
    frames = []

    # Leishmaniasis
    leish_path = data_dir / "leishmania" / "leish_selectivity.csv"
    if leish_path.exists():
        df = pd.read_csv(leish_path)
        if "disease" not in df.columns:
            df["disease"] = "Leishmaniasis"
        frames.append(df)
        print(f"  Loaded {len(df)} Leishmaniasis compounds from {leish_path.name}")

    # Trypanosomiasis
    tryp_path = data_dir / "trypanosoma" / "tryp_selectivity_expanded.csv"
    if not tryp_path.exists():
        tryp_path = data_dir / "trypanosoma" / "tryp_selectivity.csv"
    if tryp_path.exists():
        df = pd.read_csv(tryp_path)
        if "disease" not in df.columns:
            df["disease"] = "Trypanosomiasis"
        frames.append(df)
        print(f"  Loaded {len(df)} Trypanosomiasis compounds from {tryp_path.name}")

    # Schistosomiasis (from discovery_candidates.csv)
    schisto_path = data_dir / "publication" / "discovery_candidates.csv"
    if schisto_path.exists():
        dc = pd.read_csv(schisto_path)
        # Filter to targets with known orthologues
        schisto_targets = {
            "Dihydroorotate dehydrogenase (quinone), mitochondrial",
            "Histone deacetylase 8",
        }
        dc = dc[dc.best_target.isin(schisto_targets)].copy()
        dc = dc[dc.best_selectivity_ratio.notna()].copy()
        if len(dc) > 0:
            # Reshape to match other CSVs
            schisto_df = pd.DataFrame({
                "molecule_chembl_id": dc.molecule_chembl_id,
                "selectivity_ratio": dc.best_selectivity_ratio,
                "parasite_target": dc.best_target,
                "disease": "Schistosomiasis",
            })
            frames.append(schisto_df)
            print(f"  Loaded {len(schisto_df)} Schistosomiasis compounds from {schisto_path.name}")

    if not frames:
        raise FileNotFoundError(f"No selectivity CSVs found in {data_dir}")

    # Combine
    all_data = pd.concat(frames, ignore_index=True)

    # Standardize columns
    if "selectivity_ratio" not in all_data.columns:
        raise ValueError("Missing selectivity_ratio column")

    # Map target names to pair IDs
    all_data["target_pair_id"] = all_data["parasite_target"].map(TARGET_NAME_MAP)

    # Binary label: selective = ratio >= 10
    all_data["is_selective"] = (all_data["selectivity_ratio"] >= 10.0).astype(int)

    # Drop rows without a known target pair
    n_before = len(all_data)
    all_data = all_data.dropna(subset=["target_pair_id"])
    n_after = len(all_data)
    if n_before > n_after:
        print(f"  Dropped {n_before - n_after} compounds with unmapped targets")

    return all_data


def compute_all_pocket_features() -> dict[str, PocketFeatures]:
    """Compute binding-site divergence features for all target pairs.

    Returns
    -------
    dict[str, PocketFeatures]
        Map from target pair ID to computed features.
    """
    features = {}
    for pair_id, pair in TARGET_PAIRS.items():
        pf = compute_pocket_features(
            pair.parasite_pocket_sequence,
            pair.human_pocket_sequence,
            pair_name=pair_id,
            disease=pair.disease,
        )
        features[pair_id] = pf
    return features


def build_feature_matrix(
    data: pd.DataFrame,
    pocket_features: dict[str, PocketFeatures],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the feature matrix and label vector.

    For each compound, the feature vector is:
        [15 pocket divergence features]

    (Morgan fingerprints excluded for now to isolate the effect
    of pocket features alone. Adding FPs is a follow-up experiment.)

    Parameters
    ----------
    data : pd.DataFrame
        Selectivity data with target_pair_id and is_selective columns.
    pocket_features : dict[str, PocketFeatures]
        Pre-computed pocket features per target pair.

    Returns
    -------
    X : np.ndarray
        (N, 15) feature matrix.
    y : np.ndarray
        (N,) binary labels (1 = selective, 0 = not selective).
    diseases : np.ndarray
        (N,) disease labels for LODO evaluation.
    """
    X_rows = []
    y_rows = []
    disease_rows = []

    for _, row in data.iterrows():
        pair_id = row["target_pair_id"]
        if pair_id not in pocket_features:
            continue

        pf = pocket_features[pair_id]
        X_rows.append(pf.to_array())
        y_rows.append(int(row["is_selective"]))
        disease_rows.append(row["disease"])

    X = np.array(X_rows, dtype=np.float32)
    y = np.array(y_rows, dtype=np.int32)
    diseases = np.array(disease_rows)

    return X, y, diseases


def evaluate_cv(X: np.ndarray, y: np.ndarray, n_splits: int = 5) -> float:
    """5-fold stratified cross-validation AUROC.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix.
    y : np.ndarray
        Binary labels.
    n_splits : int
        Number of CV folds.

    Returns
    -------
    float
        Mean AUROC across folds.
    """
    if len(np.unique(y)) < 2:
        print("  WARNING: Only one class present — cannot compute AUROC")
        return 0.5

    # CRITICAL: Scaler must be inside the CV loop to prevent feature leakage.
    # Using sklearn Pipeline ensures fit_transform happens only on training folds.
    from sklearn.pipeline import Pipeline

    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", GradientBoostingClassifier(
            n_estimators=100,
            max_depth=3,
            random_state=42,
            min_samples_leaf=5,
        )),
    ])

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = cross_val_score(pipe, X, y, cv=cv, scoring="roc_auc")

    return float(np.mean(scores))


def evaluate_lodo(
    X: np.ndarray, y: np.ndarray, diseases: np.ndarray
) -> dict[str, float]:
    """Leave-one-disease-out evaluation.

    Train on two diseases, test on the third. This is the critical test:
    do the features TRANSFER across diseases?

    Parameters
    ----------
    X : np.ndarray
        Feature matrix.
    y : np.ndarray
        Binary labels.
    diseases : np.ndarray
        Disease labels.

    Returns
    -------
    dict[str, float]
        Disease name → AUROC when that disease is held out.
    """
    unique_diseases = sorted(set(diseases))
    results = {}

    clf = GradientBoostingClassifier(
        n_estimators=100,
        max_depth=3,
        random_state=42,
        min_samples_leaf=5,
    )

    scaler = StandardScaler()

    for held_out in unique_diseases:
        train_mask = diseases != held_out
        test_mask = diseases == held_out

        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]

        if len(np.unique(y_test)) < 2:
            print(f"  {held_out}: Only one class in test set — skipping")
            results[held_out] = float("nan")
            continue

        if len(np.unique(y_train)) < 2:
            print(f"  {held_out}: Only one class in training set — skipping")
            results[held_out] = float("nan")
            continue

        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        clf.fit(X_train_s, y_train)
        y_prob = clf.predict_proba(X_test_s)[:, 1]
        auroc = roc_auc_score(y_test, y_prob)
        results[held_out] = float(auroc)

    return results


def print_pocket_summary(pocket_features: dict[str, PocketFeatures]) -> None:
    """Print a summary of pocket features for each target pair."""
    print("\n" + "=" * 70)
    print("POCKET DIVERGENCE PROFILES")
    print("=" * 70)

    for pair_id, pf in sorted(pocket_features.items()):
        print(f"\n  {pair_id} ({pf.disease})")
        print(f"    Pocket size: {pf.pocket_size} residues")
        print(f"    Pocket identity: {pf.pocket_identity:.1%}")
        print(f"    Mean physicochemical distance: {pf.mean_physicochemical_distance:.3f}")
        print(f"    Max physicochemical distance: {pf.max_physicochemical_distance:.3f}")
        print(f"    Charge changes: {pf.n_charge_changes}")
        print(f"    Hydrophobicity flips: {pf.n_hydrophobicity_flips}")
        print(f"    Volume changes: {pf.n_volume_changes}")


def print_results(
    cv_auroc: float,
    lodo_results: dict[str, float],
    n_compounds: int,
    n_selective: int,
) -> None:
    """Print the comparison between baseline and new features."""
    print("\n" + "=" * 70)
    print("RESULTS: BINDING-SITE FEATURES vs ESM-2 BASELINE")
    print("=" * 70)

    print(f"\n  Data: {n_compounds} compounds, {n_selective} selective ({n_selective/n_compounds:.1%})")

    print(f"\n  {'Metric':<30} {'ESM-2 (Script 19)':<20} {'Pocket Features':<20} {'Delta':<10}")
    print(f"  {'-'*80}")

    # 5-fold CV
    delta_cv = cv_auroc - BASELINE_CV_AUROC
    sign_cv = "+" if delta_cv >= 0 else ""
    print(f"  {'5-fold CV AUROC':<30} {BASELINE_CV_AUROC:<20.3f} {cv_auroc:<20.3f} {sign_cv}{delta_cv:.3f}")

    # LODO
    for disease in sorted(BASELINE_LODO.keys()):
        baseline = BASELINE_LODO[disease]
        new = lodo_results.get(disease, float("nan"))
        if np.isnan(new):
            delta_str = "N/A"
            new_str = "N/A"
        else:
            delta = new - baseline
            sign = "+" if delta >= 0 else ""
            delta_str = f"{sign}{delta:.3f}"
            new_str = f"{new:.3f}"
        print(f"  {'LODO ' + disease:<30} {baseline:<20.3f} {new_str:<20} {delta_str}")

    # Interpretation
    print(f"\n  {'='*70}")
    print("  INTERPRETATION:")

    lodo_values = [v for v in lodo_results.values() if not np.isnan(v)]
    if lodo_values:
        mean_lodo = np.mean(lodo_values)
        mean_baseline = np.mean(list(BASELINE_LODO.values()))
        if mean_lodo > mean_baseline:
            print(f"  Pocket features IMPROVE cross-disease transfer "
                  f"(mean LODO: {mean_lodo:.3f} vs {mean_baseline:.3f})")
            print("  → Binding-site divergence captures transferable selectivity signal")
            print("  → ESM-2 global embeddings miss this signal")
        elif mean_lodo > 0.5:
            print(f"  Pocket features show SOME cross-disease transfer "
                  f"(mean LODO: {mean_lodo:.3f})")
            print("  → Signal exists but may need more pocket position data")
        else:
            print(f"  Pocket features do NOT improve LODO (mean: {mean_lodo:.3f})")
            print("  → Check pocket definitions — may need structural verification")

    if cv_auroc < BASELINE_CV_AUROC:
        print(f"\n  NOTE: CV AUROC is lower ({cv_auroc:.3f} vs {BASELINE_CV_AUROC:.3f}).")
        print("  This is EXPECTED — Script 19's high CV comes from compound fingerprints")
        print("  which overfit to the training distribution. Pocket features alone")
        print("  capture the TARGET-PAIR signal, not compound-specific patterns.")
        print("  The real test is LODO, not CV.")


def main():
    """Run the selectivity prediction v3 experiment."""
    print("=" * 70)
    print("SELECTIVITY PREDICTION v3")
    print("Binding-Site Divergence Features vs ESM-2 Global Embeddings")
    print("=" * 70)

    # 1. Compute pocket features for all target pairs
    print("\n[1/4] Computing pocket divergence features...")
    pocket_features = compute_all_pocket_features()
    print_pocket_summary(pocket_features)

    # 2. Load selectivity data
    print("\n[2/4] Loading selectivity data...")
    data_dir = find_data_dir()
    data = load_selectivity_data(data_dir)
    print(f"  Total: {len(data)} compounds across "
          f"{data['target_pair_id'].nunique()} target pairs")
    print(f"  Selective (≥10x): {data['is_selective'].sum()} "
          f"({data['is_selective'].mean():.1%})")
    print(f"  Diseases: {', '.join(sorted(data['disease'].unique()))}")

    # Per-target breakdown
    print("\n  Per-target breakdown:")
    for pair_id, group in data.groupby("target_pair_id"):
        n_sel = group["is_selective"].sum()
        pct = group["is_selective"].mean() * 100
        print(f"    {pair_id}: {len(group)} compounds, {n_sel} selective ({pct:.1f}%)")

    # 3. Build feature matrix
    print("\n[3/4] Building feature matrix...")
    X, y, diseases = build_feature_matrix(data, pocket_features)
    print(f"  Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"  Labels: {y.sum()} selective, {len(y) - y.sum()} non-selective")
    print(f"  Feature names: {PocketFeatures.feature_names()[:5]}... (15 total)")

    # 4. Evaluate
    print("\n[4/4] Evaluating...")

    # 5-fold CV
    print("\n  Running 5-fold stratified CV...")
    cv_auroc = evaluate_cv(X, y)
    print(f"  5-fold CV AUROC: {cv_auroc:.3f}")

    # Leave-one-disease-out
    print("\n  Running leave-one-disease-out...")
    lodo_results = evaluate_lodo(X, y, diseases)
    for disease, auroc in sorted(lodo_results.items()):
        if np.isnan(auroc):
            print(f"  LODO {disease}: N/A (single class)")
        else:
            print(f"  LODO {disease}: {auroc:.3f}")

    # Print comparison
    n_selective = int(y.sum())
    print_results(cv_auroc, lodo_results, len(y), n_selective)

    # 5. Feature importance analysis
    print("\n" + "=" * 70)
    print("FEATURE IMPORTANCE (what drives selectivity prediction)")
    print("=" * 70)

    clf_full = GradientBoostingClassifier(
        n_estimators=100, max_depth=3, random_state=42, min_samples_leaf=5,
    )
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    clf_full.fit(X_scaled, y)

    importances = clf_full.feature_importances_
    names = PocketFeatures.feature_names()
    ranked = sorted(zip(names, importances), key=lambda x: x[1], reverse=True)

    print(f"\n  {'Feature':<40} {'Importance':<12} {'Interpretation'}")
    print(f"  {'-'*80}")

    interpretations = {
        "pocket_identity": "How similar the pocket is overall",
        "pocket_divergence": "How different the pocket is (1 - identity)",
        "mean_physicochemical_distance": "Average property difference per position",
        "max_physicochemical_distance": "Worst-case single-position difference",
        "std_physicochemical_distance": "How variable the divergence is",
        "n_charge_changes": "Positions where electrostatics flip",
        "n_hydrophobicity_flips": "Positions where polarity reverses",
        "n_volume_changes": "Positions with large size changes",
        "fraction_charge_changes": "Fraction of pocket with charge flips",
        "fraction_hydro_flips": "Fraction with polarity reversals",
        "fraction_volume_changes": "Fraction with size changes",
        "total_hydrophobicity_delta": "Cumulative polarity difference",
        "total_charge_delta": "Cumulative charge difference",
        "total_volume_delta": "Cumulative size difference",
        "pocket_size": "Number of pocket residues analyzed",
    }

    for name, imp in ranked:
        interp = interpretations.get(name, "")
        bar = "█" * int(imp * 50)
        print(f"  {name:<40} {imp:<12.3f} {bar} {interp}")

    print(f"\n  Top driver: {ranked[0][0]} ({ranked[0][1]:.3f})")
    if ranked[0][1] > 0.15:
        print(f"  → This feature alone explains {ranked[0][1]*100:.0f}% of the selectivity signal")


if __name__ == "__main__":
    main()
