#!/usr/bin/env python
"""
Kira v5 Exact-Core Benchmark Runner
Compares B0 (pair-only), B1 (compound-only), B2 (compound + pair)
with scaffold-grouped cross-validation.
"""

import argparse

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupKFold


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", default="data/processed/selectivity_v5_exact_core_rows.csv")
    parser.add_argument("--output-dir", default="results/selectivity_v5_exact_core")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    trainable = df[df["is_trainable_exact_core"] == 1].copy()

    print(f"Running benchmark on {len(trainable)} trainable rows")

    # Features (simple baselines)
    # TODO: replace with proper fingerprints / embeddings later
    trainable["pair_feature"] = trainable["pair_id"].astype("category").cat.codes
    # Placeholder for compound features (e.g., Morgan fingerprint later)

    X = trainable[["pair_feature"]]  # B0 example
    y = trainable["is_selective"].astype(int)
    groups = trainable["murcko_scaffold"]  # or pair_id::murcko_scaffold

    # Grouped CV
    gkf = GroupKFold(n_splits=5)
    results = []

    for train_idx, val_idx in gkf.split(X, y, groups=groups):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight="balanced")
        model.fit(X_train, y_train)
        pred = model.predict(X_val)
        prob = model.predict_proba(X_val)[:, 1]

        results.append({
            "accuracy": accuracy_score(y_val, pred),
            "balanced_acc": balanced_accuracy_score(y_val, pred),
            "roc_auc": roc_auc_score(y_val, prob) if len(set(y_val)) > 1 else 0.5
        })

    print("Benchmark results (B0 placeholder):")
    for metric in ["accuracy", "balanced_acc", "roc_auc"]:
        print(f"{metric}: {pd.Series([r[metric] for r in results]).mean():.4f}")

if __name__ == "__main__":
    main()
