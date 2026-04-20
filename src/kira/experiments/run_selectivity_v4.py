from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from kira.experiments.selectivity_v4_features import (
    block_slices,
    featurize_dataframe,
    select_feature_blocks,
)

DEFAULT_TRAINABLE_PATH = Path("data/processed/selectivity_v4_rows_primary_trainable.csv")
DEFAULT_RESULTS_DIR = Path("results/selectivity_v4")
DEFAULT_N_BITS = 256
DEFAULT_N_SPLITS = 5
RANDOM_STATE = 42

ABLATIONS: dict[str, list[str]] = {
    "A0_pair_only": ["pair_delta"],
    "A1_compound_only": ["compound_fp", "compound_desc"],
    "A2_compound_plus_pair": ["compound_fp", "compound_desc", "pair_delta"],
    "A2b_pair_plus_side": [
        "compound_fp",
        "compound_desc",
        "pair_delta",
        "parasite_side",
        "human_side",
    ],
    "A2c_pair_plus_compat": [
        "compound_fp",
        "compound_desc",
        "pair_delta",
        "compat_diff",
    ],
    "A3_full_v4": [
        "compound_fp",
        "compound_desc",
        "pair_delta",
        "parasite_side",
        "human_side",
        "compat_diff",
    ],
}

BASE_COLUMNS = [
    "molecule_chembl_id",
    "compound_name",
    "target_pair_id",
    "murcko_scaffold",
    "selectivity_ratio",
    "is_selective",
]


def load_trainable_table(path: Path = DEFAULT_TRAINABLE_PATH) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing trainable table: {path}")

    df = pd.read_csv(path).copy()

    required = {
        "canonical_smiles",
        "target_pair_id",
        "murcko_scaffold",
        "selectivity_ratio",
        "is_selective",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Trainable table missing columns: {sorted(missing)}")

    df["is_selective"] = pd.to_numeric(df["is_selective"], errors="raise").astype(int)
    df["selectivity_ratio"] = pd.to_numeric(df["selectivity_ratio"], errors="coerce")

    return df


def make_scaffold_groups(df: pd.DataFrame) -> np.ndarray:
    return (
        df["target_pair_id"].astype(str)
        + "::"
        + df["murcko_scaffold"].fillna("NA").astype(str)
    ).to_numpy()


def build_estimator() -> Pipeline:
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    l1_ratio=0.25,
                    C=0.5,
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=RANDOM_STATE,
                ),
            ),
        ]
    )


def safe_auroc(y_true: np.ndarray, y_score: np.ndarray) -> float | None:
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)

    mask = np.isfinite(y_true) & np.isfinite(y_score)
    y_true = y_true[mask]
    y_score = y_score[mask]

    if len(y_true) == 0 or np.unique(y_true).size < 2:
        return None

    return float(roc_auc_score(y_true, y_score))


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 2:
        return None
    if np.unique(x).size < 2:
        return None
    if np.unique(y).size < 2:
        return None

    corr, _ = spearmanr(x, y)
    if corr is None or not np.isfinite(corr):
        return None
    return float(corr)


def mean_or_none(values: pd.Series) -> float | None:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if len(vals) == 0:
        return None
    return float(vals.mean())


def out_of_fold_scores(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    n_splits: int = DEFAULT_N_SPLITS,
) -> tuple[np.ndarray, pd.DataFrame]:
    splitter = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    scores = np.full(X.shape[0], np.nan, dtype=float)
    fold_rows: list[dict] = []

    for fold_idx, (train_idx, test_idx) in enumerate(
        splitter.split(X, y, groups), start=1
    ):
        model = build_estimator()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model.fit(X[train_idx], y[train_idx])

        pred = model.predict_proba(X[test_idx])[:, 1]
        scores[test_idx] = pred

        fold_rows.append(
            {
                "fold": fold_idx,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "n_train_positive": int(y[train_idx].sum()),
                "n_test_positive": int(y[test_idx].sum()),
                "n_train_groups": int(len(np.unique(groups[train_idx]))),
                "n_test_groups": int(len(np.unique(groups[test_idx]))),
            }
        )

    if np.isnan(scores).any():
        raise RuntimeError("OOF prediction vector contains NaN entries")

    return scores, pd.DataFrame(fold_rows)


def evaluate_ablation(
    df: pd.DataFrame,
    X_full: np.ndarray,
    *,
    ablation_name: str,
    block_names: list[str],
    n_bits: int = DEFAULT_N_BITS,
    n_splits: int = DEFAULT_N_SPLITS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    X = select_feature_blocks(X_full, block_names, n_bits=n_bits)
    y = df["is_selective"].astype(int).to_numpy()
    groups = make_scaffold_groups(df)

    scores, fold_df = out_of_fold_scores(X, y, groups, n_splits=n_splits)

    pred_df = df[BASE_COLUMNS].copy()
    pred_df["ablation"] = ablation_name
    pred_df["scaffold_group"] = groups
    pred_df["pred_score"] = scores
    pred_df["log10_selectivity_ratio"] = np.where(
        df["selectivity_ratio"].to_numpy(dtype=float) > 0.0,
        np.log10(df["selectivity_ratio"].to_numpy(dtype=float)),
        np.nan,
    )

    pair_rows: list[dict] = []
    for pair_id, sub in pred_df.groupby("target_pair_id", sort=True):
        pair_rows.append(
            {
                "ablation": ablation_name,
                "target_pair_id": pair_id,
                "n_rows": int(len(sub)),
                "n_positive": int(sub["is_selective"].sum()),
                "auroc": safe_auroc(
                    sub["is_selective"].to_numpy(),
                    sub["pred_score"].to_numpy(),
                ),
                "spearman_log_ratio": safe_spearman(
                    sub["pred_score"].to_numpy(),
                    sub["log10_selectivity_ratio"].to_numpy(),
                ),
            }
        )

    per_pair_df = pd.DataFrame(pair_rows)
    fold_df = fold_df.copy()
    fold_df["ablation"] = ablation_name

    summary = {
        "ablation": ablation_name,
        "blocks": block_names,
        "n_features": int(X.shape[1]),
        "n_rows": int(len(df)),
        "n_groups": int(pd.Series(groups).nunique()),
        "n_positive": int(y.sum()),
        "pooled_auroc": safe_auroc(y, scores),
        "macro_pair_auroc": mean_or_none(per_pair_df["auroc"]),
        "macro_pair_spearman": mean_or_none(per_pair_df["spearman_log_ratio"]),
    }

    return pred_df, per_pair_df, fold_df, summary


def main() -> None:
    RDLogger.DisableLog("rdApp.*")

    results_dir = DEFAULT_RESULTS_DIR
    results_dir.mkdir(parents=True, exist_ok=True)

    df = load_trainable_table(DEFAULT_TRAINABLE_PATH)
    X_full, feature_names = featurize_dataframe(df, n_bits=DEFAULT_N_BITS)

    predictions_frames: list[pd.DataFrame] = []
    per_pair_frames: list[pd.DataFrame] = []
    fold_frames: list[pd.DataFrame] = []
    summaries: list[dict] = []

    for ablation_name, block_names in ABLATIONS.items():
        pred_df, pair_df, fold_df, summary = evaluate_ablation(
            df,
            X_full,
            ablation_name=ablation_name,
            block_names=block_names,
            n_bits=DEFAULT_N_BITS,
            n_splits=DEFAULT_N_SPLITS,
        )
        predictions_frames.append(pred_df)
        per_pair_frames.append(pair_df)
        fold_frames.append(fold_df)
        summaries.append(summary)

    predictions = pd.concat(predictions_frames, ignore_index=True)
    per_pair = pd.concat(per_pair_frames, ignore_index=True)
    folds = pd.concat(fold_frames, ignore_index=True)

    summary = {
        "n_rows_trainable": int(len(df)),
        "n_feature_names_full": int(len(feature_names)),
        "n_splits": DEFAULT_N_SPLITS,
        "n_bits": DEFAULT_N_BITS,
        "block_slices": {
            name: [sl.start, sl.stop]
            for name, sl in block_slices(n_bits=DEFAULT_N_BITS).items()
        },
        "ablations": summaries,
    }

    predictions_path = results_dir / "predictions.csv"
    per_pair_path = results_dir / "per_pair_metrics.csv"
    folds_path = results_dir / "fold_metrics.csv"
    summary_path = results_dir / "summary.json"

    predictions.to_csv(predictions_path, index=False)
    per_pair.to_csv(per_pair_path, index=False)
    folds.to_csv(folds_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2))

    print("=" * 80)
    print("run_selectivity_v4: summary")
    print("=" * 80)
    print(json.dumps(summary, indent=2))
    print("\nWrote:")
    print(f"  - {predictions_path}")
    print(f"  - {per_pair_path}")
    print(f"  - {folds_path}")
    print(f"  - {summary_path}")


if __name__ == "__main__":
    main()
