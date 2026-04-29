"""Build a tiered-evidence v5 selectivity core from candidate expansion rows.

The tiered core is a data substrate, not a new benchmark claim. It keeps strict
exact evidence when available, admits one-sided bounded evidence only when the
threshold decision is forced, and excludes unresolved or conflicting evidence
from the trainable subset.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

DEFAULT_SELECTIVITY_THRESHOLD = 10.0
DECIDABLE_STATUSES = {
    "exact_matched_ratio",
    "lower_bound_ratio",
    "upper_bound_ratio",
    "interval_ratio",
}


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _float_series(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def _joined_unique(values: Iterable[Any]) -> str:
    cleaned = sorted({str(value) for value in values if _is_present(value)})
    return ";".join(cleaned)


def _first_present(values: Iterable[Any]) -> Any:
    for value in values:
        if _is_present(value):
            return value
    return None


def _count_unique_present(values: Iterable[Any]) -> int:
    return len({value for value in values if _is_present(value)})


def add_compound_key(df: pd.DataFrame) -> pd.DataFrame:
    """Prefer standard InChIKey and fall back to ChEMBL molecule id."""

    out = df.copy()
    if "standard_inchi_key" not in out.columns:
        out["standard_inchi_key"] = None
    if "molecule_chembl_id" not in out.columns:
        raise ValueError("candidate table is missing molecule_chembl_id")

    out["compound_key"] = out["standard_inchi_key"].where(
        out["standard_inchi_key"].map(_is_present),
        out["molecule_chembl_id"],
    )
    return out


def classify_tiered_candidate(
    row: pd.Series,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> tuple[str | None, str]:
    """Classify one candidate row using the v5 tiered evidence rules."""

    status = row.get("candidate_status")
    exact = pd.to_numeric(row.get("ratio_exact_human_div_parasite"), errors="coerce")
    lower = pd.to_numeric(row.get("ratio_lower_bound_human_div_parasite"), errors="coerce")
    upper = pd.to_numeric(row.get("ratio_upper_bound_human_div_parasite"), errors="coerce")

    if status == "exact_matched_ratio":
        if pd.notna(exact) and exact > 0:
            if exact >= selectivity_threshold:
                return "positive", "exact_matched_ratio_at_or_above_threshold"
            return "negative", "exact_matched_ratio_below_threshold"
        return None, "exact_matched_ratio_missing_positive_ratio"

    if status == "lower_bound_ratio":
        if pd.notna(lower) and lower >= selectivity_threshold:
            return "positive", "lower_bound_ratio_at_or_above_threshold"
        return None, "lower_bound_ratio_below_threshold_unresolved"

    if status == "upper_bound_ratio":
        if pd.notna(upper) and upper < selectivity_threshold:
            return "negative", "upper_bound_ratio_below_threshold"
        return None, "upper_bound_ratio_at_or_above_threshold_unresolved"

    if status == "interval_ratio":
        if pd.notna(lower) and lower >= selectivity_threshold:
            return "positive", "interval_lower_bound_at_or_above_threshold"
        if pd.notna(upper) and upper < selectivity_threshold:
            return "negative", "interval_upper_bound_below_threshold"
        return None, "interval_crosses_threshold_unresolved"

    return None, f"{status}_not_tiered_decidable"


def add_tiered_evidence_labels(
    candidate_rows: pd.DataFrame,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> pd.DataFrame:
    """Add compound keys, tiered labels, and tiered evidence reasons."""

    required = {
        "pair_id",
        "molecule_chembl_id",
        "candidate_status",
        "canonical_smiles",
        "murcko_scaffold",
        "ratio_exact_human_div_parasite",
        "ratio_lower_bound_human_div_parasite",
        "ratio_upper_bound_human_div_parasite",
    }
    missing = sorted(required - set(candidate_rows.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {missing}")

    out = add_compound_key(candidate_rows)
    classifications = [
        classify_tiered_candidate(row, selectivity_threshold=selectivity_threshold)
        for _, row in out.iterrows()
    ]
    out["tiered_label"] = [label for label, _reason in classifications]
    out["tiered_evidence_reason"] = [reason for _label, reason in classifications]
    out["is_tiered_decidable"] = out["tiered_label"].map(_is_present)

    for ratio_col in (
        "ratio_exact_human_div_parasite",
        "ratio_lower_bound_human_div_parasite",
        "ratio_upper_bound_human_div_parasite",
    ):
        out[ratio_col] = _float_series(out[ratio_col])

    exact_ratio = out["ratio_exact_human_div_parasite"]
    out["log10_ratio_exact_human_div_parasite"] = exact_ratio.where(exact_ratio > 0).map(
        lambda value: math.log10(value) if pd.notna(value) else math.nan
    )
    return out


def aggregate_tiered_core(
    labeled_rows: pd.DataFrame,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> pd.DataFrame:
    """Collapse decidable tiered evidence to one row per pair_id/compound_key."""

    decidable = labeled_rows[labeled_rows["is_tiered_decidable"].astype(bool)].copy()
    if decidable.empty:
        return pd.DataFrame(
            columns=[
                "tiered_core_row_id",
                "pair_id",
                "compound_key",
                "molecule_chembl_id",
                "canonical_smiles",
                "standard_inchi_key",
                "murcko_scaffold",
                "tiered_label",
                "is_selective",
                "label_status",
                "is_trainable_tiered_core",
                "tiered_evidence_count",
            ]
        )

    rows: list[dict[str, Any]] = []
    for (pair_id, compound_key), group in decidable.groupby(["pair_id", "compound_key"], sort=True, dropna=False):
        labels = set(group["tiered_label"].dropna())
        has_positive = "positive" in labels
        has_negative = "negative" in labels
        label_conflict = has_positive and has_negative
        tiered_label = "conflicting" if label_conflict else _first_present(group["tiered_label"])

        exact_ratios = _float_series(group["ratio_exact_human_div_parasite"]).dropna()
        exact_ratios = exact_ratios[exact_ratios > 0]
        exact_log_ratios = exact_ratios.map(math.log10)
        median_log_ratio = float(exact_log_ratios.median()) if not exact_log_ratios.empty else math.nan
        median_ratio = float(10**median_log_ratio) if not math.isnan(median_log_ratio) else math.nan

        canonical_smiles = _first_present(group["canonical_smiles"])
        murcko_scaffold = _first_present(group["murcko_scaffold"])
        standard_inchi_key = _first_present(group.get("standard_inchi_key", []))
        molecule_chembl_id = _first_present(group["molecule_chembl_id"])
        is_selective = 1 if tiered_label == "positive" else 0 if tiered_label == "negative" else None
        label_status = "tiered_core_conflicting" if label_conflict else "tiered_core_consistent"
        is_trainable = bool(
            not label_conflict
            and tiered_label in {"positive", "negative"}
            and _is_present(compound_key)
            and _is_present(canonical_smiles)
            and _is_present(murcko_scaffold)
        )

        lower_bounds = _float_series(group["ratio_lower_bound_human_div_parasite"]).dropna()
        upper_bounds = _float_series(group["ratio_upper_bound_human_div_parasite"]).dropna()

        rows.append(
            {
                "tiered_core_row_id": f"{pair_id}::{compound_key}",
                "pair_id": pair_id,
                "parasite_label": _first_present(group.get("parasite_label", [])),
                "human_label": _first_present(group.get("human_label", [])),
                "compound_key": compound_key,
                "molecule_chembl_id": molecule_chembl_id,
                "canonical_smiles": canonical_smiles,
                "standard_inchi_key": standard_inchi_key,
                "murcko_scaffold": murcko_scaffold,
                "tiered_label": tiered_label,
                "is_selective": is_selective,
                "selectivity_threshold": float(selectivity_threshold),
                "label_status": label_status,
                "is_trainable_tiered_core": int(is_trainable),
                "tiered_evidence_count": int(len(group)),
                "exact_evidence_count": int(group["candidate_status"].eq("exact_matched_ratio").sum()),
                "bounded_evidence_count": int(group["candidate_status"].isin({"lower_bound_ratio", "upper_bound_ratio", "interval_ratio"}).sum()),
                "n_positive_tiered_observations": int(group["tiered_label"].eq("positive").sum()),
                "n_negative_tiered_observations": int(group["tiered_label"].eq("negative").sum()),
                "selectivity_ratio_human_div_parasite": median_ratio,
                "log10_selectivity_ratio": median_log_ratio,
                "ratio_exact_min_human_div_parasite": float(exact_ratios.min()) if not exact_ratios.empty else math.nan,
                "ratio_exact_median_human_div_parasite": median_ratio,
                "ratio_exact_max_human_div_parasite": float(exact_ratios.max()) if not exact_ratios.empty else math.nan,
                "ratio_lower_bound_min_human_div_parasite": float(lower_bounds.min()) if not lower_bounds.empty else math.nan,
                "ratio_lower_bound_max_human_div_parasite": float(lower_bounds.max()) if not lower_bounds.empty else math.nan,
                "ratio_upper_bound_min_human_div_parasite": float(upper_bounds.min()) if not upper_bounds.empty else math.nan,
                "ratio_upper_bound_max_human_div_parasite": float(upper_bounds.max()) if not upper_bounds.empty else math.nan,
                "candidate_row_ids": _joined_unique(group.get("candidate_row_id", [])),
                "evidence_statuses": _joined_unique(group["candidate_status"]),
                "comparison_tiers": _joined_unique(group.get("comparison_tier", [])),
                "tiered_evidence_reasons": _joined_unique(group["tiered_evidence_reason"]),
                "parasite_standard_types": _joined_unique(group.get("parasite_standard_type", [])),
                "human_standard_types": _joined_unique(group.get("human_standard_type", [])),
                "parasite_standard_relations": _joined_unique(group.get("parasite_standard_relation", [])),
                "human_standard_relations": _joined_unique(group.get("human_standard_relation", [])),
                "parasite_measurement_interval_kinds": _joined_unique(group.get("parasite_measurement_interval_kind", [])),
                "human_measurement_interval_kinds": _joined_unique(group.get("human_measurement_interval_kind", [])),
                "parasite_assay_types": _joined_unique(group.get("parasite_assay_type", [])),
                "human_assay_types": _joined_unique(group.get("human_assay_type", [])),
                "parasite_activity_chembl_ids": _joined_unique(group.get("parasite_activity_chembl_id", [])),
                "human_activity_chembl_ids": _joined_unique(group.get("human_activity_chembl_id", [])),
                "parasite_assay_chembl_ids": _joined_unique(group.get("parasite_assay_chembl_id", [])),
                "human_assay_chembl_ids": _joined_unique(group.get("human_assay_chembl_id", [])),
                "parasite_document_chembl_ids": _joined_unique(group.get("parasite_document_chembl_id", [])),
                "human_document_chembl_ids": _joined_unique(group.get("human_document_chembl_id", [])),
                "n_parasite_documents": _count_unique_present(group.get("parasite_document_chembl_id", [])),
                "n_human_documents": _count_unique_present(group.get("human_document_chembl_id", [])),
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["pair_id", "compound_key"]).reset_index(drop=True)
    return out


def _counter_by_status(rows: pd.DataFrame) -> dict[str, int]:
    if rows.empty or "candidate_status" not in rows.columns:
        return {}
    return dict(sorted(Counter(rows["candidate_status"].dropna()).items()))


def _trainable_subset(core_rows: pd.DataFrame) -> pd.DataFrame:
    if core_rows.empty:
        return core_rows
    return core_rows[core_rows["is_trainable_tiered_core"].astype(bool)]


def summarize_tiered_core(
    candidate_rows: pd.DataFrame,
    labeled_rows: pd.DataFrame,
    core_rows: pd.DataFrame,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> dict[str, Any]:
    trainable = _trainable_subset(core_rows)
    decidable = labeled_rows[labeled_rows["is_tiered_decidable"].astype(bool)]
    unresolved = labeled_rows[~labeled_rows["is_tiered_decidable"].astype(bool)]
    conflicts = core_rows[core_rows["label_status"] == "tiered_core_conflicting"] if not core_rows.empty else core_rows

    trainable_rows_by_pair = dict(sorted(Counter(trainable["pair_id"]).items())) if not trainable.empty else {}
    trainable_positive_rows_by_pair = (
        dict(sorted({pair: int(rows["is_selective"].sum()) for pair, rows in trainable.groupby("pair_id", sort=True)}.items()))
        if not trainable.empty
        else {}
    )
    trainable_negative_rows_by_pair = (
        dict(
            sorted(
                {
                    pair: int(len(rows) - rows["is_selective"].sum())
                    for pair, rows in trainable.groupby("pair_id", sort=True)
                }.items()
            )
        )
        if not trainable.empty
        else {}
    )
    trainable_unique_scaffolds_by_pair = (
        dict(
            sorted(
                {
                    pair: int(rows["murcko_scaffold"].nunique(dropna=True))
                    for pair, rows in trainable.groupby("pair_id", sort=True)
                }.items()
            )
        )
        if not trainable.empty
        else {}
    )
    pairs_with_both_classes = (
        sorted(
            pair
            for pair, rows in trainable.groupby("pair_id", sort=True)
            if set(rows["is_selective"].dropna().astype(int)) == {0, 1}
        )
        if not trainable.empty
        else []
    )

    return {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "scientific_scope": (
            "Kira v5 tiered-core aggregation is a tiered evidence substrate after the strict exact core. "
            "It is not a final modeling claim or a new benchmark model; v4 remains the current modeling claim."
        ),
        "selectivity_threshold_human_div_parasite": float(selectivity_threshold),
        "input_candidate_rows": int(len(candidate_rows)),
        "decidable_candidate_rows": int(len(decidable)),
        "unresolved_candidate_rows": int(len(unresolved)),
        "tiered_core_rows": int(len(core_rows)),
        "trainable_tiered_core_rows": int(len(trainable)),
        "conflicting_tiered_core_rows": int(len(conflicts)),
        "candidate_rows_by_evidence_status": _counter_by_status(labeled_rows),
        "decidable_rows_by_evidence_status": _counter_by_status(decidable),
        "trainable_rows_by_pair": trainable_rows_by_pair,
        "trainable_positive_rows_by_pair": trainable_positive_rows_by_pair,
        "trainable_negative_rows_by_pair": trainable_negative_rows_by_pair,
        "trainable_unique_scaffolds_by_pair": trainable_unique_scaffolds_by_pair,
        "pairs_with_both_classes": pairs_with_both_classes,
    }


def build_tiered_core(
    candidate_csv: Path,
    output_csv: Path,
    output_summary_json: Path,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    candidate_rows = pd.read_csv(candidate_csv)
    labeled_rows = add_tiered_evidence_labels(candidate_rows, selectivity_threshold=selectivity_threshold)
    core_rows = aggregate_tiered_core(labeled_rows, selectivity_threshold=selectivity_threshold)
    summary = summarize_tiered_core(
        candidate_rows=candidate_rows,
        labeled_rows=labeled_rows,
        core_rows=core_rows,
        selectivity_threshold=selectivity_threshold,
    )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_summary_json.parent.mkdir(parents=True, exist_ok=True)
    core_rows.to_csv(output_csv, index=False)
    with output_summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return core_rows, summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-csv",
        type=Path,
        default=Path("data/processed/selectivity_v5_candidate_rows.csv"),
        help="Input v5 candidate rows CSV.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/processed/selectivity_v5_tiered_core_rows.csv"),
        help="Output tiered-core CSV, one row per target-pair/compound key.",
    )
    parser.add_argument(
        "--output-summary-json",
        type=Path,
        default=Path("data/processed/selectivity_v5_tiered_core_summary.json"),
        help="Output tiered-core summary JSON.",
    )
    parser.add_argument(
        "--selectivity-threshold",
        type=float,
        default=DEFAULT_SELECTIVITY_THRESHOLD,
        help="Human/parasite ratio threshold for selectivity classification. Default: 10.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    core_rows, summary = build_tiered_core(
        candidate_csv=args.candidate_csv,
        output_csv=args.output_csv,
        output_summary_json=args.output_summary_json,
        selectivity_threshold=args.selectivity_threshold,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote:\n  - {args.output_csv}\n  - {args.output_summary_json}")
    print(f"\nTiered-core rows: {len(core_rows)}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
