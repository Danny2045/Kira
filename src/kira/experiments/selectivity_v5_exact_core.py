"""Build a benchmark-ready exact-ratio core from v5 candidate expansion rows.

This module is intentionally conservative. The v5 expansion table is an evidence
substrate, not a final training table. This script collapses strict exact paired
activity evidence into one row per target-pair/compound key, preserving replicate
counts and flagging contradictory labels instead of hiding them.
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

DEFAULT_ALLOWED_STANDARD_TYPES = ("IC50",)
DEFAULT_SELECTIVITY_THRESHOLD = 10.0


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return False


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


def add_compound_key(df: pd.DataFrame) -> pd.DataFrame:
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


def filter_exact_core_candidates(
    candidate_rows: pd.DataFrame,
    allowed_standard_types: Sequence[str] = DEFAULT_ALLOWED_STANDARD_TYPES,
) -> pd.DataFrame:
    """Return strict exact rows suitable for compound-pair aggregation."""

    required = {
        "pair_id",
        "molecule_chembl_id",
        "candidate_status",
        "comparison_tier",
        "same_standard_type",
        "same_assay_type",
        "canonical_smiles",
        "murcko_scaffold",
        "ratio_exact_human_div_parasite",
        "parasite_standard_type",
        "human_standard_type",
        "parasite_standard_relation",
        "human_standard_relation",
    }
    missing = sorted(required - set(candidate_rows.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {missing}")

    allowed = {value.upper() for value in allowed_standard_types}
    df = add_compound_key(candidate_rows)
    ratio = _float_series(df["ratio_exact_human_div_parasite"])

    mask = (
        (df["candidate_status"] == "exact_matched_ratio")
        & (df["comparison_tier"] == "strict")
        & df["same_standard_type"].map(_boolish)
        & df["same_assay_type"].map(_boolish)
        & df["canonical_smiles"].map(_is_present)
        & df["murcko_scaffold"].map(_is_present)
        & df["compound_key"].map(_is_present)
        & (df["parasite_standard_relation"] == "=")
        & (df["human_standard_relation"] == "=")
        & df["parasite_standard_type"].astype(str).str.upper().isin(allowed)
        & df["human_standard_type"].astype(str).str.upper().isin(allowed)
        & ratio.notna()
        & (ratio > 0)
    )

    out = df.loc[mask].copy()
    out["ratio_exact_human_div_parasite"] = ratio.loc[mask].astype(float)
    out["log10_ratio_exact_human_div_parasite"] = out["ratio_exact_human_div_parasite"].map(math.log10)
    return out


def aggregate_exact_core(
    exact_rows: pd.DataFrame,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> pd.DataFrame:
    """Collapse exact evidence to one row per pair_id/compound_key."""

    if exact_rows.empty:
        columns = [
            "exact_core_row_id",
            "pair_id",
            "compound_key",
            "molecule_chembl_id",
            "canonical_smiles",
            "standard_inchi_key",
            "murcko_scaffold",
            "selectivity_ratio_human_div_parasite",
            "log10_selectivity_ratio",
            "is_selective",
            "label_status",
            "is_trainable_exact_core",
            "exact_evidence_count",
        ]
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    group_cols = ["pair_id", "compound_key"]

    for (pair_id, compound_key), group in exact_rows.groupby(group_cols, sort=True, dropna=False):
        ratios = _float_series(group["ratio_exact_human_div_parasite"]).dropna()
        ratios = ratios[ratios > 0]
        if ratios.empty:
            continue

        log_ratios = ratios.map(math.log10)
        median_log_ratio = float(log_ratios.median())
        median_ratio = float(10**median_log_ratio)
        labels = ratios >= selectivity_threshold
        label_conflict = bool(labels.nunique(dropna=True) > 1)
        is_selective = bool(median_ratio >= selectivity_threshold)

        canonical_smiles = _first_present(group["canonical_smiles"])
        murcko_scaffold = _first_present(group["murcko_scaffold"])
        standard_inchi_key = _first_present(group.get("standard_inchi_key", []))
        molecule_chembl_id = _first_present(group["molecule_chembl_id"])

        label_status = "exact_core_conflicting" if label_conflict else "exact_core_consistent"
        is_trainable = bool(
            not label_conflict
            and _is_present(canonical_smiles)
            and _is_present(murcko_scaffold)
            and _is_present(compound_key)
        )

        rows.append(
            {
                "exact_core_row_id": f"{pair_id}::{compound_key}",
                "pair_id": pair_id,
                "parasite_label": _first_present(group.get("parasite_label", [])),
                "human_label": _first_present(group.get("human_label", [])),
                "compound_key": compound_key,
                "molecule_chembl_id": molecule_chembl_id,
                "canonical_smiles": canonical_smiles,
                "standard_inchi_key": standard_inchi_key,
                "murcko_scaffold": murcko_scaffold,
                "selectivity_ratio_human_div_parasite": median_ratio,
                "log10_selectivity_ratio": median_log_ratio,
                "is_selective": int(is_selective),
                "selectivity_threshold": float(selectivity_threshold),
                "label_status": label_status,
                "is_trainable_exact_core": int(is_trainable),
                "exact_evidence_count": int(len(group)),
                "ratio_min_human_div_parasite": float(ratios.min()),
                "ratio_median_human_div_parasite": median_ratio,
                "ratio_max_human_div_parasite": float(ratios.max()),
                "log10_ratio_min": float(log_ratios.min()),
                "log10_ratio_median": median_log_ratio,
                "log10_ratio_max": float(log_ratios.max()),
                "n_unique_ratio_values": int(ratios.nunique(dropna=True)),
                "n_positive_exact_observations": int(labels.sum()),
                "n_negative_exact_observations": int((~labels).sum()),
                "parasite_standard_types": _joined_unique(group.get("parasite_standard_type", [])),
                "human_standard_types": _joined_unique(group.get("human_standard_type", [])),
                "parasite_assay_types": _joined_unique(group.get("parasite_assay_type", [])),
                "human_assay_types": _joined_unique(group.get("human_assay_type", [])),
                "parasite_activity_chembl_ids": _joined_unique(group.get("parasite_activity_chembl_id", [])),
                "human_activity_chembl_ids": _joined_unique(group.get("human_activity_chembl_id", [])),
                "parasite_assay_chembl_ids": _joined_unique(group.get("parasite_assay_chembl_id", [])),
                "human_assay_chembl_ids": _joined_unique(group.get("human_assay_chembl_id", [])),
                "parasite_document_chembl_ids": _joined_unique(group.get("parasite_document_chembl_id", [])),
                "human_document_chembl_ids": _joined_unique(group.get("human_document_chembl_id", [])),
                "n_parasite_documents": len({v for v in group.get("parasite_document_chembl_id", []) if _is_present(v)}),
                "n_human_documents": len({v for v in group.get("human_document_chembl_id", []) if _is_present(v)}),
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["pair_id", "compound_key"]).reset_index(drop=True)
    return out


def summarize_exact_core(
    candidate_rows: pd.DataFrame,
    exact_rows: pd.DataFrame,
    core_rows: pd.DataFrame,
    allowed_standard_types: Sequence[str],
    selectivity_threshold: float,
) -> dict[str, Any]:
    trainable = core_rows[core_rows["is_trainable_exact_core"].astype(bool)] if not core_rows.empty else core_rows
    conflicts = core_rows[core_rows["label_status"] == "exact_core_conflicting"] if not core_rows.empty else core_rows

    pair_summary: dict[str, Any] = {}
    pair_ids = sorted(set(candidate_rows.get("pair_id", [])) | set(core_rows.get("pair_id", [])))
    for pair_id in pair_ids:
        pair_core = core_rows[core_rows["pair_id"] == pair_id] if not core_rows.empty else core_rows
        pair_trainable = trainable[trainable["pair_id"] == pair_id] if not trainable.empty else trainable
        pair_summary[pair_id] = {
            "exact_candidate_rows": int((exact_rows["pair_id"] == pair_id).sum()) if not exact_rows.empty else 0,
            "exact_core_rows": int(len(pair_core)),
            "trainable_exact_core_rows": int(len(pair_trainable)),
            "conflicting_exact_core_rows": int((pair_core.get("label_status", pd.Series(dtype=str)) == "exact_core_conflicting").sum()),
            "trainable_positive_rows": int(pair_trainable.get("is_selective", pd.Series(dtype=int)).sum()) if not pair_trainable.empty else 0,
            "trainable_negative_rows": int(len(pair_trainable) - pair_trainable.get("is_selective", pd.Series(dtype=int)).sum()) if not pair_trainable.empty else 0,
            "trainable_unique_scaffolds": int(pair_trainable["murcko_scaffold"].nunique(dropna=True)) if not pair_trainable.empty else 0,
        }

    label_counts = Counter(core_rows["label_status"]) if not core_rows.empty else Counter()
    return {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "scientific_scope": (
            "Benchmark-ready exact-ratio v5 core collapsed to one row per target-pair/compound key; "
            "conflicting compound-pair labels are flagged and excluded from the trainable exact core."
        ),
        "input_candidate_rows": int(len(candidate_rows)),
        "strict_exact_candidate_rows_after_filter": int(len(exact_rows)),
        "exact_core_rows": int(len(core_rows)),
        "trainable_exact_core_rows": int(len(trainable)),
        "conflicting_exact_core_rows": int(len(conflicts)),
        "selectivity_threshold_human_div_parasite": float(selectivity_threshold),
        "allowed_standard_types": list(allowed_standard_types),
        "label_status_counts": dict(sorted(label_counts.items())),
        "trainable_rows_by_pair": dict(sorted(Counter(trainable["pair_id"]).items())) if not trainable.empty else {},
        "trainable_positive_rows_by_pair": dict(
            sorted(
                {
                    pair: int(rows["is_selective"].sum())
                    for pair, rows in trainable.groupby("pair_id", sort=True)
                }.items()
            )
        ) if not trainable.empty else {},
        "trainable_negative_rows_by_pair": dict(
            sorted(
                {
                    pair: int(len(rows) - rows["is_selective"].sum())
                    for pair, rows in trainable.groupby("pair_id", sort=True)
                }.items()
            )
        ) if not trainable.empty else {},
        "trainable_unique_scaffolds_by_pair": dict(
            sorted(
                {
                    pair: int(rows["murcko_scaffold"].nunique(dropna=True))
                    for pair, rows in trainable.groupby("pair_id", sort=True)
                }.items()
            )
        ) if not trainable.empty else {},
        "pair_summaries": pair_summary,
    }


def build_exact_core(
    candidate_csv: Path,
    output_csv: Path,
    output_summary_json: Path,
    allowed_standard_types: Sequence[str] = DEFAULT_ALLOWED_STANDARD_TYPES,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    candidate_rows = pd.read_csv(candidate_csv)
    exact_rows = filter_exact_core_candidates(candidate_rows, allowed_standard_types=allowed_standard_types)
    core_rows = aggregate_exact_core(exact_rows, selectivity_threshold=selectivity_threshold)
    summary = summarize_exact_core(
        candidate_rows=candidate_rows,
        exact_rows=exact_rows,
        core_rows=core_rows,
        allowed_standard_types=allowed_standard_types,
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
        default=Path("data/processed/selectivity_v5_exact_core_rows.csv"),
        help="Output exact-core CSV, one row per target-pair/compound key.",
    )
    parser.add_argument(
        "--output-summary-json",
        type=Path,
        default=Path("data/processed/selectivity_v5_exact_core_summary.json"),
        help="Output exact-core summary JSON.",
    )
    parser.add_argument(
        "--standard-type",
        dest="standard_types",
        action="append",
        help="Repeatable. Exact-core standard type filter. Default: IC50 only.",
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
    standard_types = tuple(value.upper() for value in (args.standard_types or DEFAULT_ALLOWED_STANDARD_TYPES))
    core_rows, summary = build_exact_core(
        candidate_csv=args.candidate_csv,
        output_csv=args.output_csv,
        output_summary_json=args.output_summary_json,
        allowed_standard_types=standard_types,
        selectivity_threshold=args.selectivity_threshold,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote:\n  - {args.output_csv}\n  - {args.output_summary_json}")
    print(f"\nExact-core rows: {len(core_rows)}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
