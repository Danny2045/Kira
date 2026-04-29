"""Design v6 lab-campaign tickets from the v5 selectivity evidence substrate.

This module ranks missing or weak parasite-vs-human comparator assays that would
make the selectivity evidence base more complete. It is an experiment
prioritization layer, not wet-lab validation and not a predictive model.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

DEFAULT_SELECTIVITY_THRESHOLD = 10.0
DEFAULT_TOP_N = 100
DEFAULT_CAMPAIGN_MODE = "benchmark_repair"
CAMPAIGN_MODES = ("benchmark_repair", "potency_discovery", "pair_specific")
DEFAULT_MAX_TICKETS_PER_PAIR = 25
DEFAULT_MIN_TICKETS_PER_PAIR = 5
UNKNOWN = "unknown_from_candidate_table"
SIDES = ("parasite", "human")
BOUNDED_STATUSES = {"lower_bound_ratio", "upper_bound_ratio", "interval_ratio"}
TICKET_FIELDS = [
    "ticket_id",
    "pair_id",
    "compound_key",
    "molecule_chembl_id",
    "canonical_smiles",
    "murcko_scaffold",
    "missing_side",
    "known_side",
    "known_target_label",
    "missing_target_label",
    "known_target_chembl_id",
    "missing_target_chembl_id",
    "recommended_assay_type",
    "recommended_standard_type",
    "recommended_units",
    "known_standard_type",
    "known_standard_relation",
    "known_standard_value_nM",
    "known_activity_chembl_id",
    "known_assay_chembl_id",
    "provenance_note",
    "evidence_status_source",
    "priority_score",
    "priority_reason",
    "expected_benchmark_impact",
    "data_return_schema",
]
SUMMARY_FIELDS = {
    "campaign_mode",
    "max_tickets_per_pair",
    "min_tickets_per_pair",
    "eligible_pairs_with_tickets",
    "input_candidate_rows",
    "tiered_core_rows",
    "generated_ticket_count",
    "tickets_by_pair",
    "tickets_by_missing_side",
    "top_priority_pairs",
    "target_imbalance_basis",
    "trainable_tiered_rows_by_pair",
    "trainable_positive_rows_by_pair",
    "trainable_negative_rows_by_pair",
    "pairs_with_both_classes",
    "scientific_scope",
    "top_10_ticket_ids",
    "top_10_priority_reasons",
    "tickets_with_known_assay_id",
    "tickets_with_known_activity_id",
    "tickets_missing_known_activity_id",
}
SCIENTIFIC_SCOPE = (
    "Kira v6 lab-campaign design ranks missing or weak parasite-vs-human paired "
    "selectivity measurements that would improve the evidence substrate. It is an "
    "experiment-prioritization layer, not wet-lab validation, not a predictive model, "
    "and not a new benchmark claim; v4 remains the current modeling claim."
)


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if value is pd.NA:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    if isinstance(value, str) and not value.strip():
        return False
    return True


def _clean_string(value: Any, default: str = UNKNOWN) -> str:
    if not _is_present(value):
        return default
    return str(value).strip()


def _has_known_provenance_id(value: Any) -> bool:
    return _is_present(value) and str(value).strip() != UNKNOWN


def _provenance_note(known_assay_id: Any, known_activity_id: Any) -> str:
    has_assay = _has_known_provenance_id(known_assay_id)
    has_activity = _has_known_provenance_id(known_activity_id)
    if has_assay and has_activity:
        return "known assay and activity provenance present"
    if has_assay:
        return "known assay provenance present; known activity id unavailable in v5 candidate table"
    if has_activity:
        return "known activity provenance present; known assay id unavailable in v5 candidate table"
    return "known assay and activity provenance unavailable in v5 candidate table"


def _float_or_none(value: Any) -> float | None:
    if not _is_present(value):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(numeric):
        return None
    return numeric


def _boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return False


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, (pd.Series, pd.Index)):
        return [_json_ready(item) for item in value.tolist()]
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, (pd.NA.__class__,)):
        return None
    if hasattr(value, "item"):
        try:
            return _json_ready(value.item())
        except ValueError:
            return None
    return value


def _side_col(side: str, suffix: str) -> str:
    return f"{side}_{suffix}"


def _other_side(side: str) -> str:
    if side == "parasite":
        return "human"
    if side == "human":
        return "parasite"
    raise ValueError(f"Unsupported side: {side!r}")


def _compound_key(row: pd.Series) -> str:
    inchi_key = row.get("standard_inchi_key")
    if _is_present(inchi_key):
        return str(inchi_key).strip()
    molecule_id = row.get("molecule_chembl_id")
    if _is_present(molecule_id):
        return str(molecule_id).strip()
    return UNKNOWN


def add_compound_keys(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    if "standard_inchi_key" not in out.columns:
        out["standard_inchi_key"] = None
    if "molecule_chembl_id" not in out.columns:
        out["molecule_chembl_id"] = None
    out["compound_key"] = [_compound_key(row) for _, row in out.iterrows()]
    return out


def _read_tiered_summary(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected tiered summary JSON object from {path}")
    return payload


def _counter_from_mapping(value: Any) -> dict[str, int]:
    if not isinstance(value, dict):
        return {}
    out: dict[str, int] = {}
    for key, item in value.items():
        numeric = _float_or_none(item)
        out[str(key)] = int(numeric or 0)
    return dict(sorted(out.items()))


def _trainable_core(tiered_core_rows: pd.DataFrame) -> pd.DataFrame:
    if tiered_core_rows.empty:
        return tiered_core_rows.copy()
    if "is_trainable_tiered_core" not in tiered_core_rows.columns:
        return tiered_core_rows.copy()
    return tiered_core_rows[tiered_core_rows["is_trainable_tiered_core"].map(_boolish)].copy()


def _derive_trainable_counts(
    tiered_core_rows: pd.DataFrame,
    tiered_summary: dict[str, Any],
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    rows_by_pair = _counter_from_mapping(tiered_summary.get("trainable_rows_by_pair"))
    positive_by_pair = _counter_from_mapping(tiered_summary.get("trainable_positive_rows_by_pair"))
    negative_by_pair = _counter_from_mapping(tiered_summary.get("trainable_negative_rows_by_pair"))
    if rows_by_pair and positive_by_pair:
        return rows_by_pair, positive_by_pair, negative_by_pair

    trainable = _trainable_core(tiered_core_rows)
    if trainable.empty or "pair_id" not in trainable.columns:
        return rows_by_pair, positive_by_pair, negative_by_pair

    rows_by_pair = dict(sorted(Counter(trainable["pair_id"].dropna()).items()))
    if "is_selective" in trainable.columns:
        for pair_id, group in trainable.groupby("pair_id", sort=True):
            labels = pd.to_numeric(group["is_selective"], errors="coerce")
            positive_by_pair[str(pair_id)] = int(labels.fillna(0).sum())
            negative_by_pair[str(pair_id)] = int(len(group) - labels.fillna(0).sum())
    return rows_by_pair, dict(sorted(positive_by_pair.items())), dict(sorted(negative_by_pair.items()))


def build_target_imbalance_basis(
    candidate_rows: pd.DataFrame,
    tiered_core_rows: pd.DataFrame,
    tiered_summary: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    rows_by_pair, positive_by_pair, negative_by_pair = _derive_trainable_counts(
        tiered_core_rows=tiered_core_rows,
        tiered_summary=tiered_summary,
    )
    pair_ids = set(rows_by_pair) | set(positive_by_pair) | set(negative_by_pair)
    if "pair_id" in candidate_rows.columns:
        pair_ids.update(str(pair_id) for pair_id in candidate_rows["pair_id"].dropna())
    if "pair_id" in tiered_core_rows.columns:
        pair_ids.update(str(pair_id) for pair_id in tiered_core_rows["pair_id"].dropna())

    basis: dict[str, dict[str, Any]] = {}
    for pair_id in sorted(pair_ids):
        rows = int(rows_by_pair.get(pair_id, 0))
        positives = int(positive_by_pair.get(pair_id, 0))
        negatives = int(negative_by_pair.get(pair_id, 0))
        total = rows or positives + negatives

        imbalance_ratio: float | None = None
        if positives > 0 and negatives > 0:
            imbalance_ratio = max(positives, negatives) / min(positives, negatives)

        if total == 0:
            class_status = "no_trainable_tiered_rows"
        elif positives == 0:
            class_status = "class_degenerate_zero_positive"
        elif negatives == 0:
            class_status = "class_degenerate_zero_negative"
        elif imbalance_ratio is not None and imbalance_ratio >= 10:
            class_status = "both_classes_extreme_imbalance"
        elif positives <= 2:
            class_status = "both_classes_few_positive_rows"
        else:
            class_status = "both_classes_present"

        basis[pair_id] = {
            "trainable_tiered_rows": total,
            "trainable_positive_rows": positives,
            "trainable_negative_rows": negatives,
            "class_status": class_status,
            "imbalance_ratio": imbalance_ratio,
        }
    return basis


def _trainable_compound_keys_by_pair(tiered_core_rows: pd.DataFrame) -> set[tuple[str, str]]:
    trainable = _trainable_core(add_compound_keys(tiered_core_rows))
    if trainable.empty or not {"pair_id", "compound_key"} <= set(trainable.columns):
        return set()
    return {
        (str(row["pair_id"]), str(row["compound_key"]))
        for _, row in trainable.iterrows()
        if _is_present(row.get("pair_id")) and _is_present(row.get("compound_key"))
    }


def _trainable_scaffolds_by_pair(tiered_core_rows: pd.DataFrame) -> dict[str, set[str]]:
    trainable = _trainable_core(tiered_core_rows)
    scaffolds: dict[str, set[str]] = {}
    if trainable.empty or not {"pair_id", "murcko_scaffold"} <= set(trainable.columns):
        return scaffolds
    for pair_id, group in trainable.groupby("pair_id", sort=True):
        scaffolds[str(pair_id)] = {
            str(value).strip()
            for value in group["murcko_scaffold"]
            if _is_present(value)
        }
    return scaffolds


def _exact_side(row: pd.Series, side: str) -> bool:
    relation = row.get(_side_col(side, "standard_relation"))
    interval_kind = row.get(_side_col(side, "measurement_interval_kind"))
    return relation == "=" or interval_kind == "exact"


def _ticket_seed(row: pd.Series) -> tuple[str, str] | None:
    status = row.get("candidate_status")
    if status == "single_side_only":
        known_side = str(row.get("single_side") or "").strip().lower()
        if known_side in SIDES:
            return known_side, _other_side(known_side)
        return None

    if status in BOUNDED_STATUSES:
        parasite_exact = _exact_side(row, "parasite")
        human_exact = _exact_side(row, "human")
        if parasite_exact and not human_exact:
            return "parasite", "human"
        if human_exact and not parasite_exact:
            return "human", "parasite"
    return None


def _target_label(row: pd.Series, side: str) -> str:
    return _clean_string(row.get(_side_col(side, "label")))


def _target_id(row: pd.Series, side: str) -> str:
    return _clean_string(row.get(_side_col(side, "target_chembl_id")))


def _recommended_value(row: pd.Series, known_side: str, missing_side: str, suffix: str) -> str:
    known_value = row.get(_side_col(known_side, suffix))
    if _is_present(known_value):
        return str(known_value).strip()
    missing_value = row.get(_side_col(missing_side, suffix))
    if _is_present(missing_value):
        return str(missing_value).strip()
    return UNKNOWN


def _data_return_schema(row: pd.Series, missing_side: str) -> str:
    schema = {
        "compound_key": _compound_key(row),
        "pair_id": _clean_string(row.get("pair_id")),
        "measured_side": missing_side,
        "target_chembl_id": _target_id(row, missing_side),
        "assay_type": "string",
        "standard_type": "string",
        "standard_relation": "string",
        "standard_value": "number",
        "standard_units": "nM",
        "replicate_count": "integer",
        "data_validity_comment": "string_or_null",
        "assay_chembl_id_or_external_id": "string",
        "activity_chembl_id_or_external_id": "string",
        "notes": "string",
    }
    return json.dumps(schema, sort_keys=True, separators=(",", ":"))


def _pair_score_components(pair_id: str, imbalance_basis: dict[str, dict[str, Any]]) -> tuple[float, list[str]]:
    basis = imbalance_basis.get(pair_id, {})
    total = int(basis.get("trainable_tiered_rows") or 0)
    positives = int(basis.get("trainable_positive_rows") or 0)
    negatives = int(basis.get("trainable_negative_rows") or 0)
    imbalance_ratio = _float_or_none(basis.get("imbalance_ratio"))
    class_status = str(basis.get("class_status") or "unknown")
    score = 0.0
    reasons: list[str] = []

    if total == 0:
        score += 22.0
        reasons.append("pair has no trainable tiered-core rows")
    elif positives == 0 or negatives == 0:
        score += 28.0
        reasons.append("pair is class-degenerate in the trainable tiered core")

    if positives == 0:
        score += 10.0
        reasons.append("pair has zero trainable positives")
    elif positives <= 2:
        score += 7.0
        reasons.append("pair has very few trainable positives")

    if imbalance_ratio is not None and imbalance_ratio >= 10:
        score += 16.0
        reasons.append("pair has extreme class imbalance")
    elif imbalance_ratio is not None and imbalance_ratio >= 5:
        score += 9.0
        reasons.append("pair has strong class imbalance")

    if total and total < 5:
        score += 6.0
        reasons.append("pair has very few trainable tiered-core rows")

    if class_status != "both_classes_present":
        score += 5.0
        reasons.append("missing assay could help repair benchmark imbalance")

    return score, reasons


def _potency_score(value_nM: float | None) -> tuple[float, str]:
    if value_nM is None:
        return -12.0, "known-side potency is missing"
    if value_nM <= 10:
        return 18.0, "known side has very strong potency"
    if value_nM <= 50:
        return 14.0, "known side has strong potency"
    if value_nM <= 100:
        return 10.0, "known side has sub-100 nM potency"
    if value_nM <= 1000:
        return 5.0, "known side has measurable potency"
    if value_nM <= 5000:
        return -4.0, "known-side potency is weak"
    return -10.0, "known-side potency is very weak"


def _bounded_threshold_score(row: pd.Series, threshold: float) -> tuple[float, list[str]]:
    status = row.get("candidate_status")
    if status not in BOUNDED_STATUSES:
        return 0.0, []

    lower = _float_or_none(row.get("ratio_lower_bound_human_div_parasite"))
    upper = _float_or_none(row.get("ratio_upper_bound_human_div_parasite"))
    score = 0.0
    reasons: list[str] = []

    if lower is not None and upper is not None and lower < threshold <= upper:
        score += 18.0
        reasons.append("bounded interval crosses the selectivity threshold")
    elif lower is not None and lower >= threshold:
        score += 14.0
        reasons.append("bounded evidence is already at or above the threshold")
    elif lower is not None and lower >= threshold * 0.5:
        score += 8.0
        reasons.append("bounded evidence is near the threshold")
    elif upper is not None and threshold <= upper <= threshold * 2:
        score += 8.0
        reasons.append("upper-bound evidence is near the threshold")
    elif upper is not None and upper < threshold:
        score += 4.0
        reasons.append("upper-bound evidence could confirm a negative paired ratio")
    return score, reasons


def _metadata_penalty(ticket: dict[str, Any]) -> tuple[float, list[str]]:
    score = 0.0
    reasons: list[str] = []
    important_fields = [
        "known_target_label",
        "missing_target_label",
        "known_target_chembl_id",
        "missing_target_chembl_id",
        "recommended_assay_type",
        "recommended_standard_type",
        "known_standard_type",
        "known_standard_relation",
    ]
    missing = [field for field in important_fields if ticket.get(field) == UNKNOWN]
    if missing:
        score -= min(24.0, 4.0 * len(missing))
        reasons.append(f"metadata placeholders reduce interpretability: {','.join(missing)}")
    if not _is_present(ticket.get("known_assay_chembl_id")) or ticket["known_assay_chembl_id"] == UNKNOWN:
        score -= 3.0
        reasons.append("known assay id is missing")
    return score, reasons


def _expected_benchmark_impact(pair_id: str, missing_side: str, imbalance_basis: dict[str, dict[str, Any]]) -> str:
    basis = imbalance_basis.get(pair_id, {})
    status = basis.get("class_status")
    positives = int(basis.get("trainable_positive_rows") or 0)
    negatives = int(basis.get("trainable_negative_rows") or 0)

    if status in {"class_degenerate_zero_positive", "class_degenerate_zero_negative"}:
        return (
            f"Could help close a class-degenerate tiered-core gap for {pair_id} "
            f"(positives={positives}, negatives={negatives}) by measuring the {missing_side} side."
        )
    if status == "both_classes_extreme_imbalance":
        return (
            f"Could add paired evidence to an extremely imbalanced {pair_id} tiered core "
            f"(positives={positives}, negatives={negatives})."
        )
    if status == "no_trainable_tiered_rows":
        return f"Could create an initial paired-evidence opportunity for {pair_id}."
    if positives <= 2:
        return (
            f"Could improve sparse positive coverage for {pair_id} "
            f"(positives={positives}, negatives={negatives})."
        )
    return f"Could add matched comparator evidence and scaffold coverage for {pair_id}."


def _score_ticket(
    row: pd.Series,
    ticket: dict[str, Any],
    trainable_keys: set[tuple[str, str]],
    trainable_scaffolds: dict[str, set[str]],
    imbalance_basis: dict[str, dict[str, Any]],
    selectivity_threshold: float,
) -> tuple[float, list[str]]:
    pair_id = ticket["pair_id"]
    compound_key = ticket["compound_key"]
    known_side = ticket["known_side"]
    status = row.get("candidate_status")
    score = 0.0
    reasons: list[str] = []

    if status == "single_side_only":
        if _exact_side(row, known_side):
            score += 40.0
            reasons.append("single-side row can become an exact matched ratio")
        else:
            score += 28.0
            reasons.append("single-side row can become paired comparator evidence")
    elif status in BOUNDED_STATUSES:
        score += 22.0
        reasons.append("bounded row can be resolved by a comparator-side exact measurement")

    pair_score, pair_reasons = _pair_score_components(pair_id, imbalance_basis)
    score += pair_score
    reasons.extend(pair_reasons)

    if (pair_id, compound_key) in trainable_keys:
        score -= 30.0
        reasons.append("compound is already represented in the trainable tiered core for this pair")
    else:
        score += 8.0
        reasons.append("compound is not represented in the trainable tiered core for this pair")

    if ticket["canonical_smiles"] == UNKNOWN:
        score -= 16.0
        reasons.append("canonical SMILES is missing")
    else:
        score += 3.0
        reasons.append("compound has canonical SMILES")

    scaffold = ticket["murcko_scaffold"]
    if scaffold == UNKNOWN:
        score -= 12.0
        reasons.append("Murcko scaffold is missing")
    elif scaffold not in trainable_scaffolds.get(pair_id, set()):
        score += 8.0
        reasons.append("scaffold is not represented in the trainable tiered core for this pair")
    else:
        score -= 3.0
        reasons.append("scaffold is already represented in the trainable tiered core for this pair")

    potency_score, potency_reason = _potency_score(_float_or_none(row.get(_side_col(known_side, "standard_value_nM"))))
    score += potency_score
    reasons.append(potency_reason)

    threshold_score, threshold_reasons = _bounded_threshold_score(row, selectivity_threshold)
    score += threshold_score
    reasons.extend(threshold_reasons)

    metadata_score, metadata_reasons = _metadata_penalty(ticket)
    score += metadata_score
    reasons.extend(metadata_reasons)

    return score, reasons


def _make_ticket(
    row: pd.Series,
    known_side: str,
    missing_side: str,
    trainable_keys: set[tuple[str, str]],
    trainable_scaffolds: dict[str, set[str]],
    imbalance_basis: dict[str, dict[str, Any]],
    selectivity_threshold: float,
) -> dict[str, Any]:
    pair_id = _clean_string(row.get("pair_id"))
    compound_key = _compound_key(row)
    known_standard_value = _float_or_none(row.get(_side_col(known_side, "standard_value_nM")))
    known_activity_id = _clean_string(row.get(_side_col(known_side, "activity_chembl_id")))
    known_assay_id = _clean_string(row.get(_side_col(known_side, "assay_chembl_id")))

    ticket: dict[str, Any] = {
        "ticket_id": "",
        "pair_id": pair_id,
        "compound_key": compound_key,
        "molecule_chembl_id": _clean_string(row.get("molecule_chembl_id")),
        "canonical_smiles": _clean_string(row.get("canonical_smiles")),
        "murcko_scaffold": _clean_string(row.get("murcko_scaffold")),
        "missing_side": missing_side,
        "known_side": known_side,
        "known_target_label": _target_label(row, known_side),
        "missing_target_label": _target_label(row, missing_side),
        "known_target_chembl_id": _target_id(row, known_side),
        "missing_target_chembl_id": _target_id(row, missing_side),
        "recommended_assay_type": _recommended_value(row, known_side, missing_side, "assay_type"),
        "recommended_standard_type": _recommended_value(row, known_side, missing_side, "standard_type"),
        "recommended_units": "nM",
        "known_standard_type": _clean_string(row.get(_side_col(known_side, "standard_type"))),
        "known_standard_relation": _clean_string(row.get(_side_col(known_side, "standard_relation"))),
        "known_standard_value_nM": known_standard_value,
        "known_activity_chembl_id": known_activity_id,
        "known_assay_chembl_id": known_assay_id,
        "provenance_note": _provenance_note(known_assay_id, known_activity_id),
        "evidence_status_source": _clean_string(row.get("candidate_status")),
        "priority_score": 0.0,
        "priority_reason": "",
        "expected_benchmark_impact": _expected_benchmark_impact(pair_id, missing_side, imbalance_basis),
        "data_return_schema": _data_return_schema(row, missing_side),
    }
    score, reasons = _score_ticket(
        row=row,
        ticket=ticket,
        trainable_keys=trainable_keys,
        trainable_scaffolds=trainable_scaffolds,
        imbalance_basis=imbalance_basis,
        selectivity_threshold=selectivity_threshold,
    )
    ticket["priority_score"] = round(score, 3)
    ticket["priority_reason"] = "; ".join(dict.fromkeys(reasons))
    return ticket


def _apply_scaffold_diversity(tickets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(
        tickets,
        key=lambda item: (
            -float(item["priority_score"]),
            item["pair_id"],
            item["murcko_scaffold"],
            item["compound_key"],
            item["missing_side"],
        ),
    )
    scaffold_counts: Counter[tuple[str, str]] = Counter()
    adjusted: list[dict[str, Any]] = []
    for ticket in ordered:
        item = dict(ticket)
        scaffold = item["murcko_scaffold"]
        if scaffold != UNKNOWN:
            key = (item["pair_id"], scaffold)
            scaffold_counts[key] += 1
            rank = scaffold_counts[key]
            if rank == 1:
                item["priority_score"] = round(float(item["priority_score"]) + 5.0, 3)
                item["priority_reason"] = f"{item['priority_reason']}; first ticket for this pair/scaffold in the ranked campaign"
            else:
                penalty = min(16.0, 4.0 * (rank - 1))
                item["priority_score"] = round(float(item["priority_score"]) - penalty, 3)
                item["priority_reason"] = (
                    f"{item['priority_reason']}; repeated pair/scaffold ticket penalized for campaign diversity"
                )
        adjusted.append(item)
    return adjusted


def _sort_tickets(tickets: pd.DataFrame) -> pd.DataFrame:
    if tickets.empty:
        return pd.DataFrame(columns=TICKET_FIELDS)
    return tickets.sort_values(
        ["priority_score", "pair_id", "compound_key"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)


def _assign_ticket_ids(tickets: pd.DataFrame) -> pd.DataFrame:
    out = tickets.copy().reset_index(drop=True)
    if out.empty:
        return pd.DataFrame(columns=TICKET_FIELDS)
    out["ticket_id"] = [f"V6-{index:04d}" for index in range(1, len(out) + 1)]
    return out[TICKET_FIELDS]


def _limit_top_n(tickets: pd.DataFrame, top_n: int | None) -> pd.DataFrame:
    if top_n is None:
        return tickets.copy()
    return tickets.head(max(0, int(top_n))).copy()


def _validate_campaign_args(
    campaign_mode: str,
    pair_id: str | None,
    max_tickets_per_pair: int,
    min_tickets_per_pair: int,
) -> None:
    if campaign_mode not in CAMPAIGN_MODES:
        raise ValueError(f"Unsupported campaign mode {campaign_mode!r}; expected one of {CAMPAIGN_MODES}")
    if campaign_mode == "pair_specific" and not _is_present(pair_id):
        raise ValueError("--pair-id is required when --campaign-mode pair_specific")
    if max_tickets_per_pair < 1:
        raise ValueError("max_tickets_per_pair must be at least 1")
    if min_tickets_per_pair < 0:
        raise ValueError("min_tickets_per_pair must be non-negative")


def _select_benchmark_repair_tickets(
    ranked_tickets: pd.DataFrame,
    top_n: int | None,
    max_tickets_per_pair: int,
    min_tickets_per_pair: int,
) -> pd.DataFrame:
    if ranked_tickets.empty:
        return pd.DataFrame(columns=TICKET_FIELDS)
    if top_n is not None and top_n <= 0:
        return pd.DataFrame(columns=TICKET_FIELDS)

    capped_min = min(min_tickets_per_pair, max_tickets_per_pair)
    pair_order = sorted(
        ranked_tickets["pair_id"].unique(),
        key=lambda pair: (
            -float(ranked_tickets.loc[ranked_tickets["pair_id"] == pair, "priority_score"].max()),
            str(pair),
        ),
    )
    grouped = {
        str(pair): ranked_tickets[ranked_tickets["pair_id"] == pair].reset_index(drop=False)
        for pair in pair_order
    }
    selected_indices: list[int] = []
    selected_index_set: set[int] = set()
    counts_by_pair: Counter[str] = Counter()

    def can_add(pair: str, original_index: int) -> bool:
        if top_n is not None and len(selected_indices) >= top_n:
            return False
        return counts_by_pair[pair] < max_tickets_per_pair and original_index not in selected_index_set

    def add_index(pair: str, original_index: int) -> None:
        selected_indices.append(original_index)
        selected_index_set.add(original_index)
        counts_by_pair[pair] += 1

    for rank in range(capped_min):
        for pair in pair_order:
            rows = grouped[pair]
            if len(rows) <= rank:
                continue
            original_index = int(rows.iloc[rank]["index"])
            if can_add(pair, original_index):
                add_index(pair, original_index)

    for original_index, ticket in ranked_tickets.iterrows():
        pair = str(ticket["pair_id"])
        if can_add(pair, int(original_index)):
            add_index(pair, int(original_index))

    selected = ranked_tickets.loc[selected_indices].copy()
    return _sort_tickets(selected)


def _select_campaign_tickets(
    ranked_tickets: pd.DataFrame,
    top_n: int | None,
    campaign_mode: str,
    pair_id: str | None,
    max_tickets_per_pair: int,
    min_tickets_per_pair: int,
) -> pd.DataFrame:
    _validate_campaign_args(
        campaign_mode=campaign_mode,
        pair_id=pair_id,
        max_tickets_per_pair=max_tickets_per_pair,
        min_tickets_per_pair=min_tickets_per_pair,
    )
    if ranked_tickets.empty:
        return pd.DataFrame(columns=TICKET_FIELDS)

    if campaign_mode == "potency_discovery":
        return _limit_top_n(ranked_tickets, top_n)
    if campaign_mode == "pair_specific":
        requested = ranked_tickets[ranked_tickets["pair_id"] == str(pair_id)].copy()
        return _limit_top_n(requested, top_n)
    return _select_benchmark_repair_tickets(
        ranked_tickets=ranked_tickets,
        top_n=top_n,
        max_tickets_per_pair=max_tickets_per_pair,
        min_tickets_per_pair=min_tickets_per_pair,
    )


def generate_lab_tickets(
    candidate_rows: pd.DataFrame,
    tiered_core_rows: pd.DataFrame,
    tiered_summary: dict[str, Any],
    top_n: int | None = DEFAULT_TOP_N,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
    campaign_mode: str = DEFAULT_CAMPAIGN_MODE,
    pair_id: str | None = None,
    max_tickets_per_pair: int = DEFAULT_MAX_TICKETS_PER_PAIR,
    min_tickets_per_pair: int = DEFAULT_MIN_TICKETS_PER_PAIR,
) -> pd.DataFrame:
    required = {
        "pair_id",
        "molecule_chembl_id",
        "candidate_status",
        "single_side",
        "canonical_smiles",
        "murcko_scaffold",
    }
    missing = sorted(required - set(candidate_rows.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {missing}")

    candidates = add_compound_keys(candidate_rows)
    tiered = add_compound_keys(tiered_core_rows)
    trainable_keys = _trainable_compound_keys_by_pair(tiered)
    trainable_scaffolds = _trainable_scaffolds_by_pair(tiered)
    imbalance_basis = build_target_imbalance_basis(candidates, tiered, tiered_summary)

    raw_tickets: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        seed = _ticket_seed(row)
        if seed is None:
            continue
        known_side, missing_side = seed
        raw_tickets.append(
            _make_ticket(
                row=row,
                known_side=known_side,
                missing_side=missing_side,
                trainable_keys=trainable_keys,
                trainable_scaffolds=trainable_scaffolds,
                imbalance_basis=imbalance_basis,
                selectivity_threshold=selectivity_threshold,
            )
        )

    if not raw_tickets:
        return pd.DataFrame(columns=TICKET_FIELDS)

    deduped = (
        pd.DataFrame(raw_tickets)
        .sort_values(
            ["priority_score", "pair_id", "compound_key", "missing_side", "known_activity_chembl_id"],
            ascending=[False, True, True, True, True],
            kind="mergesort",
        )
        .drop_duplicates(["pair_id", "compound_key", "missing_side"], keep="first")
    )
    adjusted = _apply_scaffold_diversity(deduped.to_dict(orient="records"))
    ranked_tickets = _sort_tickets(pd.DataFrame(adjusted))
    eligible_pairs = sorted(ranked_tickets["pair_id"].unique()) if not ranked_tickets.empty else []
    tickets = _select_campaign_tickets(
        ranked_tickets=ranked_tickets,
        top_n=top_n,
        campaign_mode=campaign_mode,
        pair_id=pair_id,
        max_tickets_per_pair=max_tickets_per_pair,
        min_tickets_per_pair=min_tickets_per_pair,
    )
    tickets = _assign_ticket_ids(tickets)
    tickets.attrs["eligible_pairs_with_tickets"] = eligible_pairs
    tickets.attrs["campaign_mode"] = campaign_mode
    tickets.attrs["max_tickets_per_pair"] = int(max_tickets_per_pair)
    tickets.attrs["min_tickets_per_pair"] = int(min_tickets_per_pair)
    return tickets


def summarize_lab_campaign(
    candidate_rows: pd.DataFrame,
    tiered_core_rows: pd.DataFrame,
    tiered_summary: dict[str, Any],
    tickets: pd.DataFrame,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
    campaign_mode: str = DEFAULT_CAMPAIGN_MODE,
    max_tickets_per_pair: int = DEFAULT_MAX_TICKETS_PER_PAIR,
    min_tickets_per_pair: int = DEFAULT_MIN_TICKETS_PER_PAIR,
    eligible_pairs_with_tickets: Sequence[str] | None = None,
) -> dict[str, Any]:
    rows_by_pair, positive_by_pair, negative_by_pair = _derive_trainable_counts(
        tiered_core_rows=tiered_core_rows,
        tiered_summary=tiered_summary,
    )
    imbalance_basis = build_target_imbalance_basis(candidate_rows, tiered_core_rows, tiered_summary)
    tickets_by_pair = dict(sorted(Counter(tickets["pair_id"]).items())) if not tickets.empty else {}
    tickets_by_missing_side = dict(sorted(Counter(tickets["missing_side"]).items())) if not tickets.empty else {}
    tickets_with_known_assay_id = (
        int(tickets["known_assay_chembl_id"].map(_has_known_provenance_id).sum()) if not tickets.empty else 0
    )
    tickets_with_known_activity_id = (
        int(tickets["known_activity_chembl_id"].map(_has_known_provenance_id).sum()) if not tickets.empty else 0
    )
    tickets_missing_known_activity_id = int(len(tickets) - tickets_with_known_activity_id)
    pairs_with_both_classes = sorted(tiered_summary.get("pairs_with_both_classes", []))
    eligible_pairs = (
        sorted(str(pair) for pair in eligible_pairs_with_tickets)
        if eligible_pairs_with_tickets is not None
        else sorted(str(pair) for pair in tickets.attrs.get("eligible_pairs_with_tickets", tickets_by_pair))
    )

    top_priority_pairs: list[dict[str, Any]] = []
    if not tickets.empty:
        for pair_id, group in tickets.groupby("pair_id", sort=True):
            top_priority_pairs.append(
                {
                    "pair_id": str(pair_id),
                    "ticket_count": int(len(group)),
                    "max_priority_score": float(group["priority_score"].max()),
                    "mean_priority_score": round(float(group["priority_score"].mean()), 3),
                }
            )
        top_priority_pairs = sorted(
            top_priority_pairs,
            key=lambda item: (-item["max_priority_score"], -item["ticket_count"], item["pair_id"]),
        )[:10]

    return {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "campaign_mode": campaign_mode,
        "max_tickets_per_pair": int(max_tickets_per_pair),
        "min_tickets_per_pair": int(min_tickets_per_pair),
        "eligible_pairs_with_tickets": eligible_pairs,
        "scientific_scope": SCIENTIFIC_SCOPE,
        "selectivity_threshold_human_div_parasite": float(selectivity_threshold),
        "input_candidate_rows": int(len(candidate_rows)),
        "tiered_core_rows": int(len(tiered_core_rows)),
        "generated_ticket_count": int(len(tickets)),
        "tickets_with_known_assay_id": tickets_with_known_assay_id,
        "tickets_with_known_activity_id": tickets_with_known_activity_id,
        "tickets_missing_known_activity_id": tickets_missing_known_activity_id,
        "tickets_by_pair": tickets_by_pair,
        "tickets_by_missing_side": tickets_by_missing_side,
        "top_priority_pairs": top_priority_pairs,
        "target_imbalance_basis": imbalance_basis,
        "trainable_tiered_rows_by_pair": dict(sorted(rows_by_pair.items())),
        "trainable_positive_rows_by_pair": dict(sorted(positive_by_pair.items())),
        "trainable_negative_rows_by_pair": dict(sorted(negative_by_pair.items())),
        "pairs_with_both_classes": pairs_with_both_classes,
        "top_10_ticket_ids": tickets["ticket_id"].head(10).tolist() if not tickets.empty else [],
        "top_10_priority_reasons": tickets["priority_reason"].head(10).tolist() if not tickets.empty else [],
    }


def _campaign_payload(
    tickets: pd.DataFrame,
    summary: dict[str, Any],
    top_n: int,
    output_csv: Path,
    output_summary_json: Path,
) -> dict[str, Any]:
    return {
        "generated_at_utc": summary["generated_at_utc"],
        "scientific_scope": SCIENTIFIC_SCOPE,
        "selectivity_threshold_human_div_parasite": summary["selectivity_threshold_human_div_parasite"],
        "campaign_mode": summary["campaign_mode"],
        "max_tickets_per_pair": summary["max_tickets_per_pair"],
        "min_tickets_per_pair": summary["min_tickets_per_pair"],
        "eligible_pairs_with_tickets": summary["eligible_pairs_with_tickets"],
        "top_n": int(top_n),
        "ticket_count": int(len(tickets)),
        "ticket_csv": str(output_csv),
        "summary_json": str(output_summary_json),
        "required_return_fields": [
            "compound_key",
            "pair_id",
            "measured_side",
            "target_chembl_id",
            "assay_type",
            "standard_type",
            "standard_relation",
            "standard_value",
            "standard_units",
            "replicate_count",
            "data_validity_comment",
            "assay_chembl_id_or_external_id",
            "activity_chembl_id_or_external_id",
            "notes",
        ],
        "tickets": tickets.to_dict(orient="records"),
    }


def build_lab_campaign(
    candidate_csv: Path,
    tiered_core_csv: Path,
    tiered_summary_json: Path,
    output_csv: Path,
    output_campaign_json: Path,
    output_summary_json: Path,
    top_n: int = DEFAULT_TOP_N,
    selectivity_threshold: float = DEFAULT_SELECTIVITY_THRESHOLD,
    campaign_mode: str = DEFAULT_CAMPAIGN_MODE,
    pair_id: str | None = None,
    max_tickets_per_pair: int = DEFAULT_MAX_TICKETS_PER_PAIR,
    min_tickets_per_pair: int = DEFAULT_MIN_TICKETS_PER_PAIR,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    candidate_rows = pd.read_csv(candidate_csv)
    tiered_core_rows = pd.read_csv(tiered_core_csv)
    tiered_summary = _read_tiered_summary(tiered_summary_json)
    tickets = generate_lab_tickets(
        candidate_rows=candidate_rows,
        tiered_core_rows=tiered_core_rows,
        tiered_summary=tiered_summary,
        top_n=top_n,
        selectivity_threshold=selectivity_threshold,
        campaign_mode=campaign_mode,
        pair_id=pair_id,
        max_tickets_per_pair=max_tickets_per_pair,
        min_tickets_per_pair=min_tickets_per_pair,
    )
    summary = summarize_lab_campaign(
        candidate_rows=candidate_rows,
        tiered_core_rows=tiered_core_rows,
        tiered_summary=tiered_summary,
        tickets=tickets,
        selectivity_threshold=selectivity_threshold,
        campaign_mode=campaign_mode,
        max_tickets_per_pair=max_tickets_per_pair,
        min_tickets_per_pair=min_tickets_per_pair,
        eligible_pairs_with_tickets=tickets.attrs.get("eligible_pairs_with_tickets"),
    )
    campaign = _campaign_payload(
        tickets=tickets,
        summary=summary,
        top_n=top_n,
        output_csv=output_csv,
        output_summary_json=output_summary_json,
    )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_campaign_json.parent.mkdir(parents=True, exist_ok=True)
    output_summary_json.parent.mkdir(parents=True, exist_ok=True)
    tickets.to_csv(output_csv, index=False)
    with output_campaign_json.open("w", encoding="utf-8") as handle:
        json.dump(_json_ready(campaign), handle, indent=2, sort_keys=True)
        handle.write("\n")
    with output_summary_json.open("w", encoding="utf-8") as handle:
        json.dump(_json_ready(summary), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return tickets, campaign, summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-csv",
        type=Path,
        default=Path("data/processed/selectivity_v5_candidate_rows.csv"),
        help="Input v5 candidate rows CSV.",
    )
    parser.add_argument(
        "--tiered-core-csv",
        type=Path,
        default=Path("data/processed/selectivity_v5_tiered_core_rows.csv"),
        help="Input v5 tiered-core rows CSV.",
    )
    parser.add_argument(
        "--tiered-summary-json",
        type=Path,
        default=Path("data/processed/selectivity_v5_tiered_core_summary.json"),
        help="Input v5 tiered-core summary JSON.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/lab_requests/v6_top_assay_tickets.csv"),
        help="Output top assay-ticket CSV.",
    )
    parser.add_argument(
        "--output-campaign-json",
        type=Path,
        default=Path("data/lab_requests/v6_gap_closure_campaign.json"),
        help="Output campaign JSON.",
    )
    parser.add_argument(
        "--output-summary-json",
        type=Path,
        default=Path("data/processed/selectivity_v6_campaign_summary.json"),
        help="Output campaign summary JSON.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_TOP_N,
        help="Maximum number of tickets to write. Default: 100.",
    )
    parser.add_argument(
        "--campaign-mode",
        choices=CAMPAIGN_MODES,
        default=DEFAULT_CAMPAIGN_MODE,
        help="Campaign portfolio mode. Default: benchmark_repair.",
    )
    parser.add_argument(
        "--pair-id",
        default=None,
        help="Target pair to emit when --campaign-mode pair_specific.",
    )
    parser.add_argument(
        "--max-tickets-per-pair",
        type=int,
        default=DEFAULT_MAX_TICKETS_PER_PAIR,
        help="Benchmark-repair cap for tickets from one target pair. Default: 25.",
    )
    parser.add_argument(
        "--min-tickets-per-pair",
        type=int,
        default=DEFAULT_MIN_TICKETS_PER_PAIR,
        help="Benchmark-repair target floor for eligible target pairs. Default: 5.",
    )
    parser.add_argument(
        "--selectivity-threshold",
        type=float,
        default=DEFAULT_SELECTIVITY_THRESHOLD,
        help="Human/parasite ratio threshold used for selectivity context. Default: 10.",
    )
    args = parser.parse_args(argv)
    if args.campaign_mode == "pair_specific" and not _is_present(args.pair_id):
        parser.error("--pair-id is required when --campaign-mode pair_specific")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    tickets, _campaign, summary = build_lab_campaign(
        candidate_csv=args.candidate_csv,
        tiered_core_csv=args.tiered_core_csv,
        tiered_summary_json=args.tiered_summary_json,
        output_csv=args.output_csv,
        output_campaign_json=args.output_campaign_json,
        output_summary_json=args.output_summary_json,
        top_n=args.top_n,
        selectivity_threshold=args.selectivity_threshold,
        campaign_mode=args.campaign_mode,
        pair_id=args.pair_id,
        max_tickets_per_pair=args.max_tickets_per_pair,
        min_tickets_per_pair=args.min_tickets_per_pair,
    )
    print(json.dumps(_json_ready(summary), indent=2, sort_keys=True))
    print(f"\nWrote:\n  - {args.output_csv}\n  - {args.output_campaign_json}\n  - {args.output_summary_json}")
    print(f"\nAssay tickets: {len(tickets)}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
