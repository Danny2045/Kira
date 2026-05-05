"""Parasite selectivity v4/v5/v6 benchmark-repair report.

The report is intentionally a reader and renderer over existing local artifacts.
It does not rerun models, regenerate selectivity data, or modify historical
result files.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any

PAIR_ORDER = ("LmDHFR", "LmPTR1", "SmDHODH", "SmHDAC8", "TbCathB", "TbPDEB1")
SIDES = ("parasite", "human")

SOURCE_PATHS = {
    "v4_model_summary": Path("results/selectivity_v4/summary.json"),
    "v4_pair_metrics": Path("results/selectivity_v4/per_pair_metrics.csv"),
    "v4_data_summary": Path("data/processed/selectivity_v4_summary.json"),
    "v5_expansion_summary": Path("data/processed/selectivity_v5_expansion_summary.json"),
    "v5_exact_core_summary": Path("data/processed/selectivity_v5_exact_core_summary.json"),
    "v5_tiered_core_summary": Path("data/processed/selectivity_v5_tiered_core_summary.json"),
    "v6_campaign_summary": Path("data/processed/selectivity_v6_campaign_summary.json"),
    "v6_top_assay_tickets": Path("data/lab_requests/v6_top_assay_tickets.csv"),
    "v6_gap_closure_campaign": Path("data/lab_requests/v6_gap_closure_campaign.json"),
    "target_pairs": Path("data/reference/selectivity_v5_target_pairs.csv"),
    "target_manifest_report": Path("data/processed/canonical_target_manifest_report.txt"),
}

REQUIRED_JSON_SOURCES = (
    "v4_model_summary",
    "v5_expansion_summary",
    "v5_exact_core_summary",
    "v5_tiered_core_summary",
    "v6_campaign_summary",
)

NON_CLAIMS = [
    "No drug-discovery claim: this dossier summarizes evidence substrate and repair tickets only.",
    "No wet-lab validation claim: v6 tickets are measurement requests, not completed experiments.",
    "No fresh model-performance claim: the report repeats existing v4/v5/v6 artifacts and does not rerun a benchmark.",
    "No clinical or therapeutic-success claim.",
    "No claim that parasite disease problems are broadly resolved.",
]

VALIDATION_COMMANDS = [
    "ruff check src/kira/selectivity tests/test_selectivity_benchmark_repair_report.py",
    "pytest -q tests/test_selectivity_benchmark_repair_report.py",
    "ruff check .",
    "pytest -q",
]


@dataclass(frozen=True)
class SelectivityRepairInputs:
    """Local artifacts used to build the benchmark-repair report."""

    source_files: list[str]
    v4_model_summary: dict[str, Any]
    v5_expansion_summary: dict[str, Any]
    v5_exact_core_summary: dict[str, Any]
    v5_tiered_core_summary: dict[str, Any]
    v6_campaign_summary: dict[str, Any]
    v6_ticket_rows: list[dict[str, str]]
    v6_gap_closure_campaign: dict[str, Any]
    target_pair_rows: list[dict[str, str]]


@dataclass(frozen=True)
class V4ModelingClaimSummary:
    """Existing v4 modeling claim distilled from the committed summary."""

    current_modeling_claim: str
    trainable_rows: int
    feature_count: int
    n_splits: int
    n_bits: int
    ablation_macro_pair_aurocs: dict[str, float | None]
    pair_only_macro_pair_auroc: float | None
    compound_only_macro_pair_auroc: float | None
    compound_plus_pair_macro_pair_auroc: float | None
    best_macro_pair_auroc_ablation: str | None
    best_macro_pair_auroc: float | None
    compound_plus_pair_is_best_by_macro_pair_auroc: bool
    interpretation: str


@dataclass(frozen=True)
class V5EvidenceSubstrateSummary:
    """v5 public-data evidence substrate summary."""

    candidate_evidence_rows: int
    curated_activity_rows: int
    exact_matched_ratio_candidate_rows: int
    exact_core_rows: int
    trainable_exact_core_rows: int
    tiered_core_rows: int
    trainable_tiered_core_rows: int
    tiered_pairs_with_both_classes: list[str]
    exact_pairs_with_both_classes: list[str]
    activity_rows_by_side: dict[str, int]
    candidate_rows_by_pair: dict[str, int]


@dataclass(frozen=True)
class V6LabCampaignSummary:
    """v6 benchmark-repair lab-campaign summary."""

    campaign_mode: str
    benchmark_repair_ticket_count: int
    tickets_by_pair: dict[str, int]
    tickets_by_missing_side: dict[str, int]
    eligible_pairs_with_tickets: list[str]
    required_return_fields: list[str]


@dataclass(frozen=True)
class TargetPairRepairSummary:
    """Per-target-pair benchmark-repair interpretation."""

    pair_id: str
    parasite_label: str
    human_label: str
    evidence_state: str
    exact_core_both_classes: bool
    tiered_core_both_classes: bool
    exact_trainable_rows: int
    exact_positive_rows: int
    exact_negative_rows: int
    tiered_trainable_rows: int
    tiered_positive_rows: int
    tiered_negative_rows: int
    missing_side_counts: dict[str, int]
    missing_side_summary: str
    ticket_count: int
    repair_priority: str
    explanation: str


@dataclass(frozen=True)
class NextActions:
    """Actionable benchmark-repair conclusions."""

    closest_to_benchmark_ready: list[str]
    parasite_side_measurement_pairs: list[str]
    human_side_comparator_pairs: list[str]
    returned_data_that_repairs_benchmark: list[str]


@dataclass(frozen=True)
class SelectivityRepairReport:
    """Full parasite selectivity benchmark-repair report."""

    title: str
    purpose: str
    source_files: list[str]
    v4_modeling_claim: V4ModelingClaimSummary
    v5_evidence_substrate: V5EvidenceSubstrateSummary
    v6_lab_campaign: V6LabCampaignSummary
    target_pair_repairs: list[TargetPairRepairSummary]
    next_actions: NextActions
    non_claims: list[str]
    validation_commands: list[str]


def load_selectivity_repair_inputs(
    repo_root: str | PathLike[str] | None = None,
) -> SelectivityRepairInputs:
    """Load all committed local artifacts needed for the report."""

    root = _resolve_repo_root(repo_root)
    json_payloads: dict[str, dict[str, Any]] = {}
    for source_name in REQUIRED_JSON_SOURCES:
        json_payloads[source_name] = _read_json_object(root / SOURCE_PATHS[source_name])

    source_files = [
        str(relative_path)
        for relative_path in SOURCE_PATHS.values()
        if (root / relative_path).exists()
    ]

    return SelectivityRepairInputs(
        source_files=source_files,
        v4_model_summary=json_payloads["v4_model_summary"],
        v5_expansion_summary=json_payloads["v5_expansion_summary"],
        v5_exact_core_summary=json_payloads["v5_exact_core_summary"],
        v5_tiered_core_summary=json_payloads["v5_tiered_core_summary"],
        v6_campaign_summary=json_payloads["v6_campaign_summary"],
        v6_ticket_rows=_read_csv_rows_if_present(root / SOURCE_PATHS["v6_top_assay_tickets"]),
        v6_gap_closure_campaign=_read_json_object_if_present(
            root / SOURCE_PATHS["v6_gap_closure_campaign"]
        ),
        target_pair_rows=_read_csv_rows_if_present(root / SOURCE_PATHS["target_pairs"]),
    )


def build_selectivity_repair_summary(
    inputs: SelectivityRepairInputs | None = None,
) -> SelectivityRepairReport:
    """Build the v4/v5/v6 benchmark-repair report from local artifacts."""

    loaded_inputs = inputs or load_selectivity_repair_inputs()
    v4_summary = _build_v4_modeling_claim(loaded_inputs.v4_model_summary)
    v5_summary = _build_v5_evidence_substrate(
        loaded_inputs.v5_expansion_summary,
        loaded_inputs.v5_exact_core_summary,
        loaded_inputs.v5_tiered_core_summary,
    )
    v6_summary = _build_v6_lab_campaign(
        loaded_inputs.v6_campaign_summary,
        loaded_inputs.v6_gap_closure_campaign,
    )
    pair_repairs = _build_target_pair_repairs(loaded_inputs, v5_summary, v6_summary)

    return SelectivityRepairReport(
        title="Parasite Selectivity Benchmark-Repair Report",
        purpose=(
            "Summarize the existing v4 modeling claim, the v5 public-data evidence "
            "substrate, and the v6 benchmark-repair ticket campaign so Kira's "
            "original empirical work is easier to inspect and act on."
        ),
        source_files=list(loaded_inputs.source_files),
        v4_modeling_claim=v4_summary,
        v5_evidence_substrate=v5_summary,
        v6_lab_campaign=v6_summary,
        target_pair_repairs=pair_repairs,
        next_actions=_build_next_actions(pair_repairs, v6_summary),
        non_claims=list(NON_CLAIMS),
        validation_commands=list(VALIDATION_COMMANDS),
    )


def write_markdown_report(path: str | PathLike[str]) -> Path:
    """Write a deterministic Markdown dossier to `path`."""

    output_path = Path(path)
    report = build_selectivity_repair_summary()
    output_path.write_text(_render_markdown_report(report), encoding="utf-8")
    return output_path


def report_to_dict(report: SelectivityRepairReport) -> dict[str, Any]:
    """Return a JSON-serializable dictionary representation of a report."""

    return _json_ready(report)


def _resolve_repo_root(repo_root: str | PathLike[str] | None) -> Path:
    if repo_root is not None:
        return Path(repo_root)
    return Path(__file__).resolve().parents[3]


def _read_json_object(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return payload


def _read_json_object_if_present(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return _read_json_object(path)


def _read_csv_rows_if_present(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [
            {str(key): "" if value is None else str(value).strip() for key, value in row.items()}
            for row in reader
        ]


def _build_v4_modeling_claim(summary: Mapping[str, Any]) -> V4ModelingClaimSummary:
    ablations = [item for item in summary.get("ablations", []) if isinstance(item, Mapping)]
    macro_pair_aurocs = {
        str(item.get("ablation")): _float_or_none(item.get("macro_pair_auroc"))
        for item in ablations
        if item.get("ablation")
    }
    best_name, best_value = _best_metric(macro_pair_aurocs)
    compound_plus_pair = macro_pair_aurocs.get("A2_compound_plus_pair")
    compound_only = macro_pair_aurocs.get("A1_compound_only")
    pair_only = macro_pair_aurocs.get("A0_pair_only")
    compound_plus_pair_is_best = (
        best_name == "A2_compound_plus_pair"
        and compound_plus_pair is not None
        and best_value is not None
    )

    if compound_only is not None and pair_only is not None and compound_only > pair_only:
        signal_clause = (
            "Compound chemistry carries most current predictive signal: compound-only "
            "features outperform pair-only features by macro pair AUROC"
        )
    else:
        signal_clause = "The existing v4 summary does not support a stronger pair-only signal"

    if compound_plus_pair_is_best:
        best_clause = "compound + pair is the best classifier by macro pair AUROC"
    else:
        best_clause = "the best macro pair AUROC ablation should be read directly from the summary file"

    return V4ModelingClaimSummary(
        current_modeling_claim="v4 remains the current modeling claim",
        trainable_rows=_int(summary.get("n_rows_trainable")),
        feature_count=_int(summary.get("n_feature_names_full")),
        n_splits=_int(summary.get("n_splits")),
        n_bits=_int(summary.get("n_bits")),
        ablation_macro_pair_aurocs=dict(sorted(macro_pair_aurocs.items())),
        pair_only_macro_pair_auroc=pair_only,
        compound_only_macro_pair_auroc=compound_only,
        compound_plus_pair_macro_pair_auroc=compound_plus_pair,
        best_macro_pair_auroc_ablation=best_name,
        best_macro_pair_auroc=best_value,
        compound_plus_pair_is_best_by_macro_pair_auroc=compound_plus_pair_is_best,
        interpretation=f"{signal_clause}; {best_clause} in the committed v4 summary.",
    )


def _build_v5_evidence_substrate(
    expansion: Mapping[str, Any],
    exact: Mapping[str, Any],
    tiered: Mapping[str, Any],
) -> V5EvidenceSubstrateSummary:
    candidate_status_counts = _int_mapping(expansion.get("candidate_rows_by_status"))
    exact_positive = _int_mapping(exact.get("trainable_positive_rows_by_pair"))
    exact_negative = _int_mapping(exact.get("trainable_negative_rows_by_pair"))
    tiered_pairs = _ordered_pairs(tiered.get("pairs_with_both_classes", []))

    return V5EvidenceSubstrateSummary(
        candidate_evidence_rows=_int(expansion.get("n_candidate_rows")),
        curated_activity_rows=_int(expansion.get("n_curated_activity_rows")),
        exact_matched_ratio_candidate_rows=_int(
            candidate_status_counts.get("exact_matched_ratio")
        ),
        exact_core_rows=_int(exact.get("exact_core_rows")),
        trainable_exact_core_rows=_int(exact.get("trainable_exact_core_rows")),
        tiered_core_rows=_int(tiered.get("tiered_core_rows")),
        trainable_tiered_core_rows=_int(tiered.get("trainable_tiered_core_rows")),
        tiered_pairs_with_both_classes=tiered_pairs,
        exact_pairs_with_both_classes=_pairs_with_both_classes(exact_positive, exact_negative),
        activity_rows_by_side=_int_mapping(expansion.get("activity_rows_by_side")),
        candidate_rows_by_pair=_ordered_int_mapping(expansion.get("candidate_rows_by_pair")),
    )


def _build_v6_lab_campaign(
    summary: Mapping[str, Any],
    campaign: Mapping[str, Any],
) -> V6LabCampaignSummary:
    return V6LabCampaignSummary(
        campaign_mode=str(summary.get("campaign_mode") or "unknown"),
        benchmark_repair_ticket_count=_int(summary.get("generated_ticket_count")),
        tickets_by_pair=_ordered_int_mapping(summary.get("tickets_by_pair")),
        tickets_by_missing_side=_ordered_side_mapping(summary.get("tickets_by_missing_side")),
        eligible_pairs_with_tickets=_ordered_pairs(summary.get("eligible_pairs_with_tickets", [])),
        required_return_fields=[
            str(field) for field in campaign.get("required_return_fields", []) if field
        ],
    )


def _build_target_pair_repairs(
    inputs: SelectivityRepairInputs,
    v5_summary: V5EvidenceSubstrateSummary,
    v6_summary: V6LabCampaignSummary,
) -> list[TargetPairRepairSummary]:
    expansion_pairs = _mapping_of_mappings(inputs.v5_expansion_summary.get("pair_summaries"))
    exact_pair_summaries = _mapping_of_mappings(inputs.v5_exact_core_summary.get("pair_summaries"))
    target_pair_rows = {row.get("pair_id", ""): row for row in inputs.target_pair_rows}
    ticket_missing_by_pair = _ticket_missing_side_counts(inputs.v6_ticket_rows)
    top_scores = _top_priority_scores(inputs.v6_campaign_summary)
    imbalance_basis = _mapping_of_mappings(inputs.v6_campaign_summary.get("target_imbalance_basis"))
    tiered_positive = _int_mapping(inputs.v5_tiered_core_summary.get("trainable_positive_rows_by_pair"))
    tiered_negative = _int_mapping(inputs.v5_tiered_core_summary.get("trainable_negative_rows_by_pair"))
    tiered_rows = _int_mapping(inputs.v5_tiered_core_summary.get("trainable_rows_by_pair"))

    pair_ids = _ordered_pairs(
        list(v6_summary.eligible_pairs_with_tickets)
        + list(v6_summary.tickets_by_pair)
        + list(v5_summary.candidate_rows_by_pair)
    )
    repairs: list[TargetPairRepairSummary] = []

    for pair_id in pair_ids:
        exact_pair = exact_pair_summaries.get(pair_id, {})
        target_pair = target_pair_rows.get(pair_id, {})
        expansion_pair = expansion_pairs.get(pair_id, {})
        exact_positive = _int(exact_pair.get("trainable_positive_rows"))
        exact_negative = _int(exact_pair.get("trainable_negative_rows"))
        tiered_pair_positive = _int(tiered_positive.get(pair_id))
        tiered_pair_negative = _int(tiered_negative.get(pair_id))
        exact_both = exact_positive > 0 and exact_negative > 0
        tiered_both = tiered_pair_positive > 0 and tiered_pair_negative > 0
        missing_side_counts = _ordered_side_mapping(ticket_missing_by_pair.get(pair_id, {}))
        ticket_count = _int(v6_summary.tickets_by_pair.get(pair_id)) or sum(
            missing_side_counts.values()
        )
        evidence_state = str(
            imbalance_basis.get(pair_id, {}).get("class_status")
            or _fallback_evidence_state(tiered_pair_positive, tiered_pair_negative)
        )

        repairs.append(
            TargetPairRepairSummary(
                pair_id=pair_id,
                parasite_label=_first_present(
                    expansion_pair.get("parasite_label"),
                    target_pair.get("parasite_label"),
                ),
                human_label=_first_present(
                    expansion_pair.get("human_label"),
                    target_pair.get("human_label"),
                ),
                evidence_state=evidence_state,
                exact_core_both_classes=exact_both,
                tiered_core_both_classes=tiered_both,
                exact_trainable_rows=_int(exact_pair.get("trainable_exact_core_rows")),
                exact_positive_rows=exact_positive,
                exact_negative_rows=exact_negative,
                tiered_trainable_rows=_int(tiered_rows.get(pair_id)),
                tiered_positive_rows=tiered_pair_positive,
                tiered_negative_rows=tiered_pair_negative,
                missing_side_counts=missing_side_counts,
                missing_side_summary=_format_side_counts(missing_side_counts),
                ticket_count=ticket_count,
                repair_priority=_repair_priority(
                    evidence_state=evidence_state,
                    ticket_count=ticket_count,
                    tiered_core_both_classes=tiered_both,
                    priority_score=top_scores.get(pair_id),
                ),
                explanation=_pair_explanation(
                    pair_id=pair_id,
                    evidence_state=evidence_state,
                    tiered_core_both_classes=tiered_both,
                    positives=tiered_pair_positive,
                    negatives=tiered_pair_negative,
                    ticket_count=ticket_count,
                    missing_side_counts=missing_side_counts,
                ),
            )
        )

    return repairs


def _build_next_actions(
    pair_repairs: list[TargetPairRepairSummary],
    v6_summary: V6LabCampaignSummary,
) -> NextActions:
    closest = [
        pair.pair_id
        for pair in pair_repairs
        if pair.tiered_core_both_classes and pair.evidence_state != "both_classes_extreme_imbalance"
    ]
    parasite_side = [
        pair.pair_id for pair in pair_repairs if pair.missing_side_counts.get("parasite", 0) > 0
    ]
    human_side = [
        pair.pair_id for pair in pair_repairs if pair.missing_side_counts.get("human", 0) > 0
    ]
    returned_data_fields = v6_summary.required_return_fields or [
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
    ]

    return NextActions(
        closest_to_benchmark_ready=closest,
        parasite_side_measurement_pairs=parasite_side,
        human_side_comparator_pairs=human_side,
        returned_data_that_repairs_benchmark=[
            "Return the missing-side comparator assay for the ticketed compound and target pair.",
            "Use comparable potency units and relation fields so a human-divided-by-parasite ratio can be reconstructed.",
            "Include returned-data fields: " + ", ".join(f"`{field}`" for field in returned_data_fields) + ".",
            "Flag data-validity concerns so repaired rows can be audited before joining any benchmark core.",
        ],
    )


def _render_markdown_report(report: SelectivityRepairReport) -> str:
    lines = [
        f"# {report.title}",
        "",
        "## Why This Report Exists",
        "",
        report.purpose,
        "",
        (
            "This dossier returns to Kira's strongest empirical public-data substrate: "
            "parasite-vs-human selectivity evidence from the existing v4/v5/v6 artifacts."
        ),
        "",
        "## Source Artifacts",
        "",
    ]
    lines.extend(f"- `{source}`" for source in report.source_files)
    lines.extend(
        [
            "",
            "## How v4, v5, and v6 Connect",
            "",
            (
                "v4 is the current modeling claim. v5 expands and cleans the assay-aware "
                "public-data evidence substrate. v6 turns gaps in that substrate into "
                "benchmark-repair assay tickets."
            ),
            "",
        ]
    )
    lines.extend(_v4_section(report.v4_modeling_claim))
    lines.extend(_v5_section(report.v5_evidence_substrate))
    lines.extend(_v6_section(report.v6_lab_campaign))
    lines.extend(_pair_repair_section(report.target_pair_repairs))
    lines.extend(_next_actions_section(report.next_actions))
    lines.extend(_direction_section())
    lines.extend(_non_claims_section(report.non_claims))
    lines.extend(_validation_section(report.validation_commands))
    return "\n".join(lines).rstrip() + "\n"


def _v4_section(v4: V4ModelingClaimSummary) -> list[str]:
    lines = [
        "## v4 Modeling Claim Summary",
        "",
        f"- {v4.current_modeling_claim}.",
        f"- Trainable rows: {v4.trainable_rows}.",
        f"- Feature count: {v4.feature_count}.",
        f"- Cross-validation splits: {v4.n_splits}.",
        f"- Morgan fingerprint bits: {v4.n_bits}.",
        f"- Interpretation: {v4.interpretation}",
        "",
        "| Ablation | Macro pair AUROC |",
        "|---|---:|",
    ]
    for ablation, auroc in v4.ablation_macro_pair_aurocs.items():
        lines.append(f"| `{ablation}` | {_format_optional_float(auroc)} |")
    lines.append("")
    return lines


def _v5_section(v5: V5EvidenceSubstrateSummary) -> list[str]:
    lines = [
        "## v5 Evidence Substrate Summary",
        "",
        "| Quantity | Count |",
        "|---|---:|",
        f"| Candidate evidence rows | {v5.candidate_evidence_rows} |",
        f"| Curated activity rows | {v5.curated_activity_rows} |",
        f"| Exact matched-ratio candidate rows | {v5.exact_matched_ratio_candidate_rows} |",
        f"| Exact-core rows | {v5.exact_core_rows} |",
        f"| Trainable exact-core rows | {v5.trainable_exact_core_rows} |",
        f"| Tiered-core rows | {v5.tiered_core_rows} |",
        f"| Trainable tiered-core rows | {v5.trainable_tiered_core_rows} |",
        "",
        f"Tiered pairs with both classes: {', '.join(v5.tiered_pairs_with_both_classes)}.",
        "",
        "| Pair | Candidate rows |",
        "|---|---:|",
    ]
    for pair_id, count in v5.candidate_rows_by_pair.items():
        lines.append(f"| `{pair_id}` | {count} |")
    lines.append("")
    return lines


def _v6_section(v6: V6LabCampaignSummary) -> list[str]:
    lines = [
        "## v6 Lab-Campaign Summary",
        "",
        f"- Campaign mode: `{v6.campaign_mode}`.",
        f"- Benchmark-repair tickets: {v6.benchmark_repair_ticket_count}.",
        "",
        "| Pair | Tickets |",
        "|---|---:|",
    ]
    for pair_id, count in v6.tickets_by_pair.items():
        lines.append(f"| `{pair_id}` | {count} |")
    lines.extend(["", "| Missing side | Tickets |", "|---|---:|"])
    for side, count in v6.tickets_by_missing_side.items():
        lines.append(f"| `{side}` | {count} |")
    lines.append("")
    return lines


def _pair_repair_section(pair_repairs: list[TargetPairRepairSummary]) -> list[str]:
    lines = [
        "## Benchmark Repair Interpretation",
        "",
        "| Pair | Evidence state | Exact both classes | Tiered both classes | Missing side | Tickets | Repair priority |",
        "|---|---|---|---|---|---:|---|",
    ]
    for pair in pair_repairs:
        lines.append(
            "| "
            f"`{pair.pair_id}` | `{pair.evidence_state}` | "
            f"{_yes_no(pair.exact_core_both_classes)} | "
            f"{_yes_no(pair.tiered_core_both_classes)} | "
            f"{pair.missing_side_summary} | {pair.ticket_count} | `{pair.repair_priority}` |"
        )
    lines.append("")
    for pair in pair_repairs:
        lines.extend(
            [
                f"### {pair.pair_id}",
                "",
                f"- Parasite target: {pair.parasite_label}.",
                f"- Human comparator: {pair.human_label}.",
                (
                    f"- Exact core: {pair.exact_trainable_rows} "
                    f"{_plural('trainable row', pair.exact_trainable_rows)} "
                    f"({pair.exact_positive_rows} positive, {pair.exact_negative_rows} negative)."
                ),
                (
                    f"- Tiered core: {pair.tiered_trainable_rows} "
                    f"{_plural('trainable row', pair.tiered_trainable_rows)} "
                    f"({pair.tiered_positive_rows} positive, {pair.tiered_negative_rows} negative)."
                ),
                f"- Repair note: {pair.explanation}",
                "",
            ]
        )
    return lines


def _next_actions_section(next_actions: NextActions) -> list[str]:
    lines = [
        "## Next Actions",
        "",
        (
            "Closest to benchmark-ready: "
            + _format_pair_list(next_actions.closest_to_benchmark_ready)
            + "."
        ),
        (
            "Need parasite-side measurements: "
            + _format_pair_list(next_actions.parasite_side_measurement_pairs)
            + "."
        ),
        (
            "Need human-side comparator measurements: "
            + _format_pair_list(next_actions.human_side_comparator_pairs)
            + "."
        ),
        "",
        "Returned data that would repair the benchmark:",
        "",
    ]
    lines.extend(f"- {action}" for action in next_actions.returned_data_that_repairs_benchmark)
    lines.append("")
    return lines


def _direction_section() -> list[str]:
    return [
        "## Contrast-Driven Inverse-Biology Direction",
        "",
        (
            "The report keeps the inverse-biology loop empirical: start from a parasite-vs-human "
            "contrast, label the evidence state, expose benchmark weakness, issue missing-measurement "
            "tickets, and wait for auditable returned data before making stronger claims."
        ),
        "",
    ]


def _non_claims_section(non_claims: list[str]) -> list[str]:
    lines = ["## Non-Claims", ""]
    lines.extend(f"- {claim}" for claim in non_claims)
    lines.append("")
    return lines


def _validation_section(commands: list[str]) -> list[str]:
    lines = ["## Validation Commands", "", "```bash"]
    lines.extend(commands)
    lines.extend(["```", ""])
    return lines


def _json_ready(value: Any) -> Any:
    if is_dataclass(value):
        return _json_ready(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def _int(value: Any) -> int:
    if value is None or value == "":
        return 0
    return int(float(value))


def _float_or_none(value: Any) -> float | None:
    if value is None or value == "":
        return None
    numeric = float(value)
    if numeric != numeric:
        return None
    return numeric


def _int_mapping(value: Any) -> dict[str, int]:
    if not isinstance(value, Mapping):
        return {}
    return {str(key): _int(item) for key, item in value.items()}


def _ordered_int_mapping(value: Any) -> dict[str, int]:
    mapping = _int_mapping(value)
    return {key: mapping[key] for key in _ordered_pairs(mapping)}


def _ordered_side_mapping(value: Any) -> dict[str, int]:
    mapping = _int_mapping(value)
    ordered = {side: mapping[side] for side in SIDES if side in mapping}
    for key in sorted(set(mapping) - set(SIDES)):
        ordered[key] = mapping[key]
    return ordered


def _mapping_of_mappings(value: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(value, Mapping):
        return {}
    return {
        str(key): dict(item)
        for key, item in value.items()
        if isinstance(item, Mapping)
    }


def _ordered_pairs(values: Any) -> list[str]:
    seen = {str(value) for value in values if value}
    ordered = [pair_id for pair_id in PAIR_ORDER if pair_id in seen]
    ordered.extend(sorted(seen - set(ordered)))
    return ordered


def _pairs_with_both_classes(positive: Mapping[str, int], negative: Mapping[str, int]) -> list[str]:
    return _ordered_pairs(
        pair_id
        for pair_id in set(positive) | set(negative)
        if positive.get(pair_id, 0) > 0 and negative.get(pair_id, 0) > 0
    )


def _best_metric(metrics: Mapping[str, float | None]) -> tuple[str | None, float | None]:
    present = [(name, value) for name, value in metrics.items() if value is not None]
    if not present:
        return None, None
    name, value = max(present, key=lambda item: (item[1], item[0]))
    return name, value


def _ticket_missing_side_counts(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        pair_id = row.get("pair_id", "")
        missing_side = row.get("missing_side", "")
        if pair_id and missing_side:
            counts[pair_id][missing_side] += 1
    return {pair_id: dict(counter) for pair_id, counter in counts.items()}


def _top_priority_scores(summary: Mapping[str, Any]) -> dict[str, float]:
    scores: dict[str, float] = {}
    rows = summary.get("top_priority_pairs", [])
    if not isinstance(rows, list):
        return scores
    for row in rows:
        if not isinstance(row, Mapping) or not row.get("pair_id"):
            continue
        score = _float_or_none(row.get("max_priority_score"))
        if score is not None:
            scores[str(row["pair_id"])] = score
    return scores


def _fallback_evidence_state(positives: int, negatives: int) -> str:
    if positives > 0 and negatives > 0:
        return "both_classes_present"
    if positives == 0 and negatives == 0:
        return "no_trainable_tiered_rows"
    if positives == 0:
        return "class_degenerate_zero_positive"
    return "class_degenerate_zero_negative"


def _repair_priority(
    *,
    evidence_state: str,
    ticket_count: int,
    tiered_core_both_classes: bool,
    priority_score: float | None,
) -> str:
    if evidence_state.startswith("class_degenerate") and ticket_count >= 20:
        return "highest"
    if evidence_state == "both_classes_extreme_imbalance":
        return "high"
    if tiered_core_both_classes and ticket_count <= 5:
        return "closest-to-ready"
    if priority_score is not None and priority_score >= 100:
        return "high"
    if ticket_count >= 10:
        return "medium"
    return "focused"


def _pair_explanation(
    *,
    pair_id: str,
    evidence_state: str,
    tiered_core_both_classes: bool,
    positives: int,
    negatives: int,
    ticket_count: int,
    missing_side_counts: Mapping[str, int],
) -> str:
    missing = _format_side_counts(missing_side_counts)
    if not tiered_core_both_classes:
        return (
            f"{pair_id} has {positives} positive and {negatives} negative trainable tiered rows, "
            f"so benchmark repair should prioritize {ticket_count} missing-side tickets "
            f"({missing}) that can create matched comparator evidence."
        )
    if evidence_state == "both_classes_extreme_imbalance":
        return (
            f"{pair_id} already has both classes but remains highly imbalanced "
            f"({positives} positive, {negatives} negative); {ticket_count} tickets "
            f"({missing}) should add comparator evidence and reduce benchmark skew."
        )
    return (
        f"{pair_id} is closest to benchmark-ready because both classes are present "
        f"({positives} positive, {negatives} negative); the {ticket_count} tickets "
        f"({missing}) mainly improve coverage and matched-evidence robustness."
    )


def _first_present(*values: Any) -> str:
    for value in values:
        if value:
            return str(value)
    return "unknown"


def _format_side_counts(counts: Mapping[str, int]) -> str:
    if not counts:
        return "none"
    parts = [f"{side}={counts[side]}" for side in SIDES if side in counts]
    parts.extend(f"{side}={counts[side]}" for side in sorted(set(counts) - set(SIDES)))
    return ", ".join(parts)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def _format_pair_list(pair_ids: list[str]) -> str:
    if not pair_ids:
        return "none"
    return ", ".join(pair_ids)


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"


def _plural(label: str, count: int) -> str:
    if count == 1:
        return label
    return f"{label}s"
