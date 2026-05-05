from __future__ import annotations

import json
from pathlib import Path

from kira.selectivity.benchmark_repair_report import (
    build_selectivity_repair_summary,
    load_selectivity_repair_inputs,
    report_to_dict,
    write_markdown_report,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_V6_PAIRS = ["LmDHFR", "LmPTR1", "SmDHODH", "SmHDAC8", "TbCathB", "TbPDEB1"]
EXPECTED_TIERED_BOTH_CLASSES = ["LmDHFR", "SmHDAC8", "TbCathB"]
FORBIDDEN_CLAIMS = (
    "discovered drugs",
    "wet-lab validated",
    "clinical cure",
    "new model-performance result",
    "solved AMR",
    "solved regeneration",
)


def _json(path: str) -> dict:
    return json.loads((REPO_ROOT / path).read_text(encoding="utf-8"))


def test_report_loads_from_existing_local_files_without_network() -> None:
    inputs = load_selectivity_repair_inputs(REPO_ROOT)
    report = build_selectivity_repair_summary(inputs)

    assert "results/selectivity_v4/summary.json" in inputs.source_files
    assert "data/processed/selectivity_v6_campaign_summary.json" in report.source_files
    assert all(not source.startswith(("http://", "https://")) for source in report.source_files)
    assert report.v6_lab_campaign.benchmark_repair_ticket_count == 100


def test_key_v5_v6_numbers_match_existing_summary_files() -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))
    expansion = _json("data/processed/selectivity_v5_expansion_summary.json")
    exact = _json("data/processed/selectivity_v5_exact_core_summary.json")
    tiered = _json("data/processed/selectivity_v5_tiered_core_summary.json")
    v6 = _json("data/processed/selectivity_v6_campaign_summary.json")

    assert report.v5_evidence_substrate.candidate_evidence_rows == expansion["n_candidate_rows"] == 12091
    assert report.v5_evidence_substrate.curated_activity_rows == expansion["n_curated_activity_rows"] == 11703
    assert (
        report.v5_evidence_substrate.exact_matched_ratio_candidate_rows
        == expansion["candidate_rows_by_status"]["exact_matched_ratio"]
        == 1227
    )
    assert report.v5_evidence_substrate.exact_core_rows == exact["exact_core_rows"] == 114
    assert (
        report.v5_evidence_substrate.trainable_exact_core_rows
        == exact["trainable_exact_core_rows"]
        == 110
    )
    assert report.v5_evidence_substrate.tiered_core_rows == tiered["tiered_core_rows"] == 135
    assert (
        report.v5_evidence_substrate.trainable_tiered_core_rows
        == tiered["trainable_tiered_core_rows"]
        == 131
    )
    assert report.v6_lab_campaign.benchmark_repair_ticket_count == v6["generated_ticket_count"] == 100
    assert report.v6_lab_campaign.tickets_by_pair == v6["tickets_by_pair"]


def test_report_includes_all_six_v6_target_pairs() -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))

    assert [pair.pair_id for pair in report.target_pair_repairs] == EXPECTED_V6_PAIRS
    assert list(report.v6_lab_campaign.tickets_by_pair) == EXPECTED_V6_PAIRS


def test_report_identifies_v5_tiered_pairs_with_both_classes() -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))

    assert report.v5_evidence_substrate.tiered_pairs_with_both_classes == EXPECTED_TIERED_BOTH_CLASSES
    tiered_both = [
        pair.pair_id for pair in report.target_pair_repairs if pair.tiered_core_both_classes
    ]
    assert tiered_both == EXPECTED_TIERED_BOTH_CLASSES


def test_report_identifies_v6_missing_side_distribution() -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))

    assert report.v6_lab_campaign.tickets_by_missing_side == {"parasite": 99, "human": 1}
    lmptr1 = next(pair for pair in report.target_pair_repairs if pair.pair_id == "LmPTR1")
    assert lmptr1.missing_side_counts == {"parasite": 24, "human": 1}


def test_markdown_report_is_deterministic(tmp_path: Path) -> None:
    first_path = write_markdown_report(tmp_path / "first.md")
    second_path = write_markdown_report(tmp_path / "second.md")

    assert first_path.read_text(encoding="utf-8") == second_path.read_text(encoding="utf-8")


def test_report_to_dict_is_json_serializable() -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))
    payload = report_to_dict(report)

    assert json.loads(json.dumps(payload, sort_keys=True)) == payload
    assert payload["v5_evidence_substrate"]["candidate_evidence_rows"] == 12091
    assert payload["v6_lab_campaign"]["tickets_by_missing_side"] == {"parasite": 99, "human": 1}


def test_generated_report_does_not_contain_forbidden_claims(tmp_path: Path) -> None:
    report = build_selectivity_repair_summary(load_selectivity_repair_inputs(REPO_ROOT))
    markdown_path = write_markdown_report(tmp_path / "report.md")
    generated = json.dumps(
        {
            "markdown": markdown_path.read_text(encoding="utf-8"),
            "report": report_to_dict(report),
        },
        sort_keys=True,
    ).lower()

    for claim in FORBIDDEN_CLAIMS:
        assert claim.lower() not in generated
