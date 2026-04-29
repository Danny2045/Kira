from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from kira.experiments.design_v6_lab_campaign import (
    SUMMARY_FIELDS,
    TICKET_FIELDS,
    UNKNOWN,
    generate_lab_tickets,
    main,
    summarize_lab_campaign,
)


def _candidate(
    *,
    pair_id: str = "PairA",
    molecule_chembl_id: str = "CHEMBL1",
    standard_inchi_key: str = "INCHI1",
    candidate_status: str = "single_side_only",
    single_side: str | None = "parasite",
    canonical_smiles: str | None = "CCO",
    murcko_scaffold: str | None = "scaffold-a",
    known_value_nM: float = 50.0,
) -> dict:
    row = {
        "candidate_row_id": f"{pair_id}::{molecule_chembl_id}::1",
        "pair_id": pair_id,
        "parasite_label": f"{pair_id} parasite target",
        "human_label": f"{pair_id} human target",
        "molecule_chembl_id": molecule_chembl_id,
        "canonical_smiles": canonical_smiles,
        "standard_inchi_key": standard_inchi_key,
        "candidate_status": candidate_status,
        "comparison_tier": candidate_status,
        "ratio_exact_human_div_parasite": None,
        "ratio_lower_bound_human_div_parasite": None,
        "ratio_upper_bound_human_div_parasite": None,
        "murcko_scaffold": murcko_scaffold,
        "structure_status": "resolved",
        "same_standard_type": None,
        "same_assay_type": None,
        "same_document": None,
        "single_side": single_side,
        "parasite_target_chembl_id": "CHEMBL_PARASITE",
        "human_target_chembl_id": "CHEMBL_HUMAN",
        "parasite_activity_chembl_id": None,
        "human_activity_chembl_id": None,
        "parasite_assay_chembl_id": None,
        "human_assay_chembl_id": None,
        "parasite_standard_type": None,
        "human_standard_type": None,
        "parasite_standard_relation": None,
        "human_standard_relation": None,
        "parasite_standard_value_nM": None,
        "human_standard_value_nM": None,
        "parasite_assay_type": None,
        "human_assay_type": None,
        "parasite_measurement_interval_kind": None,
        "human_measurement_interval_kind": None,
    }
    if single_side in {"parasite", "human"}:
        row[f"{single_side}_activity_chembl_id"] = f"ACT_{single_side}_{molecule_chembl_id}"
        row[f"{single_side}_assay_chembl_id"] = f"ASSAY_{single_side}_{molecule_chembl_id}"
        row[f"{single_side}_standard_type"] = "IC50"
        row[f"{single_side}_standard_relation"] = "="
        row[f"{single_side}_standard_value_nM"] = known_value_nM
        row[f"{single_side}_assay_type"] = "B"
        row[f"{single_side}_measurement_interval_kind"] = "exact"
    return row


def _exact_candidate(**overrides) -> dict:
    row = _candidate(candidate_status="exact_matched_ratio", single_side=None, **overrides)
    for side, value in (("parasite", 10.0), ("human", 200.0)):
        row[f"{side}_activity_chembl_id"] = f"ACT_{side}"
        row[f"{side}_assay_chembl_id"] = f"ASSAY_{side}"
        row[f"{side}_standard_type"] = "IC50"
        row[f"{side}_standard_relation"] = "="
        row[f"{side}_standard_value_nM"] = value
        row[f"{side}_assay_type"] = "B"
        row[f"{side}_measurement_interval_kind"] = "exact"
    row["ratio_exact_human_div_parasite"] = 20.0
    row["ratio_lower_bound_human_div_parasite"] = 20.0
    row["ratio_upper_bound_human_div_parasite"] = 20.0
    return row


def _summary(
    *,
    rows: dict[str, int] | None = None,
    positives: dict[str, int] | None = None,
    negatives: dict[str, int] | None = None,
    both_classes: list[str] | None = None,
) -> dict:
    return {
        "trainable_rows_by_pair": rows or {"PairA": 0},
        "trainable_positive_rows_by_pair": positives or {"PairA": 0},
        "trainable_negative_rows_by_pair": negatives or {"PairA": 0},
        "pairs_with_both_classes": both_classes or [],
    }


def _tiered_core(*rows: dict) -> pd.DataFrame:
    columns = [
        "pair_id",
        "compound_key",
        "molecule_chembl_id",
        "standard_inchi_key",
        "murcko_scaffold",
        "is_trainable_tiered_core",
        "is_selective",
    ]
    return pd.DataFrame(rows, columns=columns)


def test_single_side_only_row_becomes_ticket() -> None:
    tickets = generate_lab_tickets(
        pd.DataFrame([_candidate(single_side="parasite")]),
        _tiered_core(),
        _summary(),
        top_n=None,
    )

    assert len(tickets) == 1
    assert tickets.iloc[0]["known_side"] == "parasite"
    assert tickets.iloc[0]["missing_side"] == "human"
    assert tickets.iloc[0]["evidence_status_source"] == "single_side_only"


def test_exact_matched_ratio_row_does_not_become_missing_side_ticket() -> None:
    tickets = generate_lab_tickets(
        pd.DataFrame([_exact_candidate()]),
        _tiered_core(),
        _summary(),
        top_n=None,
    )

    assert tickets.empty


def test_priority_score_increases_for_class_degenerate_pairs() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="Degenerate", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1"),
            _candidate(pair_id="Balanced", molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2"),
        ]
    )
    tickets = generate_lab_tickets(
        candidates,
        _tiered_core(),
        _summary(
            rows={"Degenerate": 4, "Balanced": 6},
            positives={"Degenerate": 0, "Balanced": 3},
            negatives={"Degenerate": 4, "Balanced": 3},
            both_classes=["Balanced"],
        ),
        top_n=None,
    )
    scores = dict(zip(tickets["pair_id"], tickets["priority_score"], strict=True))

    assert scores["Degenerate"] > scores["Balanced"]


def test_priority_score_rewards_strong_known_side_potency() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1", known_value_nM=10.0),
            _candidate(molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2", known_value_nM=6000.0),
        ]
    )
    tickets = generate_lab_tickets(candidates, _tiered_core(), _summary(), top_n=None)
    scores = dict(zip(tickets["compound_key"], tickets["priority_score"], strict=True))

    assert scores["INCHI1"] > scores["INCHI2"]


def test_missing_known_activity_id_does_not_reduce_priority_score_when_assay_id_is_present() -> None:
    missing_activity = _candidate(
        molecule_chembl_id="CHEMBL_MISSING",
        standard_inchi_key="INCHI_MISSING",
        murcko_scaffold="scaffold-missing",
    )
    missing_activity["parasite_activity_chembl_id"] = None
    present_activity = _candidate(
        molecule_chembl_id="CHEMBL_PRESENT",
        standard_inchi_key="INCHI_PRESENT",
        murcko_scaffold="scaffold-present",
    )

    tickets = generate_lab_tickets(
        pd.DataFrame([missing_activity, present_activity]),
        _tiered_core(),
        _summary(),
        top_n=None,
    )
    by_key = tickets.set_index("compound_key")

    assert by_key.loc["INCHI_MISSING", "known_activity_chembl_id"] == UNKNOWN
    assert by_key.loc["INCHI_MISSING", "known_assay_chembl_id"] == "ASSAY_parasite_CHEMBL_MISSING"
    assert by_key.loc["INCHI_MISSING", "priority_score"] == by_key.loc["INCHI_PRESENT", "priority_score"]
    assert "known activity id is missing" not in by_key.loc["INCHI_MISSING", "priority_reason"]
    assert (
        by_key.loc["INCHI_MISSING", "provenance_note"]
        == "known assay provenance present; known activity id unavailable in v5 candidate table"
    )


def test_provenance_note_describes_available_candidate_table_ids() -> None:
    both_present = _candidate(
        molecule_chembl_id="CHEMBL_BOTH",
        standard_inchi_key="INCHI_BOTH",
        murcko_scaffold="scaffold-both",
    )
    assay_only = _candidate(
        molecule_chembl_id="CHEMBL_ASSAY_ONLY",
        standard_inchi_key="INCHI_ASSAY_ONLY",
        murcko_scaffold="scaffold-assay-only",
    )
    assay_only["parasite_activity_chembl_id"] = None
    both_missing = _candidate(
        molecule_chembl_id="CHEMBL_MISSING",
        standard_inchi_key="INCHI_MISSING",
        murcko_scaffold="scaffold-missing",
    )
    both_missing["parasite_activity_chembl_id"] = None
    both_missing["parasite_assay_chembl_id"] = None

    tickets = generate_lab_tickets(
        pd.DataFrame([both_present, assay_only, both_missing]),
        _tiered_core(),
        _summary(),
        top_n=None,
    )
    notes = dict(zip(tickets["compound_key"], tickets["provenance_note"], strict=True))

    assert notes["INCHI_BOTH"] == "known assay and activity provenance present"
    assert (
        notes["INCHI_ASSAY_ONLY"]
        == "known assay provenance present; known activity id unavailable in v5 candidate table"
    )
    assert notes["INCHI_MISSING"] == "known assay and activity provenance unavailable in v5 candidate table"


def test_compounds_already_represented_in_tiered_core_are_penalized() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI_REP", murcko_scaffold="scaffold-rep"),
            _candidate(molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI_NEW", murcko_scaffold="scaffold-new"),
        ]
    )
    tiered = _tiered_core(
        {
            "pair_id": "PairA",
            "compound_key": "INCHI_REP",
            "molecule_chembl_id": "CHEMBL1",
            "standard_inchi_key": "INCHI_REP",
            "murcko_scaffold": "scaffold-rep",
            "is_trainable_tiered_core": 1,
            "is_selective": 0,
        }
    )

    tickets = generate_lab_tickets(candidates, tiered, _summary(rows={"PairA": 1}, negatives={"PairA": 1}), top_n=None)
    scores = dict(zip(tickets["compound_key"], tickets["priority_score"], strict=True))

    assert scores["INCHI_REP"] < scores["INCHI_NEW"]


def test_scaffold_diversity_behavior_is_deterministic() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1", murcko_scaffold="scaffold-a"),
            _candidate(molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2", murcko_scaffold="scaffold-a"),
            _candidate(molecule_chembl_id="CHEMBL3", standard_inchi_key="INCHI3", murcko_scaffold="scaffold-b"),
        ]
    )

    first = generate_lab_tickets(candidates, _tiered_core(), _summary(), top_n=None)
    second = generate_lab_tickets(candidates, _tiered_core(), _summary(), top_n=None)
    same_scaffold_scores = first[first["murcko_scaffold"] == "scaffold-a"]["priority_score"].tolist()

    assert first[["compound_key", "priority_score"]].to_dict(orient="records") == second[
        ["compound_key", "priority_score"]
    ].to_dict(orient="records")
    assert same_scaffold_scores[0] > same_scaffold_scores[1]


def test_ticket_has_all_required_fields() -> None:
    tickets = generate_lab_tickets(pd.DataFrame([_candidate()]), _tiered_core(), _summary(), top_n=None)

    assert list(tickets.columns) == TICKET_FIELDS
    assert set(TICKET_FIELDS) <= set(tickets.iloc[0].index)


def test_summary_json_has_all_required_fields() -> None:
    candidates = pd.DataFrame([_candidate()])
    tiered = _tiered_core()
    summary_source = _summary()
    tickets = generate_lab_tickets(candidates, tiered, summary_source, top_n=None)
    summary = summarize_lab_campaign(candidates, tiered, summary_source, tickets)

    assert SUMMARY_FIELDS <= set(summary)
    assert summary["input_candidate_rows"] == 1
    assert summary["generated_ticket_count"] == 1
    assert "experiment-prioritization layer" in summary["scientific_scope"]


def test_summary_json_includes_known_provenance_counts() -> None:
    both_present = _candidate(
        molecule_chembl_id="CHEMBL_BOTH",
        standard_inchi_key="INCHI_BOTH",
        murcko_scaffold="scaffold-both",
    )
    assay_only = _candidate(
        molecule_chembl_id="CHEMBL_ASSAY_ONLY",
        standard_inchi_key="INCHI_ASSAY_ONLY",
        murcko_scaffold="scaffold-assay-only",
    )
    assay_only["parasite_activity_chembl_id"] = None
    both_missing = _candidate(
        molecule_chembl_id="CHEMBL_MISSING",
        standard_inchi_key="INCHI_MISSING",
        murcko_scaffold="scaffold-missing",
    )
    both_missing["parasite_activity_chembl_id"] = None
    both_missing["parasite_assay_chembl_id"] = None
    candidates = pd.DataFrame([both_present, assay_only, both_missing])
    tiered = _tiered_core()
    summary_source = _summary()

    tickets = generate_lab_tickets(candidates, tiered, summary_source, top_n=None)
    summary = summarize_lab_campaign(candidates, tiered, summary_source, tickets)

    assert summary["tickets_with_known_assay_id"] == 2
    assert summary["tickets_with_known_activity_id"] == 1
    assert summary["tickets_missing_known_activity_id"] == 2


def test_output_sorting_is_deterministic() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="PairB", molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2"),
            _candidate(pair_id="PairA", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1"),
            _candidate(pair_id="PairA", molecule_chembl_id="CHEMBL3", standard_inchi_key="INCHI3", known_value_nM=10.0),
        ]
    )
    tickets = generate_lab_tickets(candidates, _tiered_core(), _summary(), top_n=None)
    expected = tickets.sort_values(
        ["priority_score", "pair_id", "compound_key"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    assert tickets[["priority_score", "pair_id", "compound_key"]].equals(
        expected[["priority_score", "pair_id", "compound_key"]]
    )


def test_benchmark_repair_mode_does_not_collapse_to_two_pairs_when_more_exist() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="PairA", molecule_chembl_id=f"CHEMBLA{i}", standard_inchi_key=f"INCHIA{i}", known_value_nM=1.0)
            for i in range(4)
        ]
        + [
            _candidate(pair_id="PairB", molecule_chembl_id=f"CHEMBLB{i}", standard_inchi_key=f"INCHIB{i}", known_value_nM=100.0)
            for i in range(2)
        ]
        + [
            _candidate(pair_id="PairC", molecule_chembl_id=f"CHEMBLC{i}", standard_inchi_key=f"INCHIC{i}", known_value_nM=1000.0)
            for i in range(2)
        ]
    )

    tickets = generate_lab_tickets(
        candidates,
        _tiered_core(),
        _summary(),
        top_n=6,
        campaign_mode="benchmark_repair",
        max_tickets_per_pair=4,
        min_tickets_per_pair=1,
    )

    assert tickets["pair_id"].nunique() == 3


def test_potency_discovery_preserves_highest_score_ordering_behavior() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="Potent", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1", known_value_nM=1.0),
            _candidate(pair_id="Potent", molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2", known_value_nM=2.0),
            _candidate(pair_id="Weak", molecule_chembl_id="CHEMBL3", standard_inchi_key="INCHI3", known_value_nM=6000.0),
        ]
    )

    tickets = generate_lab_tickets(
        candidates,
        _tiered_core(),
        _summary(),
        top_n=2,
        campaign_mode="potency_discovery",
    )

    assert tickets["pair_id"].tolist() == ["Potent", "Potent"]
    assert tickets["priority_score"].tolist() == sorted(tickets["priority_score"], reverse=True)


def test_pair_specific_only_returns_requested_pair() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="PairA", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1"),
            _candidate(pair_id="PairB", molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2"),
        ]
    )

    tickets = generate_lab_tickets(
        candidates,
        _tiered_core(),
        _summary(),
        top_n=None,
        campaign_mode="pair_specific",
        pair_id="PairB",
    )

    assert set(tickets["pair_id"]) == {"PairB"}


def test_per_pair_cap_works_in_benchmark_repair_mode() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="PairA", molecule_chembl_id=f"CHEMBLA{i}", standard_inchi_key=f"INCHIA{i}")
            for i in range(5)
        ]
        + [
            _candidate(pair_id="PairB", molecule_chembl_id=f"CHEMBLB{i}", standard_inchi_key=f"INCHIB{i}")
            for i in range(5)
        ]
    )

    tickets = generate_lab_tickets(
        candidates,
        _tiered_core(),
        _summary(),
        top_n=10,
        campaign_mode="benchmark_repair",
        max_tickets_per_pair=2,
        min_tickets_per_pair=0,
    )

    assert tickets.groupby("pair_id").size().max() <= 2


def test_summary_includes_campaign_mode_and_eligible_pairs() -> None:
    candidates = pd.DataFrame(
        [
            _candidate(pair_id="PairA", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1"),
            _candidate(pair_id="PairB", molecule_chembl_id="CHEMBL2", standard_inchi_key="INCHI2"),
        ]
    )
    tiered = _tiered_core()
    summary_source = _summary()
    tickets = generate_lab_tickets(candidates, tiered, summary_source, top_n=1, campaign_mode="benchmark_repair")
    summary = summarize_lab_campaign(
        candidates,
        tiered,
        summary_source,
        tickets,
        campaign_mode="benchmark_repair",
        max_tickets_per_pair=25,
        min_tickets_per_pair=5,
        eligible_pairs_with_tickets=tickets.attrs["eligible_pairs_with_tickets"],
    )

    assert summary["campaign_mode"] == "benchmark_repair"
    assert summary["eligible_pairs_with_tickets"] == ["PairA", "PairB"]


def test_cli_writes_csv_and_json_outputs(tmp_path: Path) -> None:
    candidate_csv = tmp_path / "candidates.csv"
    tiered_csv = tmp_path / "tiered.csv"
    tiered_summary_json = tmp_path / "tiered_summary.json"
    output_csv = tmp_path / "tickets.csv"
    output_campaign_json = tmp_path / "campaign.json"
    output_summary_json = tmp_path / "summary.json"

    pd.DataFrame([_candidate()]).to_csv(candidate_csv, index=False)
    _tiered_core().to_csv(tiered_csv, index=False)
    tiered_summary_json.write_text(json.dumps(_summary()), encoding="utf-8")

    exit_code = main(
        [
            "--candidate-csv",
            str(candidate_csv),
            "--tiered-core-csv",
            str(tiered_csv),
            "--tiered-summary-json",
            str(tiered_summary_json),
            "--output-csv",
            str(output_csv),
            "--output-campaign-json",
            str(output_campaign_json),
            "--output-summary-json",
            str(output_summary_json),
            "--top-n",
            "5",
        ]
    )

    assert exit_code == 0
    assert pd.read_csv(output_csv).shape[0] == 1
    campaign = json.loads(output_campaign_json.read_text(encoding="utf-8"))
    summary = json.loads(output_summary_json.read_text(encoding="utf-8"))
    assert campaign["ticket_count"] == 1
    assert SUMMARY_FIELDS <= set(summary)
