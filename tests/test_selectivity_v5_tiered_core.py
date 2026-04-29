from __future__ import annotations

import math

import pandas as pd

from kira.experiments.selectivity_v5_tiered_core import (
    add_tiered_evidence_labels,
    aggregate_tiered_core,
    summarize_tiered_core,
)


def _row(**overrides):
    row = {
        "candidate_row_id": "PairA::CHEMBL1::1",
        "pair_id": "PairA",
        "parasite_label": "Parasite target",
        "human_label": "Human target",
        "molecule_chembl_id": "CHEMBL1",
        "canonical_smiles": "CCO",
        "standard_inchi_key": "INCHI1",
        "murcko_scaffold": "c1ccccc1",
        "candidate_status": "exact_matched_ratio",
        "comparison_tier": "strict",
        "ratio_exact_human_div_parasite": 12.0,
        "ratio_lower_bound_human_div_parasite": 12.0,
        "ratio_upper_bound_human_div_parasite": 12.0,
        "parasite_standard_type": "IC50",
        "human_standard_type": "IC50",
        "parasite_standard_relation": "=",
        "human_standard_relation": "=",
        "parasite_measurement_interval_kind": "exact",
        "human_measurement_interval_kind": "exact",
        "parasite_assay_type": "B",
        "human_assay_type": "B",
        "parasite_activity_chembl_id": "ACTP1",
        "human_activity_chembl_id": "ACTH1",
        "parasite_assay_chembl_id": "ASSP1",
        "human_assay_chembl_id": "ASSH1",
        "parasite_document_chembl_id": "DOCP1",
        "human_document_chembl_id": "DOCH1",
    }
    row.update(overrides)
    return row


def _labeled(*rows):
    return add_tiered_evidence_labels(pd.DataFrame(rows), selectivity_threshold=10.0)


def test_exact_positive_classification():
    labeled = _labeled(_row(ratio_exact_human_div_parasite=10.0))

    assert labeled.iloc[0]["tiered_label"] == "positive"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "exact_matched_ratio_at_or_above_threshold"


def test_exact_negative_classification():
    labeled = _labeled(_row(ratio_exact_human_div_parasite=9.99))

    assert labeled.iloc[0]["tiered_label"] == "negative"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "exact_matched_ratio_below_threshold"


def test_lower_bound_definite_positive_classification():
    labeled = _labeled(
        _row(
            candidate_status="lower_bound_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=10.0,
            ratio_upper_bound_human_div_parasite=None,
        )
    )

    assert labeled.iloc[0]["tiered_label"] == "positive"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "lower_bound_ratio_at_or_above_threshold"


def test_upper_bound_definite_negative_classification():
    labeled = _labeled(
        _row(
            candidate_status="upper_bound_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=None,
            ratio_upper_bound_human_div_parasite=9.99,
        )
    )

    assert labeled.iloc[0]["tiered_label"] == "negative"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "upper_bound_ratio_below_threshold"


def test_interval_definite_positive_classification():
    labeled = _labeled(
        _row(
            candidate_status="interval_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=10.0,
            ratio_upper_bound_human_div_parasite=100.0,
        )
    )

    assert labeled.iloc[0]["tiered_label"] == "positive"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "interval_lower_bound_at_or_above_threshold"


def test_interval_definite_negative_classification():
    labeled = _labeled(
        _row(
            candidate_status="interval_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=0.1,
            ratio_upper_bound_human_div_parasite=9.99,
        )
    )

    assert labeled.iloc[0]["tiered_label"] == "negative"
    assert labeled.iloc[0]["tiered_evidence_reason"] == "interval_upper_bound_below_threshold"


def test_interval_crossing_threshold_unresolved():
    labeled = _labeled(
        _row(
            candidate_status="interval_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=2.0,
            ratio_upper_bound_human_div_parasite=20.0,
        )
    )

    assert labeled.iloc[0]["tiered_label"] is None
    assert not bool(labeled.iloc[0]["is_tiered_decidable"])
    assert labeled.iloc[0]["tiered_evidence_reason"] == "interval_crosses_threshold_unresolved"


def test_aggregation_to_one_row_per_pair_compound():
    labeled = _labeled(
        _row(ratio_exact_human_div_parasite=20.0, parasite_activity_chembl_id="ACTP1"),
        _row(
            candidate_row_id="PairA::CHEMBL1::2",
            ratio_exact_human_div_parasite=80.0,
            ratio_lower_bound_human_div_parasite=80.0,
            ratio_upper_bound_human_div_parasite=80.0,
            parasite_activity_chembl_id="ACTP2",
        ),
    )

    core = aggregate_tiered_core(labeled, selectivity_threshold=10.0)

    assert len(core) == 1
    row = core.iloc[0]
    assert row["tiered_label"] == "positive"
    assert row["is_trainable_tiered_core"] == 1
    assert row["tiered_evidence_count"] == 2
    assert math.isclose(row["selectivity_ratio_human_div_parasite"], math.sqrt(20.0 * 80.0))
    assert row["parasite_activity_chembl_ids"] == "ACTP1;ACTP2"


def test_conflict_detection_and_exclusion_from_trainable_core():
    labeled = _labeled(
        _row(ratio_exact_human_div_parasite=20.0),
        _row(
            candidate_row_id="PairA::CHEMBL1::2",
            candidate_status="upper_bound_ratio",
            ratio_exact_human_div_parasite=None,
            ratio_lower_bound_human_div_parasite=None,
            ratio_upper_bound_human_div_parasite=5.0,
            parasite_activity_chembl_id="ACTP2",
        ),
    )

    core = aggregate_tiered_core(labeled, selectivity_threshold=10.0)

    assert len(core) == 1
    row = core.iloc[0]
    assert row["tiered_label"] == "conflicting"
    assert row["label_status"] == "tiered_core_conflicting"
    assert row["is_trainable_tiered_core"] == 0
    assert row["n_positive_tiered_observations"] == 1
    assert row["n_negative_tiered_observations"] == 1


def test_fallback_compound_key_logic():
    labeled = _labeled(_row(standard_inchi_key=""))
    core = aggregate_tiered_core(labeled, selectivity_threshold=10.0)

    assert labeled.iloc[0]["compound_key"] == "CHEMBL1"
    assert core.iloc[0]["compound_key"] == "CHEMBL1"


def test_summary_fields_are_present():
    candidate_rows = pd.DataFrame(
        [
            _row(pair_id="PairA", molecule_chembl_id="CHEMBL1", standard_inchi_key="INCHI1"),
            _row(
                candidate_row_id="PairA::CHEMBL2::1",
                pair_id="PairA",
                molecule_chembl_id="CHEMBL2",
                standard_inchi_key="INCHI2",
                ratio_exact_human_div_parasite=5.0,
                ratio_lower_bound_human_div_parasite=5.0,
                ratio_upper_bound_human_div_parasite=5.0,
            ),
            _row(
                candidate_row_id="PairB::CHEMBL3::1",
                pair_id="PairB",
                molecule_chembl_id="CHEMBL3",
                standard_inchi_key="INCHI3",
                candidate_status="interval_ratio",
                ratio_exact_human_div_parasite=None,
                ratio_lower_bound_human_div_parasite=2.0,
                ratio_upper_bound_human_div_parasite=20.0,
            ),
        ]
    )
    labeled = add_tiered_evidence_labels(candidate_rows, selectivity_threshold=10.0)
    core = aggregate_tiered_core(labeled, selectivity_threshold=10.0)

    summary = summarize_tiered_core(candidate_rows, labeled, core, selectivity_threshold=10.0)

    expected_fields = {
        "input_candidate_rows",
        "decidable_candidate_rows",
        "unresolved_candidate_rows",
        "tiered_core_rows",
        "trainable_tiered_core_rows",
        "conflicting_tiered_core_rows",
        "candidate_rows_by_evidence_status",
        "decidable_rows_by_evidence_status",
        "trainable_rows_by_pair",
        "trainable_positive_rows_by_pair",
        "trainable_negative_rows_by_pair",
        "trainable_unique_scaffolds_by_pair",
        "pairs_with_both_classes",
        "scientific_scope",
    }
    assert expected_fields <= set(summary)
    assert summary["input_candidate_rows"] == 3
    assert summary["decidable_candidate_rows"] == 2
    assert summary["unresolved_candidate_rows"] == 1
    assert summary["trainable_positive_rows_by_pair"] == {"PairA": 1}
    assert summary["trainable_negative_rows_by_pair"] == {"PairA": 1}
    assert summary["pairs_with_both_classes"] == ["PairA"]
    assert "not a final modeling claim" in summary["scientific_scope"]
