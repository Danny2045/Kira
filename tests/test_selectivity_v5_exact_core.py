from __future__ import annotations

import math

import pandas as pd

from kira.experiments.selectivity_v5_exact_core import (
    aggregate_exact_core,
    filter_exact_core_candidates,
)


def _row(**overrides):
    row = {
        "pair_id": "PairA",
        "molecule_chembl_id": "CHEMBL1",
        "canonical_smiles": "CCO",
        "standard_inchi_key": "INCHI1",
        "murcko_scaffold": "c1ccccc1",
        "candidate_status": "exact_matched_ratio",
        "comparison_tier": "strict",
        "same_standard_type": True,
        "same_assay_type": True,
        "ratio_exact_human_div_parasite": 12.0,
        "parasite_standard_type": "IC50",
        "human_standard_type": "IC50",
        "parasite_standard_relation": "=",
        "human_standard_relation": "=",
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


def test_filter_exact_core_candidates_keeps_only_strict_ic50_exact_rows():
    df = pd.DataFrame(
        [
            _row(),
            _row(molecule_chembl_id="CHEMBL2", candidate_status="approximate_ratio"),
            _row(molecule_chembl_id="CHEMBL3", comparison_tier="approximate"),
            _row(molecule_chembl_id="CHEMBL4", same_assay_type=False),
            _row(molecule_chembl_id="CHEMBL5", parasite_standard_type="EC50", human_standard_type="EC50"),
            _row(molecule_chembl_id="CHEMBL6", murcko_scaffold=""),
            _row(molecule_chembl_id="CHEMBL7", ratio_exact_human_div_parasite=-1),
        ]
    )

    exact = filter_exact_core_candidates(df)

    assert len(exact) == 1
    assert exact.iloc[0]["molecule_chembl_id"] == "CHEMBL1"
    assert exact.iloc[0]["compound_key"] == "INCHI1"


def test_filter_exact_core_can_allow_ec50_when_requested():
    df = pd.DataFrame([_row(parasite_standard_type="EC50", human_standard_type="EC50")])

    exact = filter_exact_core_candidates(df, allowed_standard_types=("IC50", "EC50"))

    assert len(exact) == 1


def test_aggregate_exact_core_uses_median_log_ratio_and_flags_conflicts():
    df = pd.DataFrame(
        [
            _row(ratio_exact_human_div_parasite=5.0),
            _row(ratio_exact_human_div_parasite=20.0, parasite_activity_chembl_id="ACTP2"),
        ]
    )
    exact = filter_exact_core_candidates(df)

    core = aggregate_exact_core(exact, selectivity_threshold=10.0)

    assert len(core) == 1
    row = core.iloc[0]
    assert row["label_status"] == "exact_core_conflicting"
    assert row["is_trainable_exact_core"] == 0
    assert row["exact_evidence_count"] == 2
    assert math.isclose(row["selectivity_ratio_human_div_parasite"], 10.0)


def test_aggregate_exact_core_keeps_consistent_labels_trainable():
    df = pd.DataFrame(
        [
            _row(ratio_exact_human_div_parasite=20.0),
            _row(ratio_exact_human_div_parasite=40.0, parasite_activity_chembl_id="ACTP2"),
        ]
    )
    exact = filter_exact_core_candidates(df)

    core = aggregate_exact_core(exact, selectivity_threshold=10.0)

    assert len(core) == 1
    row = core.iloc[0]
    assert row["label_status"] == "exact_core_consistent"
    assert row["is_trainable_exact_core"] == 1
    assert row["is_selective"] == 1
    assert math.isclose(row["selectivity_ratio_human_div_parasite"], math.sqrt(20.0 * 40.0))
