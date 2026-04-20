from __future__ import annotations

from pathlib import Path

import pytest

from kira.experiments.selectivity_v5_expand_data import (
    ExpansionPolicy,
    PairSpec,
    TargetMetadata,
    assay_passes_policy,
    build_pair_candidates,
    candidate_status,
    comparison_tier,
    data_validity_is_acceptable,
    get_activity_interval,
    ratio_bounds,
    read_pair_specs,
)


def _activity(
    *,
    pair_id: str = "TbPDEB1",
    side: str,
    molecule_chembl_id: str = "CHEMBL1",
    standard_type: str = "IC50",
    standard_family: str = "inhibition",
    standard_relation: str = "=",
    standard_value_nM: float = 50.0,
    assay_type: str = "B",
    activity_chembl_id: str | None = None,
    assay_chembl_id: str | None = None,
) -> dict:
    return {
        "pair_id": pair_id,
        "side": side,
        "target_chembl_id": "CHEMBL_PARASITE" if side == "parasite" else "CHEMBL_HUMAN",
        "target_pref_name": "parasite target" if side == "parasite" else "human target",
        "target_organism": "parasite" if side == "parasite" else "Homo sapiens",
        "target_type": "SINGLE PROTEIN",
        "target_accession": "P12345" if side == "parasite" else "Q54321",
        "molecule_chembl_id": molecule_chembl_id,
        "parent_molecule_chembl_id": molecule_chembl_id,
        "canonical_smiles": "CCO",
        "standard_inchi_key": "TEST-INCHIKEY",
        "compound_pref_name": None,
        "activity_chembl_id": activity_chembl_id or f"ACT_{side}",
        "assay_chembl_id": assay_chembl_id or f"ASSAY_{side}",
        "assay_type": assay_type,
        "assay_type_label": "binding" if assay_type == "B" else "functional",
        "assay_description": None,
        "assay_cell_type": None,
        "assay_organism": "parasite" if side == "parasite" else "Homo sapiens",
        "assay_tax_id": None,
        "document_chembl_id": f"DOC_{side}",
        "src_id": None,
        "source_description": None,
        "standard_type": standard_type,
        "standard_family": standard_family,
        "standard_relation": standard_relation,
        "standard_value_nM": standard_value_nM,
        "standard_units": "nM",
        "pchembl_value": None,
        "activity_comment": None,
        "data_validity_comment": None,
        "data_validity_description": None,
        "confidence_score": 9,
        "relationship_type": "D",
        "variant_sequence_accession": None,
        "variant_sequence_mutation": None,
        "measurement_interval_kind": "exact" if standard_relation == "=" else "bound",
    }



def test_get_activity_interval_behaves_as_expected() -> None:
    assert get_activity_interval("=", 50.0).lower_nM == 50.0
    assert get_activity_interval("=", 50.0).upper_nM == 50.0
    assert get_activity_interval("<", 50.0).lower_nM == 0.0
    assert get_activity_interval("<", 50.0).upper_nM == 50.0
    assert get_activity_interval(">", 50.0).lower_nM == 50.0
    assert get_activity_interval(">", 50.0).upper_nM is None



def test_ratio_bounds_exact_case() -> None:
    lower, upper = ratio_bounds("=", 50.0, "=", 500.0)
    assert lower == pytest.approx(10.0)
    assert upper == pytest.approx(10.0)



def test_ratio_bounds_lower_bound_case() -> None:
    lower, upper = ratio_bounds("=", 50.0, ">", 500.0)
    assert lower == pytest.approx(10.0)
    assert upper is None



def test_candidate_status_exact_strict() -> None:
    status, exact_ratio, lower, upper = candidate_status(
        match_tier="strict",
        parasite_relation="=",
        parasite_value_nM=25.0,
        human_relation="=",
        human_value_nM=250.0,
    )
    assert status == "exact_matched_ratio"
    assert exact_ratio == pytest.approx(10.0)
    assert lower == pytest.approx(10.0)
    assert upper == pytest.approx(10.0)



def test_candidate_status_exact_approximate() -> None:
    status, exact_ratio, lower, upper = candidate_status(
        match_tier="approximate",
        parasite_relation="=",
        parasite_value_nM=25.0,
        human_relation="=",
        human_value_nM=250.0,
    )
    assert status == "approximate_ratio"
    assert exact_ratio == pytest.approx(10.0)
    assert lower == pytest.approx(10.0)
    assert upper == pytest.approx(10.0)



def test_comparison_tier_strict_vs_approximate_vs_unmatched() -> None:
    parasite = _activity(side="parasite", standard_type="IC50", standard_family="inhibition", assay_type="B")
    human_strict = _activity(side="human", standard_type="IC50", standard_family="inhibition", assay_type="B")
    human_approx = _activity(side="human", standard_type="IC50", standard_family="inhibition", assay_type="F")
    human_unmatched = _activity(side="human", standard_type="EC50", standard_family="functional", assay_type="F")

    assert comparison_tier(parasite, human_strict) == "strict"
    assert comparison_tier(parasite, human_approx) == "approximate"
    assert comparison_tier(parasite, human_unmatched) == "unmatched"



def test_build_pair_candidates_yields_exact_and_single_side_rows() -> None:
    pair_spec = PairSpec(
        pair_id="TbPDEB1",
        parasite_label="Trypanosoma brucei PDEB1",
        human_label="Human PDE4",
        parasite_target="CHEMBL_PARASITE",
        human_target="CHEMBL_HUMAN",
    )
    parasite_rows = [
        _activity(side="parasite", molecule_chembl_id="CHEMBL1", standard_value_nM=50.0),
        _activity(side="parasite", molecule_chembl_id="CHEMBL2", standard_value_nM=40.0),
    ]
    human_rows = [
        _activity(side="human", molecule_chembl_id="CHEMBL1", standard_value_nM=500.0),
    ]

    candidates = build_pair_candidates(pair_spec, parasite_rows, human_rows)
    by_status = {row["molecule_chembl_id"]: row["candidate_status"] for row in candidates}

    assert by_status["CHEMBL1"] == "exact_matched_ratio"
    assert by_status["CHEMBL2"] == "single_side_only"



def test_build_pair_candidates_unmatched_comparable_for_family_mismatch() -> None:
    pair_spec = PairSpec(
        pair_id="TbPDEB1",
        parasite_label="Trypanosoma brucei PDEB1",
        human_label="Human PDE4",
        parasite_target="CHEMBL_PARASITE",
        human_target="CHEMBL_HUMAN",
    )
    parasite_rows = [_activity(side="parasite", standard_type="KI", standard_family="affinity", assay_type="B")]
    human_rows = [_activity(side="human", standard_type="EC50", standard_family="functional", assay_type="F")]

    candidates = build_pair_candidates(pair_spec, parasite_rows, human_rows)

    assert len(candidates) == 1
    assert candidates[0]["candidate_status"] == "unmatched_comparable"



def test_assay_passes_policy_requires_direct_single_protein_non_variant() -> None:
    policy = ExpansionPolicy()
    target_meta = TargetMetadata(
        target_chembl_id="CHEMBL1",
        pref_name="Target",
        organism="Homo sapiens",
        target_type="SINGLE PROTEIN",
        accession="P12345",
    )
    assay_row = {
        "assay_type": "B",
        "confidence_score": 9,
        "relationship_type": "D",
        "variant_sequence_accession": None,
        "variant_sequence_mutation": None,
        "variant_id": None,
    }
    assert assay_passes_policy(assay_row, target_meta, policy) is True

    assay_variant = dict(assay_row)
    assay_variant["variant_sequence_mutation"] = "V600E"
    assert assay_passes_policy(assay_variant, target_meta, policy) is False

    assay_indirect = dict(assay_row)
    assay_indirect["relationship_type"] = "H"
    assert assay_passes_policy(assay_indirect, target_meta, policy) is False



def test_read_pair_specs_round_trip(tmp_path: Path) -> None:
    csv_path = tmp_path / "pairs.csv"
    csv_path.write_text(
        "pair_id,parasite_label,human_label,parasite_target,human_target,notes\n"
        "TbPDEB1,Trypanosome PDEB1,Human PDE4,CHEMBL_PARASITE,CHEMBL_HUMAN,test note\n",
        encoding="utf-8",
    )

    specs = read_pair_specs(csv_path)

    assert len(specs) == 1
    assert specs[0].pair_id == "TbPDEB1"
    assert specs[0].notes == "test note"


def test_data_validity_policy_filters_flagged_values() -> None:
    policy = ExpansionPolicy()
    assert data_validity_is_acceptable(None, policy) is True
    assert data_validity_is_acceptable("", policy) is True
    assert data_validity_is_acceptable("Manually validated", policy) is True
    assert data_validity_is_acceptable("Outside typical range", policy) is False

    relaxed = ExpansionPolicy(exclude_flagged_data_validity=False)
    assert data_validity_is_acceptable("Outside typical range", relaxed) is True

def test_reference_pair_config_has_expected_target_identifiers() -> None:
    config_path = Path("data/reference/selectivity_v5_target_pairs.csv")
    specs = read_pair_specs(config_path)

    assert len(specs) == 6

    expected = {
        "SmDHODH": ("CHEMBL4523950", "chembl_id", "CHEMBL1966", "chembl_id"),
        "SmHDAC8": ("CHEMBL3797017", "chembl_id", "CHEMBL3192", "chembl_id"),
        "TbCathB": ("Q6R7Z5", "accession", "CHEMBL4072", "chembl_id"),
        "TbPDEB1": ("CHEMBL2010636", "chembl_id", "CHEMBL275", "chembl_id"),
        "LmPTR1": ("Q01782", "accession", "CHEMBL202", "chembl_id"),
        "LmDHFR": ("CHEMBL4614", "chembl_id", "CHEMBL202", "chembl_id"),
    }

    got = {
        spec.pair_id: (
            spec.parasite_target,
            spec.parasite_target_namespace,
            spec.human_target,
            spec.human_target_namespace,
        )
        for spec in specs
    }

    assert got == expected
