from __future__ import annotations

import subprocess
import sys

import pytest

from kira.contrast.schemas import (
    CONTRAST_SPEC_REQUIRED_FIELDS,
    DATA_RETURN_REQUIRED_FIELDS,
    ContrastSpec,
    EvidenceStatus,
    make_data_return_schema,
    validate_contrast_spec,
)


def _contrast_payload(**overrides):
    payload = {
        "contrast_id": "ntd-selectivity::compound-a::target-pair-1",
        "domain": "ntd_selectivity",
        "intervention_type": "small_molecule",
        "intervention_id": "compound-a",
        "desired_context": "parasite target inhibition",
        "control_context": "human target safety comparator",
        "readout_type": "IC50 ratio human_div_parasite",
        "readout_units": "fold",
        "evidence_status": EvidenceStatus.SINGLE_SIDE_ONLY,
        "known_side": "parasite",
        "missing_side": "human",
        "label": "benchmark_repair_candidate",
        "uncertainty_note": "Human comparator measurement is missing.",
        "provenance": "synthetic unit-test fixture",
    }
    payload.update(overrides)
    return payload


def test_valid_contrast_spec_passes_validation() -> None:
    spec = ContrastSpec(**_contrast_payload())

    validated = validate_contrast_spec(spec)

    assert validated == spec


def test_missing_required_fields_fail_validation() -> None:
    payload = _contrast_payload()
    del payload["contrast_id"]

    with pytest.raises(ValueError, match="contrast_id"):
        validate_contrast_spec(payload)


def test_data_return_schema_contains_required_fields() -> None:
    schema = make_data_return_schema(ContrastSpec(**_contrast_payload()))

    assert set(DATA_RETURN_REQUIRED_FIELDS).issubset(schema.required_fields)
    assert schema.contrast_id == "ntd-selectivity::compound-a::target-pair-1"
    assert schema.readout_type == "IC50 ratio human_div_parasite"


def test_ntd_selectivity_can_be_represented_without_rdkit() -> None:
    spec = validate_contrast_spec(
        _contrast_payload(
            contrast_id="ntd-selectivity::CHEMBL155771::TcCYP51-vs-human-CYP51",
            desired_context="T. cruzi CYP51 inhibition",
            control_context="human CYP51 inhibition",
            provenance="v5 selectivity candidate table",
        )
    )

    assert spec.domain == "ntd_selectivity"
    assert spec.known_side == "parasite"
    assert spec.missing_side == "human"

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import kira.contrast; raise SystemExit(1 if 'rdkit' in sys.modules else 0)",
        ],
        check=False,
    )
    assert result.returncode == 0


def test_regeneration_crispr_style_contrast_can_be_represented_without_claims() -> None:
    spec = validate_contrast_spec(
        _contrast_payload(
            contrast_id="regeneration-crispr::guide-set-a::scar-vs-regenerative-state",
            domain="regeneration_crispr",
            intervention_type="crispr_perturbation",
            intervention_id="guide-set-a",
            desired_context="regenerative cell-state marker panel",
            control_context="fibrotic or scar-like marker panel",
            readout_type="marker panel delta",
            readout_units="normalized expression score",
            evidence_status=EvidenceStatus.PROPOSED,
            known_side="none",
            missing_side="both",
            label="hypothesis_only",
            uncertainty_note=(
                "Schema fixture only; no solved-domain, efficacy, or safety claim is implied."
            ),
            provenance="contrast core unit-test fixture",
        )
    )

    assert spec.domain == "regeneration_crispr"
    assert spec.evidence_status == EvidenceStatus.PROPOSED
    assert "solved-domain" in spec.uncertainty_note
    assert set(CONTRAST_SPEC_REQUIRED_FIELDS).issubset(spec.__dataclass_fields__)
