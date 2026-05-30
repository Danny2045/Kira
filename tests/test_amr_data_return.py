from __future__ import annotations

import csv
import hashlib
import json
from io import StringIO
from pathlib import Path

import pytest

from kira.amr import (
    DATA_STATUS_REAL,
    SYNTHETIC_BANNER,
    TEMPLATE_COLUMNS,
    audit_amr_csv,
    audit_ast_completeness,
    classify_data_status,
    load_amr_csv,
    load_amr_csv_with_provenance,
    make_amr_csv_report,
    make_data_return_template_rows,
    make_markdown_report,
    report_to_dict,
    required_aggregate_columns,
    required_isolate_columns,
    write_data_return_template,
    write_markdown_report,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_CSV = REPO_ROOT / "examples" / "rwanda_amr_ast_example.csv"
TEMPLATE_CSV = REPO_ROOT / "examples" / "rwanda_amr_ast_template.csv"
EXPECTED_EXAMPLE_RECORDS = 4

FORBIDDEN_OUTPUT_TERMS = (
    "clinical recommendation",
    "prescribe",
    "solved AMR",
    "wet-lab validated",
    "clinical cure",
    "facility ranking",
    "public-health outcome claim",
    "partnership with Biohub",
    "partnership with Arc",
    "partnership with Ginkgo",
)


def test_template_columns_include_all_required_isolate_fields() -> None:
    with TEMPLATE_CSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        template_columns = tuple(reader.fieldnames or ())

    assert set(required_isolate_columns()).issubset(template_columns)
    assert tuple(TEMPLATE_COLUMNS) == template_columns
    assert set(required_aggregate_columns()).issubset(template_columns)
    assert "resistant_count" in required_aggregate_columns()
    assert "susceptible_count" in required_aggregate_columns()


def test_example_csv_loads_expected_number_of_records() -> None:
    rows = load_amr_csv(EXAMPLE_CSV)
    stream_rows = load_amr_csv(StringIO(EXAMPLE_CSV.read_text(encoding="utf-8")))

    assert len(rows) == EXPECTED_EXAMPLE_RECORDS
    assert stream_rows == rows
    assert all(row["synthetic_data_notice"].startswith("SYNTHETIC EXAMPLE") for row in rows)


def test_audit_amr_csv_matches_ast_completeness_report_shape() -> None:
    rows = load_amr_csv(EXAMPLE_CSV)
    csv_report = audit_amr_csv(EXAMPLE_CSV)
    direct_report = audit_ast_completeness(rows)

    assert report_to_dict(csv_report) == report_to_dict(direct_report)
    assert set(report_to_dict(csv_report)) == set(report_to_dict(direct_report))


def test_make_markdown_report_is_deterministic() -> None:
    report = audit_amr_csv(EXAMPLE_CSV)

    assert make_markdown_report(report) == make_markdown_report(report)


def test_markdown_report_includes_ready_and_missing_field_sections() -> None:
    markdown = make_markdown_report(audit_amr_csv(EXAMPLE_CSV))

    assert "## Benchmark-Ready Rows" in markdown
    assert "## Missing Field Counts" in markdown
    assert "## Incomplete Rows" in markdown
    assert "## Repair Actions By Missing Field" in markdown
    assert "`ast_result`" in markdown
    assert "`resistant_count_or_susceptible_count`" in markdown


def test_write_data_return_template_writes_usable_csv(tmp_path: Path) -> None:
    output_path = write_data_return_template(tmp_path / "rwanda_amr_ast_template.csv")
    rows = load_amr_csv(output_path)

    assert output_path.exists()
    assert rows == make_data_return_template_rows()
    assert set(required_isolate_columns()).issubset(rows[0])


def test_write_markdown_report_writes_report(tmp_path: Path) -> None:
    report = audit_amr_csv(EXAMPLE_CSV)
    output_path = write_markdown_report(report, tmp_path / "report.md")

    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == make_markdown_report(report)


def test_forbidden_claims_do_not_appear_in_generated_output() -> None:
    report = audit_amr_csv(EXAMPLE_CSV)
    generated_output = json.dumps(
        {
            "markdown": make_markdown_report(report),
            "report": report_to_dict(report),
            "template_rows": make_data_return_template_rows(),
        },
        sort_keys=True,
    ).lower()

    for claim in FORBIDDEN_OUTPUT_TERMS:
        assert claim.lower() not in generated_output


def test_full_report_is_json_serializable_through_report_to_dict() -> None:
    payload = report_to_dict(audit_amr_csv(EXAMPLE_CSV))

    assert json.loads(json.dumps(payload, sort_keys=True)) == payload
    assert payload["total_records"] == EXPECTED_EXAMPLE_RECORDS
    assert "repair_tickets" in payload


AMR_CSV_HEADER = (
    "synthetic_data_notice,record_type,facility_id,reporting_period,organism,"
    "specimen_source,antibiotic,ast_method,ast_result,breakpoint_version,qc_status,isolate_id"
)


def _isolate_row(notice: str, isolate_id: str) -> str:
    return (
        f"{notice},isolate_level,fake-facility-001,2026-04,Escherichia coli,urine,"
        f"ceftriaxone,disk diffusion,resistant,CLSI M100 2026,pass,{isolate_id}"
    )


def _csv_text(*rows: str) -> str:
    return "\n".join((AMR_CSV_HEADER, *rows)) + "\n"


def test_synthetic_csv_report_has_banner_and_provenance_block_with_unchanged_counts() -> None:
    report = audit_amr_csv(EXAMPLE_CSV)
    base = make_markdown_report(report)  # no provenance: numeric body only
    full = make_amr_csv_report(EXAMPLE_CSV)  # provenance-stamped

    # synthetic banner is the very first line; provenance block follows the title
    assert full.splitlines()[0] == SYNTHETIC_BANNER
    assert "## Data Provenance" in full
    assert "| Data status | `SYNTHETIC` |" in full

    # the no-provenance path carries neither banner nor provenance block
    assert SYNTHETIC_BANNER not in base
    assert "## Data Provenance" not in base

    # every numeric section from '## Summary' onward is byte-identical
    assert base.split("## Summary", 1)[1] == full.split("## Summary", 1)[1]


def test_synthetic_csv_report_is_byte_identical_for_same_input() -> None:
    assert make_amr_csv_report(EXAMPLE_CSV) == make_amr_csv_report(EXAMPLE_CSV)


def test_blank_notice_file_classifies_real_with_no_banner() -> None:
    csv_text = _csv_text(_isolate_row("", "iso-1"), _isolate_row("", "iso-2"))

    assert classify_data_status(load_amr_csv(StringIO(csv_text))) == DATA_STATUS_REAL

    report = make_amr_csv_report(StringIO(csv_text))
    assert SYNTHETIC_BANNER not in report
    assert "## Data Provenance" in report
    assert "| Data status | `REAL` |" in report


def test_mixed_provenance_file_is_refused_by_report_and_audit() -> None:
    csv_text = _csv_text(_isolate_row("SYNTHETIC EXAMPLE", "iso-1"), _isolate_row("", "iso-2"))

    with pytest.raises(
        ValueError, match="mixed-provenance file: 1 rows marked synthetic, 1 rows unmarked"
    ):
        make_amr_csv_report(StringIO(csv_text))

    with pytest.raises(ValueError, match="mixed-provenance file"):
        audit_amr_csv(StringIO(csv_text))


def test_inconsistent_notice_file_is_refused() -> None:
    csv_text = _csv_text(
        _isolate_row("SYNTHETIC EXAMPLE A", "iso-1"),
        _isolate_row("SYNTHETIC EXAMPLE B", "iso-2"),
    )

    with pytest.raises(ValueError, match="inconsistent synthetic_data_notice column"):
        make_amr_csv_report(StringIO(csv_text))


def test_empty_file_reports_no_data_audited_without_banner() -> None:
    report = make_amr_csv_report(StringIO(AMR_CSV_HEADER + "\n"))

    assert SYNTHETIC_BANNER not in report
    assert "| Data status | `EMPTY` |" in report
    assert "No data audited" in report


def test_provenance_block_carries_source_and_stable_content_hash() -> None:
    _rows_first, prov_first = load_amr_csv_with_provenance(EXAMPLE_CSV)
    _rows_second, prov_second = load_amr_csv_with_provenance(EXAMPLE_CSV)

    # same input -> identical source and content hash (deterministic provenance)
    assert prov_first.source == prov_second.source == str(EXAMPLE_CSV)
    assert prov_first.content_sha256 == prov_second.content_sha256
    assert prov_first.content_sha256 == hashlib.sha256(EXAMPLE_CSV.read_bytes()).hexdigest()
    assert prov_first.n_rows == EXPECTED_EXAMPLE_RECORDS

    report = make_amr_csv_report(EXAMPLE_CSV)
    assert prov_first.source in report
    assert prov_first.content_sha256 in report
    assert f"| Rows audited | {prov_first.n_rows} |" in report
