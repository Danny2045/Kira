from __future__ import annotations

import csv
import json
from io import StringIO
from pathlib import Path

from kira.amr import (
    TEMPLATE_COLUMNS,
    audit_amr_csv,
    audit_ast_completeness,
    load_amr_csv,
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
