"""Tests for the canonical parasite target manifest pipeline."""

import sys
import unittest
from pathlib import Path

import pytest

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kira.data.target_manifest import (  # noqa: E402
    STATUS_RESOLVED,
    STATUS_UNIQUE,
    VALIDATION_ERROR,
    build_canonical_manifest,
    build_validation_table,
    load_target_datasets,
    resolve_mapping_rule,
    run_target_manifest_pipeline,
)


@pytest.mark.skipif(
    not (Path(__file__).resolve().parents[1] / "data" / "processed" / "schisto_parasite_targets.csv").exists(),
    reason="Requires processed data files (data/processed/*.csv) that are not tracked in git. Run the pipeline locally to regenerate.",
)
class TargetManifestPipelineTest(unittest.TestCase):
    """Verify the manifest pipeline because cross-disease mappings need regression coverage."""

    def test_exact_rule_precedence_uses_canonical_pdeb4_mapping(self) -> None:
        """Ensure exact-ID rules win because older scripts disagree on PDE-family human mappings."""

        rule, match_source = resolve_mapping_rule("CHEMBL2010636", "class 1 phosphodiesterase pdeb1")

        self.assertEqual(match_source, "exact_target_id")
        self.assertEqual(rule.human_target_chembl_id, "CHEMBL275")
        self.assertEqual(rule.human_target_name, "Phosphodiesterase 4B")

    def test_manifest_covers_all_observed_targets(self) -> None:
        """Ensure every observed parasite target is represented because omissions would hide mapping gaps."""

        targets = load_target_datasets()
        manifest = build_canonical_manifest(targets)

        self.assertEqual(len(manifest), len(targets))
        self.assertEqual(
            {row["parasite_target_chembl_id"] for row in manifest},
            {row["parasite_target_chembl_id"] for row in targets},
        )

    def test_manifest_contains_expected_resolved_and_unique_targets(self) -> None:
        """Check key curated targets because they anchor all downstream selectivity analyses."""

        targets = load_target_datasets()
        manifest = {row["parasite_target_chembl_id"]: row for row in build_canonical_manifest(targets)}

        self.assertEqual(manifest["CHEMBL3797017"]["mapping_status"], STATUS_RESOLVED)
        self.assertEqual(manifest["CHEMBL3797017"]["human_target_chembl_id"], "CHEMBL3192")
        self.assertEqual(manifest["CHEMBL1837"]["mapping_status"], STATUS_UNIQUE)
        self.assertEqual(manifest["CHEMBL1837"]["human_target_chembl_id"], "")

    def test_validation_has_no_errors_on_current_repo_data(self) -> None:
        """Validate repo datasets end to end because current artifacts should agree with the canonical manifest."""

        targets = load_target_datasets()
        manifest = build_canonical_manifest(targets)
        validation = build_validation_table(manifest, targets)

        error_count = sum(1 for row in validation if row["validation_status"] == VALIDATION_ERROR)
        self.assertEqual(error_count, 0)

    def test_pipeline_writes_expected_outputs(self) -> None:
        """Ensure the end-to-end runner emits artifacts because the manifest must be reproducible from one command."""

        tmp_path = Path(self._testMethodName)
        tmp_root = Path.cwd() / ".tmp_target_manifest_test"
        tmp_root.mkdir(exist_ok=True)
        output_dir = tmp_root / tmp_path
        output_dir.mkdir(exist_ok=True)

        results = run_target_manifest_pipeline(output_dir=output_dir)

        self.assertTrue(results["manifest_path"].exists())
        self.assertTrue(results["validation_path"].exists())
        self.assertTrue(results["report_path"].exists())


if __name__ == "__main__":
    unittest.main()
