"""Build the canonical parasite target manifest and validate mapping datasets."""

from pathlib import Path
import sys

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kira.data.target_manifest import run_target_manifest_pipeline


def main() -> None:
    """Run the pipeline because this repo prefers reproducible single-command workflows."""

    results = run_target_manifest_pipeline()
    print("Canonical target-manifest pipeline complete.")
    print(f"Manifest: {results['manifest_path']}")
    print(f"Validation: {results['validation_path']}")
    print(f"Report: {results['report_path']}")


if __name__ == "__main__":
    main()
