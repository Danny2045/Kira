"""Canonical parasite target manifest and validation pipeline.

This module centralizes parasite-to-human target mapping rules so the repo
can build one reproducible manifest instead of relying on divergent script-
local dictionaries. It also validates existing selectivity datasets against
that manifest so mapping assumptions stay explicit and testable.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Iterable

BASE_DIR = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_DIR = BASE_DIR / "data" / "processed"
DEFAULT_REPORT_PATH = DEFAULT_OUTPUT_DIR / "canonical_target_manifest_report.txt"
DEFAULT_MANIFEST_PATH = DEFAULT_OUTPUT_DIR / "canonical_target_manifest.csv"
DEFAULT_VALIDATION_PATH = DEFAULT_OUTPUT_DIR / "canonical_target_manifest_validation.csv"

STATUS_RESOLVED = "resolved"
STATUS_PROPOSED = "proposed"
STATUS_UNIQUE = "parasite_unique"
STATUS_UNRESOLVED = "unresolved"

TYPE_DIRECT = "direct_orthologue"
TYPE_FUNCTIONAL = "functional_equivalent"
TYPE_UNIQUE = "parasite_unique"
TYPE_NONE = "no_mapping"

MATCH_EXACT = "exact_target_id"
MATCH_PATTERN = "name_pattern"
MATCH_NONE = "no_match"

VALIDATION_OK = "ok"
VALIDATION_WARNING = "warning"
VALIDATION_ERROR = "error"

TARGET_DATASETS = (
    {
        "disease": "Schistosomiasis",
        "dataset_key": "schisto_targets",
        "path": BASE_DIR / "data" / "processed" / "schisto_parasite_targets.csv",
    },
    {
        "disease": "Trypanosomiasis",
        "dataset_key": "tryp_targets",
        "path": BASE_DIR / "data" / "trypanosoma" / "tryp_targets.csv",
    },
    {
        "disease": "Leishmaniasis",
        "dataset_key": "leish_targets",
        "path": BASE_DIR / "data" / "leishmania" / "leish_targets.csv",
    },
)

SELECTIVITY_DATASETS = (
    {
        "disease": "Schistosomiasis",
        "dataset_key": "schisto_selectivity",
        "path": BASE_DIR / "data" / "processed" / "kira_selectivity_analysis.csv",
        "parasite_target_column": "parasite_target",
        "human_target_id_column": None,
        "human_target_name_column": None,
    },
    {
        "disease": "Trypanosomiasis",
        "dataset_key": "tryp_selectivity",
        "path": BASE_DIR / "data" / "trypanosoma" / "tryp_selectivity.csv",
        "parasite_target_column": "parasite_target",
        "human_target_id_column": "human_target_id",
        "human_target_name_column": None,
    },
    {
        "disease": "Leishmaniasis",
        "dataset_key": "leish_selectivity",
        "path": BASE_DIR / "data" / "leishmania" / "leish_selectivity.csv",
        "parasite_target_column": "parasite_target",
        "human_target_id_column": None,
        "human_target_name_column": "human_target",
    },
)

MANIFEST_FIELDNAMES = [
    "disease",
    "dataset_key",
    "parasite_organism",
    "parasite_target_name",
    "parasite_target_chembl_id",
    "normalized_target_name",
    "mapping_status",
    "mapping_type",
    "human_target_chembl_id",
    "human_target_name",
    "human_gene_symbol",
    "evidence_level",
    "notes",
    "match_source",
]

VALIDATION_FIELDNAMES = [
    "validation_scope",
    "dataset_key",
    "disease",
    "parasite_target_chembl_id",
    "parasite_target_name",
    "validation_status",
    "issue_code",
    "details",
    "expected_human_target_chembl_id",
    "observed_human_target_chembl_id",
    "expected_human_target_name",
    "observed_human_target_name",
]


@dataclass(frozen=True)
class MappingRule:
    """Describe HOW a parasite target should map and WHY the rule exists."""

    mapping_status: str
    mapping_type: str
    human_target_chembl_id: str | None
    human_target_name: str | None
    human_gene_symbol: str | None
    evidence_level: str
    notes: str


@dataclass(frozen=True)
class ValidationRecord:
    """Capture WHAT was validated and WHY any mismatch matters."""

    validation_scope: str
    dataset_key: str
    disease: str
    parasite_target_chembl_id: str
    parasite_target_name: str
    validation_status: str
    issue_code: str
    details: str
    expected_human_target_chembl_id: str
    observed_human_target_chembl_id: str
    expected_human_target_name: str
    observed_human_target_name: str


EXACT_TARGET_RULES = {
    "CHEMBL6110": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_FUNCTIONAL,
        human_target_chembl_id="CHEMBL3952",
        human_target_name="Thioredoxin reductase 1",
        human_gene_symbol="TXNRD1",
        evidence_level="legacy_curated",
        notes="SmTGR is a parasite fusion enzyme; human TXNRD1 is the closest tractable comparator.",
    ),
    "CHEMBL3797017": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_DIRECT,
        human_target_chembl_id="CHEMBL3192",
        human_target_name="Histone deacetylase 8",
        human_gene_symbol="HDAC8",
        evidence_level="legacy_curated",
        notes="Existing schistosomiasis selectivity analysis uses human HDAC8 as the orthologue.",
    ),
    "CHEMBL4523950": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_DIRECT,
        human_target_chembl_id="CHEMBL1966",
        human_target_name="Dihydroorotate dehydrogenase",
        human_gene_symbol="DHODH",
        evidence_level="legacy_curated",
        notes="Existing schistosomiasis selectivity analysis uses human DHODH as the orthologue.",
    ),
    "CHEMBL5832": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_DIRECT,
        human_target_chembl_id="CHEMBL3837",
        human_target_name="Cathepsin B",
        human_gene_symbol="CTSB",
        evidence_level="validated_by_dataset",
        notes="Trypanosoma selectivity dataset already records human cathepsin B measurements.",
    ),
    "CHEMBL4188": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_DIRECT,
        human_target_chembl_id="CHEMBL3837",
        human_target_name="Cathepsin B",
        human_gene_symbol="CTSB",
        evidence_level="legacy_curated",
        notes="Rhodesain is a cathepsin-L-like cysteine protease; current platform compares against CTSB-family human proteases.",
    ),
    "CHEMBL2010636": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_FUNCTIONAL,
        human_target_chembl_id="CHEMBL275",
        human_target_name="Phosphodiesterase 4B",
        human_gene_symbol="PDE4B",
        evidence_level="canonicalized_from_cross_disease_model",
        notes="Canonical choice resolves the older CHEMBL2094253 family-level mapping in favor of explicit PDE4B.",
    ),
    "CHEMBL2169733": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_FUNCTIONAL,
        human_target_chembl_id="CHEMBL275",
        human_target_name="Phosphodiesterase 4B",
        human_gene_symbol="PDE4B",
        evidence_level="validated_by_dataset",
        notes="Leishmania selectivity dataset uses PDE4B as the human comparator.",
    ),
    "CHEMBL1837": MappingRule(
        mapping_status=STATUS_UNIQUE,
        mapping_type=TYPE_UNIQUE,
        human_target_chembl_id=None,
        human_target_name=None,
        human_gene_symbol=None,
        evidence_level="legacy_curated",
        notes="Trypanothione reductase is unique to trypanosomatids and intentionally has no direct human orthologue mapping.",
    ),
    "CHEMBL2176840": MappingRule(
        mapping_status=STATUS_UNIQUE,
        mapping_type=TYPE_UNIQUE,
        human_target_chembl_id=None,
        human_target_name=None,
        human_gene_symbol=None,
        evidence_level="legacy_curated",
        notes="Trypanothione reductase is unique to trypanosomatids and intentionally has no direct human orthologue mapping.",
    ),
    "CHEMBL1944501": MappingRule(
        mapping_status=STATUS_UNIQUE,
        mapping_type=TYPE_UNIQUE,
        human_target_chembl_id=None,
        human_target_name=None,
        human_gene_symbol=None,
        evidence_level="legacy_curated",
        notes="Trypanothione reductase is unique to trypanosomatids and intentionally has no direct human orthologue mapping.",
    ),
    "CHEMBL1944500": MappingRule(
        mapping_status=STATUS_UNIQUE,
        mapping_type=TYPE_UNIQUE,
        human_target_chembl_id=None,
        human_target_name=None,
        human_gene_symbol=None,
        evidence_level="legacy_curated",
        notes="Trypanothione reductase is unique to trypanosomatids and intentionally has no direct human orthologue mapping.",
    ),
    "CHEMBL4614": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_FUNCTIONAL,
        human_target_chembl_id="CHEMBL202",
        human_target_name="Dihydrofolate reductase",
        human_gene_symbol="DHFR",
        evidence_level="validated_by_dataset",
        notes="Leishmania DHFR-TS is bifunctional; human DHFR is the established functional comparator in the dataset.",
    ),
    "CHEMBL6194": MappingRule(
        mapping_status=STATUS_RESOLVED,
        mapping_type=TYPE_FUNCTIONAL,
        human_target_chembl_id="CHEMBL202",
        human_target_name="Dihydrofolate reductase",
        human_gene_symbol="DHFR",
        evidence_level="validated_by_dataset",
        notes="PTR1 lacks a direct human orthologue; the platform uses human DHFR as the closest functional comparator.",
    ),
}

NAME_PATTERN_RULES = (
    (
        "trypanothione reductase",
        MappingRule(
            mapping_status=STATUS_UNIQUE,
            mapping_type=TYPE_UNIQUE,
            human_target_chembl_id=None,
            human_target_name=None,
            human_gene_symbol=None,
            evidence_level="pattern_rule",
            notes="Trypanothione reductase is parasite-specific across kinetoplastids.",
        ),
    ),
    (
        "histone deacetylase 8",
        MappingRule(
            mapping_status=STATUS_RESOLVED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL3192",
            human_target_name="Histone deacetylase 8",
            human_gene_symbol="HDAC8",
            evidence_level="pattern_rule",
            notes="HDAC8 naming is stable enough to use an exact-name fallback.",
        ),
    ),
    (
        "dihydroorotate dehydrogenase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL1966",
            human_target_name="Dihydroorotate dehydrogenase",
            human_gene_symbol="DHODH",
            evidence_level="pattern_rule",
            notes="Pattern fallback maps DHODH-family targets to human DHODH when no exact-ID curation exists.",
        ),
    ),
    (
        "cathepsin b",
        MappingRule(
            mapping_status=STATUS_RESOLVED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL3837",
            human_target_name="Cathepsin B",
            human_gene_symbol="CTSB",
            evidence_level="pattern_rule",
            notes="Cathepsin-B-like parasite proteases are compared to human CTSB in current platform analyses.",
        ),
    ),
    (
        "cysteine protease",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL3524",
            human_target_name="Cathepsin L",
            human_gene_symbol="CTSL",
            evidence_level="pattern_rule",
            notes="Generic cysteine proteases are mapped to human CTSL-family comparators pending tighter curation.",
        ),
    ),
    (
        "ornithine decarboxylase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL3920",
            human_target_name="Ornithine decarboxylase",
            human_gene_symbol="ODC1",
            evidence_level="pattern_rule",
            notes="Legacy trypanosome platform maps ODC targets to human ODC1.",
        ),
    ),
    (
        "methionyl trna synthetase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL4523",
            human_target_name="Methionine--tRNA ligase",
            human_gene_symbol="MARS1",
            evidence_level="pattern_rule",
            notes="Legacy platform uses human MetRS as the comparator for parasite MetRS targets.",
        ),
    ),
    (
        "glycylpeptide n tetradecanoyltransferase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL2146302",
            human_target_name="N-myristoyltransferase 1",
            human_gene_symbol="NMT1",
            evidence_level="pattern_rule",
            notes="Leishmania NMT target names differ but map to the same NMT1 comparator.",
        ),
    ),
    (
        "n myristoyltransferase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL2146302",
            human_target_name="N-myristoyltransferase 1",
            human_gene_symbol="NMT1",
            evidence_level="pattern_rule",
            notes="Legacy trypanosome and leishmania platform maps NMT targets to human NMT1.",
        ),
    ),
    (
        "phosphodiesterase pdeb1",
        MappingRule(
            mapping_status=STATUS_RESOLVED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL275",
            human_target_name="Phosphodiesterase 4B",
            human_gene_symbol="PDE4B",
            evidence_level="pattern_rule",
            notes="Cross-disease model and leishmania dataset both converge on PDE4B as the canonical human comparator.",
        ),
    ),
    (
        "phosphodiesterase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL275",
            human_target_name="Phosphodiesterase 4B",
            human_gene_symbol="PDE4B",
            evidence_level="pattern_rule",
            notes="Generic PDE-family parasite targets default to PDE4B until isoform-specific curation is added.",
        ),
    ),
    (
        "sterol 14",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_DIRECT,
            human_target_chembl_id="CHEMBL1978",
            human_target_name="Lanosterol 14-alpha demethylase",
            human_gene_symbol="CYP51A1",
            evidence_level="pattern_rule",
            notes="Legacy leishmania platform maps sterol 14-alpha demethylase to human CYP51A1.",
        ),
    ),
    (
        "topoisomerase",
        MappingRule(
            mapping_status=STATUS_PROPOSED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL1806",
            human_target_name="DNA topoisomerase II alpha",
            human_gene_symbol="TOP2A",
            evidence_level="pattern_rule",
            notes="Legacy platform compares parasite topoisomerase assays to human TOP2A-family activity.",
        ),
    ),
    (
        "pteridine reductase",
        MappingRule(
            mapping_status=STATUS_RESOLVED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL202",
            human_target_name="Dihydrofolate reductase",
            human_gene_symbol="DHFR",
            evidence_level="pattern_rule",
            notes="Existing leishmania dataset already uses human DHFR as the functional comparator for PTR1.",
        ),
    ),
    (
        "dihydrofolate reductase",
        MappingRule(
            mapping_status=STATUS_RESOLVED,
            mapping_type=TYPE_FUNCTIONAL,
            human_target_chembl_id="CHEMBL202",
            human_target_name="Dihydrofolate reductase",
            human_gene_symbol="DHFR",
            evidence_level="pattern_rule",
            notes="DHFR-family parasite targets are paired to human DHFR in existing analyses.",
        ),
    ),
)


def normalize_text(value: str | None) -> str:
    """Normalize free text so fuzzy target-name matching is reproducible."""

    if value is None:
        return ""

    normalized = "".join(character.lower() if character.isalnum() else " " for character in str(value))
    return " ".join(normalized.split())


def log_step(step_name: str, input_desc: str, output_desc: str, start_time: float) -> str:
    """Return a stable pipeline log line so every step records input, output, and runtime."""

    elapsed_seconds = perf_counter() - start_time
    return (
        f"{step_name}: input={input_desc} | output={output_desc} | "
        f"duration_seconds={elapsed_seconds:.3f}"
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows because the pipeline should run in minimal Python environments."""

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    """Write CSV rows because manifest outputs must be easy to inspect and version."""

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def load_target_datasets() -> list[dict[str, str]]:
    """Load all parasite target catalogs because the manifest must cover every observed target."""

    unique_rows: dict[tuple[str, str], dict[str, str]] = {}
    for dataset in TARGET_DATASETS:
        for row in read_csv_rows(dataset["path"]):
            manifest_seed = {
                "disease": dataset["disease"],
                "dataset_key": dataset["dataset_key"],
                "parasite_organism": row.get("organism", ""),
                "parasite_target_name": row.get("pref_name", ""),
                "parasite_target_chembl_id": row.get("target_chembl_id", ""),
                "normalized_target_name": normalize_text(row.get("pref_name")),
            }
            unique_rows[(manifest_seed["disease"], manifest_seed["parasite_target_chembl_id"])] = manifest_seed

    return sorted(
        unique_rows.values(),
        key=lambda row: (row["disease"], row["parasite_target_name"], row["parasite_target_chembl_id"]),
    )


def resolve_mapping_rule(target_chembl_id: str, normalized_target_name: str) -> tuple[MappingRule, str]:
    """Resolve the canonical mapping rule because exact IDs should override name heuristics."""

    if target_chembl_id in EXACT_TARGET_RULES:
        return EXACT_TARGET_RULES[target_chembl_id], MATCH_EXACT

    for pattern, rule in NAME_PATTERN_RULES:
        if pattern in normalized_target_name:
            return rule, MATCH_PATTERN

    return (
        MappingRule(
            mapping_status=STATUS_UNRESOLVED,
            mapping_type=TYPE_NONE,
            human_target_chembl_id=None,
            human_target_name=None,
            human_gene_symbol=None,
            evidence_level="unmapped",
            notes="No canonical parasite-human mapping has been curated for this target yet.",
        ),
        MATCH_NONE,
    )


def build_canonical_manifest(targets: list[dict[str, str]]) -> list[dict[str, str]]:
    """Build the manifest from observed targets so the output stays anchored to repo data."""

    manifest_rows: list[dict[str, str]] = []
    for target in targets:
        rule, match_source = resolve_mapping_rule(
            target_chembl_id=target["parasite_target_chembl_id"],
            normalized_target_name=target["normalized_target_name"],
        )
        manifest_rows.append(
            {
                **target,
                **{
                    key: value if value is not None else ""
                    for key, value in asdict(rule).items()
                },
                "match_source": match_source,
            }
        )

    return sorted(
        manifest_rows,
        key=lambda row: (row["disease"], row["parasite_target_name"], row["parasite_target_chembl_id"]),
    )


def validate_manifest_integrity(manifest: list[dict[str, str]]) -> list[ValidationRecord]:
    """Validate manifest shape because duplicate or missing mappings would break downstream interpretation."""

    records: list[ValidationRecord] = []
    counts = Counter((row["disease"], row["parasite_target_chembl_id"]) for row in manifest)
    for row in manifest:
        if counts[(row["disease"], row["parasite_target_chembl_id"])] > 1:
            records.append(
                ValidationRecord(
                    validation_scope="manifest_integrity",
                    dataset_key=row["dataset_key"],
                    disease=row["disease"],
                    parasite_target_chembl_id=row["parasite_target_chembl_id"],
                    parasite_target_name=row["parasite_target_name"],
                    validation_status=VALIDATION_ERROR,
                    issue_code="duplicate_parasite_target_id",
                    details="Manifest contains duplicate parasite target IDs within the same disease catalog.",
                    expected_human_target_chembl_id="",
                    observed_human_target_chembl_id="",
                    expected_human_target_name="",
                    observed_human_target_name="",
                )
            )
        if not row["mapping_status"]:
            records.append(
                ValidationRecord(
                    validation_scope="manifest_integrity",
                    dataset_key=row["dataset_key"],
                    disease=row["disease"],
                    parasite_target_chembl_id=row["parasite_target_chembl_id"],
                    parasite_target_name=row["parasite_target_name"],
                    validation_status=VALIDATION_ERROR,
                    issue_code="missing_mapping_status",
                    details="Manifest row is missing a mapping status classification.",
                    expected_human_target_chembl_id="",
                    observed_human_target_chembl_id="",
                    expected_human_target_name="",
                    observed_human_target_name="",
                )
            )

    return records


def validate_target_coverage(manifest: list[dict[str, str]], targets: list[dict[str, str]]) -> list[ValidationRecord]:
    """Validate coverage because every observed parasite target must appear in the canonical manifest."""

    records: list[ValidationRecord] = []
    manifest_index = {(row["disease"], row["parasite_target_chembl_id"]) for row in manifest}

    for row in targets:
        key = (row["disease"], row["parasite_target_chembl_id"])
        if key not in manifest_index:
            records.append(
                ValidationRecord(
                    validation_scope="target_coverage",
                    dataset_key=row["dataset_key"],
                    disease=row["disease"],
                    parasite_target_chembl_id=row["parasite_target_chembl_id"],
                    parasite_target_name=row["parasite_target_name"],
                    validation_status=VALIDATION_ERROR,
                    issue_code="target_missing_from_manifest",
                    details="Observed target row is missing from the canonical manifest.",
                    expected_human_target_chembl_id="",
                    observed_human_target_chembl_id="",
                    expected_human_target_name="",
                    observed_human_target_name="",
                )
            )

    return records


def select_manifest_row(manifest: list[dict[str, str]], disease: str, parasite_target_name: str) -> dict[str, str] | None:
    """Locate the canonical manifest row because selectivity files key off parasite target names."""

    normalized_target_name = normalize_text(parasite_target_name)
    matches = [
        row for row in manifest
        if row["disease"] == disease and row["normalized_target_name"] == normalized_target_name
    ]
    if not matches:
        return None
    if len(matches) == 1:
        return matches[0]

    prioritized_statuses = {STATUS_RESOLVED: 0, STATUS_PROPOSED: 1, STATUS_UNIQUE: 2, STATUS_UNRESOLVED: 3}
    return sorted(
        matches,
        key=lambda row: (prioritized_statuses.get(row["mapping_status"], 99), row["parasite_target_chembl_id"]),
    )[0]


def collect_unique_values(rows: list[dict[str, str]], parasite_target_column: str, target_name: str, value_column: str | None) -> str:
    """Collect joined unique values because selectivity datasets may carry repeated metadata rows."""

    if value_column is None:
        return ""

    values = sorted(
        {
            row.get(value_column, "")
            for row in rows
            if row.get(parasite_target_column, "") == target_name and row.get(value_column, "")
        }
    )
    return "|".join(values)


def human_target_name_matches(expected_name: str, expected_gene_symbol: str, observed_name: str) -> bool:
    """Compare expected and observed human labels because historical datasets use both full names and gene symbols."""

    normalized_expected = normalize_text(expected_name)
    normalized_symbol = normalize_text(expected_gene_symbol)
    normalized_observed = normalize_text(observed_name.split("|")[0].split("(")[0].strip())

    if not normalized_observed:
        return True
    if normalized_observed == normalized_expected:
        return True
    if normalized_symbol and normalized_observed == normalized_symbol:
        return True
    if normalized_symbol and normalized_symbol in normalized_observed:
        return True
    if normalized_expected and normalized_observed in normalized_expected:
        return True
    return False


def validate_selectivity_dataset(manifest: list[dict[str, str]], dataset: dict[str, str | Path | None]) -> list[ValidationRecord]:
    """Validate selectivity files because historical outputs should agree with the canonical mapping table."""

    rows = read_csv_rows(dataset["path"])
    parasite_target_column = str(dataset["parasite_target_column"])
    human_target_id_column = dataset["human_target_id_column"]
    human_target_name_column = dataset["human_target_name_column"]
    disease = str(dataset["disease"])
    dataset_key = str(dataset["dataset_key"])
    records: list[ValidationRecord] = []

    parasite_targets = sorted({row.get(parasite_target_column, "") for row in rows if row.get(parasite_target_column, "")})
    for parasite_target_name in parasite_targets:
        manifest_row = select_manifest_row(manifest, disease=disease, parasite_target_name=parasite_target_name)
        if manifest_row is None:
            records.append(
                ValidationRecord(
                    validation_scope="selectivity_dataset",
                    dataset_key=dataset_key,
                    disease=disease,
                    parasite_target_chembl_id="",
                    parasite_target_name=parasite_target_name,
                    validation_status=VALIDATION_ERROR,
                    issue_code="selectivity_target_missing",
                    details="Selectivity dataset references a parasite target missing from the canonical manifest.",
                    expected_human_target_chembl_id="",
                    observed_human_target_chembl_id="",
                    expected_human_target_name="",
                    observed_human_target_name="",
                )
            )
            continue

        observed_human_target_id = collect_unique_values(rows, parasite_target_column, parasite_target_name, human_target_id_column)
        observed_human_target_name = collect_unique_values(rows, parasite_target_column, parasite_target_name, human_target_name_column)
        expected_human_target_id = manifest_row.get("human_target_chembl_id", "")
        expected_human_target_name = manifest_row.get("human_target_name", "")

        status = VALIDATION_OK
        issue_code = "mapping_validated"
        details = "Selectivity dataset is compatible with the canonical mapping."

        if human_target_id_column is not None and observed_human_target_id and observed_human_target_id != expected_human_target_id:
            status = VALIDATION_ERROR
            issue_code = "human_target_id_mismatch"
            details = "Observed human target ID does not match the canonical manifest."

        observed_name_token = observed_human_target_name.split("|")[0] if observed_human_target_name else ""
        if (
            status == VALIDATION_OK
            and human_target_name_column is not None
            and observed_name_token
            and not human_target_name_matches(
                expected_name=expected_human_target_name,
                expected_gene_symbol=manifest_row.get("human_gene_symbol", ""),
                observed_name=observed_name_token,
            )
        ):
            status = VALIDATION_ERROR
            issue_code = "human_target_name_mismatch"
            details = "Observed human target name does not match the canonical manifest."

        if (
            status == VALIDATION_OK
            and not observed_human_target_id
            and not observed_human_target_name
            and manifest_row["mapping_status"] in (STATUS_RESOLVED, STATUS_PROPOSED)
        ):
            status = VALIDATION_WARNING
            issue_code = "human_target_not_recorded_in_dataset"
            details = "Dataset omits human target metadata even though the canonical manifest resolves a comparator."

        records.append(
            ValidationRecord(
                validation_scope="selectivity_dataset",
                dataset_key=dataset_key,
                disease=disease,
                parasite_target_chembl_id=manifest_row["parasite_target_chembl_id"],
                parasite_target_name=manifest_row["parasite_target_name"],
                validation_status=status,
                issue_code=issue_code,
                details=details,
                expected_human_target_chembl_id=expected_human_target_id,
                observed_human_target_chembl_id=observed_human_target_id,
                expected_human_target_name=expected_human_target_name,
                observed_human_target_name=observed_human_target_name,
            )
        )

    return records


def build_validation_table(manifest: list[dict[str, str]], targets: list[dict[str, str]]) -> list[dict[str, str]]:
    """Assemble all validation findings so downstream checks can consume one table."""

    validation_records: list[ValidationRecord] = []
    validation_records.extend(validate_manifest_integrity(manifest))
    validation_records.extend(validate_target_coverage(manifest, targets))

    for dataset in SELECTIVITY_DATASETS:
        validation_records.extend(validate_selectivity_dataset(manifest, dataset))

    return sorted(
        [asdict(record) for record in validation_records],
        key=lambda row: (row["validation_status"], row["disease"], row["dataset_key"], row["parasite_target_name"]),
    )


def format_report_lines(manifest: list[dict[str, str]], validation: list[dict[str, str]], logs: Iterable[str]) -> list[str]:
    """Build a compact report because end-to-end validation should be easy to scan and audit."""

    status_counts = Counter(row["mapping_status"] for row in manifest)
    disease_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in manifest:
        disease_counts[row["disease"]][row["mapping_status"]] += 1

    lines = [
        "Kira canonical target-manifest report",
        "=" * 40,
        "",
        "Pipeline step logs:",
        *[f"  - {line}" for line in logs],
        "",
        "Manifest summary:",
        f"  - Total parasite targets: {len(manifest)}",
        f"  - Resolved mappings: {status_counts.get(STATUS_RESOLVED, 0)}",
        f"  - Proposed mappings: {status_counts.get(STATUS_PROPOSED, 0)}",
        f"  - Parasite-unique targets: {status_counts.get(STATUS_UNIQUE, 0)}",
        f"  - Unresolved targets: {status_counts.get(STATUS_UNRESOLVED, 0)}",
        "",
        "Per-disease coverage:",
    ]

    for disease in sorted(disease_counts):
        counts = disease_counts[disease]
        lines.append(
            "  - "
            f"{disease}: resolved={counts.get(STATUS_RESOLVED, 0)}, "
            f"proposed={counts.get(STATUS_PROPOSED, 0)}, "
            f"unique={counts.get(STATUS_UNIQUE, 0)}, "
            f"unresolved={counts.get(STATUS_UNRESOLVED, 0)}"
        )

    validation_counts = Counter(row["validation_status"] for row in validation)
    lines.extend(["", "Validation summary:"])
    if not validation:
        lines.append("  - No validation findings. All checks passed cleanly.")
    else:
        for status in [VALIDATION_ERROR, VALIDATION_WARNING, VALIDATION_OK]:
            lines.append(f"  - {status}: {validation_counts.get(status, 0)}")

        important_findings = [row for row in validation if row["validation_status"] != VALIDATION_OK]
        if not important_findings:
            lines.append("  - No warnings or errors detected.")
        else:
            lines.extend(["", "Findings:"])
            for row in important_findings:
                lines.append(
                    "  - "
                    f"{row['validation_status'].upper()} [{row['dataset_key']}] "
                    f"{row['parasite_target_name'] or row['parasite_target_chembl_id']}: "
                    f"{row['details']}"
                )

    return lines


def run_target_manifest_pipeline(output_dir: Path | None = None) -> dict[str, object]:
    """Run the manifest pipeline end to end so one function can build and validate all mappings."""

    resolved_output_dir = output_dir or DEFAULT_OUTPUT_DIR
    resolved_output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = resolved_output_dir / DEFAULT_MANIFEST_PATH.name
    validation_path = resolved_output_dir / DEFAULT_VALIDATION_PATH.name
    report_path = resolved_output_dir / DEFAULT_REPORT_PATH.name

    logs: list[str] = []

    start_time = perf_counter()
    targets = load_target_datasets()
    logs.append(
        log_step(
            step_name="load_target_datasets",
            input_desc="3 parasite target catalogs",
            output_desc=f"{len(targets)} unique parasite targets",
            start_time=start_time,
        )
    )

    start_time = perf_counter()
    manifest = build_canonical_manifest(targets)
    logs.append(
        log_step(
            step_name="build_canonical_manifest",
            input_desc=f"{len(targets)} unique parasite targets",
            output_desc=f"{len(manifest)} manifest rows",
            start_time=start_time,
        )
    )

    start_time = perf_counter()
    validation = build_validation_table(manifest, targets)
    logs.append(
        log_step(
            step_name="build_validation_table",
            input_desc=f"{len(manifest)} manifest rows + 3 selectivity datasets",
            output_desc=f"{len(validation)} validation records",
            start_time=start_time,
        )
    )

    start_time = perf_counter()
    write_csv_rows(manifest_path, MANIFEST_FIELDNAMES, manifest)
    write_csv_rows(validation_path, VALIDATION_FIELDNAMES, validation)
    report_lines = format_report_lines(manifest, validation, logs)
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    logs.append(
        log_step(
            step_name="write_outputs",
            input_desc="manifest + validation tables",
            output_desc=f"{manifest_path.name}, {validation_path.name}, {report_path.name}",
            start_time=start_time,
        )
    )

    return {
        "targets": targets,
        "manifest": manifest,
        "validation": validation,
        "logs": logs,
        "manifest_path": manifest_path,
        "validation_path": validation_path,
        "report_path": report_path,
    }
