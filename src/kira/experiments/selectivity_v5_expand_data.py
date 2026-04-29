"""Assay-aware v5 selectivity data expansion for curated parasite-vs-human pairs.

This module is intentionally narrow and conservative.

Scientific stance
-----------------
- It expands a *curated* set of parasite-vs-human comparator pairs.
- It does not claim to solve orthology discovery automatically.
- It preserves assay provenance and exposes uncertainty explicitly.
- It distinguishes exact ratios, approximate ratios, bounds, and weaker evidence.

Operational stance
------------------
- Pair definitions are provided explicitly in a CSV file.
- ChEMBL is used as the primary structured source.
- Filtering is strict by default but policy-controlled.
- Output is a candidate table plus a summary JSON, not a final training table.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence
from urllib.parse import parse_qs, urlparse

import requests

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
except Exception:  # pragma: no cover - optional dependency in some environments
    Chem = None
    MurckoScaffold = None


DEFAULT_CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_ALLOWED_STANDARD_TYPES = ("IC50", "KI", "KD", "EC50", "AC50")
DEFAULT_ALLOWED_RELATIONS = ("=", "<", "<=", ">", ">=")
DEFAULT_ALLOWED_ASSAY_TYPES = ("B", "F")
ASSAY_TYPE_LABEL = {
    "A": "adme",
    "B": "binding",
    "F": "functional",
    "P": "physicochemical",
    "T": "toxicity",
    "U": "unassigned",
}
STANDARD_TYPE_FAMILY = {
    "IC50": "inhibition",
    "KI": "affinity",
    "KD": "affinity",
    "EC50": "functional",
    "AC50": "functional",
}
STRICT_MATCH_STATUS_RANK = {
    "exact_matched_ratio": 0,
    "lower_bound_ratio": 1,
    "upper_bound_ratio": 2,
    "interval_ratio": 3,
    "approximate_ratio": 4,
    "unmatched_comparable": 5,
    "single_side_only": 6,
}


@dataclass(frozen=True)
class PairSpec:
    pair_id: str
    parasite_label: str
    human_label: str
    parasite_target: str
    human_target: str
    parasite_target_namespace: str = "chembl_id"
    human_target_namespace: str = "chembl_id"
    notes: str = ""


@dataclass(frozen=True)
class ExpansionPolicy:
    allowed_standard_types: tuple[str, ...] = DEFAULT_ALLOWED_STANDARD_TYPES
    allowed_relations: tuple[str, ...] = DEFAULT_ALLOWED_RELATIONS
    allowed_assay_types: tuple[str, ...] = DEFAULT_ALLOWED_ASSAY_TYPES
    standard_units: str = "nM"
    min_confidence_score: int = 9
    exclude_flagged_data_validity: bool = True
    require_direct_relationship: bool = True
    require_single_protein_target: bool = True
    exclude_variants: bool = True
    page_limit: int = 1000
    set_chunk_size: int = 50
    request_timeout_seconds: int = 60
    max_retries: int = 3
    retry_backoff_seconds: float = 1.0

    def to_json_ready(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class TargetMetadata:
    target_chembl_id: str
    pref_name: str | None
    organism: str | None
    target_type: str | None
    accession: str | None


@dataclass(frozen=True)
class ActivityInterval:
    lower_nM: float | None
    upper_nM: float | None


def _log(message: str) -> None:
    print(f"[selectivity_v5] {message}", file=sys.stderr, flush=True)


class ChEMBLClientError(RuntimeError):
    pass


class ChEMBLClient:
    """Thin REST client for the ChEMBL web services API."""

    def __init__(
        self,
        base_url: str = DEFAULT_CHEMBL_BASE_URL,
        timeout_seconds: int = 60,
        max_retries: int = 3,
        session: requests.Session | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.max_retries = max_retries
        self.session = session or requests.Session()

    def _request_json(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        url = f"{self.base_url}/{path.lstrip('/')}"
        merged_params = {"format": "json"}
        if params:
            merged_params.update({k: v for k, v in params.items() if v is not None})

        last_error: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            try:
                response = self.session.get(url, params=merged_params, timeout=self.timeout_seconds)
                response.raise_for_status()
                payload = response.json()
                if not isinstance(payload, dict):
                    raise ChEMBLClientError(f"Expected JSON object from {url}, got {type(payload)!r}")
                return payload
            except Exception as exc:  # pragma: no cover - exercised via monkeypatch in integration use
                last_error = exc
                if attempt == self.max_retries:
                    break
        raise ChEMBLClientError(f"Failed to fetch {url}: {last_error}")

    def _infer_list_key(self, payload: dict[str, Any]) -> str:
        candidate_keys = [key for key, value in payload.items() if isinstance(value, list)]
        if len(candidate_keys) == 1:
            return candidate_keys[0]
        if not candidate_keys:
            raise ChEMBLClientError(f"Could not infer list key from payload keys: {sorted(payload.keys())}")
        # Prefer non-metadata plural keys if multiple list fields exist.
        for key in candidate_keys:
            if key != "page_meta":
                return key
        return candidate_keys[0]

    def iter_resource(self, resource: str, params: dict[str, Any] | None = None) -> Iterator[dict[str, Any]]:
        page_limit = params.get("limit") if params else None
        if page_limit is None:
            page_limit = 1000
        offset = params.get("offset") if params else 0
        base_params = dict(params or {})
        base_params["limit"] = page_limit
        base_params["offset"] = offset

        while True:
            payload = self._request_json(f"{resource}.json", base_params)
            list_key = self._infer_list_key(payload)
            for row in payload.get(list_key, []):
                if isinstance(row, dict):
                    yield row
            next_value = payload.get("page_meta", {}).get("next")
            if not next_value:
                break
            next_query = parse_qs(urlparse(next_value).query)
            next_offset = next_query.get("offset", [None])[0]
            if next_offset is None:
                break
            base_params["offset"] = int(next_offset)

    def get_by_id(self, resource: str, identifier: str) -> dict[str, Any]:
        return self._request_json(f"{resource}/{identifier}.json")

    def get_set(self, resource: str, identifiers: Sequence[str]) -> list[dict[str, Any]]:
        if not identifiers:
            return []
        joined = ";".join(identifiers)
        payload = self._request_json(f"{resource}/set/{joined}.json")
        list_key = self._infer_list_key(payload)
        return [row for row in payload.get(list_key, []) if isinstance(row, dict)]

    def resolve_target(self, identifier: str, namespace: str) -> dict[str, Any]:
        namespace_normalized = namespace.strip().lower()
        if namespace_normalized == "chembl_id":
            return self.get_by_id("target", identifier)
        if namespace_normalized in {"accession", "uniprot", "uniprot_accession"}:
            matches = list(
                self.iter_resource(
                    "target",
                    {
                        "target_components__accession": identifier,
                        "target_type": "SINGLE PROTEIN",
                        "limit": 1000,
                    },
                )
            )
            if len(matches) != 1:
                raise ChEMBLClientError(
                    f"Expected exactly one SINGLE PROTEIN target for accession {identifier!r}, found {len(matches)}"
                )
            return matches[0]
        raise ValueError(f"Unsupported target namespace: {namespace!r}")


def normalize_standard_type(value: str | None) -> str | None:
    if value is None:
        return None
    cleaned = value.strip().upper()
    return cleaned or None



def normalize_relation(value: str | None) -> str | None:
    if value is None:
        return None
    cleaned = value.strip()
    return cleaned or None



def safe_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(numeric):
        return None
    return numeric



def read_pair_specs(csv_path: Path) -> list[PairSpec]:
    allowed_namespaces = {"chembl_id", "accession", "uniprot", "uniprot_accession"}

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "pair_id",
            "parasite_label",
            "human_label",
            "parasite_target",
            "human_target",
        }
        fieldnames = set(reader.fieldnames or [])
        missing = sorted(required - fieldnames)
        if missing:
            raise ValueError(f"Pair config is missing required columns: {missing}")

        specs: list[PairSpec] = []
        for line_number, row in enumerate(reader, start=2):
            extra_columns = row.get(None)
            if extra_columns:
                raise ValueError(
                    f"Pair config row {line_number} has extra CSV columns {extra_columns!r}. "
                    "This usually means a label containing a comma was not quoted."
                )

            pair_id = row["pair_id"].strip()
            parasite_label = row["parasite_label"].strip()
            human_label = row["human_label"].strip()
            parasite_target = row["parasite_target"].strip()
            human_target = row["human_target"].strip()
            parasite_namespace = (row.get("parasite_target_namespace") or "chembl_id").strip() or "chembl_id"
            human_namespace = (row.get("human_target_namespace") or "chembl_id").strip() or "chembl_id"

            parasite_namespace_norm = parasite_namespace.lower()
            human_namespace_norm = human_namespace.lower()

            if parasite_namespace_norm not in allowed_namespaces:
                raise ValueError(
                    f"Pair config row {line_number} has unsupported parasite target namespace "
                    f"{parasite_namespace!r} for pair {pair_id!r}."
                )
            if human_namespace_norm not in allowed_namespaces:
                raise ValueError(
                    f"Pair config row {line_number} has unsupported human target namespace "
                    f"{human_namespace!r} for pair {pair_id!r}."
                )

            if parasite_namespace_norm == "chembl_id" and not parasite_target.upper().startswith("CHEMBL"):
                raise ValueError(
                    f"Pair config row {line_number} has parasite_target={parasite_target!r}, "
                    "but namespace is chembl_id."
                )
            if human_namespace_norm == "chembl_id" and not human_target.upper().startswith("CHEMBL"):
                raise ValueError(
                    f"Pair config row {line_number} has human_target={human_target!r}, "
                    "but namespace is chembl_id."
                )

            specs.append(
                PairSpec(
                    pair_id=pair_id,
                    parasite_label=parasite_label,
                    human_label=human_label,
                    parasite_target=parasite_target,
                    human_target=human_target,
                    parasite_target_namespace=parasite_namespace,
                    human_target_namespace=human_namespace,
                    notes=(row.get("notes") or "").strip(),
                )
            )

    if not specs:
        raise ValueError(f"Pair config is empty: {csv_path}")
    return specs



def extract_target_metadata(target_payload: dict[str, Any]) -> TargetMetadata:
    accession = None
    components = target_payload.get("target_components") or []
    if components and isinstance(components[0], dict):
        accession = components[0].get("accession")
    return TargetMetadata(
        target_chembl_id=str(target_payload.get("target_chembl_id") or ""),
        pref_name=target_payload.get("pref_name"),
        organism=target_payload.get("organism"),
        target_type=target_payload.get("target_type"),
        accession=accession,
    )



def get_activity_interval(relation: str, value_nM: float) -> ActivityInterval:
    if value_nM <= 0:
        raise ValueError(f"Activity values must be positive, got {value_nM}")
    if relation == "=":
        return ActivityInterval(lower_nM=value_nM, upper_nM=value_nM)
    if relation in {"<", "<="}:
        return ActivityInterval(lower_nM=0.0, upper_nM=value_nM)
    if relation in {">", ">="}:
        return ActivityInterval(lower_nM=value_nM, upper_nM=None)
    raise ValueError(f"Unsupported relation: {relation!r}")



def classify_standard_family(standard_type: str | None) -> str | None:
    if standard_type is None:
        return None
    return STANDARD_TYPE_FAMILY.get(standard_type)



def is_variant_assay(assay_row: dict[str, Any]) -> bool:
    return bool(
        assay_row.get("variant_sequence_accession")
        or assay_row.get("variant_sequence_mutation")
        or assay_row.get("variant_id") not in (None, "", 0, "0")
    )



def data_validity_is_acceptable(comment: Any, policy: ExpansionPolicy) -> bool:
    if not policy.exclude_flagged_data_validity:
        return True
    if comment in (None, "", "NaN"):
        return True
    return str(comment).strip() == "Manually validated"



def extract_molecule_fields(molecule_row: dict[str, Any]) -> dict[str, str | None]:
    structures = molecule_row.get("molecule_structures") or {}
    hierarchy = molecule_row.get("molecule_hierarchy") or {}
    return {
        "canonical_smiles": structures.get("canonical_smiles") or molecule_row.get("canonical_smiles"),
        "standard_inchi_key": structures.get("standard_inchi_key")
        or molecule_row.get("standard_inchi_key")
        or molecule_row.get("molecule_structures__standard_inchi_key"),
        "parent_molecule_chembl_id": hierarchy.get("parent_chembl_id") or molecule_row.get("parent_molecule_chembl_id"),
        "pref_name": molecule_row.get("pref_name"),
    }



def canonicalize_structure(smiles: str | None, inchi_key: str | None) -> tuple[str | None, str | None, str | None, str]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None, inchi_key, None, "missing_smiles"
    if Chem is None or MurckoScaffold is None:
        return smiles, inchi_key, None, "chembl_raw"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, inchi_key, None, "bad_smiles"

    canonical = Chem.MolToSmiles(mol, canonical=True)
    resolved_inchi_key = inchi_key or Chem.MolToInchiKey(mol)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or None
    return canonical, resolved_inchi_key, scaffold, "resolved"



def assay_passes_policy(
    assay_row: dict[str, Any],
    target_meta: TargetMetadata,
    policy: ExpansionPolicy,
) -> bool:
    assay_type = assay_row.get("assay_type")
    if assay_type not in policy.allowed_assay_types:
        return False

    confidence_score = safe_float(assay_row.get("confidence_score"))
    if confidence_score is None or confidence_score < policy.min_confidence_score:
        return False

    relationship_type = assay_row.get("relationship_type")
    if policy.require_direct_relationship and relationship_type != "D":
        return False

    if policy.exclude_variants and is_variant_assay(assay_row):
        return False

    if policy.require_single_protein_target and target_meta.target_type != "SINGLE PROTEIN":
        return False

    return True



def curate_activity_rows(
    raw_activities: Iterable[dict[str, Any]],
    assays_by_id: dict[str, dict[str, Any]],
    molecules_by_id: dict[str, dict[str, Any]],
    target_meta: TargetMetadata,
    pair_id: str,
    side: str,
    policy: ExpansionPolicy,
) -> list[dict[str, Any]]:
    curated: list[dict[str, Any]] = []
    for row in raw_activities:
        standard_type = normalize_standard_type(row.get("standard_type"))
        standard_relation = normalize_relation(row.get("standard_relation"))
        standard_units = row.get("standard_units")
        standard_value_nM = safe_float(row.get("standard_value"))
        molecule_chembl_id = row.get("molecule_chembl_id")
        assay_chembl_id = row.get("assay_chembl_id")

        if standard_type not in policy.allowed_standard_types:
            continue
        if standard_relation not in policy.allowed_relations:
            continue
        if standard_units != policy.standard_units:
            continue
        if standard_value_nM is None or standard_value_nM <= 0:
            continue
        if not molecule_chembl_id or not assay_chembl_id:
            continue

        assay_row = assays_by_id.get(str(assay_chembl_id))
        if not assay_row:
            continue
        if not assay_passes_policy(assay_row, target_meta, policy):
            continue
        if not data_validity_is_acceptable(row.get("data_validity_comment"), policy):
            continue

        molecule_row = molecules_by_id.get(str(molecule_chembl_id), {})
        molecule_fields = extract_molecule_fields(molecule_row)
        canonical_smiles, standard_inchi_key, murcko_scaffold, structure_status = canonicalize_structure(
            molecule_fields["canonical_smiles"],
            molecule_fields["standard_inchi_key"],
        )
        assay_type = assay_row.get("assay_type")

        curated.append(
            {
                "pair_id": pair_id,
                "side": side,
                "target_chembl_id": target_meta.target_chembl_id,
                "target_pref_name": target_meta.pref_name,
                "target_organism": target_meta.organism,
                "target_type": target_meta.target_type,
                "target_accession": target_meta.accession,
                "molecule_chembl_id": str(molecule_chembl_id),
                "parent_molecule_chembl_id": molecule_fields["parent_molecule_chembl_id"],
                "canonical_smiles": canonical_smiles,
                "standard_inchi_key": standard_inchi_key,
                "murcko_scaffold": murcko_scaffold,
                "structure_status": structure_status,
                "compound_pref_name": molecule_fields["pref_name"],
                "activity_chembl_id": row.get("activity_chembl_id"),
                "assay_chembl_id": str(assay_chembl_id),
                "assay_type": assay_type,
                "assay_type_label": ASSAY_TYPE_LABEL.get(assay_type),
                "assay_description": assay_row.get("description") or assay_row.get("assay_description"),
                "assay_cell_type": assay_row.get("cell_type"),
                "assay_organism": assay_row.get("assay_organism"),
                "assay_tax_id": assay_row.get("assay_tax_id"),
                "document_chembl_id": row.get("document_chembl_id"),
                "src_id": row.get("src_id"),
                "source_description": row.get("source_description"),
                "standard_type": standard_type,
                "standard_family": classify_standard_family(standard_type),
                "standard_relation": standard_relation,
                "standard_value_nM": standard_value_nM,
                "standard_units": standard_units,
                "pchembl_value": safe_float(row.get("pchembl_value")),
                "activity_comment": row.get("activity_comment"),
                "data_validity_comment": row.get("data_validity_comment"),
                "data_validity_description": row.get("data_validity_description"),
                "potential_duplicate": row.get("potential_duplicate"),
                "confidence_score": int(float(assay_row.get("confidence_score"))),
                "relationship_type": assay_row.get("relationship_type"),
                "variant_sequence_accession": assay_row.get("variant_sequence_accession"),
                "variant_sequence_mutation": assay_row.get("variant_sequence_mutation"),
                "measurement_interval_kind": relation_to_interval_kind(standard_relation),
            }
        )
    return curated



def relation_to_interval_kind(relation: str) -> str:
    if relation == "=":
        return "exact"
    if relation in {"<", "<="}:
        return "upper_bound"
    if relation in {">", ">="}:
        return "lower_bound"
    raise ValueError(f"Unsupported relation: {relation!r}")



def comparison_tier(parasite_activity: dict[str, Any], human_activity: dict[str, Any]) -> str:
    same_standard_type = parasite_activity["standard_type"] == human_activity["standard_type"]
    same_assay_type = parasite_activity["assay_type"] == human_activity["assay_type"]
    same_standard_family = parasite_activity["standard_family"] == human_activity["standard_family"]

    if same_standard_type and same_assay_type:
        return "strict"
    if same_standard_type:
        return "approximate"
    if same_standard_family and parasite_activity["standard_family"] is not None:
        return "approximate"
    return "unmatched"



def ratio_bounds(
    parasite_relation: str,
    parasite_value_nM: float,
    human_relation: str,
    human_value_nM: float,
) -> tuple[float | None, float | None]:
    parasite_interval = get_activity_interval(parasite_relation, parasite_value_nM)
    human_interval = get_activity_interval(human_relation, human_value_nM)

    lower_bound = None
    upper_bound = None

    if human_interval.lower_nM is not None and parasite_interval.upper_nM is not None and parasite_interval.upper_nM > 0:
        lower_bound = human_interval.lower_nM / parasite_interval.upper_nM

    if human_interval.upper_nM is not None and parasite_interval.lower_nM is not None and parasite_interval.lower_nM > 0:
        upper_bound = human_interval.upper_nM / parasite_interval.lower_nM

    return lower_bound, upper_bound



def candidate_status(
    match_tier: str,
    parasite_relation: str,
    parasite_value_nM: float,
    human_relation: str,
    human_value_nM: float,
) -> tuple[str, float | None, float | None, float | None]:
    if match_tier == "unmatched":
        return "unmatched_comparable", None, None, None

    lower_bound, upper_bound = ratio_bounds(
        parasite_relation=parasite_relation,
        parasite_value_nM=parasite_value_nM,
        human_relation=human_relation,
        human_value_nM=human_value_nM,
    )

    exact_parasite = parasite_relation == "="
    exact_human = human_relation == "="
    if exact_parasite and exact_human and lower_bound is not None and upper_bound is not None:
        exact_ratio = lower_bound
        if match_tier == "strict":
            return "exact_matched_ratio", exact_ratio, lower_bound, upper_bound
        return "approximate_ratio", exact_ratio, lower_bound, upper_bound

    if lower_bound is not None and upper_bound is not None:
        return "interval_ratio", None, lower_bound, upper_bound
    if lower_bound is not None:
        return "lower_bound_ratio", None, lower_bound, None
    if upper_bound is not None:
        return "upper_bound_ratio", None, None, upper_bound
    return "unmatched_comparable", None, None, None



def build_pair_candidates(
    pair_spec: PairSpec,
    parasite_activities: Sequence[dict[str, Any]],
    human_activities: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    parasite_by_molecule: dict[str, list[dict[str, Any]]] = defaultdict(list)
    human_by_molecule: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in parasite_activities:
        parasite_by_molecule[row["molecule_chembl_id"]].append(row)
    for row in human_activities:
        human_by_molecule[row["molecule_chembl_id"]].append(row)

    candidate_rows: list[dict[str, Any]] = []
    all_molecules = sorted(set(parasite_by_molecule) | set(human_by_molecule))
    candidate_index = 0

    for molecule_chembl_id in all_molecules:
        p_rows = parasite_by_molecule.get(molecule_chembl_id, [])
        h_rows = human_by_molecule.get(molecule_chembl_id, [])

        if not p_rows or not h_rows:
            lone_side = "parasite" if p_rows else "human"
            lone_rows = p_rows or h_rows
            for lone in lone_rows:
                candidate_index += 1
                candidate_rows.append(
                    {
                        "candidate_row_id": f"{pair_spec.pair_id}::{molecule_chembl_id}::{candidate_index}",
                        "pair_id": pair_spec.pair_id,
                        "parasite_label": pair_spec.parasite_label,
                        "human_label": pair_spec.human_label,
                        "molecule_chembl_id": molecule_chembl_id,
                        "canonical_smiles": lone.get("canonical_smiles"),
                        "standard_inchi_key": lone.get("standard_inchi_key"),
                        "murcko_scaffold": lone.get("murcko_scaffold"),
                        "structure_status": lone.get("structure_status"),
                        "candidate_status": "single_side_only",
                        "comparison_tier": "single_side_only",
                        "same_standard_type": None,
                        "same_assay_type": None,
                        "same_document": None,
                        "ratio_exact_human_div_parasite": None,
                        "ratio_lower_bound_human_div_parasite": None,
                        "ratio_upper_bound_human_div_parasite": None,
                        "single_side": lone_side,
                        "parasite_target_chembl_id": p_rows[0]["target_chembl_id"] if p_rows else pair_spec.parasite_target,
                        "human_target_chembl_id": h_rows[0]["target_chembl_id"] if h_rows else pair_spec.human_target,
                        "parasite_activity_chembl_id": lone.get("activity_chembl_id") if lone_side == "parasite" else None,
                        "human_activity_chembl_id": lone.get("activity_chembl_id") if lone_side == "human" else None,
                        "parasite_assay_chembl_id": lone.get("assay_chembl_id") if lone_side == "parasite" else None,
                        "human_assay_chembl_id": lone.get("assay_chembl_id") if lone_side == "human" else None,
                        "parasite_standard_type": lone.get("standard_type") if lone_side == "parasite" else None,
                        "human_standard_type": lone.get("standard_type") if lone_side == "human" else None,
                        "parasite_standard_relation": lone.get("standard_relation") if lone_side == "parasite" else None,
                        "human_standard_relation": lone.get("standard_relation") if lone_side == "human" else None,
                        "parasite_standard_value_nM": lone.get("standard_value_nM") if lone_side == "parasite" else None,
                        "human_standard_value_nM": lone.get("standard_value_nM") if lone_side == "human" else None,
                        "parasite_pchembl_value": lone.get("pchembl_value") if lone_side == "parasite" else None,
                        "human_pchembl_value": lone.get("pchembl_value") if lone_side == "human" else None,
                        "parasite_confidence_score": lone.get("confidence_score") if lone_side == "parasite" else None,
                        "human_confidence_score": lone.get("confidence_score") if lone_side == "human" else None,
                        "parasite_relationship_type": lone.get("relationship_type") if lone_side == "parasite" else None,
                        "human_relationship_type": lone.get("relationship_type") if lone_side == "human" else None,
                        "parasite_data_validity_comment": lone.get("data_validity_comment") if lone_side == "parasite" else None,
                        "human_data_validity_comment": lone.get("data_validity_comment") if lone_side == "human" else None,
                        "parasite_potential_duplicate": lone.get("potential_duplicate") if lone_side == "parasite" else None,
                        "human_potential_duplicate": lone.get("potential_duplicate") if lone_side == "human" else None,
                        "parasite_assay_type": lone.get("assay_type") if lone_side == "parasite" else None,
                        "human_assay_type": lone.get("assay_type") if lone_side == "human" else None,
                        "parasite_document_chembl_id": lone.get("document_chembl_id") if lone_side == "parasite" else None,
                        "human_document_chembl_id": lone.get("document_chembl_id") if lone_side == "human" else None,
                        "parasite_measurement_interval_kind": lone.get("measurement_interval_kind") if lone_side == "parasite" else None,
                        "human_measurement_interval_kind": lone.get("measurement_interval_kind") if lone_side == "human" else None,
                    }
                )
            continue

        for p_row in p_rows:
            for h_row in h_rows:
                candidate_index += 1
                match_tier = comparison_tier(p_row, h_row)
                status, exact_ratio, lower_bound, upper_bound = candidate_status(
                    match_tier=match_tier,
                    parasite_relation=p_row["standard_relation"],
                    parasite_value_nM=p_row["standard_value_nM"],
                    human_relation=h_row["standard_relation"],
                    human_value_nM=h_row["standard_value_nM"],
                )
                candidate_rows.append(
                    {
                        "candidate_row_id": f"{pair_spec.pair_id}::{molecule_chembl_id}::{candidate_index}",
                        "pair_id": pair_spec.pair_id,
                        "parasite_label": pair_spec.parasite_label,
                        "human_label": pair_spec.human_label,
                        "molecule_chembl_id": molecule_chembl_id,
                        "canonical_smiles": p_row.get("canonical_smiles") or h_row.get("canonical_smiles"),
                        "standard_inchi_key": p_row.get("standard_inchi_key") or h_row.get("standard_inchi_key"),
                        "murcko_scaffold": p_row.get("murcko_scaffold") or h_row.get("murcko_scaffold"),
                        "structure_status": p_row.get("structure_status") or h_row.get("structure_status"),
                        "candidate_status": status,
                        "comparison_tier": match_tier,
                        "same_standard_type": p_row.get("standard_type") == h_row.get("standard_type"),
                        "same_assay_type": p_row.get("assay_type") == h_row.get("assay_type"),
                        "same_document": p_row.get("document_chembl_id") == h_row.get("document_chembl_id"),
                        "ratio_exact_human_div_parasite": exact_ratio,
                        "ratio_lower_bound_human_div_parasite": lower_bound,
                        "ratio_upper_bound_human_div_parasite": upper_bound,
                        "single_side": None,
                        "parasite_target_chembl_id": p_row["target_chembl_id"],
                        "human_target_chembl_id": h_row["target_chembl_id"],
                        "parasite_activity_chembl_id": p_row.get("activity_chembl_id"),
                        "human_activity_chembl_id": h_row.get("activity_chembl_id"),
                        "parasite_assay_chembl_id": p_row.get("assay_chembl_id"),
                        "human_assay_chembl_id": h_row.get("assay_chembl_id"),
                        "parasite_standard_type": p_row.get("standard_type"),
                        "human_standard_type": h_row.get("standard_type"),
                        "parasite_standard_family": p_row.get("standard_family"),
                        "human_standard_family": h_row.get("standard_family"),
                        "parasite_standard_relation": p_row.get("standard_relation"),
                        "human_standard_relation": h_row.get("standard_relation"),
                        "parasite_standard_value_nM": p_row.get("standard_value_nM"),
                        "human_standard_value_nM": h_row.get("standard_value_nM"),
                        "parasite_pchembl_value": p_row.get("pchembl_value"),
                        "human_pchembl_value": h_row.get("pchembl_value"),
                        "parasite_confidence_score": p_row.get("confidence_score"),
                        "human_confidence_score": h_row.get("confidence_score"),
                        "parasite_relationship_type": p_row.get("relationship_type"),
                        "human_relationship_type": h_row.get("relationship_type"),
                        "parasite_data_validity_comment": p_row.get("data_validity_comment"),
                        "human_data_validity_comment": h_row.get("data_validity_comment"),
                        "parasite_potential_duplicate": p_row.get("potential_duplicate"),
                        "human_potential_duplicate": h_row.get("potential_duplicate"),
                        "parasite_assay_type": p_row.get("assay_type"),
                        "human_assay_type": h_row.get("assay_type"),
                        "parasite_document_chembl_id": p_row.get("document_chembl_id"),
                        "human_document_chembl_id": h_row.get("document_chembl_id"),
                        "parasite_measurement_interval_kind": p_row.get("measurement_interval_kind"),
                        "human_measurement_interval_kind": h_row.get("measurement_interval_kind"),
                    }
                )
    return candidate_rows



def fetch_target_activity_rows(
    client: ChEMBLClient,
    target_chembl_id: str,
    policy: ExpansionPolicy,
) -> list[dict[str, Any]]:
    params = {
        "target_chembl_id": target_chembl_id,
        "standard_type__in": ",".join(policy.allowed_standard_types),
        "standard_units": policy.standard_units,
        "limit": policy.page_limit,
    }
    return list(client.iter_resource("activity", params=params))



def fetch_assays_by_id(
    client: ChEMBLClient,
    assay_ids: Iterable[str],
    chunk_size: int,
) -> dict[str, dict[str, Any]]:
    assay_id_list = sorted({str(assay_id) for assay_id in assay_ids if assay_id})
    assays_by_id: dict[str, dict[str, Any]] = {}
    if not assay_id_list:
        return assays_by_id

    effective_chunk_size = max(1, chunk_size)
    total_chunks = math.ceil(len(assay_id_list) / effective_chunk_size)
    _log(f"fetching {len(assay_id_list)} assay records in {total_chunks} chunks")

    for chunk_index, start in enumerate(range(0, len(assay_id_list), effective_chunk_size), start=1):
        chunk = assay_id_list[start : start + effective_chunk_size]
        _log(f"assay chunk {chunk_index}/{total_chunks} ({len(chunk)} ids)")
        try:
            rows = client.get_set("assay", chunk)
        except ChEMBLClientError as exc:
            _log(f"assay set lookup failed for chunk {chunk_index}; falling back to per-id lookup: {exc}")
            rows = []
            for assay_id in chunk:
                try:
                    rows.append(client.get_by_id("assay", assay_id))
                except ChEMBLClientError as single_exc:
                    _log(f"skipping assay {assay_id}: {single_exc}")
        for row in rows:
            assay_chembl_id = row.get("assay_chembl_id")
            if assay_chembl_id:
                assays_by_id[str(assay_chembl_id)] = row
    return assays_by_id


def activity_row_passes_nonstructure_policy(
    row: dict[str, Any],
    assays_by_id: dict[str, dict[str, Any]],
    target_meta: TargetMetadata,
    policy: ExpansionPolicy,
) -> bool:
    standard_type = normalize_standard_type(row.get("standard_type"))
    standard_relation = normalize_relation(row.get("standard_relation"))
    standard_units = row.get("standard_units")
    standard_value_nM = safe_float(row.get("standard_value"))
    molecule_chembl_id = row.get("molecule_chembl_id")
    assay_chembl_id = row.get("assay_chembl_id")

    if standard_type not in policy.allowed_standard_types:
        return False
    if standard_relation not in policy.allowed_relations:
        return False
    if standard_units != policy.standard_units:
        return False
    if standard_value_nM is None or standard_value_nM <= 0:
        return False
    if not molecule_chembl_id or not assay_chembl_id:
        return False

    assay_row = assays_by_id.get(str(assay_chembl_id))
    if not assay_row:
        return False
    if not assay_passes_policy(assay_row, target_meta, policy):
        return False
    if not data_validity_is_acceptable(row.get("data_validity_comment"), policy):
        return False
    return True


def fetch_molecules_by_id(
    client: ChEMBLClient,
    molecule_ids: Iterable[str],
    chunk_size: int,
) -> dict[str, dict[str, Any]]:
    molecule_id_list = sorted({str(molecule_id) for molecule_id in molecule_ids if molecule_id})
    molecules_by_id: dict[str, dict[str, Any]] = {}
    if not molecule_id_list:
        return molecules_by_id

    effective_chunk_size = max(1, chunk_size)
    total_chunks = math.ceil(len(molecule_id_list) / effective_chunk_size)
    _log(f"fetching {len(molecule_id_list)} molecule records in {total_chunks} chunks")

    for chunk_index, start in enumerate(range(0, len(molecule_id_list), effective_chunk_size), start=1):
        chunk = molecule_id_list[start : start + effective_chunk_size]
        _log(f"molecule chunk {chunk_index}/{total_chunks} ({len(chunk)} ids)")
        try:
            rows = client.get_set("molecule", chunk)
        except ChEMBLClientError as exc:
            _log(f"molecule set lookup failed for chunk {chunk_index}; falling back to per-id lookup: {exc}")
            rows = []
            for molecule_id in chunk:
                try:
                    rows.append(client.get_by_id("molecule", molecule_id))
                except ChEMBLClientError as single_exc:
                    _log(f"skipping molecule {molecule_id}: {single_exc}")
        for row in rows:
            molecule_chembl_id = row.get("molecule_chembl_id")
            if molecule_chembl_id:
                molecules_by_id[str(molecule_chembl_id)] = row
    return molecules_by_id



def summarize_expansion(
    pair_specs: Sequence[PairSpec],
    target_metadata_by_id: dict[str, TargetMetadata],
    curated_activities: Sequence[dict[str, Any]],
    candidate_rows: Sequence[dict[str, Any]],
    policy: ExpansionPolicy,
) -> dict[str, Any]:
    activities_by_side = Counter(row["side"] for row in curated_activities)
    activities_by_pair_and_side: dict[str, Counter[str]] = defaultdict(Counter)
    unique_molecules_by_pair_and_side: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"parasite": set(), "human": set()})
    standard_type_counts = Counter()
    assay_type_counts = Counter()
    status_counts = Counter()
    rows_by_pair = Counter()
    unique_candidate_molecules_by_pair: dict[str, set[str]] = defaultdict(set)

    for row in curated_activities:
        activities_by_pair_and_side[row["pair_id"]][row["side"]] += 1
        unique_molecules_by_pair_and_side[row["pair_id"]][row["side"]].add(row["molecule_chembl_id"])
        standard_type_counts[row["standard_type"]] += 1
        assay_type_counts[row["assay_type"]] += 1

    for row in candidate_rows:
        status_counts[row["candidate_status"]] += 1
        rows_by_pair[row["pair_id"]] += 1
        unique_candidate_molecules_by_pair[row["pair_id"]].add(row["molecule_chembl_id"])

    pair_summaries: dict[str, Any] = {}
    for spec in pair_specs:
        pair_candidates = [row for row in candidate_rows if row["pair_id"] == spec.pair_id]
        pair_status_counts = Counter(row["candidate_status"] for row in pair_candidates)
        pair_summaries[spec.pair_id] = {
            "parasite_label": spec.parasite_label,
            "human_label": spec.human_label,
            "parasite_target_chembl_id": spec.parasite_target,
            "human_target_chembl_id": spec.human_target,
            "parasite_activity_rows": activities_by_pair_and_side[spec.pair_id]["parasite"],
            "human_activity_rows": activities_by_pair_and_side[spec.pair_id]["human"],
            "unique_parasite_molecules": len(unique_molecules_by_pair_and_side[spec.pair_id]["parasite"]),
            "unique_human_molecules": len(unique_molecules_by_pair_and_side[spec.pair_id]["human"]),
            "candidate_rows": len(pair_candidates),
            "unique_candidate_molecules": len(unique_candidate_molecules_by_pair[spec.pair_id]),
            "candidate_status_counts": dict(sorted(pair_status_counts.items())),
        }

    target_metadata_json = {
        target_id: asdict(metadata)
        for target_id, metadata in sorted(target_metadata_by_id.items())
    }

    return {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "scientific_scope": (
            "Curated parasite-vs-human comparator expansion from ChEMBL with assay-aware provenance; "
            "candidate rows are not final trainable labels."
        ),
        "policy": policy.to_json_ready(),
        "n_pair_specs": len(pair_specs),
        "n_curated_activity_rows": len(curated_activities),
        "n_candidate_rows": len(candidate_rows),
        "activity_rows_by_side": dict(sorted(activities_by_side.items())),
        "activity_rows_by_standard_type": dict(sorted(standard_type_counts.items())),
        "activity_rows_by_assay_type": dict(sorted(assay_type_counts.items())),
        "candidate_rows_by_status": dict(sorted(status_counts.items())),
        "candidate_rows_by_pair": dict(sorted(rows_by_pair.items())),
        "pair_summaries": pair_summaries,
        "target_metadata": target_metadata_json,
    }



def write_candidate_csv(rows: Sequence[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    preferred_prefix = [
        "candidate_row_id",
        "pair_id",
        "parasite_label",
        "human_label",
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_inchi_key",
        "candidate_status",
        "comparison_tier",
        "ratio_exact_human_div_parasite",
        "ratio_lower_bound_human_div_parasite",
        "ratio_upper_bound_human_div_parasite",
    ]
    fieldnames: list[str] = []
    seen: set[str] = set()
    for key in preferred_prefix:
        if key not in seen:
            fieldnames.append(key)
            seen.add(key)
    for row in rows:
        for key in row.keys():
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)
    if not rows:
        fieldnames = preferred_prefix
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})



def write_summary_json(summary: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")



def expand_selectivity_dataset(
    pair_specs: Sequence[PairSpec],
    policy: ExpansionPolicy,
    client: ChEMBLClient,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    target_metadata_by_id: dict[str, TargetMetadata] = {}
    resolved_specs: list[PairSpec] = []

    for spec in pair_specs:
        parasite_payload = client.resolve_target(spec.parasite_target, spec.parasite_target_namespace)
        human_payload = client.resolve_target(spec.human_target, spec.human_target_namespace)
        parasite_meta = extract_target_metadata(parasite_payload)
        human_meta = extract_target_metadata(human_payload)
        target_metadata_by_id[parasite_meta.target_chembl_id] = parasite_meta
        target_metadata_by_id[human_meta.target_chembl_id] = human_meta
        resolved_specs.append(
            PairSpec(
                pair_id=spec.pair_id,
                parasite_label=spec.parasite_label,
                human_label=spec.human_label,
                parasite_target=parasite_meta.target_chembl_id,
                human_target=human_meta.target_chembl_id,
                parasite_target_namespace="chembl_id",
                human_target_namespace="chembl_id",
                notes=spec.notes,
            )
        )

    curated_activities: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []

    for pair_index, spec in enumerate(resolved_specs, start=1):
        parasite_meta = target_metadata_by_id[spec.parasite_target]
        human_meta = target_metadata_by_id[spec.human_target]

        _log(
            f"[{pair_index}/{len(resolved_specs)}] {spec.pair_id}: fetching activity rows "
            f"for {parasite_meta.target_chembl_id} and {human_meta.target_chembl_id}"
        )
        parasite_raw = fetch_target_activity_rows(client, parasite_meta.target_chembl_id, policy)
        human_raw = fetch_target_activity_rows(client, human_meta.target_chembl_id, policy)
        _log(f"{spec.pair_id}: raw activity rows parasite={len(parasite_raw)} human={len(human_raw)}")

        all_assay_ids = [row.get("assay_chembl_id") for row in [*parasite_raw, *human_raw]]
        assays_by_id = fetch_assays_by_id(client, all_assay_ids, chunk_size=policy.set_chunk_size)

        parasite_prefiltered = [
            row
            for row in parasite_raw
            if activity_row_passes_nonstructure_policy(row, assays_by_id, parasite_meta, policy)
        ]
        human_prefiltered = [
            row
            for row in human_raw
            if activity_row_passes_nonstructure_policy(row, assays_by_id, human_meta, policy)
        ]
        _log(
            f"{spec.pair_id}: assay-policy survivors before structure fetch "
            f"parasite={len(parasite_prefiltered)} human={len(human_prefiltered)}"
        )

        all_molecule_ids = [row.get("molecule_chembl_id") for row in [*parasite_prefiltered, *human_prefiltered]]
        molecules_by_id = fetch_molecules_by_id(client, all_molecule_ids, chunk_size=policy.set_chunk_size)

        parasite_curated = curate_activity_rows(
            raw_activities=parasite_prefiltered,
            assays_by_id=assays_by_id,
            molecules_by_id=molecules_by_id,
            target_meta=parasite_meta,
            pair_id=spec.pair_id,
            side="parasite",
            policy=policy,
        )
        human_curated = curate_activity_rows(
            raw_activities=human_prefiltered,
            assays_by_id=assays_by_id,
            molecules_by_id=molecules_by_id,
            target_meta=human_meta,
            pair_id=spec.pair_id,
            side="human",
            policy=policy,
        )

        pair_candidates = build_pair_candidates(spec, parasite_curated, human_curated)
        curated_activities.extend(parasite_curated)
        curated_activities.extend(human_curated)
        candidate_rows.extend(pair_candidates)
        _log(
            f"{spec.pair_id}: curated parasite={len(parasite_curated)} human={len(human_curated)} "
            f"candidate_rows={len(pair_candidates)}"
        )

    summary = summarize_expansion(
        pair_specs=resolved_specs,
        target_metadata_by_id=target_metadata_by_id,
        curated_activities=curated_activities,
        candidate_rows=candidate_rows,
        policy=policy,
    )
    return candidate_rows, summary



def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pair-config",
        type=Path,
        required=True,
        help="CSV file defining curated parasite-vs-human target pairs.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/processed/selectivity_v5_candidate_rows.csv"),
        help="Output CSV path for candidate expansion rows.",
    )
    parser.add_argument(
        "--output-summary-json",
        type=Path,
        default=Path("data/processed/selectivity_v5_expansion_summary.json"),
        help="Output JSON path for expansion summary.",
    )
    parser.add_argument(
        "--chembl-base-url",
        default=DEFAULT_CHEMBL_BASE_URL,
        help="Base URL for the ChEMBL REST API.",
    )
    parser.add_argument(
        "--standard-type",
        dest="standard_types",
        action="append",
        help="Repeatable. Restrict to one or more standard types (default: IC50, KI, KD, EC50, AC50).",
    )
    parser.add_argument(
        "--standard-units",
        default="nM",
        help="Standard activity units to require (default: nM).",
    )
    parser.add_argument(
        "--min-confidence-score",
        type=int,
        default=9,
        help="Minimum ChEMBL assay confidence score (default: 9).",
    )
    parser.add_argument(
        "--allow-assay-type",
        dest="assay_types",
        action="append",
        help="Repeatable. Allowed ChEMBL assay types (default: B and F).",
    )
    parser.add_argument(
        "--keep-indirect-relationships",
        action="store_true",
        help="Keep assays whose relationship_type is not D.",
    )
    parser.add_argument(
        "--include-variants",
        action="store_true",
        help="Keep variant assays instead of excluding them.",
    )
    parser.add_argument(
        "--keep-flagged-data-validity",
        action="store_true",
        help="Keep activities with a non-empty data_validity_comment other than 'Manually validated'.",
    )
    parser.add_argument(
        "--allow-non-single-protein-targets",
        action="store_true",
        help="Do not require target_type=SINGLE PROTEIN.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=60,
        help="HTTP timeout in seconds for each API call.",
    )
    parser.add_argument(
        "--set-chunk-size",
        type=int,
        default=50,
        help="Chunk size for ChEMBL set endpoint lookups.",
    )
    return parser.parse_args(argv)



def build_policy_from_args(args: argparse.Namespace) -> ExpansionPolicy:
    standard_types = tuple(normalize_standard_type(value) for value in (args.standard_types or DEFAULT_ALLOWED_STANDARD_TYPES))
    assay_types = tuple((value or "").strip().upper() for value in (args.assay_types or DEFAULT_ALLOWED_ASSAY_TYPES))
    return ExpansionPolicy(
        allowed_standard_types=tuple(value for value in standard_types if value),
        allowed_assay_types=tuple(value for value in assay_types if value),
        standard_units=args.standard_units,
        min_confidence_score=args.min_confidence_score,
        require_direct_relationship=not args.keep_indirect_relationships,
        require_single_protein_target=not args.allow_non_single_protein_targets,
        exclude_variants=not args.include_variants,
        exclude_flagged_data_validity=not args.keep_flagged_data_validity,
        page_limit=1000,
        set_chunk_size=args.set_chunk_size,
        request_timeout_seconds=args.timeout_seconds,
        max_retries=3,
    )



def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    pair_specs = read_pair_specs(args.pair_config)
    policy = build_policy_from_args(args)
    client = ChEMBLClient(
        base_url=args.chembl_base_url,
        timeout_seconds=policy.request_timeout_seconds,
        max_retries=policy.max_retries,
    )
    candidate_rows, summary = expand_selectivity_dataset(
        pair_specs=pair_specs,
        policy=policy,
        client=client,
    )
    write_candidate_csv(candidate_rows, args.output_csv)
    write_summary_json(summary, args.output_summary_json)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
