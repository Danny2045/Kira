"""Verify UniProt accession mappings for the eleven Kira target proteins.

Queries the UniProt REST API for each candidate accession and checks
organism, recommended name, and canonical sequence length against
expected values. Writes a structured provenance record to
data/models/uniprot_mappings_provenance.json and prints a human-readable
report to stdout. Read-only: does not modify uniprot_ids_v2.json or
sequences_v2.json.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
USER_AGENT = "kira-uniprot-verifier/1.0 (github.com/Danny2045/Kira)"
REQUEST_DELAY_SECONDS = 0.5
PROVENANCE_PATH = Path("data/models/uniprot_mappings_provenance.json")

# (key, candidate_accession, organism_substring, name_keyword, min_length, max_length)
TARGETS: list[tuple[str, str, str, str, int, int]] = [
    ("SmHDAC8", "A5H660", "Schistosoma mansoni", "histone deacetylase", 350, 450),
    ("SmDHODH", "G4VFD7", "Schistosoma mansoni", "dihydroorotate dehydrogenase", 300, 400),
    ("LmDHFR", "P07382", "Leishmania major", "dihydrofolate reductase", 500, 540),
    ("LmPTR1", "Q01782", "Leishmania major", "pteridine reductase", 250, 300),
    ("TbCatB", "Q6R7Z5", "Trypanosoma brucei", "cathepsin", 300, 360),
    ("TbPDEB1", "Q8WQX9", "Trypanosoma brucei", "phosphodiesterase", 850, 950),
    ("HsHDAC8", "Q9BY41", "Homo sapiens", "histone deacetylase", 350, 400),
    ("HsDHODH", "Q02127", "Homo sapiens", "dihydroorotate dehydrogenase", 380, 410),
    ("HsDHFR", "P00374", "Homo sapiens", "dihydrofolate reductase", 180, 220),
    ("HsCatB", "P07858", "Homo sapiens", "cathepsin B", 300, 360),
    ("HsPDE4B", "Q07343", "Homo sapiens", "phosphodiesterase", 700, 780),
]


def fetch_uniprot(accession: str) -> dict | None:
    """Fetch a UniProt entry. Returns parsed JSON, or None for HTTP 404."""
    url = UNIPROT_URL.format(accession=accession)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def extract_recommended_name(entry: dict) -> str:
    """Pull the protein name out of a UniProt JSON entry.

    Swiss-Prot entries use recommendedName; TrEMBL entries use submittedName.
    """
    description = entry.get("proteinDescription", {})
    rec = description.get("recommendedName")
    if rec and rec.get("fullName", {}).get("value"):
        return rec["fullName"]["value"]
    for field in ("submissionNames", "submittedNames"):
        submitted = description.get(field) or []
        if submitted and submitted[0].get("fullName", {}).get("value"):
            return submitted[0]["fullName"]["value"]
    return ""


def verify_one(key: str, accession: str, organism_sub: str, name_kw: str,
               min_len: int, max_len: int) -> dict:
    """Verify a single candidate accession. Returns a provenance record."""
    retrieved_at = datetime.now(timezone.utc).isoformat()
    entry = fetch_uniprot(accession)
    if entry is None:
        return {
            "key": key,
            "accession": accession,
            "status": "FAIL",
            "reason": "HTTP 404: accession not found",
            "retrieved_at": retrieved_at,
        }

    organism = entry.get("organism", {}).get("scientificName", "")
    name = extract_recommended_name(entry)
    sequence = entry.get("sequence", {}).get("value", "")
    length = entry.get("sequence", {}).get("length", len(sequence))
    audit = entry.get("entryAudit", {})

    failures: list[str] = []
    if organism_sub.lower() not in organism.lower():
        failures.append(
            f"organism mismatch: expected substring '{organism_sub}', got '{organism}'"
        )
    if name_kw.lower() not in name.lower():
        failures.append(
            f"name mismatch: expected keyword '{name_kw}', got '{name}'"
        )
    if not (min_len <= length <= max_len):
        failures.append(
            f"length out of range: expected {min_len}-{max_len}, got {length}"
        )

    return {
        "key": key,
        "accession": accession,
        "organism": organism,
        "recommended_name": name,
        "sequence_length": length,
        "sequence": sequence,
        "entry_version": audit.get("entryVersion"),
        "sequence_version": audit.get("sequenceVersion"),
        "last_sequence_update_date": audit.get("lastSequenceUpdateDate"),
        "retrieved_at": retrieved_at,
        "status": "PASS" if not failures else "FAIL",
        "reason": "; ".join(failures) if failures else "",
        "expected": {
            "organism_substring": organism_sub,
            "name_keyword": name_kw,
            "length_min": min_len,
            "length_max": max_len,
        },
    }


def main() -> int:
    print("Verifying eleven Kira target UniProt mappings.")
    print("Querying:", UNIPROT_URL.format(accession="{accession}"))
    print()

    records: dict[str, dict] = {}
    failures: list[str] = []

    for i, target in enumerate(TARGETS):
        key = target[0]
        if i > 0:
            time.sleep(REQUEST_DELAY_SECONDS)
        record = verify_one(*target)
        records[key] = record

        if record["status"] == "PASS":
            print(
                f"  PASS  {key:<8s}  {record['accession']:<8s}  "
                f"{record['organism']:<25s}  "
                f"len={record['sequence_length']:<4d}  "
                f"{record['recommended_name']}"
            )
        else:
            failures.append(key)
            print(f"  FAIL  {key:<8s}  {record['accession']:<8s}  {record['reason']}")
            if "organism" in record:
                print(
                    f"        organism={record.get('organism')!r}  "
                    f"name={record.get('recommended_name')!r}  "
                    f"length={record.get('sequence_length')}"
                )

    PROVENANCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with PROVENANCE_PATH.open("w", encoding="utf-8") as fh:
        json.dump(records, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print()
    print(f"Wrote provenance to {PROVENANCE_PATH}")
    print(f"Result: {len(TARGETS) - len(failures)}/{len(TARGETS)} PASS")

    if failures:
        print(f"FAILED: {', '.join(failures)}")
        print("Stopping. Do not proceed to correction.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
