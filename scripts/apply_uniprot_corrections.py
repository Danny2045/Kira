"""Apply verified UniProt mappings to the repository data files.

Reads data/models/uniprot_mappings_provenance.json (produced by
verify_uniprot_mappings.py) and writes:

  - data/models/uniprot_ids_v2.json: {protein_key: accession}
  - data/models/sequences_v2.json:   {accession: {sequence, header,
                                                  length, label}}

Both files use sorted keys and 2-space indentation. The script is
idempotent: running it twice produces identical files.

Aborts if any provenance record is not status=PASS, on the principle
that only verified mappings should be written into the reference data.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROVENANCE_PATH = Path("data/models/uniprot_mappings_provenance.json")
IDS_PATH = Path("data/models/uniprot_ids_v2.json")
SEQS_PATH = Path("data/models/sequences_v2.json")


def build_header(record: dict) -> str:
    accession = record["accession"]
    name = record.get("recommended_name", "")
    organism = record.get("organism", "")
    sv = record.get("sequence_version")
    parts = [f">{accession}", name, f"OS={organism}"]
    if sv is not None:
        parts.append(f"SV={sv}")
    return " ".join(p for p in parts if p)


def main() -> int:
    if not PROVENANCE_PATH.exists():
        print(f"ERROR: {PROVENANCE_PATH} not found. Run verify_uniprot_mappings.py first.")
        return 1

    with PROVENANCE_PATH.open("r", encoding="utf-8") as fh:
        provenance = json.load(fh)

    failed = [k for k, r in provenance.items() if r.get("status") != "PASS"]
    if failed:
        print(f"ERROR: provenance contains non-PASS entries: {', '.join(failed)}")
        print("Refusing to write unverified mappings into reference data.")
        return 1

    uniprot_ids: dict[str, str] = {}
    sequences: dict[str, dict] = {}
    for key, record in provenance.items():
        accession = record["accession"]
        uniprot_ids[key] = accession
        sequences[accession] = {
            "sequence": record["sequence"],
            "header": build_header(record),
            "length": record["sequence_length"],
            "label": key,
        }

    IDS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with IDS_PATH.open("w", encoding="utf-8") as fh:
        json.dump(uniprot_ids, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with SEQS_PATH.open("w", encoding="utf-8") as fh:
        json.dump(sequences, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"Wrote {IDS_PATH} ({len(uniprot_ids)} mappings)")
    print(f"Wrote {SEQS_PATH} ({len(sequences)} sequences)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
