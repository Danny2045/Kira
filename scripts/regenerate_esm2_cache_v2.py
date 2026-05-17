"""Regenerate the ESM-2 v2 embeddings cache against the corrected UniProt
mappings committed in PR #17.

Reads the canonical, verified sequences from
``data/models/uniprot_mappings_provenance.json``. Writes the per-protein
mean and standard-deviation embeddings to
``data/models/esm2_embeddings_v2.npz`` and a JSON sidecar with the
generation provenance to
``data/models/esm2_embeddings_v2_provenance.json``.

Model recipe matches archived/scripts/20_selectivity_model_v2.py exactly:

  * Model      : ESM-2 ``esm2_t33_650M_UR50D`` (fair-esm)
  * Pooling    : mean (and std) over residue tokens excluding BOS/EOS,
                 i.e. ``representations[33][0, 1:len(seq)+1]``
  * Layer      : 33 (final transformer layer)
  * Truncation : first 1022 residues (ESM-2 token limit)
  * Output dim : 1280, dtype float32

Closes audit Finding 2 (contaminated protein embeddings cache).

Usage:

    conda activate bio-builder
    KMP_DUPLICATE_LIB_OK=TRUE python scripts/regenerate_esm2_cache_v2.py
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# macOS OpenMP duplicate-symbol workaround used by Script 20.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

REPO_ROOT = Path(__file__).resolve().parents[1]
PROVENANCE_IN = REPO_ROOT / "data" / "models" / "uniprot_mappings_provenance.json"
EMB_OUT = REPO_ROOT / "data" / "models" / "esm2_embeddings_v2.npz"
PROV_OUT = REPO_ROOT / "data" / "models" / "esm2_embeddings_v2_provenance.json"

MODEL_ID = "esm2_t33_650M_UR50D"
LAYER = 33
EMBED_DIM = 1280
MAX_TOKENS = 1022

# Old cache keys that the new cache must NOT contain (contamination from
# pre-PR-17 mappings).
SUPERSEDED_REMOVED = {
    "O00444": ("Human serine/threonine-protein kinase PLK4 — was incorrectly "
               "mapped to SmHDAC8 in the pre-PR-17 cache."),
    "P07711": ("Human cathepsin L — was incorrectly mapped to HsCatB in the "
               "pre-PR-17 cache."),
}
# New cache keys added by this regeneration (correct accessions for
# SmHDAC8, SmDHODH, HsCatB).
SUPERSEDED_ADDED = {
    "A5H660": "SmHDAC8 (Schistosoma mansoni histone deacetylase) — correct mapping per PR #17.",
    "G4VFD7": ("SmDHODH (Schistosoma mansoni dihydroorotate dehydrogenase) — "
               "missing entirely from the pre-PR-17 cache; downstream code "
               "silently zeroed its ESM-2 features."),
    "P07858": "HsCatB (human cathepsin B) — correct human orthologue per PR #17.",
}


def sha256_seq(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def get_git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def get_esm_version() -> str:
    try:
        import esm
        version = getattr(esm, "__version__", None)
        if version:
            return version
        from importlib.metadata import version as pkg_version
        try:
            return pkg_version("fair-esm")
        except Exception:
            return "unknown"
    except Exception:
        return "unknown"


def main() -> None:
    with PROVENANCE_IN.open() as f:
        prov = json.load(f)

    # Order the proteins deterministically by accession.
    entries = sorted(prov.values(), key=lambda v: v["accession"])

    print(f"Loaded {len(entries)} verified UniProt mappings from "
          f"{PROVENANCE_IN.relative_to(REPO_ROOT)}")
    for e in entries:
        print(f"  {e['key']:<10} {e['accession']:<8} "
              f"{e['sequence_length']:>4} aa  {e['organism']}")

    # Lazy-import esm + torch so the script can be linted without them.
    import esm
    import torch  # noqa: F401

    print(f"\nLoading {MODEL_ID} (first call downloads ~2.5 GB) ...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    print("Model ready. Embedding sequences (CPU, no grad) ...")

    embeddings: dict[str, np.ndarray] = {}
    protein_records: list[dict] = []

    for e in entries:
        accession = e["accession"]
        full_seq = e["sequence"]
        seq = full_seq[:MAX_TOKENS]
        truncated = len(full_seq) > MAX_TOKENS

        suffix = f" (truncated to {MAX_TOKENS} for ESM-2)" if truncated else ""
        print(f"\n  {e['key']} ({accession}) — {len(full_seq)} aa{suffix}")

        _, _, batch_tokens = batch_converter([(accession, seq)])
        with torch.no_grad():
            results = model(
                batch_tokens, repr_layers=[LAYER], return_contacts=False
            )
        # Tokens are [BOS, ...residues, EOS]; slice [1:len(seq)+1] excludes both.
        token_emb = results["representations"][LAYER][0, 1:len(seq) + 1]
        mean_emb = token_emb.mean(0).numpy().astype(np.float32)
        std_emb = token_emb.std(0).numpy().astype(np.float32)

        assert mean_emb.shape == (EMBED_DIM,), mean_emb.shape
        assert std_emb.shape == (EMBED_DIM,), std_emb.shape

        embeddings[accession] = mean_emb
        embeddings[f"{accession}_std"] = std_emb

        protein_records.append({
            "key": e["key"],
            "accession": accession,
            "organism": e["organism"],
            "recommended_name": e.get("recommended_name", ""),
            "sequence_length": len(full_seq),
            "sequence_length_embedded": len(seq),
            "truncated": truncated,
            "sequence_sha256": sha256_seq(full_seq),
            "uniprot_sequence_version": e.get("sequence_version"),
            "uniprot_last_sequence_update_date": e.get(
                "last_sequence_update_date"
            ),
            "embedding_mean_norm": float(np.linalg.norm(mean_emb)),
            "embedding_std_mean": float(std_emb.mean()),
        })

        print(f"    mean dim={mean_emb.shape[0]} "
              f"‖mean‖={np.linalg.norm(mean_emb):.3f} "
              f"⟨std⟩={std_emb.mean():.3f}")

    EMB_OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(EMB_OUT, **embeddings)
    print(f"\nWrote {EMB_OUT.relative_to(REPO_ROOT)} "
          f"({len(embeddings)} arrays = {len(entries)} mean + "
          f"{len(entries)} std)")

    provenance_doc = {
        "purpose": (
            "Regenerated ESM-2 embeddings cache against corrected UniProt "
            "mappings from PR #17. Closes audit Finding 2 (contaminated "
            "protein embeddings cache)."
        ),
        "model_id": MODEL_ID,
        "layer": LAYER,
        "embedding_dim": EMBED_DIM,
        "dtype": "float32",
        "pooling": (
            "mean and std over residue tokens [1:len(seq)+1], excluding "
            "BOS/EOS, of representations[33]"
        ),
        "max_sequence_length": MAX_TOKENS,
        "fair_esm_version": get_esm_version(),
        "code_git_sha": get_git_sha(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_provenance_file": str(PROVENANCE_IN.relative_to(REPO_ROOT)),
        "output_npz_file": str(EMB_OUT.relative_to(REPO_ROOT)),
        "supersedes": {
            "removed_keys": sorted(SUPERSEDED_REMOVED.keys()),
            "removed_keys_explanation": SUPERSEDED_REMOVED,
            "added_keys": sorted(SUPERSEDED_ADDED.keys()),
            "added_keys_explanation": SUPERSEDED_ADDED,
        },
        "proteins": protein_records,
    }

    with PROV_OUT.open("w") as f:
        json.dump(provenance_doc, f, indent=2)
    print(f"Wrote {PROV_OUT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
