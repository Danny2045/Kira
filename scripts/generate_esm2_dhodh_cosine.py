"""Recompute the SmDHODH / HsDHODH ESM-2 cosine against the corrected sequences.

The original case study reported a 0.9897 ESM-2 cosine similarity between
SmDHODH and HsDHODH. That value was computed against C1L5Z2 mislabeled as
SmDHODH (C1L5Z2 is not S. mansoni DHODH). PR #17 corrected the UniProt
mapping to G4VFD7 (379 aa) and committed the canonical sequence to
``data/models/sequences_v2.json``. This script regenerates the ESM-2 cosine,
euclidean L2 distance, and 3-mer Jaccard similarity against the correct
sequences and writes a structured provenance artifact.

Outputs
-------
``data/models/esm2_embeddings_dhodh.npz``
    Mean-pooled ESM-2 embeddings for G4VFD7 and Q02127 (dim 1280 each).
``data/models/esm2_dhodh_cosine.json``
    Provenance record: model id, both UniProt IDs, sequence lengths and
    SHA-256 hashes, embedding dim, cosine similarity, euclidean L2,
    3-mer Jaccard, UTC timestamp, git SHA, and a ``supersedes`` block
    recording the prior incorrect numbers.

Run from the repo root with the ``bio-builder`` conda env active::

    python scripts/generate_esm2_dhodh_cosine.py

The first invocation downloads the esm2_t33_650M_UR50D weights (~2.5 GB)
into the local torch hub cache.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import esm
import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parents[1]
SEQUENCES_PATH = REPO_ROOT / "data" / "models" / "sequences_v2.json"
EMBED_OUT = REPO_ROOT / "data" / "models" / "esm2_embeddings_dhodh.npz"
JSON_OUT = REPO_ROOT / "data" / "models" / "esm2_dhodh_cosine.json"

PARASITE_UID = "G4VFD7"  # SmDHODH, Schistosoma mansoni
HUMAN_UID = "Q02127"     # HsDHODH, Homo sapiens

MODEL_ID = "esm2_t33_650M_UR50D"
REPR_LAYER = 33
ESM2_MAX_TOKENS = 1022

# Previous (incorrect) values from src/kira/experiments/case_study_chembl155771.py
PREVIOUS_COSINE = 0.9897
PREVIOUS_DISTANCE = 6.09
PREVIOUS_JACCARD = 0.0326


def load_sequences() -> dict[str, str]:
    with SEQUENCES_PATH.open() as fh:
        data = json.load(fh)
    out: dict[str, str] = {}
    for uid in (PARASITE_UID, HUMAN_UID):
        if uid not in data:
            raise KeyError(f"{uid} missing from {SEQUENCES_PATH}")
        out[uid] = data[uid]["sequence"]
    return out


def embed(model, batch_converter, uid: str, seq: str) -> np.ndarray:
    truncated = seq[:ESM2_MAX_TOKENS]
    _, _, batch_tokens = batch_converter([(uid, truncated)])
    with torch.no_grad():
        results = model(
            batch_tokens, repr_layers=[REPR_LAYER], return_contacts=False
        )
    token_emb = results["representations"][REPR_LAYER][0, 1 : len(truncated) + 1]
    return token_emb.mean(0).cpu().numpy().astype(np.float32)


def kmer_jaccard(a: str, b: str, k: int = 3) -> float:
    ka = {a[i : i + k] for i in range(len(a) - k + 1)}
    kb = {b[i : i + k] for i in range(len(b) - k + 1)}
    if not ka or not kb:
        return 0.0
    return len(ka & kb) / len(ka | kb)


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def git_sha() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def main() -> None:
    print(f"Loading sequences from {SEQUENCES_PATH.relative_to(REPO_ROOT)}")
    sequences = load_sequences()
    p_seq = sequences[PARASITE_UID]
    h_seq = sequences[HUMAN_UID]
    print(f"  {PARASITE_UID} (SmDHODH): {len(p_seq)} aa")
    print(f"  {HUMAN_UID} (HsDHODH): {len(h_seq)} aa")

    print(f"\nLoading ESM-2 model: {MODEL_ID}")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    print("Computing mean-pooled embeddings (layer 33, excluding BOS/EOS)...")
    p_emb = embed(model, batch_converter, PARASITE_UID, p_seq)
    h_emb = embed(model, batch_converter, HUMAN_UID, h_seq)
    print(f"  {PARASITE_UID}: dim={p_emb.shape[0]}")
    print(f"  {HUMAN_UID}: dim={h_emb.shape[0]}")

    cosine = float(
        np.dot(p_emb, h_emb) / (np.linalg.norm(p_emb) * np.linalg.norm(h_emb))
    )
    euclidean = float(np.linalg.norm(p_emb - h_emb))
    jaccard = float(kmer_jaccard(p_seq, h_seq, k=3))

    artifact = {
        "model_id": MODEL_ID,
        "model_source": "fair-esm",
        "repr_layer": REPR_LAYER,
        "embedding_dim": int(p_emb.shape[0]),
        "pooling": (
            "mean over token positions [1:len(seq)+1], excluding BOS and EOS"
        ),
        "max_tokens": ESM2_MAX_TOKENS,
        "cosine_formula": "dot(p, h) / (norm(p) * norm(h))",
        "parasite": {
            "uniprot_id": PARASITE_UID,
            "label": "SmDHODH",
            "organism": "Schistosoma mansoni",
            "sequence_length": len(p_seq),
            "sequence_sha256": sha256(p_seq),
        },
        "human": {
            "uniprot_id": HUMAN_UID,
            "label": "HsDHODH",
            "organism": "Homo sapiens",
            "sequence_length": len(h_seq),
            "sequence_sha256": sha256(h_seq),
        },
        "cosine_similarity": round(cosine, 6),
        "embedding_distance": round(euclidean, 6),
        "kmer3_jaccard": round(jaccard, 6),
        "computed_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"
        ),
        "code_git_sha": git_sha(),
        "supersedes": {
            "previous_cosine_similarity": PREVIOUS_COSINE,
            "previous_embedding_distance": PREVIOUS_DISTANCE,
            "previous_kmer3_jaccard": PREVIOUS_JACCARD,
            "note": (
                "Previous values were computed against UniProt C1L5Z2 "
                "mislabeled as SmDHODH. C1L5Z2 is not Schistosoma mansoni "
                "dihydroorotate dehydrogenase; the correct SmDHODH "
                "accession is G4VFD7 (379 aa, mitochondrial DHODH). See "
                "data/models/uniprot_mappings_provenance.json and PR #17 "
                "for the verification audit."
            ),
        },
    }

    EMBED_OUT.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting embeddings -> {EMBED_OUT.relative_to(REPO_ROOT)}")
    np.savez(EMBED_OUT, **{PARASITE_UID: p_emb, HUMAN_UID: h_emb})

    print(f"Writing artifact   -> {JSON_OUT.relative_to(REPO_ROOT)}")
    JSON_OUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print()
    print("=" * 60)
    print("RESULTS (SmDHODH G4VFD7 vs HsDHODH Q02127)")
    print("=" * 60)
    print(f"  ESM-2 cosine similarity : {cosine:.4f}")
    print(f"  ESM-2 L2 distance       : {euclidean:.4f}")
    print(f"  k-mer 3 Jaccard         : {jaccard:.4f}")
    print()
    print(
        f"  Previous (incorrect)    : cos={PREVIOUS_COSINE}, "
        f"L2={PREVIOUS_DISTANCE}, J={PREVIOUS_JACCARD}"
    )
    print(f"  Δ cosine vs old         : {cosine - PREVIOUS_COSINE:+.4f}")


if __name__ == "__main__":
    main()
