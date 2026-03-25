#!/bin/bash
cat > 19_selectivity_prediction.py << 'PYTHONSCRIPT'
"""
Kira - Script 19: Selectivity Prediction from Protein Embeddings
===================================================================

THE QUESTION:
  Can we predict whether a compound will be selective for a parasite
  enzyme over its human orthologue, using only the protein sequences
  of the parasite and human targets?

WHY THIS MATTERS:
  If selectivity is partially predictable from sequence embeddings,
  it means the structural divergence that enables selective inhibition
  leaves a detectable signature in protein language model representations.
  That would let us triage targets for selectivity BEFORE any compounds
  are tested — a massive acceleration of the repurposing workflow.

APPROACH:
  1. Retrieve protein sequences for all parasite targets and human orthologues
  2. Compute ESM-2 embeddings for each protein
  3. For each target-orthologue pair, compute embedding distance metrics
  4. Correlate embedding distances with experimental selectivity outcomes
  5. Train a simple classifier: given (parasite_embedding, human_embedding,
     compound_fingerprint), predict selectivity class

DATA:
  311 selectivity labels from Kira (Scripts 10, 15, 16, 18)
  7 target-orthologue pairs across 3 diseases

REQUIREMENTS:
  pip install fair-esm torch
  (ESM-2 650M runs on CPU, ~2.5 GB RAM)

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 19_selectivity_prediction.py
"""

import os
import sys
import time
import json
import requests
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

BASE_DIR = os.path.dirname(__file__)
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
LEISH_DIR = os.path.join(BASE_DIR, "data", "leishmania")
MODEL_DIR = os.path.join(BASE_DIR, "data", "models")
PUB_DIR = os.path.join(BASE_DIR, "data", "publication")

# Target-orthologue pairs with UniProt IDs for sequence retrieval
TARGET_PAIRS = [
    {
        "disease": "Schistosomiasis",
        "parasite_target": "Histone deacetylase 8",
        "human_target": "HDAC8",
        "parasite_uniprot": "A0A094ZS53",  # SmHDAC8
        "human_uniprot": "Q9BY41",          # Human HDAC8
        "selectivity_source": "schisto",
    },
    {
        "disease": "Schistosomiasis",
        "parasite_target": "Dihydroorotate dehydrogenase (quinone), mitochondrial",
        "human_target": "DHODH",
        "parasite_uniprot": "G4LZI2",      # SmDHODH
        "human_uniprot": "Q02127",          # Human DHODH
        "selectivity_source": "schisto",
    },
    {
        "disease": "Trypanosomiasis",
        "parasite_target": "Cathepsin B-like cysteine protease",
        "human_target": "Cathepsin B",
        "parasite_uniprot": "Q7YWB5",      # T. brucei cathepsin B
        "human_uniprot": "P07858",          # Human cathepsin B
        "selectivity_source": "tryp",
    },
    {
        "disease": "Trypanosomiasis",
        "parasite_target": "Class 1 phosphodiesterase PDEB1",
        "human_target": "PDE4B",
        "parasite_uniprot": "Q57UX2",      # T. brucei PDEB1
        "human_uniprot": "Q07343",          # Human PDE4B
        "selectivity_source": "tryp",
    },
    {
        "disease": "Leishmaniasis",
        "parasite_target": "Pteridine reductase 1",
        "human_target": "DHFR",
        "parasite_uniprot": "Q01782",      # L. major PTR1
        "human_uniprot": "P00374",          # Human DHFR
        "selectivity_source": "leish",
    },
    {
        "disease": "Leishmaniasis",
        "parasite_target": "Bifunctional dihydrofolate reductase-thymidylate synthase",
        "human_target": "DHFR",
        "parasite_uniprot": "Q01781",      # L. major DHFR-TS
        "human_uniprot": "P00374",          # Human DHFR
        "selectivity_source": "leish",
    },
]


# ---------------------------------------------------------------------------
# STEP 1: Retrieve protein sequences from UniProt
# ---------------------------------------------------------------------------

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            lines = response.text.strip().split("\n")
            header = lines[0]
            sequence = "".join(lines[1:])
            return sequence, header
        else:
            print(f"    UniProt {uniprot_id}: HTTP {response.status_code}")
            return None, None
    except Exception as e:
        print(f"    UniProt {uniprot_id}: {e}")
        return None, None


def retrieve_sequences():
    """Get all protein sequences for target-orthologue pairs."""

    print("\n" + "=" * 60)
    print("STEP 1: RETRIEVING PROTEIN SEQUENCES")
    print("=" * 60)

    os.makedirs(MODEL_DIR, exist_ok=True)
    seq_cache = os.path.join(MODEL_DIR, "sequences.json")

    if os.path.exists(seq_cache):
        with open(seq_cache) as f:
            sequences = json.load(f)
        print(f"  Loaded cached sequences: {len(sequences)} proteins")
        return sequences

    sequences = {}
    all_ids = set()
    for pair in TARGET_PAIRS:
        all_ids.add(pair["parasite_uniprot"])
        all_ids.add(pair["human_uniprot"])

    for uid in sorted(all_ids):
        seq, header = fetch_uniprot_sequence(uid)
        if seq:
            sequences[uid] = {"sequence": seq, "header": header, "length": len(seq)}
            print(f"  {uid}: {len(seq)} aa — {header[:60]}")
        else:
            print(f"  {uid}: FAILED")
        time.sleep(0.5)

    with open(seq_cache, "w") as f:
        json.dump(sequences, f, indent=2)

    print(f"\n  Retrieved {len(sequences)} sequences")
    return sequences


# ---------------------------------------------------------------------------
# STEP 2: Compute ESM-2 embeddings
# ---------------------------------------------------------------------------

def compute_esm_embeddings(sequences):
    """Compute ESM-2 per-protein embeddings."""

    print("\n" + "=" * 60)
    print("STEP 2: COMPUTING ESM-2 EMBEDDINGS")
    print("=" * 60)

    embed_cache = os.path.join(MODEL_DIR, "esm2_embeddings.npz")
    if os.path.exists(embed_cache):
        data = np.load(embed_cache, allow_pickle=True)
        embeddings = {k: data[k] for k in data.files}
        print(f"  Loaded cached embeddings: {len(embeddings)} proteins")
        return embeddings

    try:
        import torch
        import esm
    except ImportError:
        print("  ESM-2 not installed. Run: pip install fair-esm torch")
        print("  Falling back to sequence-based features only.")
        return None

    print("  Loading ESM-2 model (esm2_t33_650M_UR50D)...")
    print("  This may take 1-2 minutes on first run.")

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    embeddings = {}

    for uid, info in sequences.items():
        seq = info["sequence"]
        # Truncate to 1022 tokens (ESM-2 max)
        if len(seq) > 1022:
            seq = seq[:1022]

        print(f"  Embedding {uid} ({len(seq)} aa)...", end=" ", flush=True)

        data = [(uid, seq)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)

        # Mean pooling over sequence length (excluding BOS/EOS)
        token_embeddings = results["representations"][33]
        seq_embedding = token_embeddings[0, 1:len(seq)+1].mean(0).numpy()

        embeddings[uid] = seq_embedding
        print(f"dim={seq_embedding.shape[0]}")

    # Save
    np.savez(embed_cache, **embeddings)
    print(f"\n  Saved embeddings: {embed_cache}")

    return embeddings


# ---------------------------------------------------------------------------
# STEP 3: Compute embedding distances between target pairs
# ---------------------------------------------------------------------------

def compute_pair_distances(embeddings, sequences):
    """Compute distance metrics between parasite and human protein embeddings."""

    print("\n" + "=" * 60)
    print("STEP 3: TARGET PAIR EMBEDDING DISTANCES")
    print("=" * 60)

    pair_features = []

    for pair in TARGET_PAIRS:
        p_uid = pair["parasite_uniprot"]
        h_uid = pair["human_uniprot"]

        record = {
            "disease": pair["disease"],
            "parasite_target": pair["parasite_target"],
            "human_target": pair["human_target"],
            "parasite_uniprot": p_uid,
            "human_uniprot": h_uid,
        }

        # Sequence-based features (always available)
        p_seq = sequences.get(p_uid, {}).get("sequence", "")
        h_seq = sequences.get(h_uid, {}).get("sequence", "")

        record["parasite_length"] = len(p_seq)
        record["human_length"] = len(h_seq)
        record["length_ratio"] = len(p_seq) / max(1, len(h_seq))

        # Simple sequence identity (global alignment would be better but this is fast)
        # Use k-mer overlap as proxy
        def kmer_set(seq, k=3):
            return set(seq[i:i+k] for i in range(len(seq) - k + 1))

        p_kmers = kmer_set(p_seq)
        h_kmers = kmer_set(h_seq)
        if p_kmers and h_kmers:
            jaccard = len(p_kmers & h_kmers) / len(p_kmers | h_kmers)
            record["kmer3_jaccard"] = round(jaccard, 4)
        else:
            record["kmer3_jaccard"] = None

        # Embedding-based features
        if embeddings and p_uid in embeddings and h_uid in embeddings:
            p_emb = embeddings[p_uid]
            h_emb = embeddings[h_uid]

            # Cosine similarity
            cos_sim = np.dot(p_emb, h_emb) / (np.linalg.norm(p_emb) * np.linalg.norm(h_emb))
            record["esm2_cosine_similarity"] = round(float(cos_sim), 4)

            # Euclidean distance
            euclid = np.linalg.norm(p_emb - h_emb)
            record["esm2_euclidean_distance"] = round(float(euclid), 2)

            # L1 distance
            l1 = np.sum(np.abs(p_emb - h_emb))
            record["esm2_l1_distance"] = round(float(l1), 2)

            # Embedding difference vector (for downstream classifier)
            record["embedding_diff_norm"] = round(float(np.linalg.norm(p_emb - h_emb)), 4)
        else:
            record["esm2_cosine_similarity"] = None
            record["esm2_euclidean_distance"] = None

        pair_features.append(record)

    pair_df = pd.DataFrame(pair_features)

    # Print
    print(f"\n  {'Disease':15s} {'Parasite target':30s} {'kmer3':>6s} {'ESM cos':>8s} {'ESM dist':>9s}")
    print(f"  {'-'*15} {'-'*30} {'-'*6} {'-'*8} {'-'*9}")
    for _, r in pair_df.iterrows():
        kmer = f"{r['kmer3_jaccard']:.3f}" if pd.notna(r.get('kmer3_jaccard')) else "N/A"
        cos = f"{r['esm2_cosine_similarity']:.4f}" if pd.notna(r.get('esm2_cosine_similarity')) else "N/A"
        dist = f"{r['esm2_euclidean_distance']:.1f}" if pd.notna(r.get('esm2_euclidean_distance')) else "N/A"
        print(f"  {r['disease']:15s} {r['parasite_target'][:30]:30s} {kmer:>6s} {cos:>8s} {dist:>9s}")

    return pair_df


# ---------------------------------------------------------------------------
# STEP 4: Load all selectivity data and merge with pair features
# ---------------------------------------------------------------------------

def load_selectivity_data():
    """Load all 311 selectivity labels from all three diseases."""

    print("\n" + "=" * 60)
    print("STEP 4: LOADING ALL SELECTIVITY DATA")
    print("=" * 60)

    all_sel = []

    # Schistosomiasis
    path = os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
        df["disease"] = "Schistosomiasis"
        all_sel.append(df)
        print(f"  Schistosomiasis: {len(df)} records")

    # Trypanosomiasis
    path = os.path.join(TRYP_DIR, "tryp_selectivity_expanded.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
        df["disease"] = "Trypanosomiasis"
        all_sel.append(df)
        print(f"  Trypanosomiasis: {len(df)} records")

    # Leishmaniasis
    path = os.path.join(LEISH_DIR, "leish_selectivity.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
        if "disease" not in df.columns:
            df["disease"] = "Leishmaniasis"
        all_sel.append(df)
        print(f"  Leishmaniasis: {len(df)} records")

    if not all_sel:
        print("  No selectivity data found!")
        return pd.DataFrame()

    combined = pd.concat(all_sel, ignore_index=True)
    print(f"\n  Total selectivity records: {len(combined)}")
    print(f"  Class distribution:")
    for cls, n in combined["selectivity_class"].value_counts().items():
        print(f"    {cls}: {n}")

    return combined


# ---------------------------------------------------------------------------
# STEP 5: Correlate embedding distances with selectivity
# ---------------------------------------------------------------------------

def correlate_embeddings_selectivity(pair_df, selectivity_df):
    """Test whether embedding distance predicts selectivity outcome."""

    print("\n" + "=" * 60)
    print("STEP 5: EMBEDDING DISTANCE vs SELECTIVITY")
    print("=" * 60)

    if pair_df is None or "esm2_cosine_similarity" not in pair_df.columns:
        print("  No embedding data. Using sequence features only.")
        return

    # For each target pair, get the aggregate selectivity outcome
    results = []
    for _, pair in pair_df.iterrows():
        ptarget = pair["parasite_target"]

        # Find matching selectivity data
        sel_data = selectivity_df[
            selectivity_df["parasite_target"].str.contains(ptarget[:20], na=False)
        ]

        if len(sel_data) == 0:
            continue

        n_total = len(sel_data)
        n_selective = len(sel_data[sel_data["selectivity_class"] == "SELECTIVE"])
        n_counter = len(sel_data[sel_data["selectivity_class"] == "COUNTER-SELECTIVE"])
        pct_selective = round(100 * n_selective / n_total, 1)
        median_ratio = sel_data["selectivity_ratio"].median()

        results.append({
            "target": ptarget[:35],
            "disease": pair["disease"],
            "n_compounds": n_total,
            "pct_selective": pct_selective,
            "median_ratio": median_ratio,
            "kmer3_jaccard": pair.get("kmer3_jaccard"),
            "esm2_cosine": pair.get("esm2_cosine_similarity"),
            "esm2_distance": pair.get("esm2_euclidean_distance"),
        })

    if not results:
        print("  Could not match targets.")
        return

    res_df = pd.DataFrame(results)

    print(f"\n  {'Target':35s} {'%Sel':>6s} {'MedR':>6s} {'kmer3':>6s} {'ESM cos':>8s} {'ESM dist':>9s}")
    print(f"  {'-'*35} {'-'*6} {'-'*6} {'-'*6} {'-'*8} {'-'*9}")
    for _, r in res_df.iterrows():
        cos = f"{r['esm2_cosine']:.4f}" if pd.notna(r.get('esm2_cosine')) else "N/A"
        dist = f"{r['esm2_distance']:.1f}" if pd.notna(r.get('esm2_distance')) else "N/A"
        kmer = f"{r['kmer3_jaccard']:.3f}" if pd.notna(r.get('kmer3_jaccard')) else "N/A"
        print(f"  {r['target']:35s} {r['pct_selective']:5.1f}% {r['median_ratio']:6.1f} "
              f"{kmer:>6s} {cos:>8s} {dist:>9s}")

    # Correlation analysis
    from scipy import stats

    for feature in ["kmer3_jaccard", "esm2_cosine", "esm2_distance"]:
        valid = res_df.dropna(subset=[feature, "median_ratio"])
        if len(valid) >= 4:
            corr, pval = stats.spearmanr(valid[feature], valid["median_ratio"])
            direction = "lower similarity → higher selectivity" if corr < 0 else "higher similarity → higher selectivity"
            sig = "SIGNIFICANT" if pval < 0.1 else "not significant"
            print(f"\n  {feature} vs median selectivity ratio:")
            print(f"    Spearman rho = {corr:.3f}, p = {pval:.3f} ({sig})")
            print(f"    Interpretation: {direction}")

    # Key finding
    print(f"\n  KEY HYPOTHESIS:")
    print(f"  If ESM-2 cosine similarity between parasite and human proteins")
    print(f"  is LOW, the proteins are structurally/functionally divergent,")
    print(f"  and selective inhibition should be EASIER.")
    print(f"  If this correlation holds, ESM-2 embeddings could predict")
    print(f"  selectivity potential from sequence alone — before any")
    print(f"  compounds are synthesized or tested.")

    return res_df


# ---------------------------------------------------------------------------
# STEP 6: Build compound-level selectivity classifier
# ---------------------------------------------------------------------------

def build_selectivity_classifier(selectivity_df, pair_df, embeddings, sequences):
    """
    Train a simple classifier to predict selectivity class from
    compound fingerprint + target pair features.
    """

    print("\n" + "=" * 60)
    print("STEP 6: SELECTIVITY CLASSIFIER")
    print("=" * 60)

    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.metrics import classification_report, roc_auc_score
    from sklearn.preprocessing import LabelEncoder

    # Binary task: SELECTIVE+MODERATE vs POOR+COUNTER-SELECTIVE
    sel_df = selectivity_df.copy()
    sel_df["binary_label"] = sel_df["selectivity_class"].map({
        "SELECTIVE": 1, "MODERATE": 1,
        "POOR": 0, "COUNTER-SELECTIVE": 0,
    })
    sel_df = sel_df.dropna(subset=["binary_label"])

    print(f"  Binary classification task:")
    print(f"    Selective/Moderate (positive): {(sel_df['binary_label'] == 1).sum()}")
    print(f"    Poor/Counter-selective (negative): {(sel_df['binary_label'] == 0).sum()}")

    # Build features for each compound
    features = []
    labels = []
    compounds_used = []

    for _, row in sel_df.iterrows():
        # Compound fingerprint
        smiles = row.get("smiles") or row.get("canonical_smiles")
        if pd.isna(smiles):
            continue

        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
        fp_array = np.array(fp)

        # Target pair features
        ptarget = row["parasite_target"]
        pair_match = None
        for _, pr in pair_df.iterrows():
            if pr["parasite_target"][:20] in ptarget or ptarget[:20] in str(pr["parasite_target"]):
                pair_match = pr
                break

        if pair_match is None:
            # Use zeros for unknown pairs
            pair_feats = np.zeros(5)
        else:
            pair_feats = np.array([
                pair_match.get("kmer3_jaccard", 0) or 0,
                pair_match.get("esm2_cosine_similarity", 0) or 0,
                pair_match.get("esm2_euclidean_distance", 0) or 0,
                pair_match.get("length_ratio", 1) or 1,
                pair_match.get("parasite_length", 300) / 1000,
            ])

        # Selectivity ratio as feature (log-transformed)
        # NO — this would be data leakage. The ratio IS the label.
        # Use only: fingerprint + target pair features

        feature_vec = np.concatenate([fp_array, pair_feats])
        features.append(feature_vec)
        labels.append(int(row["binary_label"]))
        compounds_used.append(row["molecule_chembl_id"])

    if len(features) < 20:
        print(f"  Only {len(features)} samples with complete features. Too few for ML.")
        return None

    X = np.array(features)
    y = np.array(labels)

    print(f"\n  Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    print(f"    256 Morgan fingerprint bits + 5 target pair features")
    print(f"  Labels: {sum(y)} positive, {len(y) - sum(y)} negative")

    # Stratified 5-fold cross-validation
    print(f"\n  Training Gradient Boosting Classifier (5-fold CV)...")

    clf = GradientBoostingClassifier(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        min_samples_leaf=5,
        random_state=42,
    )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # AUROC
    auroc_scores = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc")
    print(f"  AUROC: {auroc_scores.mean():.3f} +/- {auroc_scores.std():.3f}")
    print(f"  Per-fold: {[f'{s:.3f}' for s in auroc_scores]}")

    # Accuracy
    acc_scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")
    print(f"  Accuracy: {acc_scores.mean():.3f} +/- {acc_scores.std():.3f}")

    # Baseline
    majority_class = max(sum(y), len(y) - sum(y)) / len(y)
    print(f"  Majority class baseline: {majority_class:.3f}")

    # Feature importance (train on full data for interpretation)
    clf.fit(X, y)
    importances = clf.feature_importances_

    # Top features
    top_idx = np.argsort(importances)[::-1][:15]
    print(f"\n  Top 15 feature importances:")
    for i, idx in enumerate(top_idx):
        if idx < 256:
            name = f"FP_bit_{idx}"
        else:
            pair_names = ["kmer3_jaccard", "esm2_cosine", "esm2_distance",
                          "length_ratio", "parasite_length"]
            name = pair_names[idx - 256] if idx - 256 < len(pair_names) else f"feat_{idx}"
        print(f"    {i+1:2d}. {name:25s} {importances[idx]:.4f}")

    # Check if target pair features are in top features
    pair_importance = importances[256:].sum()
    fp_importance = importances[:256].sum()
    print(f"\n  Total importance:")
    print(f"    Compound fingerprint: {fp_importance:.3f} ({100*fp_importance:.1f}%)")
    print(f"    Target pair features: {pair_importance:.3f} ({100*pair_importance:.1f}%)")

    if pair_importance > 0.05:
        print(f"\n  Target pair features contribute meaningfully ({100*pair_importance:.1f}%).")
        print(f"  Selectivity is partially predictable from protein-level features.")
    else:
        print(f"\n  Selectivity prediction is dominated by compound structure.")
        print(f"  Target pair features contribute only {100*pair_importance:.1f}%.")

    return {
        "auroc_mean": auroc_scores.mean(),
        "auroc_std": auroc_scores.std(),
        "accuracy_mean": acc_scores.mean(),
        "baseline": majority_class,
        "n_samples": len(y),
        "pair_importance": pair_importance,
        "fp_importance": fp_importance,
    }


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Script 19: Selectivity Prediction")
    print("  From Protein Embeddings to Selectivity Classifier")
    print("=" * 60)

    start_time = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)

    # Step 1: Sequences
    sequences = retrieve_sequences()

    # Step 2: ESM-2 embeddings
    embeddings = compute_esm_embeddings(sequences)

    # Step 3: Pair distances
    pair_df = compute_pair_distances(embeddings, sequences)
    pair_df.to_csv(os.path.join(MODEL_DIR, "target_pair_features.csv"), index=False)

    # Step 4: Load selectivity data
    selectivity_df = load_selectivity_data()

    # Step 5: Correlate
    correlation_df = correlate_embeddings_selectivity(pair_df, selectivity_df)
    if correlation_df is not None:
        correlation_df.to_csv(os.path.join(MODEL_DIR, "embedding_selectivity_correlation.csv"),
                              index=False)

    # Step 6: Classifier
    clf_results = build_selectivity_classifier(selectivity_df, pair_df, embeddings, sequences)

    elapsed = time.time() - start_time

    print(f"\n{'=' * 60}")
    print(f"  SELECTIVITY PREDICTION COMPLETE")
    print(f"  Time: {elapsed:.0f} seconds")
    print(f"")
    if clf_results:
        print(f"  Classifier AUROC: {clf_results['auroc_mean']:.3f} +/- {clf_results['auroc_std']:.3f}")
        print(f"  Baseline (majority): {clf_results['baseline']:.3f}")
        print(f"  Improvement over baseline: {clf_results['auroc_mean'] - 0.5:.3f}")
        print(f"")
        print(f"  Feature importance split:")
        print(f"    Compound structure: {100*clf_results['fp_importance']:.1f}%")
        print(f"    Target pair (protein-level): {100*clf_results['pair_importance']:.1f}%")
    print(f"")
    print(f"  This is the first attempt to predict antiparasitic selectivity")
    print(f"  from protein language model embeddings + compound fingerprints.")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 19_selectivity_prediction.py successfully."
echo ""
echo "BEFORE RUNNING, install dependencies:"
echo "  pip install fair-esm torch --break-system-packages"
echo ""
echo "Then run:"
echo "  python 19_selectivity_prediction.py"
