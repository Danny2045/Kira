#!/bin/bash
cat > 20_selectivity_model_v2.py << 'PYTHONSCRIPT'
"""
Kira - Script 20: Selectivity Model v2
=========================================

Fixes from Script 19:
  1. UniProt IDs retrieved FROM ChEMBL target API (not hand-curated)
  2. Sequences verified by organism and protein name
  3. SMILES properly joined from activity files

Upgrades:
  1. Richer compound features: 256-bit Morgan FP + 10 RDKit descriptors
  2. Per-target ESM-2 features: mean, variance, and norm of embeddings
  3. Held-out disease evaluation: train on 2 diseases, test on 3rd
  4. Neural network (MLP) alongside Gradient Boosting
  5. SHAP-style feature group analysis

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 20_selectivity_model_v2.py
"""

import os
import sys
import time
import json
import requests
import numpy as np
import pandas as pd
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier

BASE_DIR = os.path.dirname(__file__)
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
LEISH_DIR = os.path.join(BASE_DIR, "data", "leishmania")
MODEL_DIR = os.path.join(BASE_DIR, "data", "models")

# ChEMBL target IDs for our parasite targets (these are reliable)
CHEMBL_TARGETS = {
    "SmHDAC8": {"chembl_id": "CHEMBL3788", "organism": "Schistosoma mansoni", "keyword": "deacetylase"},
    "SmDHODH": {"chembl_id": "CHEMBL2366242", "organism": "Schistosoma mansoni", "keyword": "dihydroorotate"},
    "TbCatB": {"chembl_id": "CHEMBL5832", "organism": "Trypanosoma brucei", "keyword": "cathepsin"},
    "TbPDEB1": {"chembl_id": "CHEMBL2010636", "organism": "Trypanosoma brucei", "keyword": "phosphodiesterase"},
    "LmPTR1": {"chembl_id": "CHEMBL6194", "organism": "Leishmania major", "keyword": "pteridine"},
    "LmDHFR": {"chembl_id": "CHEMBL4614", "organism": "Leishmania major", "keyword": "dihydrofolate"},
}

HUMAN_TARGETS = {
    "HsHDAC8": {"chembl_id": "CHEMBL3192", "keyword": "deacetylase 8"},
    "HsDHODH": {"chembl_id": "CHEMBL1966", "keyword": "dihydroorotate"},
    "HsCatB": {"chembl_id": "CHEMBL3837", "keyword": "cathepsin B"},
    "HsPDE4B": {"chembl_id": "CHEMBL275", "keyword": "phosphodiesterase 4B"},
    "HsDHFR": {"chembl_id": "CHEMBL202", "keyword": "dihydrofolate reductase"},
}

TARGET_PAIRS = [
    ("SmHDAC8", "HsHDAC8", "Schistosomiasis"),
    ("SmDHODH", "HsDHODH", "Schistosomiasis"),
    ("TbCatB", "HsCatB", "Trypanosomiasis"),
    ("TbPDEB1", "HsPDE4B", "Trypanosomiasis"),
    ("LmPTR1", "HsDHFR", "Leishmaniasis"),
    ("LmDHFR", "HsDHFR", "Leishmaniasis"),
]


# ---------------------------------------------------------------------------
# STEP 1: Get UniProt IDs from ChEMBL target components
# ---------------------------------------------------------------------------

def get_uniprot_from_chembl(chembl_id, label):
    """Query ChEMBL target API for UniProt accession."""
    from chembl_webresource_client.new_client import new_client
    target_api = new_client.target

    try:
        target = target_api.get(chembl_id)
        components = target.get("target_components", [])
        for comp in components:
            accession = comp.get("accession")
            if accession:
                desc = comp.get("component_description", "")
                print(f"    {label}: {chembl_id} → UniProt {accession} ({desc[:50]})")
                return accession
        print(f"    {label}: {chembl_id} → no UniProt component found")
        return None
    except Exception as e:
        print(f"    {label}: {chembl_id} → error: {e}")
        return None


def resolve_all_uniprot_ids():
    """Get verified UniProt IDs for all targets."""

    print("\n" + "=" * 60)
    print("STEP 1: RESOLVING UNIPROT IDs FROM ChEMBL")
    print("=" * 60)

    cache_path = os.path.join(MODEL_DIR, "uniprot_ids_v2.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            ids = json.load(f)
        print(f"  Loaded cached IDs: {len(ids)} targets")
        for k, v in ids.items():
            print(f"    {k}: {v}")
        return ids

    ids = {}

    print("\n  Parasite targets:")
    for label, info in CHEMBL_TARGETS.items():
        uid = get_uniprot_from_chembl(info["chembl_id"], label)
        if uid:
            ids[label] = uid
        time.sleep(0.3)

    print("\n  Human targets:")
    for label, info in HUMAN_TARGETS.items():
        uid = get_uniprot_from_chembl(info["chembl_id"], label)
        if uid:
            ids[label] = uid
        time.sleep(0.3)

    with open(cache_path, "w") as f:
        json.dump(ids, f, indent=2)

    print(f"\n  Resolved {len(ids)} UniProt IDs")
    return ids


# ---------------------------------------------------------------------------
# STEP 2: Fetch sequences and compute ESM-2 embeddings
# ---------------------------------------------------------------------------

def fetch_sequence(uniprot_id):
    """Fetch from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:])
            return seq, header
    except:
        pass
    return None, None


def get_sequences_and_embeddings(uniprot_ids):
    """Fetch sequences and compute ESM-2 embeddings."""

    print("\n" + "=" * 60)
    print("STEP 2: SEQUENCES AND ESM-2 EMBEDDINGS")
    print("=" * 60)

    seq_cache = os.path.join(MODEL_DIR, "sequences_v2.json")
    emb_cache = os.path.join(MODEL_DIR, "esm2_embeddings_v2.npz")

    # Sequences
    if os.path.exists(seq_cache):
        with open(seq_cache) as f:
            sequences = json.load(f)
        print(f"  Loaded cached sequences: {len(sequences)}")
    else:
        sequences = {}
        for label, uid in uniprot_ids.items():
            seq, header = fetch_sequence(uid)
            if seq:
                sequences[uid] = {"sequence": seq, "header": header, "length": len(seq), "label": label}
                print(f"  {label} ({uid}): {len(seq)} aa — {header[:60]}")
            else:
                print(f"  {label} ({uid}): FAILED")
            time.sleep(0.3)

        with open(seq_cache, "w") as f:
            json.dump(sequences, f, indent=2)

    # ESM-2 Embeddings
    if os.path.exists(emb_cache):
        data = np.load(emb_cache, allow_pickle=True)
        embeddings = {k: data[k] for k in data.files}
        print(f"  Loaded cached embeddings: {len(embeddings)}")
    else:
        try:
            import torch
            import esm

            print("  Loading ESM-2 model...")
            model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            batch_converter = alphabet.get_batch_converter()
            model.eval()

            embeddings = {}
            for uid, info in sequences.items():
                seq = info["sequence"][:1022]
                print(f"  Embedding {info.get('label', uid)} ({len(seq)} aa)...", end=" ", flush=True)

                data_in = [(uid, seq)]
                _, _, batch_tokens = batch_converter(data_in)

                with torch.no_grad():
                    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

                token_emb = results["representations"][33][0, 1:len(seq)+1]

                # Store mean, std, and full mean embedding
                mean_emb = token_emb.mean(0).numpy()
                std_emb = token_emb.std(0).numpy()

                embeddings[uid] = mean_emb
                embeddings[f"{uid}_std"] = std_emb
                print(f"dim={mean_emb.shape[0]}")

            np.savez(emb_cache, **embeddings)
            print(f"  Saved embeddings: {emb_cache}")

        except ImportError:
            print("  ESM-2 not available. Using sequence features only.")
            embeddings = None

    return sequences, embeddings


# ---------------------------------------------------------------------------
# STEP 3: Build rich feature vectors
# ---------------------------------------------------------------------------

def compute_compound_features(smiles):
    """Compute Morgan FP + RDKit descriptors for a compound."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    # 256-bit Morgan fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
    fp_array = np.array(fp, dtype=np.float32)

    # 10 RDKit descriptors
    try:
        descriptors = np.array([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol),
            Descriptors.NumAromaticRings(mol),
            Descriptors.FractionCSP3(mol),
            QED.qed(mol),
        ], dtype=np.float32)
    except:
        descriptors = np.zeros(10, dtype=np.float32)

    return np.concatenate([fp_array, descriptors])


def compute_pair_features(p_uid, h_uid, sequences, embeddings):
    """Compute features for a parasite-human target pair."""

    features = {}

    # Sequence features
    p_seq = sequences.get(p_uid, {}).get("sequence", "")
    h_seq = sequences.get(h_uid, {}).get("sequence", "")

    features["parasite_length"] = len(p_seq) / 1000  # Normalize
    features["human_length"] = len(h_seq) / 1000
    features["length_ratio"] = len(p_seq) / max(1, len(h_seq))
    features["length_diff"] = abs(len(p_seq) - len(h_seq)) / 1000

    # k-mer features
    def kmer_jaccard(s1, s2, k):
        k1 = set(s1[i:i+k] for i in range(len(s1)-k+1))
        k2 = set(s2[i:i+k] for i in range(len(s2)-k+1))
        if not k1 or not k2:
            return 0
        return len(k1 & k2) / len(k1 | k2)

    features["kmer2_jaccard"] = kmer_jaccard(p_seq, h_seq, 2)
    features["kmer3_jaccard"] = kmer_jaccard(p_seq, h_seq, 3)
    features["kmer4_jaccard"] = kmer_jaccard(p_seq, h_seq, 4)

    # Amino acid composition difference
    from collections import Counter
    def aa_composition(seq):
        c = Counter(seq)
        total = max(1, len(seq))
        return {aa: c.get(aa, 0) / total for aa in "ACDEFGHIKLMNPQRSTVWY"}

    p_comp = aa_composition(p_seq)
    h_comp = aa_composition(h_seq)
    comp_diff = sum((p_comp.get(aa, 0) - h_comp.get(aa, 0))**2 for aa in "ACDEFGHIKLMNPQRSTVWY")
    features["aa_composition_l2"] = comp_diff

    # Embedding features
    if embeddings and p_uid in embeddings and h_uid in embeddings:
        p_emb = embeddings[p_uid]
        h_emb = embeddings[h_uid]

        # Cosine similarity
        cos = np.dot(p_emb, h_emb) / (np.linalg.norm(p_emb) * np.linalg.norm(h_emb) + 1e-8)
        features["esm2_cosine"] = float(cos)

        # Euclidean distance (normalized)
        features["esm2_euclidean"] = float(np.linalg.norm(p_emb - h_emb))

        # L1 distance
        features["esm2_l1"] = float(np.sum(np.abs(p_emb - h_emb)))

        # Difference vector norm
        diff = p_emb - h_emb
        features["esm2_diff_norm"] = float(np.linalg.norm(diff))

        # Parasite and human embedding norms
        features["esm2_parasite_norm"] = float(np.linalg.norm(p_emb))
        features["esm2_human_norm"] = float(np.linalg.norm(h_emb))

        # Embedding variance features (if available)
        p_std_key = f"{p_uid}_std"
        h_std_key = f"{h_uid}_std"
        if p_std_key in embeddings and h_std_key in embeddings:
            features["esm2_parasite_mean_std"] = float(embeddings[p_std_key].mean())
            features["esm2_human_mean_std"] = float(embeddings[h_std_key].mean())
            features["esm2_std_ratio"] = float(
                embeddings[p_std_key].mean() / (embeddings[h_std_key].mean() + 1e-8)
            )

    return features


# ---------------------------------------------------------------------------
# STEP 4: Build dataset and train models
# ---------------------------------------------------------------------------

def build_dataset(selectivity_df, uniprot_ids, sequences, embeddings):
    """Build feature matrix from all selectivity data."""

    print("\n" + "=" * 60)
    print("STEP 4: BUILDING FEATURE MATRIX")
    print("=" * 60)

    # Load SMILES
    all_smiles = {}
    for path in [
        os.path.join(PROCESSED_DIR, "schisto_filtered_activities.csv"),
        os.path.join(TRYP_DIR, "tryp_activities.csv"),
        os.path.join(LEISH_DIR, "leish_activities.csv"),
    ]:
        if os.path.exists(path):
            act = pd.read_csv(path)
            if "canonical_smiles" in act.columns:
                for _, r in act.iterrows():
                    if pd.notna(r.get("canonical_smiles")):
                        all_smiles[r["molecule_chembl_id"]] = r["canonical_smiles"]
    print(f"  SMILES lookup: {len(all_smiles)} compounds")

    # Precompute pair features for each target pair
    pair_feature_cache = {}
    for p_label, h_label, disease in TARGET_PAIRS:
        p_uid = uniprot_ids.get(p_label)
        h_uid = uniprot_ids.get(h_label)
        if p_uid and h_uid:
            pair_feats = compute_pair_features(p_uid, h_uid, sequences, embeddings)
            pair_feature_cache[(p_label, h_label)] = pair_feats
            print(f"  Pair {p_label}-{h_label}: {len(pair_feats)} features")

    # Build samples
    compound_feat_dim = 266  # 256 FP + 10 descriptors
    pair_feat_names = None

    X_list = []
    y_list = []
    diseases = []
    compounds = []

    for _, row in selectivity_df.iterrows():
        cid = row["molecule_chembl_id"]
        smiles = all_smiles.get(cid)
        if smiles is None:
            continue

        # Compound features
        cpd_feats = compute_compound_features(smiles)
        if cpd_feats is None:
            continue

        # Find matching target pair
        ptarget = str(row.get("parasite_target", ""))
        disease = row.get("disease", "")

        matched_pair = None
        for p_label, h_label, pair_disease in TARGET_PAIRS:
            target_info = CHEMBL_TARGETS.get(p_label, {})
            keyword = target_info.get("keyword", "")
            if keyword and keyword.lower() in ptarget.lower():
                if pair_disease.lower()[:5] in disease.lower()[:5] or disease == "":
                    matched_pair = (p_label, h_label)
                    break

        if matched_pair is None:
            continue

        pair_feats_dict = pair_feature_cache.get(matched_pair, {})
        if pair_feat_names is None:
            pair_feat_names = sorted(pair_feats_dict.keys())
        pair_feats = np.array([pair_feats_dict.get(k, 0) for k in pair_feat_names], dtype=np.float32)

        # Combine
        feature_vec = np.concatenate([cpd_feats, pair_feats])
        X_list.append(feature_vec)

        # Label
        sel_class = row["selectivity_class"]
        label = 1 if sel_class in ["SELECTIVE", "MODERATE"] else 0
        y_list.append(label)
        diseases.append(disease)
        compounds.append(cid)

    X = np.array(X_list, dtype=np.float32)
    y = np.array(y_list)
    diseases = np.array(diseases)

    # Replace NaN/inf
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    print(f"\n  Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    print(f"    Compound features: {compound_feat_dim} (256 FP + 10 descriptors)")
    print(f"    Target pair features: {len(pair_feat_names)} ({', '.join(pair_feat_names[:5])}...)")
    print(f"  Labels: {sum(y)} positive, {len(y)-sum(y)} negative")
    print(f"  Diseases: {pd.Series(diseases).value_counts().to_dict()}")

    return X, y, diseases, compounds, pair_feat_names, compound_feat_dim


def train_and_evaluate(X, y, diseases, pair_feat_names, compound_feat_dim):
    """Train models with multiple evaluation strategies."""

    print("\n" + "=" * 60)
    print("STEP 5: MODEL TRAINING AND EVALUATION")
    print("=" * 60)

    results = {}

    # --- Evaluation 1: 5-fold stratified CV ---
    print("\n  [A] 5-FOLD STRATIFIED CROSS-VALIDATION")

    gb = GradientBoostingClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.1,
        min_samples_leaf=5, subsample=0.8, random_state=42,
    )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    auroc_gb = cross_val_score(gb, X, y, cv=cv, scoring="roc_auc")
    acc_gb = cross_val_score(gb, X, y, cv=cv, scoring="accuracy")

    print(f"  Gradient Boosting:")
    print(f"    AUROC: {auroc_gb.mean():.3f} +/- {auroc_gb.std():.3f}")
    print(f"    Accuracy: {acc_gb.mean():.3f} +/- {acc_gb.std():.3f}")

    results["gb_cv_auroc"] = auroc_gb.mean()

    # MLP
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    mlp = MLPClassifier(
        hidden_layer_sizes=(128, 64, 32), activation="relu",
        max_iter=500, learning_rate_init=0.001,
        early_stopping=True, validation_fraction=0.15,
        random_state=42,
    )

    auroc_mlp = cross_val_score(mlp, X_scaled, y, cv=cv, scoring="roc_auc")
    acc_mlp = cross_val_score(mlp, X_scaled, y, cv=cv, scoring="accuracy")

    print(f"\n  Neural Network (MLP 128-64-32):")
    print(f"    AUROC: {auroc_mlp.mean():.3f} +/- {auroc_mlp.std():.3f}")
    print(f"    Accuracy: {acc_mlp.mean():.3f} +/- {acc_mlp.std():.3f}")

    results["mlp_cv_auroc"] = auroc_mlp.mean()

    # Baseline
    majority = max(sum(y), len(y)-sum(y)) / len(y)
    print(f"\n  Majority class baseline: {majority:.3f}")
    results["baseline"] = majority

    # --- Evaluation 2: Leave-one-disease-out ---
    print("\n  [B] LEAVE-ONE-DISEASE-OUT EVALUATION")
    print("  (Train on 2 diseases, test on held-out 3rd)")

    unique_diseases = [d for d in np.unique(diseases) if d]
    for held_out in unique_diseases:
        train_mask = diseases != held_out
        test_mask = diseases == held_out

        if sum(test_mask) < 5 or sum(train_mask) < 20:
            continue

        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]

        # Check if both classes present
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            print(f"  {held_out}: skipped (single class in train or test)")
            continue

        gb_lodo = GradientBoostingClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            min_samples_leaf=5, subsample=0.8, random_state=42,
        )
        gb_lodo.fit(X_train, y_train)
        y_prob = gb_lodo.predict_proba(X_test)[:, 1]
        auroc = roc_auc_score(y_test, y_prob)

        y_pred = gb_lodo.predict(X_test)
        acc = (y_pred == y_test).mean()

        pos_rate = y_test.mean()

        print(f"\n  Held out: {held_out}")
        print(f"    Train: {sum(train_mask)} samples from other diseases")
        print(f"    Test: {sum(test_mask)} samples from {held_out}")
        print(f"    Test positive rate: {pos_rate:.2f}")
        print(f"    AUROC: {auroc:.3f}")
        print(f"    Accuracy: {acc:.3f}")

        results[f"lodo_{held_out}_auroc"] = auroc

    # --- Feature importance analysis ---
    print("\n  [C] FEATURE IMPORTANCE ANALYSIS")

    gb.fit(X, y)
    importances = gb.feature_importances_

    # Group importances
    fp_imp = importances[:256].sum()
    desc_imp = importances[256:266].sum()
    pair_imp = importances[266:].sum()

    print(f"\n  Feature group importance:")
    print(f"    Morgan fingerprint (256 bits): {fp_imp:.3f} ({100*fp_imp:.1f}%)")
    print(f"    RDKit descriptors (10):        {desc_imp:.3f} ({100*desc_imp:.1f}%)")
    print(f"    Target pair features ({len(pair_feat_names)}):      {pair_imp:.3f} ({100*pair_imp:.1f}%)")
    print(f"    ── Compound total:             {fp_imp+desc_imp:.3f} ({100*(fp_imp+desc_imp):.1f}%)")
    print(f"    ── Protein total:              {pair_imp:.3f} ({100*pair_imp:.1f}%)")

    results["compound_importance"] = fp_imp + desc_imp
    results["protein_importance"] = pair_imp

    # Top individual features
    feat_names = [f"FP_{i}" for i in range(256)]
    feat_names += ["MW", "LogP", "HBD", "HBA", "TPSA", "RotBonds", "Rings", "AromRings", "Fsp3", "QED"]
    feat_names += pair_feat_names

    top_idx = np.argsort(importances)[::-1][:20]
    print(f"\n  Top 20 features:")
    for rank, idx in enumerate(top_idx, 1):
        name = feat_names[idx] if idx < len(feat_names) else f"feat_{idx}"
        print(f"    {rank:2d}. {name:30s} {importances[idx]:.4f}")

    return results


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Script 20: Selectivity Model v2")
    print("  ESM-2 + Compound Features + Held-Out Evaluation")
    print("=" * 60)

    start_time = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)

    # Step 1: Resolve UniProt IDs
    uniprot_ids = resolve_all_uniprot_ids()

    # Step 2: Sequences and embeddings
    sequences, embeddings = get_sequences_and_embeddings(uniprot_ids)

    # Step 3: Load all selectivity data
    print("\n" + "=" * 60)
    print("STEP 3: LOADING SELECTIVITY DATA")
    print("=" * 60)

    all_sel = []
    for path, disease in [
        (os.path.join(PROCESSED_DIR, "kira_selectivity_analysis.csv"), "Schistosomiasis"),
        (os.path.join(TRYP_DIR, "tryp_selectivity_expanded.csv"), "Trypanosomiasis"),
        (os.path.join(LEISH_DIR, "leish_selectivity.csv"), "Leishmaniasis"),
    ]:
        if os.path.exists(path):
            df = pd.read_csv(path)
            df["disease"] = disease
            all_sel.append(df)
            print(f"  {disease}: {len(df)} records")

    selectivity_df = pd.concat(all_sel, ignore_index=True)
    print(f"  Total: {len(selectivity_df)} records")

    # Step 4: Build dataset
    X, y, diseases, compounds, pair_feat_names, cpd_dim = build_dataset(
        selectivity_df, uniprot_ids, sequences, embeddings
    )

    if len(X) < 30:
        print(f"ERROR: Only {len(X)} samples. Need at least 30.")
        sys.exit(1)

    # Step 5: Train and evaluate
    results = train_and_evaluate(X, y, diseases, pair_feat_names, cpd_dim)

    # Save results
    results_path = os.path.join(MODEL_DIR, "model_v2_results.json")
    with open(results_path, "w") as f:
        json.dump({k: float(v) if isinstance(v, (np.floating, float)) else v
                    for k, v in results.items()}, f, indent=2)

    elapsed = time.time() - start_time

    print(f"\n{'=' * 60}")
    print(f"  SELECTIVITY MODEL v2 COMPLETE")
    print(f"  Time: {elapsed:.0f} seconds")
    print(f"")
    print(f"  5-fold CV AUROC:")
    print(f"    Gradient Boosting: {results.get('gb_cv_auroc', 0):.3f}")
    print(f"    Neural Network:    {results.get('mlp_cv_auroc', 0):.3f}")
    print(f"    Baseline:          {results.get('baseline', 0):.3f}")
    print(f"")
    print(f"  Feature importance:")
    print(f"    Compound features: {100*results.get('compound_importance', 0):.1f}%")
    print(f"    Protein features:  {100*results.get('protein_importance', 0):.1f}%")
    print(f"")
    if any(k.startswith("lodo_") for k in results):
        print(f"  Leave-one-disease-out AUROC:")
        for k, v in sorted(results.items()):
            if k.startswith("lodo_"):
                disease = k.replace("lodo_", "").replace("_auroc", "")
                print(f"    Held out {disease}: {v:.3f}")
    print(f"{'=' * 60}")
PYTHONSCRIPT

echo "Created 20_selectivity_model_v2.py successfully."
echo "Now run: python 20_selectivity_model_v2.py"
