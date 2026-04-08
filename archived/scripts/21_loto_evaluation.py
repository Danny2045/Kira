"""
Kira - Script 21: Leave-One-Target-pair-Out (LOTO) Evaluation
==============================================================

Evaluates how well the selectivity model generalises to unseen target pairs.
For each of the 6 target pairs, train on the remaining 5 and test on the
held-out pair. Reports AUROC per fold.

Uses esm2_embeddings_v2.npz (clean embeddings from Script 20).

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 21_loto_evaluation.py
"""

import os, json, numpy as np, pandas as pd
from collections import Counter

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR     = os.path.dirname(__file__)
PROCESSED    = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR     = os.path.join(BASE_DIR, "data", "trypanosoma")
LEISH_DIR    = os.path.join(BASE_DIR, "data", "leishmania")
MODEL_DIR    = os.path.join(BASE_DIR, "data", "models")

EMB_PATH     = os.path.join(MODEL_DIR, "esm2_embeddings_v2.npz")
UNIPROT_PATH = os.path.join(MODEL_DIR, "uniprot_ids_v2.json")
SEQ_PATH     = os.path.join(MODEL_DIR, "sequences_v2.json")

# ---------------------------------------------------------------------------
# Target pair registry
# ---------------------------------------------------------------------------
# Each tuple: (parasite_label, human_label, disease, keyword_in_parasite_target)
TARGET_PAIRS = [
    ("SmHDAC8",  "HsHDAC8",  "Schistosomiasis",  "deacetylase"),
    ("SmDHODH",  "HsDHODH",  "Schistosomiasis",  "dihydroorotate"),
    ("TbCatB",   "HsCatB",   "Trypanosomiasis",  "cathepsin"),
    ("TbPDEB1",  "HsPDE4B",  "Trypanosomiasis",  "phosphodiesterase"),
    ("LmPTR1",   "HsDHFR",   "Leishmaniasis",    "pteridine"),
    ("LmDHFR",   "HsDHFR",   "Leishmaniasis",    "dihydrofolate"),
]

# ---------------------------------------------------------------------------
# Feature helpers (same as Script 20)
# ---------------------------------------------------------------------------

def compound_features(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256), dtype=np.float32)
    try:
        desc = np.array([
            Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol), Descriptors.NumAromaticRings(mol),
            Descriptors.FractionCSP3(mol), QED.qed(mol),
        ], dtype=np.float32)
    except Exception:
        desc = np.zeros(10, dtype=np.float32)
    return np.concatenate([fp, desc])  # dim 266


def kmer_jaccard(s1, s2, k):
    k1 = set(s1[i:i+k] for i in range(len(s1)-k+1))
    k2 = set(s2[i:i+k] for i in range(len(s2)-k+1))
    if not k1 or not k2:
        return 0.0
    return len(k1 & k2) / len(k1 | k2)


def aa_l2(s1, s2):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    c1, c2 = Counter(s1), Counter(s2)
    t1, t2 = max(1, len(s1)), max(1, len(s2))
    return sum((c1.get(a,0)/t1 - c2.get(a,0)/t2)**2 for a in aas)


def pair_features(p_uid, h_uid, sequences, embeddings):
    """Return a fixed-length feature vector for a target pair."""
    p_seq = sequences.get(p_uid, {}).get("sequence", "")
    h_seq = sequences.get(h_uid, {}).get("sequence", "")

    feats = [
        len(p_seq) / 1000,
        len(h_seq) / 1000,
        len(p_seq) / max(1, len(h_seq)),
        abs(len(p_seq) - len(h_seq)) / 1000,
        kmer_jaccard(p_seq, h_seq, 2),
        kmer_jaccard(p_seq, h_seq, 3),
        kmer_jaccard(p_seq, h_seq, 4),
        aa_l2(p_seq, h_seq),
    ]

    # ESM-2 features (zeros when embedding is unavailable, e.g. SmDHODH)
    if p_uid in embeddings and h_uid in embeddings:
        pe, he = embeddings[p_uid], embeddings[h_uid]
        cos = float(np.dot(pe, he) / (np.linalg.norm(pe) * np.linalg.norm(he) + 1e-8))
        feats += [
            cos,
            float(np.linalg.norm(pe - he)),
            float(np.sum(np.abs(pe - he))),
            float(np.linalg.norm(pe)),
            float(np.linalg.norm(he)),
        ]
        p_std_key, h_std_key = f"{p_uid}_std", f"{h_uid}_std"
        if p_std_key in embeddings and h_std_key in embeddings:
            feats += [
                float(embeddings[p_std_key].mean()),
                float(embeddings[h_std_key].mean()),
                float(embeddings[p_std_key].mean() / (embeddings[h_std_key].mean() + 1e-8)),
            ]
        else:
            feats += [0.0, 0.0, 0.0]
    else:
        feats += [0.0] * 8  # placeholders (cosine, euclidean, l1, norms x2, std x3)

    return np.array(feats, dtype=np.float32)  # dim 16


# ---------------------------------------------------------------------------
# Load everything
# ---------------------------------------------------------------------------

def load_data():
    print("Loading embeddings and sequences...")
    embeddings = dict(np.load(EMB_PATH))
    with open(UNIPROT_PATH) as f:
        uniprot = json.load(f)
    with open(SEQ_PATH) as f:
        seqs_raw = json.load(f)
    # sequences keyed by UniProt ID
    sequences = {v["label"]: v for v in seqs_raw.values() if "label" in v}
    # also make available by UniProt ID
    for uid, v in seqs_raw.items():
        sequences[uid] = v
    print(f"  Embeddings: {len([k for k in embeddings if '_std' not in k])} proteins")
    print(f"  UniProt IDs cached: {list(uniprot.keys())}")

    print("\nLoading selectivity data...")
    frames = []
    for path, disease in [
        (os.path.join(PROCESSED, "kira_selectivity_analysis.csv"), "Schistosomiasis"),
        (os.path.join(TRYP_DIR,  "tryp_selectivity_expanded.csv"), "Trypanosomiasis"),
        (os.path.join(LEISH_DIR, "leish_selectivity.csv"),         "Leishmaniasis"),
    ]:
        df = pd.read_csv(path)
        df["disease"] = disease
        frames.append(df)
        print(f"  {disease}: {len(df)} records")
    sel_df = pd.concat(frames, ignore_index=True)
    print(f"  Total: {len(sel_df)} records")

    print("\nLoading SMILES...")
    smiles_map = {}
    for path in [
        os.path.join(PROCESSED, "schisto_filtered_activities.csv"),
        os.path.join(TRYP_DIR,  "tryp_activities.csv"),
        os.path.join(LEISH_DIR, "leish_activities.csv"),
    ]:
        if os.path.exists(path):
            act = pd.read_csv(path)
            if "canonical_smiles" in act.columns:
                for _, r in act.iterrows():
                    if pd.notna(r.get("canonical_smiles")):
                        smiles_map[r["molecule_chembl_id"]] = r["canonical_smiles"]
    print(f"  SMILES: {len(smiles_map)} compounds")

    return embeddings, uniprot, sequences, sel_df, smiles_map


# ---------------------------------------------------------------------------
# Assign target pairs to rows
# ---------------------------------------------------------------------------

def assign_pairs(sel_df):
    """Tag each row with its target pair index."""
    pair_col = []
    for _, row in sel_df.iterrows():
        ptarget = str(row.get("parasite_target", "")).lower()
        disease = str(row.get("disease", ""))
        assigned = None
        for idx, (p_label, h_label, pair_disease, keyword) in enumerate(TARGET_PAIRS):
            if keyword.lower() in ptarget:
                if pair_disease[:5].lower() in disease[:5].lower() or disease == "":
                    assigned = idx
                    break
        pair_col.append(assigned)
    return pair_col


# ---------------------------------------------------------------------------
# Build full feature matrix
# ---------------------------------------------------------------------------

def build_feature_matrix(sel_df, pair_col, embeddings, uniprot, sequences, smiles_map):
    # Precompute pair feature vectors for each of the 6 pairs
    pair_feat_cache = {}
    for idx, (p_label, h_label, disease, keyword) in enumerate(TARGET_PAIRS):
        p_uid = uniprot.get(p_label)
        h_uid = uniprot.get(h_label)
        pf = pair_features(p_uid or "", h_uid or "", sequences, embeddings)
        pair_feat_cache[idx] = pf
        has_emb = (p_uid in embeddings) and (h_uid in embeddings)
        print(f"  [{idx}] {p_label}/{h_label}: pair_feats dim={len(pf)}, ESM-2={'yes' if has_emb else 'NO'}")

    X_list, y_list, pair_idx_list = [], [], []
    skipped_no_smiles, skipped_bad_mol, skipped_no_pair = 0, 0, 0

    for i, (_, row) in enumerate(sel_df.iterrows()):
        pidx = pair_col[i]
        if pidx is None:
            skipped_no_pair += 1
            continue

        cid = row["molecule_chembl_id"]
        smiles = smiles_map.get(cid)
        if smiles is None:
            skipped_no_smiles += 1
            continue

        cpd = compound_features(smiles)
        if cpd is None:
            skipped_bad_mol += 1
            continue

        pf = pair_feat_cache[pidx]
        X_list.append(np.concatenate([cpd, pf]))

        sel_class = row["selectivity_class"]
        y_list.append(1 if sel_class in ("SELECTIVE", "MODERATE") else 0)
        pair_idx_list.append(pidx)

    X = np.nan_to_num(np.array(X_list, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    y = np.array(y_list)
    pair_idx = np.array(pair_idx_list)

    print(f"\n  Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    print(f"    Compound features: 266 (256-bit Morgan FP + 10 RDKit descriptors)")
    print(f"    Pair features:     16 (sequence stats + ESM-2 distances)")
    print(f"  Labels: {y.sum()} positive ({100*y.mean():.1f}%), {(1-y).sum()} negative")
    print(f"  Skipped: {skipped_no_pair} no-pair, {skipped_no_smiles} no-SMILES, {skipped_bad_mol} bad-mol")

    return X, y, pair_idx


# ---------------------------------------------------------------------------
# LOTO evaluation
# ---------------------------------------------------------------------------

def loto_evaluation(X, y, pair_idx):
    print("\n" + "="*60)
    print("LEAVE-ONE-TARGET-PAIR-OUT EVALUATION")
    print("Train on 5 pairs → test on held-out 6th pair")
    print("="*60)

    results = []

    for held_idx, (p_label, h_label, disease, keyword) in enumerate(TARGET_PAIRS):
        test_mask  = pair_idx == held_idx
        train_mask = pair_idx != held_idx

        n_test  = test_mask.sum()
        n_train = train_mask.sum()

        pair_name = f"{p_label}/{h_label}"

        if n_test == 0:
            print(f"\n[{held_idx}] {pair_name}: NO TEST DATA — skipping")
            results.append({"pair": pair_name, "disease": disease,
                            "n_train": n_train, "n_test": 0,
                            "auroc": None, "note": "no test data"})
            continue

        y_test = y[test_mask]

        if len(np.unique(y_test)) < 2:
            only = "all positive" if y_test.sum() == len(y_test) else "all negative"
            print(f"\n[{held_idx}] {pair_name}: SINGLE CLASS in test ({only}) — AUROC undefined")
            results.append({"pair": pair_name, "disease": disease,
                            "n_train": n_train, "n_test": int(n_test),
                            "auroc": None, "note": f"single class: {only}"})
            continue

        if len(np.unique(y[train_mask])) < 2:
            print(f"\n[{held_idx}] {pair_name}: SINGLE CLASS in train — skipping")
            results.append({"pair": pair_name, "disease": disease,
                            "n_train": n_train, "n_test": int(n_test),
                            "auroc": None, "note": "single class in train"})
            continue

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_mask])
        X_test  = scaler.transform(X[test_mask])

        gb = GradientBoostingClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            min_samples_leaf=5, subsample=0.8, random_state=42,
        )
        gb.fit(X_train, y[train_mask])
        y_prob = gb.predict_proba(X_test)[:, 1]
        auroc  = roc_auc_score(y_test, y_prob)

        pos_rate = y_test.mean()
        baseline = max(pos_rate, 1 - pos_rate)

        print(f"\n[{held_idx}] Held-out: {pair_name} ({disease})")
        print(f"    Train: {n_train} samples from other 5 pairs")
        print(f"    Test:  {n_test} samples  |  pos rate: {pos_rate:.2f}  |  baseline: {baseline:.3f}")
        print(f"    AUROC: {auroc:.3f}")

        results.append({
            "pair": pair_name, "disease": disease,
            "n_train": int(n_train), "n_test": int(n_test),
            "pos_rate": float(pos_rate), "baseline": float(baseline),
            "auroc": float(auroc),
        })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("="*60)
    print("  KIRA -- Script 21: Leave-One-Target-pair-Out Evaluation")
    print("="*60 + "\n")

    embeddings, uniprot, sequences, sel_df, smiles_map = load_data()

    print("\nAssigning target pairs to rows...")
    pair_col = assign_pairs(sel_df)
    assigned = [p for p in pair_col if p is not None]
    print(f"  Assigned: {len(assigned)}/{len(pair_col)}")
    counts = Counter(assigned)
    for idx, (p_label, h_label, disease, kw) in enumerate(TARGET_PAIRS):
        print(f"    [{idx}] {p_label}/{h_label}: {counts.get(idx, 0)} rows")

    print("\nBuilding feature matrix...")
    X, y, pair_idx = build_feature_matrix(sel_df, pair_col, embeddings, uniprot, sequences, smiles_map)

    results = loto_evaluation(X, y, pair_idx)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'Pair':<22} {'Disease':<18} {'N_test':>7} {'Baseline':>9} {'AUROC':>7}  Note")
    print("-"*75)
    for r in results:
        auroc_str    = f"{r['auroc']:.3f}"    if r["auroc"]    is not None else "  N/A "
        baseline_str = f"{r['baseline']:.3f}" if r.get("baseline") is not None else "   N/A"
        note = r.get("note", "")
        print(f"{r['pair']:<22} {r['disease']:<18} {r['n_test']:>7} {baseline_str:>9} {auroc_str:>7}  {note}")

    valid = [r for r in results if r["auroc"] is not None]
    if valid:
        mean_auroc = np.mean([r["auroc"] for r in valid])
        print(f"\nMean AUROC across {len(valid)} evaluable pairs: {mean_auroc:.3f}")
        print(f"(Random baseline ≈ 0.500)")

    # Save
    out_path = os.path.join(MODEL_DIR, "loto_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")
    print("="*60)
