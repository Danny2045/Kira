#!/bin/bash
cd ~/kira

python3 << 'PATCH'
"""
Script 21: Clean selectivity model rerun.

Two bugs fixed:
1. Pair features concatenation: the compound_feat_dim was set to 266
   (256+10) but then used as the boundary between compound and pair
   features. The pair features were appended after index 266 but
   compound_feat_dim was 266, so feature_vec had length 266+17=283,
   but the importance split used compound_feat_dim=266 as the cutoff.
   ACTUAL BUG: the matching logic failed to find pairs because it
   matched disease name prefixes incorrectly. The disease column in
   leish selectivity was set differently.

2. Target mapping: SmHDAC8 and SmDHODH ChEMBL IDs were wrong.
   Fix: use correct ChEMBL target IDs found from the activity files.

Approach: Instead of relying on ChEMBL API for UniProt IDs (which
gave PLK4 for SmHDAC8), directly use known-good UniProt accessions
verified by organism + protein name from UniProt search.
"""

import os, json, sys
import numpy as np
import pandas as pd
from datetime import datetime
from collections import Counter
import time
import requests

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else os.getcwd()
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
TRYP_DIR = os.path.join(BASE_DIR, "data", "trypanosoma")
LEISH_DIR = os.path.join(BASE_DIR, "data", "leishmania")
MODEL_DIR = os.path.join(BASE_DIR, "data", "models")
os.makedirs(MODEL_DIR, exist_ok=True)

###############################################################################
# VERIFIED UNIPROT ACCESSIONS
# These were verified by searching UniProt for organism + protein name.
# Each has been checked: correct organism, correct protein, correct length.
###############################################################################

VERIFIED_TARGETS = {
    # Parasite targets
    "SmHDAC8": {
        "uniprot": "A0A3Q0KU27",   # S. mansoni HDAC8, 340 aa
        "organism": "Schistosoma mansoni",
        "name": "Histone deacetylase 8",
        "expected_length_range": (300, 400),
    },
    "SmDHODH": {
        "uniprot": "G4VQR5",       # S. mansoni DHODH, ~400 aa
        "organism": "Schistosoma mansoni",
        "name": "Dihydroorotate dehydrogenase",
        "expected_length_range": (300, 500),
    },
    "TbCatB": {
        "uniprot": "Q6R7Z5",       # T. brucei cathepsin B, 340 aa (verified in Script 20)
        "organism": "Trypanosoma brucei",
        "name": "Cathepsin B-like cysteine protease",
        "expected_length_range": (250, 400),
    },
    "TbPDEB1": {
        "uniprot": "Q8WQX9",       # T. brucei PDEB1, 930 aa (verified in Script 20)
        "organism": "Trypanosoma brucei",
        "name": "Phosphodiesterase PDEB1",
        "expected_length_range": (800, 1000),
    },
    "LmPTR1": {
        "uniprot": "Q01782",       # L. major PTR1, 288 aa (verified in Scripts 19+20)
        "organism": "Leishmania major",
        "name": "Pteridine reductase 1",
        "expected_length_range": (250, 320),
    },
    "LmDHFR": {
        "uniprot": "P07382",       # L. major DHFR-TS, 520 aa (verified in Script 20)
        "organism": "Leishmania major",
        "name": "Bifunctional DHFR-TS",
        "expected_length_range": (450, 600),
    },
    # Human targets (all verified in Script 20)
    "HsHDAC8": {
        "uniprot": "Q9BY41",       # Human HDAC8, 377 aa
        "organism": "Homo sapiens",
        "name": "Histone deacetylase 8",
        "expected_length_range": (350, 400),
    },
    "HsDHODH": {
        "uniprot": "Q02127",       # Human DHODH, 395 aa
        "organism": "Homo sapiens",
        "name": "DHODH",
        "expected_length_range": (350, 420),
    },
    "HsCatL": {
        "uniprot": "P07711",       # Human cathepsin L, 333 aa
        "organism": "Homo sapiens",
        "name": "Cathepsin L",
        "expected_length_range": (300, 370),
    },
    "HsPDE4B": {
        "uniprot": "Q07343",       # Human PDE4B, 736 aa
        "organism": "Homo sapiens",
        "name": "PDE4B",
        "expected_length_range": (700, 800),
    },
    "HsDHFR": {
        "uniprot": "P00374",       # Human DHFR, 187 aa
        "organism": "Homo sapiens",
        "name": "DHFR",
        "expected_length_range": (150, 210),
    },
}

TARGET_PAIRS = [
    ("SmHDAC8", "HsHDAC8", "Schistosomiasis", "deacetylase"),
    ("SmDHODH", "HsDHODH", "Schistosomiasis", "dihydroorotate"),
    ("TbCatB", "HsCatL", "Trypanosomiasis", "cathepsin"),
    ("TbPDEB1", "HsPDE4B", "Trypanosomiasis", "phosphodiesterase"),
    ("LmPTR1", "HsDHFR", "Leishmaniasis", "pteridine"),
    ("LmDHFR", "HsDHFR", "Leishmaniasis", "dihydrofolate"),
]

###############################################################################
# STEP 1: Fetch and verify sequences
###############################################################################

print("=" * 60)
print("  KIRA -- Script 21: Clean Selectivity Model")
print("  Fixed pair features + verified target mapping")
print("=" * 60)

seq_cache = os.path.join(MODEL_DIR, "sequences_v3.json")
emb_cache = os.path.join(MODEL_DIR, "esm2_embeddings_v3.npz")

# Fetch sequences
if os.path.exists(seq_cache):
    with open(seq_cache) as f:
        sequences = json.load(f)
    print(f"\n  Loaded cached sequences: {len(sequences)}")
else:
    print("\n  Fetching sequences from UniProt...")
    sequences = {}
    for label, info in VERIFIED_TARGETS.items():
        uid = info["uniprot"]
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200:
                lines = resp.text.strip().split("\n")
                header = lines[0]
                seq = "".join(lines[1:])
                lo, hi = info["expected_length_range"]
                length_ok = lo <= len(seq) <= hi
                status = "OK" if length_ok else f"WARNING: {len(seq)} aa outside [{lo},{hi}]"
                sequences[uid] = {"sequence": seq, "header": header, "length": len(seq), "label": label}
                print(f"  {label:12s} ({uid}): {len(seq):4d} aa — {status}")
            else:
                print(f"  {label:12s} ({uid}): HTTP {resp.status_code}")
        except Exception as e:
            print(f"  {label:12s} ({uid}): {e}")
        time.sleep(0.3)

    with open(seq_cache, "w") as f:
        json.dump(sequences, f, indent=2)

# Verify all targets resolved
for label, info in VERIFIED_TARGETS.items():
    uid = info["uniprot"]
    if uid not in sequences:
        print(f"  MISSING: {label} ({uid})")

###############################################################################
# STEP 2: ESM-2 embeddings
###############################################################################

if os.path.exists(emb_cache):
    data = np.load(emb_cache, allow_pickle=True)
    embeddings = {k: data[k] for k in data.files}
    print(f"\n  Loaded cached embeddings: {len(embeddings)}")
else:
    print("\n  Computing ESM-2 embeddings...")
    try:
        import torch, esm
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        batch_converter = alphabet.get_batch_converter()
        model.eval()

        embeddings = {}
        for uid, info in sequences.items():
            seq = info["sequence"][:1022]
            label = info.get("label", uid)
            print(f"  Embedding {label} ({len(seq)} aa)...", end=" ", flush=True)
            batch_data = [(uid, seq)]
            _, _, batch_tokens = batch_converter(batch_data)
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            emb = results["representations"][33][0, 1:len(seq)+1].mean(0).numpy()
            embeddings[uid] = emb
            print(f"dim={emb.shape[0]}")

        np.savez(emb_cache, **embeddings)
    except ImportError:
        print("  ESM-2 not available.")
        embeddings = {}

###############################################################################
# STEP 3: Compute pair features (FIXED)
###############################################################################

print("\n" + "=" * 60)
print("STEP 3: COMPUTING PAIR FEATURES")
print("=" * 60)

def compute_pair_features(p_uid, h_uid):
    """Compute all features for a target pair."""
    feats = {}
    p_info = sequences.get(p_uid, {})
    h_info = sequences.get(h_uid, {})
    p_seq = p_info.get("sequence", "")
    h_seq = h_info.get("sequence", "")

    # Sequence features
    feats["para_len"] = len(p_seq) / 1000
    feats["human_len"] = len(h_seq) / 1000
    feats["len_ratio"] = len(p_seq) / max(1, len(h_seq))
    feats["len_diff"] = abs(len(p_seq) - len(h_seq)) / 1000

    # k-mer overlap
    for k in [2, 3, 4]:
        s1 = set(p_seq[i:i+k] for i in range(len(p_seq)-k+1))
        s2 = set(h_seq[i:i+k] for i in range(len(h_seq)-k+1))
        feats[f"kmer{k}"] = len(s1 & s2) / max(1, len(s1 | s2))

    # AA composition distance
    def aa_comp(seq):
        c = Counter(seq)
        t = max(1, len(seq))
        return {aa: c.get(aa, 0)/t for aa in "ACDEFGHIKLMNPQRSTVWY"}
    pc, hc = aa_comp(p_seq), aa_comp(h_seq)
    feats["aa_l2"] = sum((pc.get(a,0)-hc.get(a,0))**2 for a in "ACDEFGHIKLMNPQRSTVWY")

    # ESM-2 features
    p_emb = embeddings.get(p_uid)
    h_emb = embeddings.get(h_uid)
    if p_emb is not None and h_emb is not None:
        cos = np.dot(p_emb, h_emb) / (np.linalg.norm(p_emb) * np.linalg.norm(h_emb) + 1e-8)
        feats["esm_cos"] = float(cos)
        feats["esm_euclid"] = float(np.linalg.norm(p_emb - h_emb))
        feats["esm_l1"] = float(np.sum(np.abs(p_emb - h_emb)))
        feats["esm_pnorm"] = float(np.linalg.norm(p_emb))
        feats["esm_hnorm"] = float(np.linalg.norm(h_emb))
    else:
        feats["esm_cos"] = 0
        feats["esm_euclid"] = 0
        feats["esm_l1"] = 0
        feats["esm_pnorm"] = 0
        feats["esm_hnorm"] = 0

    return feats

# Precompute for all pairs
pair_features = {}
PAIR_FEAT_KEYS = None

for p_label, h_label, disease, keyword in TARGET_PAIRS:
    p_uid = VERIFIED_TARGETS[p_label]["uniprot"]
    h_uid = VERIFIED_TARGETS[h_label]["uniprot"]
    feats = compute_pair_features(p_uid, h_uid)
    pair_features[(keyword, disease)] = feats
    if PAIR_FEAT_KEYS is None:
        PAIR_FEAT_KEYS = sorted(feats.keys())
    cos = feats.get("esm_cos", 0)
    print(f"  {p_label:12s} - {h_label:8s}: ESM cos={cos:.4f}, {len(feats)} features")

N_PAIR_FEATS = len(PAIR_FEAT_KEYS)
print(f"\n  Pair feature names ({N_PAIR_FEATS}): {PAIR_FEAT_KEYS}")

###############################################################################
# STEP 4: Build dataset (FIXED concatenation)
###############################################################################

print("\n" + "=" * 60)
print("STEP 4: BUILDING FEATURE MATRIX (FIXED)")
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
        for _, r in act.iterrows():
            smi = r.get("canonical_smiles")
            if pd.notna(smi):
                all_smiles[r["molecule_chembl_id"]] = smi
print(f"  SMILES lookup: {len(all_smiles)} compounds")

# Load selectivity
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
selectivity_df = pd.concat(all_sel, ignore_index=True)
print(f"  Selectivity records: {len(selectivity_df)}")

# Build feature vectors
X_list, y_list, disease_list, cpd_list = [], [], [], []
matched_count = 0
unmatched_count = 0

for _, row in selectivity_df.iterrows():
    cid = row["molecule_chembl_id"]
    smi = all_smiles.get(cid)
    if smi is None:
        continue

    mol = Chem.MolFromSmiles(str(smi))
    if mol is None:
        continue

    # Compound features: 256 FP + 10 descriptors = 266
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
    fp_arr = np.array(fp, dtype=np.float32)

    try:
        desc = np.array([
            Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol), Descriptors.NumAromaticRings(mol),
            Descriptors.FractionCSP3(mol), QED.qed(mol),
        ], dtype=np.float32)
    except:
        desc = np.zeros(10, dtype=np.float32)

    cpd_feats = np.concatenate([fp_arr, desc])  # 266

    # Match to target pair using keyword + disease
    ptarget = str(row.get("parasite_target", "")).lower()
    disease = str(row.get("disease", ""))

    matched_pair = None
    for p_label, h_label, pair_disease, keyword in TARGET_PAIRS:
        if keyword in ptarget and disease == pair_disease:
            matched_pair = (keyword, pair_disease)
            break

    if matched_pair is None:
        # Fallback: match by keyword only
        for p_label, h_label, pair_disease, keyword in TARGET_PAIRS:
            if keyword in ptarget:
                matched_pair = (keyword, pair_disease)
                break

    if matched_pair is not None and matched_pair in pair_features:
        pf = pair_features[matched_pair]
        pair_arr = np.array([pf.get(k, 0) for k in PAIR_FEAT_KEYS], dtype=np.float32)
        matched_count += 1
    else:
        pair_arr = np.zeros(N_PAIR_FEATS, dtype=np.float32)
        unmatched_count += 1

    # CONCATENATE: compound (266) + pair (N_PAIR_FEATS)
    feature_vec = np.concatenate([cpd_feats, pair_arr])
    X_list.append(feature_vec)

    label = 1 if row["selectivity_class"] in ["SELECTIVE", "MODERATE"] else 0
    y_list.append(label)
    disease_list.append(disease)
    cpd_list.append(cid)

X = np.array(X_list, dtype=np.float32)
y = np.array(y_list)
diseases = np.array(disease_list)
X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

CPD_DIM = 266
TOTAL_DIM = CPD_DIM + N_PAIR_FEATS

print(f"\n  Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
print(f"    Compound features: {CPD_DIM} (256 FP + 10 desc)")
print(f"    Pair features: {N_PAIR_FEATS}")
print(f"    Total: {TOTAL_DIM}")
print(f"  Pair matching: {matched_count} matched, {unmatched_count} unmatched")
print(f"  Labels: {sum(y)} positive, {len(y)-sum(y)} negative")
print(f"  Diseases: {pd.Series(diseases).value_counts().to_dict()}")

# VERIFY pair features are non-zero
pair_cols = X[:, CPD_DIM:]
pair_nonzero = (pair_cols != 0).any(axis=0).sum()
print(f"  Pair feature columns with non-zero values: {pair_nonzero}/{N_PAIR_FEATS}")
if pair_nonzero == 0:
    print("  ERROR: Pair features are all zero! Bug not fixed.")
    sys.exit(1)
else:
    print(f"  CONFIRMED: Pair features successfully included in feature matrix.")

###############################################################################
# STEP 5: Train and evaluate
###############################################################################

print("\n" + "=" * 60)
print("STEP 5: MODEL TRAINING AND EVALUATION")
print("=" * 60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# [A] Gradient Boosting - compound only
print("\n  [A] COMPOUND-ONLY vs COMPOUND+PAIR (5-fold CV)")

gb_cpd = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                      min_samples_leaf=5, subsample=0.8, random_state=42)
auroc_cpd = cross_val_score(gb_cpd, X[:, :CPD_DIM], y, cv=cv, scoring="roc_auc")

gb_full = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                       min_samples_leaf=5, subsample=0.8, random_state=42)
auroc_full = cross_val_score(gb_full, X, y, cv=cv, scoring="roc_auc")

print(f"  Compound only:  AUROC {auroc_cpd.mean():.3f} +/- {auroc_cpd.std():.3f}")
print(f"  Compound+Pair:  AUROC {auroc_full.mean():.3f} +/- {auroc_full.std():.3f}")
print(f"  Improvement:    {auroc_full.mean() - auroc_cpd.mean():+.3f}")
print(f"  Baseline:       {max(sum(y), len(y)-sum(y))/len(y):.3f}")

# MLP
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

mlp = MLPClassifier(hidden_layer_sizes=(128, 64, 32), activation="relu",
                     max_iter=500, learning_rate_init=0.001,
                     early_stopping=True, validation_fraction=0.15, random_state=42)
auroc_mlp = cross_val_score(mlp, X_scaled, y, cv=cv, scoring="roc_auc")
print(f"\n  MLP (compound+pair): AUROC {auroc_mlp.mean():.3f} +/- {auroc_mlp.std():.3f}")

# [B] Leave-one-disease-out
print("\n  [B] LEAVE-ONE-DISEASE-OUT (the real test)")

for held_out in ["Schistosomiasis", "Trypanosomiasis", "Leishmaniasis"]:
    train_mask = diseases != held_out
    test_mask = diseases == held_out

    if sum(test_mask) < 5 or len(np.unique(y[test_mask])) < 2:
        print(f"\n  {held_out}: skipped")
        continue

    X_tr, y_tr = X[train_mask], y[train_mask]
    X_te, y_te = X[test_mask], y[test_mask]

    # Compound only
    gb1 = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                       min_samples_leaf=5, subsample=0.8, random_state=42)
    gb1.fit(X_tr[:, :CPD_DIM], y_tr)
    auroc1 = roc_auc_score(y_te, gb1.predict_proba(X_te[:, :CPD_DIM])[:, 1])

    # Compound + pair
    gb2 = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                       min_samples_leaf=5, subsample=0.8, random_state=42)
    gb2.fit(X_tr, y_tr)
    auroc2 = roc_auc_score(y_te, gb2.predict_proba(X_te)[:, 1])

    print(f"\n  Held out: {held_out} ({sum(test_mask)} samples, {y_te.mean():.0%} positive)")
    print(f"    Compound only:  AUROC {auroc1:.3f}")
    print(f"    Compound+Pair:  AUROC {auroc2:.3f}")
    print(f"    Pair improvement: {auroc2 - auroc1:+.3f}")

# [C] Feature importance
print("\n  [C] FEATURE IMPORTANCE")
gb_full.fit(X, y)
imp = gb_full.feature_importances_

fp_imp = imp[:256].sum()
desc_imp = imp[256:266].sum()
pair_imp = imp[266:].sum()

print(f"\n  Feature group importance:")
print(f"    Morgan FP (256):    {fp_imp:.3f} ({100*fp_imp:.1f}%)")
print(f"    RDKit desc (10):    {desc_imp:.3f} ({100*desc_imp:.1f}%)")
print(f"    Pair features ({N_PAIR_FEATS}):  {pair_imp:.3f} ({100*pair_imp:.1f}%)")
print(f"    ── Compound total:  {fp_imp+desc_imp:.3f} ({100*(fp_imp+desc_imp):.1f}%)")
print(f"    ── Protein total:   {pair_imp:.3f} ({100*pair_imp:.1f}%)")

# Top pair features
pair_feat_imp = [(PAIR_FEAT_KEYS[i], imp[CPD_DIM + i]) for i in range(N_PAIR_FEATS)]
pair_feat_imp.sort(key=lambda x: x[1], reverse=True)
print(f"\n  Pair feature ranking:")
for name, importance in pair_feat_imp:
    print(f"    {name:20s} {importance:.4f}")

# Top 15 overall
all_names = [f"FP_{i}" for i in range(256)]
all_names += ["MW", "LogP", "HBD", "HBA", "TPSA", "RotBonds", "Rings", "AromRings", "Fsp3", "QED"]
all_names += PAIR_FEAT_KEYS
top_idx = np.argsort(imp)[::-1][:15]
print(f"\n  Top 15 features overall:")
for rank, idx in enumerate(top_idx, 1):
    name = all_names[idx] if idx < len(all_names) else f"feat_{idx}"
    group = "PAIR" if idx >= CPD_DIM else ("DESC" if idx >= 256 else "FP")
    print(f"    {rank:2d}. [{group:4s}] {name:20s} {imp[idx]:.4f}")

print(f"\n{'=' * 60}")
print(f"  CLEAN RERUN COMPLETE")
print(f"{'=' * 60}")
PATCH
