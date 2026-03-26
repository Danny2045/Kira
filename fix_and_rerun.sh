#!/bin/bash
cd ~/kira

# Delete old caches
rm -f data/models/sequences_v3.json data/models/esm2_embeddings_v3.npz
echo "Cleared sequence and embedding caches."

# Patch Script 21's inline python
# The SmDHODH UniProt ID G4VQR5 returned 0 aa. Correct ID is G4VFD7 (379 aa).
python3 -c "
# Read the script that was created by create_script21.sh
# It's an inline python block, so we need to check what's on disk
import os

# The script 21 code lives as inline python executed by bash
# We need to find where SmDHODH UniProt is defined and fix it

# Check if there's a standalone .py file
candidates = ['20_selectivity_model_v2.py', '21_clean_model.py']
for f in candidates:
    if os.path.exists(f):
        with open(f, 'r') as fh:
            code = fh.read()
        if 'G4VQR5' in code:
            code = code.replace('G4VQR5', 'G4VFD7')
            with open(f, 'w') as fh:
                fh.write(code)
            print(f'  Fixed G4VQR5 → G4VFD7 in {f}')
        elif 'G4VFD7' in code:
            print(f'  {f} already has correct ID')
        else:
            print(f'  {f} does not contain SmDHODH UniProt ID')
print('Done checking .py files.')
"

# Now re-run Script 21 as a fresh inline python with the correct ID
python3 << 'CLEANRUN'
import os, json, sys, time, requests
import numpy as np
import pandas as pd
from collections import Counter

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier

BASE = os.getcwd()
PROCESSED = os.path.join(BASE, "data", "processed")
TRYP = os.path.join(BASE, "data", "trypanosoma")
LEISH = os.path.join(BASE, "data", "leishmania")
MDIR = os.path.join(BASE, "data", "models")
os.makedirs(MDIR, exist_ok=True)

TARGETS = {
    "SmHDAC8":  {"uid": "A0A3Q0KU27", "expect": (300, 450)},
    "SmDHODH":  {"uid": "G4VFD7",     "expect": (350, 420)},  # FIXED
    "TbCatB":   {"uid": "Q6R7Z5",     "expect": (300, 400)},
    "TbPDEB1":  {"uid": "Q8WQX9",     "expect": (850, 1000)},
    "LmPTR1":   {"uid": "Q01782",     "expect": (250, 320)},
    "LmDHFR":   {"uid": "P07382",     "expect": (450, 600)},
    "HsHDAC8":  {"uid": "Q9BY41",     "expect": (350, 400)},
    "HsDHODH":  {"uid": "Q02127",     "expect": (350, 420)},
    "HsCatL":   {"uid": "P07711",     "expect": (300, 370)},
    "HsPDE4B":  {"uid": "Q07343",     "expect": (700, 800)},
    "HsDHFR":   {"uid": "P00374",     "expect": (150, 210)},
}

PAIRS = [
    ("SmHDAC8", "HsHDAC8", "Schistosomiasis", "deacetylase"),
    ("SmDHODH", "HsDHODH", "Schistosomiasis", "dihydroorotate"),
    ("TbCatB",  "HsCatL",  "Trypanosomiasis", "cathepsin"),
    ("TbPDEB1", "HsPDE4B", "Trypanosomiasis", "phosphodiesterase"),
    ("LmPTR1",  "HsDHFR",  "Leishmaniasis",   "pteridine"),
    ("LmDHFR",  "HsDHFR",  "Leishmaniasis",   "dihydrofolate"),
]

print("=" * 60)
print("  KIRA -- Script 21b: Clean Model (SmDHODH FIXED)")
print("=" * 60)

# ---- STEP 1: Sequences ----
print("\n  STEP 1: Fetching sequences...")
sequences = {}
all_ok = True
for label, info in TARGETS.items():
    uid = info["uid"]
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        lines = resp.text.strip().split("\n")
        header = lines[0]
        seq = "".join(lines[1:])
        lo, hi = info["expect"]
        ok = lo <= len(seq) <= hi
        status = "OK" if ok else f"WARNING ({len(seq)} outside [{lo},{hi}])"
        if not ok: all_ok = False
        sequences[uid] = {"sequence": seq, "header": header, "length": len(seq), "label": label}
        print(f"  {label:12s} ({uid:12s}): {len(seq):4d} aa  {status}")
    except Exception as e:
        print(f"  {label:12s} ({uid:12s}): FAILED — {e}")
        all_ok = False
    time.sleep(0.3)

if not all_ok:
    print("\n  WARNING: Some sequences outside expected range. Proceeding anyway.")

# ---- STEP 2: ESM-2 Embeddings ----
print("\n  STEP 2: Computing ESM-2 embeddings...")
try:
    import torch, esm
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    bc = alphabet.get_batch_converter()
    model.eval()

    embeddings = {}
    for uid, info in sequences.items():
        seq = info["sequence"][:1022]
        label = info.get("label", uid)
        print(f"  Embedding {label:12s} ({len(seq)} aa)...", end=" ", flush=True)
        _, _, tokens = bc([(uid, seq)])
        with torch.no_grad():
            res = model(tokens, repr_layers=[33], return_contacts=False)
        emb = res["representations"][33][0, 1:len(seq)+1].mean(0).numpy()
        embeddings[uid] = emb
        print(f"dim={emb.shape[0]}")
except ImportError:
    print("  ESM-2 not available!")
    sys.exit(1)

# ---- STEP 3: Pair features ----
print("\n  STEP 3: Computing pair features...")

def pair_feats(p_uid, h_uid):
    f = {}
    ps = sequences.get(p_uid, {}).get("sequence", "")
    hs = sequences.get(h_uid, {}).get("sequence", "")
    f["para_len"] = len(ps)/1000
    f["human_len"] = len(hs)/1000
    f["len_ratio"] = len(ps)/max(1,len(hs))
    f["len_diff"] = abs(len(ps)-len(hs))/1000
    for k in [2,3,4]:
        s1 = set(ps[i:i+k] for i in range(len(ps)-k+1))
        s2 = set(hs[i:i+k] for i in range(len(hs)-k+1))
        f[f"kmer{k}"] = len(s1&s2)/max(1,len(s1|s2))
    pc = Counter(ps); hc = Counter(hs)
    tl1 = max(1,len(ps)); tl2 = max(1,len(hs))
    f["aa_l2"] = sum((pc.get(a,0)/tl1 - hc.get(a,0)/tl2)**2 for a in "ACDEFGHIKLMNPQRSTVWY")
    pe = embeddings.get(p_uid)
    he = embeddings.get(h_uid)
    if pe is not None and he is not None:
        f["esm_cos"] = float(np.dot(pe,he)/(np.linalg.norm(pe)*np.linalg.norm(he)+1e-8))
        f["esm_euclid"] = float(np.linalg.norm(pe-he))
        f["esm_l1"] = float(np.sum(np.abs(pe-he)))
        f["esm_pnorm"] = float(np.linalg.norm(pe))
        f["esm_hnorm"] = float(np.linalg.norm(he))
    return f

pair_cache = {}
KEYS = None
for pl, hl, dis, kw in PAIRS:
    pu = TARGETS[pl]["uid"]; hu = TARGETS[hl]["uid"]
    pf = pair_feats(pu, hu)
    pair_cache[(kw, dis)] = pf
    if KEYS is None: KEYS = sorted(pf.keys())
    cos = pf.get("esm_cos", 0)
    print(f"  {pl:12s} - {hl:8s}: ESM cos={cos:.4f}, len_diff={pf['len_diff']:.3f}")

NP = len(KEYS)
print(f"  Pair features ({NP}): {KEYS}")

# ---- STEP 4: Build dataset ----
print("\n  STEP 4: Building feature matrix...")

smiles_lookup = {}
for path in [
    os.path.join(PROCESSED, "schisto_filtered_activities.csv"),
    os.path.join(TRYP, "tryp_activities.csv"),
    os.path.join(LEISH, "leish_activities.csv"),
]:
    if os.path.exists(path):
        for _, r in pd.read_csv(path).iterrows():
            s = r.get("canonical_smiles")
            if pd.notna(s): smiles_lookup[r["molecule_chembl_id"]] = s

all_sel = []
for path, dis in [
    (os.path.join(PROCESSED, "kira_selectivity_analysis.csv"), "Schistosomiasis"),
    (os.path.join(TRYP, "tryp_selectivity_expanded.csv"), "Trypanosomiasis"),
    (os.path.join(LEISH, "leish_selectivity.csv"), "Leishmaniasis"),
]:
    if os.path.exists(path):
        df = pd.read_csv(path); df["disease"] = dis; all_sel.append(df)
sel_df = pd.concat(all_sel, ignore_index=True)
print(f"  Selectivity: {len(sel_df)}, SMILES: {len(smiles_lookup)}")

X_list, y_list, d_list = [], [], []
matched = 0

for _, row in sel_df.iterrows():
    smi = smiles_lookup.get(row["molecule_chembl_id"])
    if smi is None: continue
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None: continue

    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=256), dtype=np.float32)
    try:
        desc = np.array([Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol), Descriptors.NumAromaticRings(mol),
            Descriptors.FractionCSP3(mol), QED.qed(mol)], dtype=np.float32)
    except: desc = np.zeros(10, dtype=np.float32)

    cpd = np.concatenate([fp, desc])

    pt = str(row.get("parasite_target","")).lower()
    dis = str(row.get("disease",""))
    mp = None
    for _,_,pd2,kw in PAIRS:
        if kw in pt and dis == pd2: mp = (kw, pd2); break
    if mp is None:
        for _,_,pd2,kw in PAIRS:
            if kw in pt: mp = (kw, pd2); break

    if mp and mp in pair_cache:
        pf = pair_cache[mp]
        pa = np.array([pf.get(k,0) for k in KEYS], dtype=np.float32)
        matched += 1
    else:
        pa = np.zeros(NP, dtype=np.float32)

    X_list.append(np.concatenate([cpd, pa]))
    y_list.append(1 if row["selectivity_class"] in ["SELECTIVE","MODERATE"] else 0)
    d_list.append(dis)

X = np.nan_to_num(np.array(X_list, dtype=np.float32))
y = np.array(y_list)
D = np.array(d_list)
CPD = 266

pnz = (X[:,CPD:] != 0).any(axis=0).sum()
print(f"\n  Matrix: {X.shape[0]} x {X.shape[1]}")
print(f"  Matched pairs: {matched}/{len(X)}")
print(f"  Non-zero pair columns: {pnz}/{NP}")
print(f"  Labels: {sum(y)} pos, {len(y)-sum(y)} neg")
print(f"  Diseases: {dict(pd.Series(D).value_counts())}")

if pnz == 0:
    print("  ERROR: pair features still zero!")
    sys.exit(1)

# ---- STEP 5: Evaluate ----
print("\n" + "=" * 60)
print("EVALUATION (CLEAN, SmDHODH FIXED)")
print("=" * 60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
baseline = max(sum(y), len(y)-sum(y)) / len(y)

# [A] 5-fold CV
gb_cpd = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1, min_samples_leaf=5, subsample=0.8, random_state=42)
gb_full = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1, min_samples_leaf=5, subsample=0.8, random_state=42)

a_cpd = cross_val_score(gb_cpd, X[:,:CPD], y, cv=cv, scoring="roc_auc")
a_full = cross_val_score(gb_full, X, y, cv=cv, scoring="roc_auc")

print(f"\n  [A] 5-FOLD CV")
print(f"  Compound only:  AUROC {a_cpd.mean():.3f} +/- {a_cpd.std():.3f}")
print(f"  Compound+Pair:  AUROC {a_full.mean():.3f} +/- {a_full.std():.3f}")
print(f"  Pair effect:    {a_full.mean()-a_cpd.mean():+.3f}")
print(f"  Baseline:       {baseline:.3f}")

# MLP
sc = StandardScaler()
Xs = sc.fit_transform(X)
mlp = MLPClassifier(hidden_layer_sizes=(128,64,32), max_iter=500, learning_rate_init=0.001,
                     early_stopping=True, validation_fraction=0.15, random_state=42)
a_mlp = cross_val_score(mlp, Xs, y, cv=cv, scoring="roc_auc")
print(f"  MLP full:       AUROC {a_mlp.mean():.3f} +/- {a_mlp.std():.3f}")

# [B] LODO
print(f"\n  [B] LEAVE-ONE-DISEASE-OUT")
for held in ["Schistosomiasis", "Trypanosomiasis", "Leishmaniasis"]:
    tr = D != held; te = D == held
    if sum(te) < 5 or len(np.unique(y[te])) < 2:
        print(f"\n  {held}: skipped"); continue

    g1 = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                      min_samples_leaf=5, subsample=0.8, random_state=42)
    g1.fit(X[tr,:CPD], y[tr])
    a1 = roc_auc_score(y[te], g1.predict_proba(X[te,:CPD])[:,1])

    g2 = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                      min_samples_leaf=5, subsample=0.8, random_state=42)
    g2.fit(X[tr], y[tr])
    a2 = roc_auc_score(y[te], g2.predict_proba(X[te])[:,1])

    print(f"\n  Held out: {held} ({sum(te)} samples, {y[te].mean():.0%} pos)")
    print(f"    Compound only:  AUROC {a1:.3f}")
    print(f"    Compound+Pair:  AUROC {a2:.3f}")
    print(f"    Pair effect:    {a2-a1:+.3f}")

# [C] Importance
print(f"\n  [C] FEATURE IMPORTANCE")
gb_full.fit(X, y)
imp = gb_full.feature_importances_
fp_i = imp[:256].sum()
desc_i = imp[256:266].sum()
pair_i = imp[266:].sum()
print(f"    Morgan FP:      {fp_i:.3f} ({100*fp_i:.1f}%)")
print(f"    RDKit desc:     {desc_i:.3f} ({100*desc_i:.1f}%)")
print(f"    Pair features:  {pair_i:.3f} ({100*pair_i:.1f}%)")
print(f"    Compound total: {fp_i+desc_i:.3f} ({100*(fp_i+desc_i):.1f}%)")
print(f"    Protein total:  {pair_i:.3f} ({100*pair_i:.1f}%)")

# Top pair features
pfi = [(KEYS[i], imp[CPD+i]) for i in range(NP)]
pfi.sort(key=lambda x: x[1], reverse=True)
print(f"\n  Pair feature ranking:")
for name, importance in pfi:
    print(f"    {name:20s} {importance:.4f}")

# Top 15 overall
names = [f"FP_{i}" for i in range(256)]
names += ["MW","LogP","HBD","HBA","TPSA","RotBonds","Rings","AromRings","Fsp3","QED"]
names += KEYS
top = np.argsort(imp)[::-1][:15]
print(f"\n  Top 15 features:")
for r, i in enumerate(top, 1):
    n = names[i] if i < len(names) else f"f_{i}"
    g = "PAIR" if i >= CPD else ("DESC" if i >= 256 else "FP")
    print(f"    {r:2d}. [{g:4s}] {n:20s} {imp[i]:.4f}")

print(f"\n{'=' * 60}")
print(f"  CLEAN RERUN COMPLETE (SmDHODH FIXED)")
print(f"{'=' * 60}")

# Save results
results = {
    "cv_compound_only": float(a_cpd.mean()),
    "cv_compound_pair": float(a_full.mean()),
    "cv_mlp": float(a_mlp.mean()),
    "baseline": float(baseline),
    "compound_importance": float(fp_i + desc_i),
    "protein_importance": float(pair_i),
}
with open(os.path.join(MDIR, "clean_model_results.json"), "w") as f:
    json.dump(results, f, indent=2)
CLEANRUN
