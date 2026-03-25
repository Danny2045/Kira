#!/bin/bash
cd ~/kira

# Delete cached sequences and embeddings (they used wrong UniProt IDs)
rm -f data/models/sequences.json data/models/esm2_embeddings.npz

python3 << 'PATCH'
import json

# Fix the UniProt IDs in the script
with open("19_selectivity_prediction.py", "r") as f:
    code = f.read()

# The wrong IDs and their corrections:
# SmHDAC8: A0A094ZS53 (wrong, 101aa fragment) → C1LV40 (correct SmHDAC8, 340aa)
# SmDHODH: G4LZI2 (not found) → C1L5Z2 (S. mansoni DHODH)
# TbCatB: Q7YWB5 (tick vWF!) → Q8WQ44 (T. brucei cathepsin B-like)
# LmDHFR-TS: Q01781 (parsley!) → P07382 (L. major DHFR-TS bifunctional)

replacements = [
    ('"A0A094ZS53"', '"C1LV40"'),       # SmHDAC8
    ('"G4LZI2"', '"C1L5Z2"'),           # SmDHODH  
    ('"Q7YWB5"', '"Q8WQ44"'),           # T. brucei cathepsin B
    ('"Q01781"', '"P07382"'),            # L. major DHFR-TS
]

for old, new in replacements:
    if old in code:
        code = code.replace(old, new)
        print(f"  Fixed: {old} → {new}")
    else:
        print(f"  Not found: {old}")

# Fix the classifier SMILES lookup:
# The selectivity CSVs don't have SMILES. We need to join from activities.

old_smiles = '''        smiles = row.get("smiles") or row.get("canonical_smiles")
        if pd.isna(smiles):
            continue'''

new_smiles = '''        smiles = row.get("smiles") or row.get("canonical_smiles")
        if pd.isna(smiles) or smiles is None:
            # Look up from activity data
            cid = row["molecule_chembl_id"]
            smiles_lookup = all_smiles.get(cid)
            if smiles_lookup is None:
                continue
            smiles = smiles_lookup'''

if old_smiles in code:
    code = code.replace(old_smiles, new_smiles)
    print("  Fixed: SMILES lookup fallback added")

# Add SMILES loading before the classifier loop
old_features = '''    features = []
    labels = []
    compounds_used = []'''

new_features = '''    # Build SMILES lookup from all activity files
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

    features = []
    labels = []
    compounds_used = []'''

if old_features in code:
    code = code.replace(old_features, new_features)
    print("  Fixed: SMILES lookup from activity files added")

with open("19_selectivity_prediction.py", "w") as f:
    f.write(code)

print("\nAll patches applied. Re-run: python 19_selectivity_prediction.py")
PATCH
