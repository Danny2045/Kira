"""
Kira — Script 01: Download and Explore PrimeKG
================================================

This is the first real piece of Kira. It does three things:

1. Downloads PrimeKG (a biomedical knowledge graph) if not already present
2. Loads it into memory and shows you what's inside
3. Finds everything connected to schistosomiasis

WHAT IS PrimeKG?
PrimeKG is a knowledge graph built by Harvard's Zitnik Lab. A knowledge graph
is a dataset of relationships between biomedical entities. Each row says:

    "Entity A has relationship R with Entity B"

For example:
    "Schistosomiasis -- associates_with --> SmTGR protein"
    "Praziquantel   -- treats           --> Schistosomiasis"
    "SmTGR          -- interacts_with   --> some other protein"

By traversing these relationships, we can find:
- Which proteins are associated with schistosomiasis (drug targets)
- Which drugs are known to interact with those targets
- Which biological pathways are involved
- What other diseases share mechanisms (potential cross-disease insights)

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 01_explore_primekg.py

This will take a few minutes the first time (downloading ~250MB).
"""

import os
import sys
import time
import pandas as pd
import networkx as nx

# ---------------------------------------------------------------------------
# CONFIGURATION — things you might change
# ---------------------------------------------------------------------------

# Where to store the downloaded data
DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "raw")

# PrimeKG download URL (Harvard Dataverse)
PRIMEKG_URL = "https://dataverse.harvard.edu/api/access/datafile/6180620"
PRIMEKG_FILE = os.path.join(DATA_DIR, "primekg.csv")

# Our target disease
DISEASE_NAME = "schistosomiasis"

# ---------------------------------------------------------------------------
# STEP 1: Download PrimeKG
# ---------------------------------------------------------------------------

def download_primekg():
    """Download PrimeKG if we don't already have it."""
    
    if os.path.exists(PRIMEKG_FILE):
        size_mb = os.path.getsize(PRIMEKG_FILE) / (1024 * 1024)
        print(f"PrimeKG already downloaded ({size_mb:.0f} MB)")
        return
    
    print("Downloading PrimeKG from Harvard Dataverse...")
    print("This is ~250 MB and may take a few minutes.\n")
    
    import requests
    from tqdm import tqdm
    
    os.makedirs(DATA_DIR, exist_ok=True)
    
    # Stream download with progress bar
    response = requests.get(PRIMEKG_URL, stream=True)
    response.raise_for_status()
    
    total_size = int(response.headers.get("content-length", 0))
    
    with open(PRIMEKG_FILE, "wb") as f:
        with tqdm(total=total_size, unit="B", unit_scale=True, desc="Downloading") as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))
    
    size_mb = os.path.getsize(PRIMEKG_FILE) / (1024 * 1024)
    print(f"\nDownload complete: {size_mb:.0f} MB saved to {PRIMEKG_FILE}")


# ---------------------------------------------------------------------------
# STEP 2: Load and understand PrimeKG
# ---------------------------------------------------------------------------

def load_primekg():
    """Load PrimeKG into a pandas DataFrame and explain what's inside."""
    
    print("\n" + "=" * 60)
    print("LOADING PrimeKG")
    print("=" * 60)
    
    start = time.time()
    df = pd.read_csv(PRIMEKG_FILE, low_memory=False)
    elapsed = time.time() - start
    
    print(f"\nLoaded {len(df):,} edges in {elapsed:.1f} seconds")
    print(f"Columns: {list(df.columns)}")
    
    # --- Understand the structure ---
    print("\n" + "-" * 40)
    print("WHAT'S IN THE GRAPH?")
    print("-" * 40)
    
    # Show relationship types
    print(f"\nRelationship types ({df['relation'].nunique()} unique):")
    for rel, count in df["relation"].value_counts().head(15).items():
        print(f"  {rel}: {count:,} edges")
    
    # Show node types
    # PrimeKG has x_type (source node type) and y_type (target node type)
    print(f"\nSource node types:")
    for ntype, count in df["x_type"].value_counts().items():
        print(f"  {ntype}: {count:,}")
    
    print(f"\nTarget node types:")
    for ntype, count in df["y_type"].value_counts().items():
        print(f"  {ntype}: {count:,}")
    
    return df


# ---------------------------------------------------------------------------
# STEP 3: Find schistosomiasis and its neighborhood
# ---------------------------------------------------------------------------

def explore_disease(df, disease_name):
    """Find everything connected to our target disease in the knowledge graph."""
    
    print("\n" + "=" * 60)
    print(f"EXPLORING: {disease_name.upper()}")
    print("=" * 60)
    
    # --- Find disease nodes ---
    # Search both x_name and y_name columns (disease can appear on either side)
    mask_x = df["x_name"].str.contains(disease_name, case=False, na=False)
    mask_y = df["y_name"].str.contains(disease_name, case=False, na=False)
    
    disease_edges = df[mask_x | mask_y].copy()
    
    if len(disease_edges) == 0:
        print(f"\nNo edges found containing '{disease_name}'.")
        print("Try checking the exact disease name in PrimeKG.")
        return None
    
    print(f"\nFound {len(disease_edges):,} edges involving '{disease_name}'")
    
    # --- Show exact disease node names ---
    disease_names_x = df.loc[mask_x, "x_name"].unique()
    disease_names_y = df.loc[mask_y, "y_name"].unique()
    all_disease_names = set(list(disease_names_x) + list(disease_names_y))
    
    print(f"\nExact disease node names found:")
    for name in sorted(all_disease_names):
        print(f"  - {name}")
    
    # --- Categorize connections by relationship type ---
    print(f"\nRelationship breakdown:")
    for rel, count in disease_edges["relation"].value_counts().items():
        print(f"  {rel}: {count}")
    
    # --- Find connected genes/proteins (potential drug targets) ---
    print("\n" + "-" * 40)
    print("CONNECTED GENES/PROTEINS (potential drug targets)")
    print("-" * 40)
    
    # Edges where disease is source, target is gene/protein
    gene_edges_1 = disease_edges[
        (mask_x) & (df["y_type"].isin(["gene/protein"]))
    ][["x_name", "relation", "y_name", "y_type"]]
    
    # Edges where disease is target, source is gene/protein
    gene_edges_2 = disease_edges[
        (mask_y) & (df["x_type"].isin(["gene/protein"]))
    ][["x_name", "x_type", "relation", "y_name"]]
    
    # Collect all connected genes
    genes_from_1 = set(gene_edges_1["y_name"].tolist()) if len(gene_edges_1) > 0 else set()
    genes_from_2 = set(gene_edges_2["x_name"].tolist()) if len(gene_edges_2) > 0 else set()
    all_genes = genes_from_1 | genes_from_2
    
    if all_genes:
        print(f"\nFound {len(all_genes)} connected genes/proteins:")
        for gene in sorted(all_genes)[:30]:  # Show first 30
            print(f"  - {gene}")
        if len(all_genes) > 30:
            print(f"  ... and {len(all_genes) - 30} more")
    else:
        print("\nNo direct gene/protein connections found.")
    
    # --- Find connected drugs ---
    print("\n" + "-" * 40)
    print("CONNECTED DRUGS")
    print("-" * 40)
    
    drug_edges_1 = disease_edges[
        (mask_x) & (df["y_type"] == "drug")
    ][["x_name", "relation", "y_name"]]
    
    drug_edges_2 = disease_edges[
        (mask_y) & (df["x_type"] == "drug")
    ][["x_name", "relation", "y_name"]]
    
    drugs_from_1 = set(drug_edges_1["y_name"].tolist()) if len(drug_edges_1) > 0 else set()
    drugs_from_2 = set(drug_edges_2["x_name"].tolist()) if len(drug_edges_2) > 0 else set()
    all_drugs = drugs_from_1 | drugs_from_2
    
    if all_drugs:
        print(f"\nFound {len(all_drugs)} connected drugs:")
        for drug in sorted(all_drugs):
            print(f"  - {drug}")
    else:
        print("\nNo direct drug connections found.")
    
    # --- Find connected pathways ---
    print("\n" + "-" * 40)
    print("CONNECTED BIOLOGICAL PATHWAYS")
    print("-" * 40)
    
    pathway_edges_1 = disease_edges[
        (mask_x) & (df["y_type"] == "biological_process")
    ][["x_name", "relation", "y_name"]]
    
    pathway_edges_2 = disease_edges[
        (mask_y) & (df["x_type"] == "biological_process")
    ][["x_name", "relation", "y_name"]]
    
    paths_from_1 = set(pathway_edges_1["y_name"].tolist()) if len(pathway_edges_1) > 0 else set()
    paths_from_2 = set(pathway_edges_2["x_name"].tolist()) if len(pathway_edges_2) > 0 else set()
    all_pathways = paths_from_1 | paths_from_2
    
    if all_pathways:
        print(f"\nFound {len(all_pathways)} connected biological processes:")
        for pathway in sorted(all_pathways)[:20]:
            print(f"  - {pathway}")
        if len(all_pathways) > 20:
            print(f"  ... and {len(all_pathways) - 20} more")
    else:
        print("\nNo direct biological process connections found.")
    
    # --- Summary ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Disease:             {disease_name}")
    print(f"  Total edges:         {len(disease_edges):,}")
    print(f"  Genes/proteins:      {len(all_genes)}")
    print(f"  Drugs:               {len(all_drugs)}")
    print(f"  Biological processes: {len(all_pathways)}")
    print(f"\n  These genes/proteins are your candidate drug targets.")
    print(f"  These drugs are your starting candidates for repurposing.")
    print(f"  Next step: investigate which targets have 3D structures")
    print(f"  available for molecular docking.")
    
    # --- Save results for next script ---
    results_dir = os.path.join(os.path.dirname(__file__), "data", "processed")
    os.makedirs(results_dir, exist_ok=True)
    
    disease_edges.to_csv(
        os.path.join(results_dir, f"{disease_name}_subgraph.csv"),
        index=False
    )
    print(f"\n  Subgraph saved to data/processed/{disease_name}_subgraph.csv")
    
    return {
        "edges": disease_edges,
        "genes": all_genes,
        "drugs": all_drugs,
        "pathways": all_pathways,
    }


# ---------------------------------------------------------------------------
# STEP 4: Two-hop exploration — drugs that target the same proteins
# ---------------------------------------------------------------------------

def find_repurposing_candidates(df, disease_genes):
    """
    This is where repurposing logic begins.
    
    We have: genes/proteins connected to schistosomiasis.
    We want: drugs that interact with those same genes/proteins,
             even if those drugs aren't directly linked to schistosomiasis.
    
    This is a 2-hop query:
        schistosomiasis → gene/protein → drug
    
    The drugs we find in hop 2 are repurposing candidates.
    """
    
    if not disease_genes:
        print("No disease genes to search from.")
        return None
    
    print("\n" + "=" * 60)
    print("TWO-HOP SEARCH: Finding drugs that target disease-linked proteins")
    print("=" * 60)
    
    # Find all drug-gene/protein edges in the entire graph
    drug_gene_mask = (
        ((df["x_type"] == "drug") & (df["y_type"] == "gene/protein")) |
        ((df["x_type"] == "gene/protein") & (df["y_type"] == "drug"))
    )
    drug_gene_edges = df[drug_gene_mask]
    
    print(f"\nTotal drug-protein edges in PrimeKG: {len(drug_gene_edges):,}")
    
    # Filter to edges involving our disease-linked genes
    relevant_mask = (
        drug_gene_edges["x_name"].isin(disease_genes) |
        drug_gene_edges["y_name"].isin(disease_genes)
    )
    relevant_edges = drug_gene_edges[relevant_mask]
    
    print(f"Drug-protein edges involving our {len(disease_genes)} disease genes: {len(relevant_edges):,}")
    
    # Extract the drugs
    drugs_x = set(relevant_edges.loc[relevant_edges["x_type"] == "drug", "x_name"])
    drugs_y = set(relevant_edges.loc[relevant_edges["y_type"] == "drug", "y_name"])
    candidate_drugs = drugs_x | drugs_y
    
    print(f"\nRepurposing candidates found: {len(candidate_drugs)}")
    
    if candidate_drugs:
        # Show which target each drug hits
        print(f"\nTop candidates (showing target connections):\n")
        
        drug_target_pairs = []
        
        for _, row in relevant_edges.iterrows():
            if row["x_type"] == "drug":
                drug_target_pairs.append({
                    "drug": row["x_name"],
                    "relation": row["relation"],
                    "target": row["y_name"],
                })
            else:
                drug_target_pairs.append({
                    "drug": row["y_name"],
                    "relation": row["relation"],
                    "target": row["x_name"],
                })
        
        pairs_df = pd.DataFrame(drug_target_pairs)
        
        # Count how many disease targets each drug hits
        drug_target_counts = pairs_df.groupby("drug")["target"].nunique().sort_values(ascending=False)
        
        print("Drugs ranked by number of disease-linked targets hit:")
        print("-" * 50)
        for drug, count in drug_target_counts.head(25).items():
            targets = pairs_df[pairs_df["drug"] == drug]["target"].unique()
            targets_str = ", ".join(sorted(targets)[:3])
            if len(targets) > 3:
                targets_str += f", +{len(targets)-3} more"
            print(f"  {drug}")
            print(f"    Targets ({count}): {targets_str}")
            print()
        
        # Save candidate list
        results_dir = os.path.join(os.path.dirname(__file__), "data", "processed")
        pairs_df.to_csv(
            os.path.join(results_dir, "repurposing_candidates_raw.csv"),
            index=False
        )
        drug_target_counts.to_csv(
            os.path.join(results_dir, "candidate_target_counts.csv")
        )
        print(f"\nResults saved to data/processed/")
        
        return pairs_df
    
    return None


# ---------------------------------------------------------------------------
# MAIN — Run everything
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    
    print("=" * 60)
    print("  KIRA — Phase 1, Script 01")
    print("  PrimeKG Exploration for Schistosomiasis")
    print("=" * 60)
    print()
    
    # Step 1: Get the data
    download_primekg()
    
    # Step 2: Load and understand it
    df = load_primekg()
    
    # Step 3: Find our disease neighborhood
    results = explore_disease(df, DISEASE_NAME)
    
    # Step 4: Two-hop search for repurposing candidates
    if results and results["genes"]:
        candidates = find_repurposing_candidates(df, results["genes"])
    
    print("\n" + "=" * 60)
    print("  DONE")
    print("  Next: 02_evaluate_targets.py — assess which targets")
    print("  have 3D structures available for molecular docking")
    print("=" * 60)
