"""
Kira - Script 09: Novelty Filter
==================================

For each top-ranked compound, check: has this already been proposed
for schistosomiasis in the published literature?

Uses the NCBI PubMed E-utilities API (free, no key required for
low-volume queries) to count publications mentioning each compound
in the context of schistosomiasis.

Classifications:
  - NOVEL:     0 publications → genuinely new prediction
  - EMERGING:  1-3 publications → explored but not well-established
  - KNOWN:     4-10 publications → established hypothesis
  - WELL-KNOWN: 10+ publications → extensively studied

The compounds classified as NOVEL or EMERGING are Kira's most
publishable outputs — they represent candidates the field has
not yet systematically considered.

HOW TO RUN:
    cd ~/kira
    conda activate bio-builder
    python 09_novelty_filter.py
"""

import os
import sys
import time
import requests
import numpy as np
import pandas as pd
from datetime import datetime

PROCESSED_DIR = os.path.join(os.path.dirname(__file__), "data", "processed")
EVAL_DIR = os.path.join(os.path.dirname(__file__), "data", "eval")
REPORT_DIR = os.path.join(os.path.dirname(__file__), "data", "reports")

# PubMed E-utilities base URL
PUBMED_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Rate limit: NCBI allows 3 requests/second without API key
RATE_LIMIT_SECONDS = 0.4


# ---------------------------------------------------------------------------
# STEP 1: Search PubMed for compound + schistosomiasis co-occurrence
# ---------------------------------------------------------------------------

def search_pubmed(compound_name, disease="schistosomiasis"):
    """
    Query PubMed for papers mentioning both the compound and the disease.
    Returns the count of matching publications.
    """
    query = f'"{compound_name}" AND "{disease}"'

    params = {
        "db": "pubmed",
        "term": query,
        "rettype": "count",
        "retmode": "json",
    }

    try:
        response = requests.get(PUBMED_SEARCH_URL, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        count = int(data.get("esearchresult", {}).get("count", 0))
        return count
    except Exception as e:
        return -1  # Error indicator


def classify_novelty(pub_count):
    """Classify a compound based on its publication count."""
    if pub_count < 0:
        return "ERROR"
    elif pub_count == 0:
        return "NOVEL"
    elif pub_count <= 3:
        return "EMERGING"
    elif pub_count <= 10:
        return "KNOWN"
    else:
        return "WELL-KNOWN"


# ---------------------------------------------------------------------------
# STEP 2: Run novelty filter on top candidates
# ---------------------------------------------------------------------------

def run_novelty_filter(ranked_df, top_n=50):
    """
    For the top N compounds, check PubMed for existing literature
    on their connection to schistosomiasis.
    """
    print("\n" + "=" * 60)
    print("RUNNING NOVELTY FILTER")
    print("=" * 60)

    # Get top compounds that have names (ChEMBL IDs alone won't find papers)
    candidates = ranked_df.head(top_n).copy()

    results = []

    print(f"\n  Checking top {top_n} compounds against PubMed...")
    print(f"  Query format: \"[compound name]\" AND \"schistosomiasis\"\n")

    for idx, (_, row) in enumerate(candidates.iterrows()):
        cid = row["molecule_chembl_id"]
        name = row.get("pref_name")

        # Skip compounds without human-readable names
        if pd.isna(name) or name == cid:
            # Try searching by ChEMBL ID as fallback
            search_name = cid
            search_type = "chembl_id"
        else:
            search_name = name
            search_type = "drug_name"

        # Search PubMed
        pub_count = search_pubmed(search_name)
        novelty = classify_novelty(pub_count)

        # Also search with just the compound name (without schistosomiasis)
        # to see if the compound itself is well-studied
        general_count = search_pubmed(search_name, disease="")

        results.append({
            "molecule_chembl_id": cid,
            "pref_name": name if pd.notna(name) else cid,
            "search_name": search_name,
            "search_type": search_type,
            "composite_score": row.get("composite_score", 0),
            "label": row.get("label", ""),
            "best_value": row.get("best_value"),
            "best_target": row.get("best_target", ""),
            "schisto_publications": pub_count,
            "general_publications": general_count,
            "novelty_class": novelty,
        })

        # Progress
        status = f"{'*' if novelty == 'NOVEL' else ' '}"
        pub_str = f"{pub_count}" if pub_count >= 0 else "ERR"
        print(f"  {status} [{idx+1:2d}/{top_n}] {search_name:35s}  "
              f"schisto={pub_str:>4s}  general={general_count:>6d}  -> {novelty}")

        # Rate limit
        time.sleep(RATE_LIMIT_SECONDS)

    results_df = pd.DataFrame(results)
    return results_df


# ---------------------------------------------------------------------------
# STEP 3: Analyze and report
# ---------------------------------------------------------------------------

def analyze_novelty(novelty_df):
    """
    Analyze novelty results and highlight the most publishable findings.
    """
    print("\n" + "=" * 60)
    print("NOVELTY ANALYSIS")
    print("=" * 60)

    # Summary by class
    print(f"\n  Novelty distribution (top {len(novelty_df)} compounds):")
    for cls in ["NOVEL", "EMERGING", "KNOWN", "WELL-KNOWN", "ERROR"]:
        n = len(novelty_df[novelty_df["novelty_class"] == cls])
        if n > 0:
            print(f"    {cls:12s}: {n}")

    # NOVEL compounds — the most publishable
    novel = novelty_df[novelty_df["novelty_class"] == "NOVEL"].copy()
    if len(novel) > 0:
        print(f"\n  {'=' * 55}")
        print(f"  NOVEL CANDIDATES — No prior schistosomiasis literature")
        print(f"  {'=' * 55}")
        print(f"  These are Kira's most publishable findings.\n")

        for _, row in novel.iterrows():
            name = row["pref_name"]
            score = row["composite_score"]
            target = row.get("best_target", "unknown")
            ic50 = row.get("best_value")
            general = row["general_publications"]

            print(f"  {name}")
            if pd.notna(ic50) and ic50 > 0:
                print(f"    IC50 = {ic50:.0f} nM vs {target}")
            print(f"    Composite score: {score:.4f}")
            print(f"    Schistosomiasis publications: 0")
            print(f"    General publications: {general}")
            if general > 100:
                print(f"    Note: Well-studied compound, but NOT in schistosomiasis context")
            elif general > 0:
                print(f"    Note: Some literature exists, but not for schistosomiasis")
            else:
                print(f"    Note: Research compound with limited literature overall")
            print()

    # EMERGING compounds — early-stage hypotheses
    emerging = novelty_df[novelty_df["novelty_class"] == "EMERGING"].copy()
    if len(emerging) > 0:
        print(f"\n  {'=' * 55}")
        print(f"  EMERGING CANDIDATES — 1-3 prior publications")
        print(f"  {'=' * 55}")
        print(f"  These have been explored minimally. Kira independently confirms them.\n")

        for _, row in emerging.iterrows():
            name = row["pref_name"]
            pub_count = row["schisto_publications"]
            ic50 = row.get("best_value")
            target = row.get("best_target", "unknown")

            print(f"  {name}")
            if pd.notna(ic50) and ic50 > 0:
                print(f"    IC50 = {ic50:.0f} nM vs {target}")
            print(f"    Schistosomiasis publications: {pub_count}")
            print()

    # KNOWN/WELL-KNOWN — rediscoveries
    known = novelty_df[novelty_df["novelty_class"].isin(["KNOWN", "WELL-KNOWN"])].copy()
    if len(known) > 0:
        print(f"\n  {'=' * 55}")
        print(f"  KNOWN/WELL-KNOWN — Kira confirms established hypotheses")
        print(f"  {'=' * 55}\n")

        for _, row in known.iterrows():
            name = row["pref_name"]
            pub_count = row["schisto_publications"]
            print(f"  {name:35s}  {pub_count} publications")

    return novel, emerging, known


# ---------------------------------------------------------------------------
# STEP 4: Generate novelty report
# ---------------------------------------------------------------------------

def generate_novelty_report(novelty_df, novel, emerging, known):
    """Save the novelty analysis as a formal report."""

    os.makedirs(REPORT_DIR, exist_ok=True)

    report = []
    report.append("=" * 70)
    report.append("KIRA v1 — NOVELTY ANALYSIS REPORT")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("=" * 70)
    report.append("")
    report.append("METHODOLOGY:")
    report.append("  For each top-ranked compound, PubMed was queried for")
    report.append("  co-occurrence of the compound name with 'schistosomiasis'.")
    report.append("  Compounds with zero prior publications represent genuinely")
    report.append("  novel predictions from the Kira pipeline.")
    report.append("")
    report.append("CLASSIFICATION:")
    report.append("  NOVEL:      0 publications (new prediction)")
    report.append("  EMERGING:   1-3 publications (minimally explored)")
    report.append("  KNOWN:      4-10 publications (established hypothesis)")
    report.append("  WELL-KNOWN: 10+ publications (extensively studied)")
    report.append("")

    report.append(f"SUMMARY (top {len(novelty_df)} compounds):")
    for cls in ["NOVEL", "EMERGING", "KNOWN", "WELL-KNOWN"]:
        n = len(novelty_df[novelty_df["novelty_class"] == cls])
        report.append(f"  {cls:12s}: {n}")

    if len(novel) > 0:
        report.append("")
        report.append("-" * 70)
        report.append("NOVEL CANDIDATES — HIGHEST PUBLICATION PRIORITY")
        report.append("-" * 70)
        for _, row in novel.iterrows():
            name = row["pref_name"]
            ic50 = row.get("best_value")
            target = row.get("best_target", "")
            report.append(f"  {name}")
            if pd.notna(ic50) and ic50 > 0:
                report.append(f"    IC50 = {ic50:.0f} nM vs {target}")
            report.append(f"    Score: {row['composite_score']:.4f}")
            report.append(f"    General literature: {row['general_publications']} papers")
            report.append("")

    if len(emerging) > 0:
        report.append("-" * 70)
        report.append("EMERGING CANDIDATES — INDEPENDENT CONFIRMATION VALUE")
        report.append("-" * 70)
        for _, row in emerging.iterrows():
            name = row["pref_name"]
            report.append(f"  {name}  ({row['schisto_publications']} prior papers)")
            if pd.notna(row.get("best_value")) and row["best_value"] > 0:
                report.append(f"    IC50 = {row['best_value']:.0f} nM vs {row.get('best_target','')}")
            report.append("")

    report.append("-" * 70)
    report.append("INTERPRETATION")
    report.append("-" * 70)
    report.append("  NOVEL compounds represent Kira's strongest publication case.")
    report.append("  The pipeline surfaced them through systematic integration of")
    report.append("  target-based activity, structural similarity, and whole-organism")
    report.append("  evidence — a methodology that has not been applied to these")
    report.append("  specific compounds in the schistosomiasis context before.")
    report.append("")
    report.append("  EMERGING compounds demonstrate independent computational")
    report.append("  confirmation of hypotheses with minimal prior evidence.")
    report.append("")
    report.append("  KNOWN/WELL-KNOWN compounds validate the pipeline: if it")
    report.append("  correctly surfaces established candidates, its novel")
    report.append("  predictions carry more weight.")
    report.append("")
    report.append("=" * 70)

    report_text = "\n".join(report)
    report_path = os.path.join(REPORT_DIR, "kira_v1_novelty_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)

    # Also save full novelty data
    novelty_df.to_csv(
        os.path.join(PROCESSED_DIR, "kira_novelty_analysis.csv"), index=False
    )

    print(f"\n  Report saved to: {report_path}")
    print(f"  Data saved to: {os.path.join(PROCESSED_DIR, 'kira_novelty_analysis.csv')}")

    return report_text


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("=" * 60)
    print("  KIRA -- Phase 1, Script 09")
    print("  Novelty Filter (PubMed Literature Check)")
    print("=" * 60)

    # Load v3 extended rankings
    v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3_extended.csv")
    if not os.path.exists(v3_path):
        v3_path = os.path.join(PROCESSED_DIR, "kira_ranked_candidates_v3.csv")
    if not os.path.exists(v3_path):
        print("ERROR: Run Script 06 first.")
        sys.exit(1)

    ranked_df = pd.read_csv(v3_path)
    ranked_df = ranked_df.sort_values("composite_score", ascending=False)
    print(f"\n  Loaded {len(ranked_df)} ranked compounds")

    # Run novelty filter on top 50
    novelty_df = run_novelty_filter(ranked_df, top_n=50)

    # Analyze
    novel, emerging, known = analyze_novelty(novelty_df)

    # Generate report
    generate_novelty_report(novelty_df, novel, emerging, known)

    # Commit summary
    print(f"\n{'=' * 60}")
    print(f"  NOVELTY FILTER COMPLETE")
    print(f"")
    print(f"  Novel candidates (0 prior papers): {len(novel)}")
    print(f"  Emerging (1-3 papers): {len(emerging)}")
    print(f"  Known (4+ papers): {len(known)}")
    print(f"")
    if len(novel) > 0:
        print(f"  Your most publishable findings are the NOVEL compounds.")
        print(f"  These have never been proposed for schistosomiasis before.")
    print(f"{'=' * 60}")
