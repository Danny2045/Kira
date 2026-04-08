# Kira — Architecture Decisions

## Decision 1: PrimeKG as primary knowledge graph

**Date:** March 2026
**Status:** Active

**Context:** We need a biomedical knowledge graph that links diseases, genes/proteins,
drugs, pathways, and biological processes. Options include PrimeKG, Hetionet,
DRKG, and building our own from source databases.

**Decision:** Use PrimeKG as the primary knowledge graph for Phase 1.

**Rationale:**
- PrimeKG integrates 20+ biomedical databases into a single unified graph
- It contains disease-gene, drug-gene, drug-disease, and gene-gene edges
- It is publicly available and well-documented
- It is large enough to provide meaningful connectivity for neglected diseases
- Harvard's Zitnik Lab maintains it with academic rigor

**Consequences:**
- We inherit PrimeKG's biases and coverage gaps
- Schistosomiasis-specific coverage may be thin (NTDs are under-represented)
- We will likely need to augment with targeted data from ChEMBL, UniProt,
  and PDB as the pipeline matures

## Decision 2: Python + pandas + NetworkX for Phase 1

**Date:** March 2026
**Status:** Active

**Context:** The pipeline needs to load, query, and traverse a knowledge graph.
Options range from Python with pandas/NetworkX to Neo4j to JAX-based
graph operations.

**Decision:** Use pandas for data loading/filtering and NetworkX for graph
traversal in Phase 1. Migrate performance-critical paths to JAX in Phase 2.

**Rationale:**
- Simplicity. The graph fits in memory on our 48GB M4 Pro.
- NetworkX is well-documented and sufficient for the queries we need now.
- Premature optimization into JAX/GPU before the pipeline logic is stable
  would slow us down.
- We can always profile and migrate later.

**Consequences:**
- Some operations will be slower than necessary (acceptable for Phase 1)
- When we scale to batched docking or large-scale scoring, we will need
  to migrate those specific bottlenecks to JAX

## Decision 3: Evaluation-first pipeline design

**Date:** March 2026
**Status:** Active

**Context:** We could build the pipeline feature-first (add capabilities, evaluate
later) or evaluation-first (define what "correct" means before building features).

**Decision:** Evaluation-first. The retrospective evaluation set is built before
the ranking algorithm is refined.

**Rationale:**
- Without ground truth, we cannot distinguish a working pipeline from a
  plausible-looking demo.
- Known schistosomiasis drug repurposing outcomes exist in literature.
- Building the eval set first forces us to understand the domain deeply.

**Consequences:**
- Slower start (we invest in evaluation infrastructure before features)
- Much higher confidence in results when we do produce them
- Publication-quality methodology from the beginning
