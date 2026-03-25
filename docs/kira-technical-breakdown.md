# Kira: Technical Breakdown
## How a drug repurposing pipeline works, from data layer to finding

*Written to explain the system to someone who is technically brilliant but not necessarily a domain expert in computational pharmacology.*

---

## The problem in one paragraph

There is a parasitic worm called Schistosoma mansoni that lives inside the blood vessels of about 200 million people in sub-Saharan Africa. The world treats it with one drug (praziquantel) and has for 40 years. If resistance emerges, there is no backup. New drug development costs $1B+ and takes 15 years, and the market won't fund it for a disease of poverty. Drug repurposing — finding that an existing approved drug also kills this worm — is the fastest viable path. The question is: how do you systematically search the space of approved drugs × parasite vulnerabilities and identify candidates that are potent, selective, safe, and deployable in the actual clinical context?

---

## The architecture

The pipeline is 14 sequential scripts (13 Phase 1 + 1 Phase 2). No ML models are trained. No GPU is used. The intelligence is in the data integration logic and the filtering cascade. Here is the system diagram:

```
INPUT LAYER
  PrimeKG (8M-edge biomedical knowledge graph, Harvard)
  ChEMBL (experimental bioactivity database, EBI)
  PubMed (36M+ publication index, NIH)
  WHO Essential Medicines List (manually curated)
  RCSB PDB (protein crystal structures: 2X99, 2ZZC)

EVIDENCE PATHWAYS (3)
  1. Target-based IC50     → "this compound inhibits this parasite protein at X nM"
  2. Structural similarity → "this compound looks like known active compounds"
  3. Whole-organism assay   → "this compound kills the worm at X nM"

SCORING (7 signals → weighted sum → 0-1 composite)
  potency | target_essentiality | confidence | drug_stage |
  multi_target | similarity | whole_organism

FILTERING CASCADE
  ADMET (Lipinski, TPSA, QED) → drug-likeness gate
  Selectivity (human orthologue IC50 / parasite IC50) → safety gate
  Novelty (PubMed co-occurrence count) → literature check
  Supply chain (WHO EML, African availability, cost, route) → deployment gate

STRUCTURAL VALIDATION (Phase 2)
  Molecular docking (AutoDock Vina) → SmTGR vs human TrxR1
  Binding energy comparison → computational selectivity estimate
  Correlation with experimental IC50 → docking setup validation

OUTPUT
  Two-tier ranked shortlist:
    Translational tier: approved drugs surviving all filters
    Discovery tier: selective research compounds for wet-lab follow-up
```

---

## The data model

Every compound in the system is characterized by a vector of heterogeneous features, most of which are sparse. For a given compound C:

```
C = {
  chembl_id:          string        // unique identifier
  smiles:             string        // molecular graph as text
  fingerprint:        bit[2048]     // Morgan/ECFP2 encoding of substructures
  
  // Per-target activity (sparse — most compounds tested against 0-2 targets)
  activities: [{
    target:           string        // e.g. "SmHDAC8"
    ic50_nM:          float         // inhibitory concentration, lower = better
    n_measurements:   int           // replication count
    assay_type:       string        // "biochemical", "cell-based", etc.
  }]
  
  // Whole-organism (very sparse — 6 out of 194 compounds)
  whole_org_ic50_nM:  float | null
  
  // Structural similarity to reference set
  max_tanimoto:       float [0,1]   // max over all known potent actives
  nearest_active:     string        // chembl_id of most similar active
  
  // Drug-likeness (computed from SMILES via RDKit)
  MW:                 float         // molecular weight, daltons
  LogP:               float         // octanol-water partition coefficient
  HBD:                int           // hydrogen bond donors
  HBA:                int           // hydrogen bond acceptors
  TPSA:               float         // topological polar surface area, Å²
  QED:                float [0,1]   // quantitative drug-likeness estimate
  lipinski_violations: int [0-4]
  
  // Selectivity (sparse — 91 out of 176 compounds)
  selectivity: [{
    parasite_target:  string
    human_orthologue: string
    parasite_ic50:    float
    human_ic50:       float
    ratio:            float         // human/parasite; >1 means parasite-selective
    class:            enum          // SELECTIVE|MODERATE|POOR|COUNTER-SELECTIVE
  }]
  
  // Drug development status
  max_phase:          float [0-4]   // 4 = approved, 0 = research
  first_approval:     int | null    // year
  
  // Supply chain (only for named drugs, manually curated)
  who_eml:            bool
  africa_availability: enum
  cost_tier:          enum
  route:              string        // "oral", "iv", etc.
  
  // Literature
  schisto_pub_count:  int           // PubMed hits for "compound AND schistosomiasis"
}
```

The fundamental challenge: this feature vector is extremely sparse. Most compounds have data for 1-2 targets, zero whole-organism measurements, zero selectivity data, and no supply chain information. The scoring function must handle missing data gracefully and not reward compounds simply for having data.

---

## The scoring function

The composite score is not learned. It is a hand-designed weighted sum with explicit handling of missing data. Design rationale: with 176 compounds and 7 signals (most sparse), there is not enough data to learn weights without overfitting. Domain knowledge is the better prior.

```python
composite = (
    0.25 * potency_score        # pIC50 → [0,1] via linear scaling
  + 0.15 * target_essentiality  # domain-knowledge weight per target
  + 0.05 * confidence_score     # log2(n_measurements) → [0,1], saturating
  + 0.10 * drug_stage_score     # max_phase → {0.1, 0.3, 0.5, 0.7, 1.0}
  + 0.05 * multitarget_score    # n_distinct_targets → {0.3, 0.7, 1.0}
  + 0.15 * similarity_score     # max Tanimoto to reference set → [0,1]
  + 0.25 * whole_org_score      # pIC50 of whole-organism assay → [0,1]
)

# Post-hoc multipliers (not part of the weighted sum):
final = composite * selectivity_mult * admet_mult

# Selectivity multiplier:
#   COUNTER-SELECTIVE: 0.10  (near-elimination)
#   POOR:              0.50
#   MODERATE:          0.85
#   SELECTIVE:         1.10  (bonus)
#   UNKNOWN (SmTGR):   0.80  (favorable biology, no data)
#   UNKNOWN (other):   0.70

# ADMET multiplier:
#   0 Lipinski violations: 1.00
#   1 violation:           0.85
#   2 violations:          0.60
#   3+ violations:         0.35
#   TPSA > 140:            additional 0.75x
```

**Design decision: selectivity as multiplier, not signal.** Selectivity is not one signal among seven. It is a gatekeeper. A compound with 30 nM IC50 but 0.02x selectivity is worse than useless — it is dangerous. Representing this as a multiplicative penalty (not an additive signal) means counter-selective compounds get their score collapsed to ~10% regardless of how high their other signals are. This is the correct behavior: you cannot compensate for toxicity with potency.

**Design decision: weights are not symmetric.** Potency and whole-organism evidence each get 0.25 (together 50%). This reflects a deliberate choice: direct experimental evidence of compound-target or compound-organism interaction is worth more than indirect signals like structural similarity or drug development stage. A compound with 10 nM IC50 and 77 nM whole-organism activity should dominate, all else being equal.

---

## The key engineering decisions and their rationale

### Why no ML model?

With 176 compounds, 7 features (most sparse), and a benchmark of 228 entries (110 positives, 50 negatives), there is not enough data for supervised learning. Any model complex enough to capture non-linear interactions between signals would overfit. The weighted sum with domain-knowledge weights is the appropriate model complexity for this data regime. The weights encode scientific priors (target essentiality, evidence type hierarchy) that a learned model would need thousands of labeled examples to discover.

If the dataset grows to 10,000+ compounds with dense feature coverage, a learned model becomes appropriate. The scoring function's signal structure (7 named features per compound) is designed to be directly usable as input features for a future learned ranker.

### Why three evidence pathways instead of one?

Praziquantel — the actual standard of care — has zero target-based IC50 data. Its molecular mechanism is still not characterized. A pipeline with only target-based evidence ranks it last (193/194). A pipeline with target-based + structural similarity ranks it 132nd. A pipeline with all three pathways ranks it 44th.

This is not a failure mode specific to praziquantel. It is a general property of drug discovery databases: different drugs have evidence in different modalities. Target-based IC50 data exists for drugs discovered through rational design. Whole-organism data exists for drugs discovered through phenotypic screening. Structural similarity bridges both. Any single-pathway pipeline will systematically miss drugs whose evidence lives in a different modality.

The mathematical structure is simple: the composite score is a weighted sum over all pathways. If a pathway has no data for a compound, that signal is zero and the compound's score is driven by the remaining pathways. The weights are chosen so that no single pathway dominates — a compound needs evidence from at least two pathways to rank highly.

### Why is selectivity the most important filter?

This was discovered empirically, not assumed. Before Script 10, the pipeline's top candidates included quisinostat (28 nM IC50, Phase II drug, 0 schistosomiasis publications). It looked like the best finding. After selectivity analysis: quisinostat has IC50 <1 nM against human HDAC8 and 28 nM against parasite HDAC8. Selectivity ratio: 0.02. It is 50x more potent against the patient's protein than the worm's protein. It would kill the patient before affecting the worm.

The broader pattern: 40 out of 91 compounds with dual-species data are counter-selective. This is not a rare edge case. It is the dominant failure mode for compounds tested against parasite targets that have human orthologues. Any repurposing pipeline that skips selectivity analysis is systematically promoting toxic candidates.

### Why hand-curated supply chain data?

There is no API for "is this drug available in a rural East African health center." The WHO Essential Medicines List is a PDF document. African formulary data is fragmented across national drug authority websites, procurement records, and institutional knowledge. The supply chain layer in Kira is manually curated from the WHO EML 2023 edition and from direct knowledge of the East African health system.

This is a deliberate design choice: encode deployment reality as structured data even when it cannot be programmatically queried. A drug repurposing pipeline that identifies atovaquone as a candidate but cannot tell you whether it is available in the target country is solving only half the problem.

---

## The evaluation setup

### Benchmark construction

228 compounds. Three classes:
- **Active (110):** IC50 < 1000 nM against any S. mansoni target in ChEMBL, plus 3 curated known drugs (praziquantel, oxamniquine, mefloquine) based on clinical evidence.
- **Weak (68):** IC50 1000-10000 nM.
- **Inactive (50):** 5 experimental negatives (tested in whole-worm assays, <30% max activity) + 45 assumed negatives (approved drugs from unrelated therapeutic areas: diabetes, depression, hypertension, etc.).

80/20 stratified split. Fixed seed. Curated drugs in both sets. Inclusion criteria documented and frozen.

### Metrics

AUROC = 1.000 (bootstrap 95% CI: [1.000, 1.000]). This is the correct number and the correct interpretation is: **the pipeline perfectly separates compounds with anti-schistosomal evidence from unrelated drugs.** It does NOT mean the pipeline is a perfect schistosomiasis drug predictor. The discrimination task is easy because the positive and negative classes occupy completely different regions of chemical and pharmacological space.

The hard discrimination — separating potent-and-selective compounds from potent-but-counter-selective compounds — is what the selectivity analysis addresses. There, the picture is much less clean: 27 counter-selective compounds have Tanimoto > 0.3 to a selective compound, meaning they are structurally similar but translationally opposite. The pipeline's composite score + selectivity multiplier handles this correctly (counter-selectives get 0.1x), but this is by design, not by learned discrimination.

### Sanity anchors

Two drugs function as integrity checks:

**Praziquantel** must rank well among translational candidates despite having no target-based IC50. Its journey through pipeline versions (rank 193 → 132 → 44 → translational #1) traces the progressive addition of evidence pathways. It reaches #1 among translational candidates because its supply chain score is 1.000 (WHO EML, donated, widely available, oral).

**Atovaquone** must rank well among all candidates. It has target-based IC50 (430 nM), confirmed selectivity (6.0x), WHO EML status, and African availability. It ranks #2 among translational candidates.

If either anchor fails, the pipeline is considered broken regardless of aggregate metrics.

---

## The central finding

### Claim

Systematic cross-species selectivity analysis of 91 anti-schistosomal compounds reveals that SmHDAC8 chemical matter is overwhelmingly non-selective (90.4%), while SmDHODH harbors compounds with >10-fold parasite selectivity.

### Evidence

**SmHDAC8:** 73/100 compounds have dual-species data. 66/73 (90.4%) are poor or counter-selective. Median ratio: 1.0x. All three clinical HDAC inhibitors (vorinostat, panobinostat, quisinostat) are severely counter-selective (ratios ~0.00-0.01). Only 1 compound exceeds 10x selectivity.

**SmDHODH:** 18/21 compounds have dual-species data. 7/18 (38.9%) are poor or counter-selective. Median ratio: 4.7x. Three compounds exceed 10x: CHEMBL155771 (30.8x, 23 nM, QED 0.89), CHEMBL4474026 (20.3x, 227 nM), CHEMBL4452960 (10.4x, 78 nM). Atovaquone: 6.0x, 430 nM, approved, WHO EML.

**SmTGR:** 43 compounds, zero experimental dual-species data. Assessed computationally via molecular docking (AutoDock Vina) against SmTGR (PDB 2X99, GSH binding site) and human TrxR1 (PDB 2ZZC, FAD binding site). 38/43 (88%) showed equivalent or preferential binding to the human enzyme. Mean docking energy: -6.4 kcal/mol (SmTGR) vs -7.2 kcal/mol (human TrxR1). Most potent compound (CHEMBL3322287, 10 nM experimental IC50) had delta -1.3 kcal/mol favoring the human target. Docking validated by significant correlation with experimental IC50 (Spearman rho = 0.367, p = 0.016). Caveats: rigid docking, non-equivalent binding sites compared, fusion interface not targeted.

### Why this matters

SmHDAC8 has received the most medicinal chemistry investment of any schistosomiasis target. Over 100 compounds have been designed and tested. The implicit assumption in this body of work is that SmHDAC8 inhibitors can be developed into anti-schistosomal drugs. The selectivity data contradicts this assumption for 90% of the chemical matter. The field has been optimizing potency against a target where the selectivity barrier is the binding constraint, not potency.

SmTGR, despite being the most biologically compelling target (single point of failure in parasite redox defense, unique fusion architecture), shows a similar pattern computationally: 88% of compounds bind the human enzyme at least as well. The current SmTGR chemical matter does not inherently exploit the fusion architecture for selectivity. Future selective SmTGR inhibitors would need to target the TrxR-GR domain interface specifically.

SmDHODH has received less attention but has a fundamentally better selectivity landscape. The parasite and human enzymes are structurally different enough that selective inhibition is achievable without heroic medicinal chemistry.

This reframing — from "which compounds are most potent?" to "which targets allow selective compounds?" — is the analytical contribution. The answer is now complete across all three major targets: SmHDAC8 fails experimentally, SmTGR fails computationally, SmDHODH passes both.

---

## The failure modes

### Known failure modes the pipeline handles correctly

1. **Missing evidence type.** Praziquantel has no target-based IC50. The pipeline recovers it through whole-organism evidence (77 nM) and structural similarity (0.524). Final rank: translational #1.

2. **Impressive but toxic compounds.** Quisinostat: 28 nM IC50, Phase II, 0 schistosomiasis publications. Looks like the best finding. Selectivity ratio: 0.02. Pipeline correctly eliminates it via 0.10x multiplier.

3. **Structurally similar compounds with opposite selectivity.** CHEMBL3798937 (Tanimoto 0.755 to selective compound CHEMBL4855490) but counter-selective at 0.31x. Pipeline separates them correctly via the selectivity multiplier.

4. **Docking setup with wrong active site coordinates.** First SmTGR docking run used estimated coordinates for human TrxR1 — 85 Angstroms from the actual FAD site. All 43 compounds scored 0.0 kcal/mol against the human receptor, producing false selectivity ratios of 500,000x. Diagnosed by inspecting HETATM records in PDB files. Fixed by centering on co-crystallized ligand positions. Second run produced real energies and the honest result: 88% non-selective.

### Known failure modes the pipeline does NOT handle

1. **Novel scaffolds with no reference.** A drug with a completely new chemical scaffold will have low Tanimoto to all known actives. If it also lacks target-based IC50 and whole-organism data, it scores near zero regardless of true potential. The pipeline is biased toward the known chemical space.

2. **SmTGR selectivity.** Phase 2 docking (Script 14) revealed that 88% of SmTGR compounds computationally prefer the human TrxR1 receptor. Docking scores correlated with experimental IC50 (Spearman rho=0.37, p=0.016), validating the setup. The non-selectivity finding is genuine but comes with caveats: rigid docking, non-equivalent binding sites (GSH vs FAD), and the unexplored fusion interface.

3. **Assay condition heterogeneity.** Selectivity ratios are computed from ChEMBL data generated across different labs, different assay formats, different buffer conditions, different temperatures. A selectivity ratio of 6x computed from two measurements by different labs in different years is less reliable than a 6x ratio from matched assays in the same lab. The pipeline treats all ratios equally.

4. **Negative data scarcity.** 5 experimental negatives out of 50 total inactives. The other 45 are assumed negatives — drugs for diabetes, depression, etc. that we assume are inactive against the worm. This is a reasonable assumption but not experimentally verified. The benchmark is easy to pass.

5. **Deployment data incompleteness.** Supply chain data is curated for a handful of named drugs. The vast majority of the 137-compound shortlist has no supply chain information because they are unnamed research compounds that have never been manufactured for clinical use.

---

## What would make it stronger

Listed in order of impact per unit of effort:

1. **Experimental selectivity data for SmTGR.** Test the top 10 SmTGR inhibitors against human thioredoxin reductase 1 in a matched biochemical assay. If any show >10x selectivity, SmTGR becomes an immediately actionable target class. If none do, the 43-compound SmTGR series is deprioritized. One experiment unlocks an entire target class.

2. **Whole-worm validation of top SmDHODH compounds.** CHEMBL155771 (30.8x selective, 23 nM) has never been tested in a living worm. Target-based potency does not guarantee whole-organism killing (the compound might not penetrate the worm's tegument). One assay converts a computational prediction into an experimental finding.

3. **Hard negative expansion.** Mine ChEMBL for compounds tested against parasite targets with IC50 > 100 μM (confirmed failures). These are experimentally verified negatives that test the pipeline's ability to discriminate within the bioactive chemical space, not just between bioactive and unrelated drugs.

4. **AlphaFold 3 docking for data-desert targets.** Five parasite targets have zero screening data: SmTPx, calcium channels (×2), SmNTPDase, SmCB1. AF3 can predict their structures. Virtual screening (docking approved drugs against predicted structures) generates binding hypotheses where no experimental data exists. This is the natural integration point for structural AI.

5. **Learned scoring function.** With enough data (expanded evaluation set, additional targets, more selectivity measurements), replace the hand-designed weighted sum with a gradient-boosted tree or similar low-complexity learned model. The current 7-signal feature vector is directly usable as input features.

---

## The implementation reality

**Runtime:** ~40 minutes total (network I/O to ChEMBL/PubMed APIs + ~4 minutes docking). All cached after first run — subsequent runs take <10 minutes.

**Compute:** Apple M4 Pro, CPU only. 48 GB RAM (uses <4 GB). No GPU, no cloud, no cluster.

**Dependencies:** Python 3.11, pandas, RDKit (cheminformatics), chembl_webresource_client, scikit-learn (metrics only), requests, NumPy, AutoDock Vina (molecular docking), meeko (ligand preparation), gemmi (macromolecular structures). Standard scientific Python stack.

**Code:** 14 scripts, each 200-500 lines. Linear pipeline — each script reads the previous script's output CSVs. No configuration files, no orchestration framework, no database. Deliberate simplicity: the pipeline should be readable end-to-end by a single person in one sitting.

**Data volume:** ~250 MB raw (PrimeKG), ~50 MB processed (all CSVs combined). Fits in memory on any modern laptop.

**Reproducibility:** Fixed random seeds. Deterministic pipeline. Every intermediate result saved as CSV. Git history records every step.

---

## The output

Two-tier shortlist:

**Translational tier (3 approved drugs):**
1. Praziquantel — standard of care, WHO EML, deployment score 0.459
2. Atovaquone — SmDHODH, 430 nM, 6.0x selective, WHO EML, deployment score 0.443
3. Idebenone — sirtuin, 1900 nM, weak, not available in Africa, score 0.210

**Discovery tier (top 4, all ≥10x selective):**
1. CHEMBL155771 — SmDHODH, 23 nM, 30.8x selective, QED 0.89
2. CHEMBL4855490 — SmHDAC8, 100 nM, 11.0x selective (1 of 73 — the needle)
3. CHEMBL4452960 — SmDHODH, 78 nM, 10.4x selective, QED 0.89
4. CHEMBL4474026 — SmDHODH, 227 nM, 20.3x selective, QED 0.88

Plus 101 discovery-tier compounds and 33 deprioritized (poor selectivity).

**SmTGR docking selectivity (Phase 2):**
- 43 compounds docked against SmTGR (PDB 2X99) and human TrxR1 (PDB 2ZZC)
- 88% (38/43) computationally non-selective
- Spearman rho = 0.37 (p = 0.016) validating setup
- SmTGR current chemical matter does NOT show inherent selectivity

---

## The one-paragraph summary

Kira is a 14-script computational pipeline that integrates target-based bioactivity, whole-organism phenotypic evidence, structural similarity, ADMET drug-likeness, cross-species selectivity analysis, literature novelty screening, and supply chain feasibility assessment to rank drug repurposing candidates for schistosomiasis. Its central finding is that two of three major targets present severe selectivity barriers: SmHDAC8 experimentally (90.4% non-selective) and SmTGR computationally via molecular docking (88% non-selective), while SmDHODH harbors compounds with >10-fold parasite selectivity including a 23 nM inhibitor with 30.8x selectivity and near-perfect drug-likeness. The pipeline runs in 35 minutes on a laptop (plus ~5 minutes for docking), uses no ML training, and produces a deployment-adjusted shortlist that encodes East African clinical reality. Built by one person in Nashville, TN.
