# Kira: Ground Level to God Level
## The Complete Scientific and Technical Breakdown
### Everything we built. Everything it means. From atoms to institutions.

*Daniel Ngabonziza — Nashville, TN — March 2026*

---

# PART ONE: THE SCIENCE BENEATH THE PIPELINE

Everything Kira does rests on a chain of scientific concepts. Each concept builds on the one below it. If you understand this chain from the bottom up, you understand not just what the pipeline does but why every decision in it was made.

---

## LEVEL 0: WHAT IS LIFE AT THE MOLECULAR LEVEL?

Your body is made of approximately 37 trillion cells. Each cell is not a blob of jelly — it is a factory. Inside every cell are thousands of molecular machines doing specific jobs. Some machines read DNA. Some machines build other machines. Some machines generate energy. Some machines transport cargo across membranes. Some machines detect signals from the outside world. Some machines destroy invaders.

These machines are called proteins.

A protein is a long chain of smaller units called amino acids, folded into a precise three-dimensional shape. There are 20 different amino acids — think of them as 20 different Lego brick types. A typical protein is a chain of 200 to 1000 of these bricks, folded into a unique shape. The sequence of bricks determines the shape, and the shape determines the function.

This is the most fundamental fact in molecular biology: shape determines function. A wrench works because its jaw is shaped to grip a bolt. An enzyme works because its active site is shaped to grab a specific molecule and catalyze a specific chemical reaction. An antibody works because its binding region is shaped to latch onto a specific pathogen surface. Change the shape, change the function.

DNA is the blueprint that specifies which amino acids go in which order. When a cell needs to build a protein, it reads the relevant section of DNA (a gene), transcribes it into a messenger molecule (mRNA), and then a molecular factory called the ribosome reads the mRNA and chains together the amino acids in the specified order. The chain then folds into its functional shape — a process that takes milliseconds to seconds and that, until recently, was extraordinarily difficult to predict computationally. AlphaFold changed that.

Every living organism — from bacteria to humans to parasitic worms — runs on proteins. The specific set of proteins an organism has determines what it can do: what environments it can survive in, what food it can eat, what threats it can defend against, how it reproduces.

This is why drugs work. A drug is a small molecule designed to fit into a protein's active site and block it from working. If you block the right protein in a pathogen, the pathogen dies. If you block the right protein in a cancer cell, the cell stops dividing. If you block the right protein in your own body, you can reduce inflammation, lower blood pressure, or suppress an immune reaction.

The challenge — and it is the central challenge of all pharmacology — is specificity. Your proteins and the pathogen's proteins evolved from common ancestors. They are similar. A drug designed to block a worm protein might also block the equivalent human protein. That is toxicity. That is the selectivity problem. And it turned out to be the single most important finding of everything we built.


## LEVEL 1: WHAT IS SCHISTOSOMIASIS?

Schistosomiasis is a disease caused by parasitic flatworms of the genus Schistosoma. Schistosoma mansoni — the species Kira targets — is a blood fluke. Not a bacterium. Not a virus. A multicellular animal with its own digestive system, nervous system, and reproductive system, roughly 1-2 centimeters long, living inside your blood vessels.

The life cycle is grotesque and elegant. The worm's eggs are shed in the feces of an infected person into freshwater. In water, the eggs hatch into free-swimming larvae called miracidia. These penetrate a specific freshwater snail species (Biomphalaria), multiply inside the snail, and emerge as a different larval form called cercariae. When a human wades, swims, or washes clothes in that water, cercariae penetrate the skin — you can be infected in under a minute. Inside the human, cercariae transform into schistosomula, migrate through the lungs and liver, mature into adult worms, pair up (male and female), and settle in the mesenteric and portal blood vessels around the intestine. There they feed on blood and produce eggs — hundreds per day — for years.

The eggs are what make you sick. Some eggs pass into the intestine and out of the body (completing the cycle). But many get swept by blood flow into the liver and other organs, where they become trapped. Your immune system recognizes them as foreign and surrounds each egg with a granuloma — a ball of immune cells. Over years, thousands of granulomas in your liver cause fibrosis (scarring), portal hypertension (increased blood pressure in liver vessels), and eventually organ failure. In the bladder (S. haematobium, a related species), chronic egg deposition causes inflammation that can progress to bladder cancer.

Approximately 200 million people are infected, overwhelmingly in sub-Saharan Africa. Children in rural areas near freshwater are most affected. The disease doesn't kill quickly — it disables slowly. Chronic fatigue, abdominal pain, bloody stool, anemia, impaired growth and cognitive development in children. It is a disease of poverty: it persists where sanitation is poor, where freshwater contact is unavoidable, and where health systems are thin.

The entire global treatment strategy depends on a single drug: praziquantel. One drug. For 200 million people. For over 40 years. The WHO has flagged this single-drug dependency as a resistance vulnerability. If praziquantel resistance emerges at scale — and there are early signals of reduced efficacy in some settings — there is no adequate backup.

That is the problem Kira addresses.


## LEVEL 2: WHAT IS A DRUG?

A drug is a small molecule — typically 150 to 500 daltons in molecular weight, containing 20 to 50 atoms — that interacts with a specific protein in a specific way to produce a therapeutic effect.

"Small" is relative. Proteins are large: tens of thousands of daltons. A drug molecule is like a key fitting into a lock. The lock is the protein's active site — a pocket or groove on its surface where the important chemistry happens. When the drug molecule enters the pocket and binds, it can block the protein from doing its job. This is inhibition.

How tightly the drug binds determines how potent it is. This is measured as IC50: the concentration of drug needed to reduce the protein's activity by half.

IC50 is measured in moles per liter (molarity). Because the numbers are very small, we use nanomolar (nM = 10^-9 moles per liter) or micromolar (μM = 10^-6 moles per liter). The scale:

- 1 nM IC50: extraordinarily potent. One billionth of a mole per liter shuts down half the enzyme.
- 10 nM: very potent.
- 100 nM: potent. Most successful drugs operate in this range.
- 1,000 nM (1 μM): moderate. A useful starting point for optimization.
- 10,000 nM (10 μM): weak. Marginally active. Might not work at achievable concentrations in the body.
- 100,000 nM (100 μM): essentially inactive.

Medicinal chemists use pIC50 — the negative log base 10 of the IC50 in molar — so that higher numbers mean more potent:

    IC50 = 10 nM → pIC50 = 8.0
    IC50 = 100 nM → pIC50 = 7.0
    IC50 = 1,000 nM → pIC50 = 6.0

Kira converts every IC50 to pIC50, then normalizes to a 0-1 score. This is Signal 1 in the composite ranking: potency.

But potency alone is not enough. A drug must also be absorbed by the body, distributed to the site of infection, not destroyed too quickly by the liver, safely excreted, and non-toxic. And it must hit the parasite's protein without hitting the patient's protein. That last requirement — selectivity — is where Kira made its most important discovery.


## LEVEL 3: WHAT IS A DRUG TARGET?

A drug target is a specific protein in the disease-causing organism that you aim to inhibit with a drug. Not every protein is a good target. Three properties make a target worth pursuing.

**Essentiality.** If you knock the protein out, does the organism suffer or die? SmTGR (thioredoxin glutathione reductase) is essential because it is the worm's only defense against oxidative attack from your immune system. It is a single point of failure — the worm has no backup antioxidant system. Block SmTGR, and the worm burns.

**Druggability.** Does the protein have a pocket that a small molecule can fit into? Some proteins have wide, flat surfaces with no crevices — impossible for a small molecule to grip. SmHDAC8 (histone deacetylase 8) has a deep, narrow active site channel — ideal for small molecules. That is why it has the most chemical data of any schistosomiasis target: 100 compounds tested, 146 activity measurements.

**Selectivity.** Can you inhibit the parasite version without inhibiting the human version? This is the requirement that Kira found most targets fail. Humans have their own HDAC8, their own DHODH, their own thioredoxin reductase. If a drug hits both versions equally, it's toxic. If it hits the parasite version 10 times more potently than the human version, you have a therapeutic window: at a dose that fully blocks the parasite enzyme, only 10% of the human enzyme is affected. That window is what makes a drug candidate viable.

Kira identified 10 Schistosoma mansoni protein targets. Here is each one, what it does for the worm, and what we learned about it:

**SmTGR (Thioredoxin Glutathione Reductase)** — Redox defense. When your immune cells attack the worm, they fire reactive oxygen species — superoxide, hydrogen peroxide — molecular fire. Most organisms defend against this with two separate enzyme systems: thioredoxin reductase and glutathione reductase. The schistosome uniquely fused both into a single enzyme. This means SmTGR is a single point of failure with no backup. It is arguably the most important target in schistosomiasis drug development. 43 compounds tested. Best hit: 10 nM (extraordinarily potent). Selectivity status: UNKNOWN. Not a single compound has been tested against human thioredoxin reductase 1 in ChEMBL. This is the most critical experimental gap Kira identified.

**SmHDAC8 (Histone Deacetylase 8)** — Gene regulation. DNA wraps around spool-like proteins called histones. When acetyl chemical groups are attached to histones (acetylation), DNA is loose and genes are active. HDAC enzymes remove these acetyl groups, tightening the DNA and silencing genes. The worm needs SmHDAC8 to switch genes on and off across its complex life cycle — larva, schistosomulum, adult. 100 compounds tested. Best hit: 28 nM. Selectivity: DEVASTATING. 90.4% of compounds with human data are non-selective or counter-selective. Median selectivity ratio: 1.0x. The typical SmHDAC8 compound hits the human enzyme exactly as hard as the parasite enzyme. This is Kira's central finding.

**SmDHODH (Dihydroorotate Dehydrogenase)** — Pyrimidine synthesis. Pyrimidines are half the building blocks of DNA and RNA (cytosine, thymine, uracil). DHODH catalyzes the fourth step in making them. Block it, and the worm cannot replicate DNA or make RNA. 21 compounds tested. Best hit: 19 nM. Selectivity: FAVORABLE. Median ratio 4.7x. Three compounds exceed 10x selectivity. The best — CHEMBL155771 — achieves 30.8x. Atovaquone, an approved antimalarial, shows 6x selectivity. SmDHODH is Kira's recommended priority target.

**SmTPx (Thioredoxin Peroxidase)** — Downstream of SmTGR. Directly neutralizes hydrogen peroxide using SmTGR's output. Zero screening data in ChEMBL.

**Calcium Channel Beta Subunits 1 and 2** — Muscle contraction. These are the likely targets of praziquantel, which causes massive calcium influx into the worm's muscle, paralyzing it. Zero screening data. The mechanism of the world's only schistosomiasis drug is still not characterized at the molecular level.

**Sirtuin (NAD-dependent Deacetylase)** — Stress response and metabolic regulation. A different class of deacetylase from HDACs. Only 5 compounds tested, all weak (best: 1,900 nM). Barely explored.

**SmNTPDase (ATP-diphosphohydrolase 1)** — Immune evasion. Sits on the worm's outer surface and degrades extracellular ATP, which is an immune danger signal. By clearing ATP, the worm dampens the local immune response. Minimal data.

**SmCB1 (Cathepsin B1)** — Gut protease. The worm feeds on blood and uses SmCB1 to digest hemoglobin. Block it, and the worm starves. Zero ChEMBL data, but published literature on cathepsin inhibitors exists.

**SmVKR2 (Venus Kinase Receptor 2)** — Reproduction. No human equivalent exists (invertebrate-specific). Controls egg production. Since eggs cause the disease symptoms, blocking reproduction could be therapeutic without killing the worm. 7 compounds tested. Best: 440 nM. Perfect selectivity opportunity because there's no human version to accidentally hit.


## LEVEL 4: WHAT IS DRUG REPURPOSING?

Traditional drug development starts from scratch. Discover a target. Screen millions of compounds. Optimize hits through years of medicinal chemistry. Test in cells, then animals, then Phase I (safety in healthy volunteers), Phase II (does it work in patients?), Phase III (does it work better than alternatives, across thousands of patients?). Timeline: 10-15 years. Cost: over $1 billion. Failure rate: over 90%.

For neglected tropical diseases like schistosomiasis, this model is broken. The populations most affected cannot pay prices that recoup billion-dollar investments. Pharmaceutical companies rationally — if unjustly — decline to develop drugs for diseases of poverty. This is the market failure at the heart of the NTD crisis.

Drug repurposing offers a shortcut. If a drug is already approved for one disease, its safety profile is established. Manufacturing processes exist. Formulations are optimized. Supply chains are in place. If that drug also happens to inhibit a parasite protein, you can potentially skip most of the development timeline and go nearly straight to Phase II trials for the new indication.

The question is: how do you systematically search the space of "existing approved drugs × parasite targets" to find promising matches?

That is what Kira does. And it does it not by asking only "does this drug hit the parasite target?" but by asking the harder question: "does this drug hit the parasite target WITHOUT also hitting the human version?"


## LEVEL 5: WHAT IS SELECTIVITY AND WHY DID IT CHANGE EVERYTHING?

Selectivity is the single most important concept in translational drug development. It is the ratio of a compound's potency against the human version of a protein versus the parasite version.

    Selectivity Ratio = Human IC50 / Parasite IC50

If a compound has IC50 = 100 nM against SmHDAC8 and IC50 = 1,000 nM against human HDAC8, the ratio is 10. The compound is 10 times more potent against the parasite enzyme. At a dose that fully blocks the parasite enzyme, only one-tenth of the human enzyme is affected. That is a therapeutic window.

If the ratio is 1, the compound hits both equally. Every dose that helps the patient also poisons them.

If the ratio is less than 1, the compound actually prefers the human target. It is a human enzyme inhibitor that weakly cross-reacts with the parasite version. This is counter-selective. It is dangerous.

Kira's Script 10 queried ChEMBL for activity data against human orthologues of each parasite target. The results:

**SmHDAC8: 73 compounds with dual-species data.**
- 1 selective (≥10x): CHEMBL4855490, ratio 11.0x
- 6 moderate (3-10x)
- 30 poor (1-3x)
- 36 counter-selective (<1x)
- Non-selective rate: **90.4%**
- Median ratio: **1.0x**

The three most famous HDAC inhibitors in clinical development — vorinostat (IC50 590 nM vs SmHDAC8, ~1 nM vs human HDAC8), panobinostat (~450 nM vs parasite, ~5 nM vs human), and quisinostat (28 nM vs parasite, <1 nM vs human) — are among the most severely counter-selective compounds in the dataset. They are potent human HDAC8 inhibitors that weakly cross-react with the parasite. They are not schistosomiasis drug candidates. They are cancer drugs.

**SmDHODH: 18 compounds with dual-species data.**
- 3 selective (≥10x)
- 8 moderate (3-10x)
- 3 poor (1-3x)
- 4 counter-selective (<1x)
- Non-selective rate: **38.9%**
- Median ratio: **4.7x**

The best: CHEMBL155771, 23 nM against parasite DHODH, 709 nM against human DHODH. 30.8x selective. MW=244, QED=0.89. Zero Lipinski violations. Zero schistosomiasis publications. This is Kira's strongest novel finding.

Atovaquone: 430 nM parasite, 2,600 nM human. 6.0x selective. Already approved. Already on the WHO Essential Medicines List. Already available in sub-Saharan Africa as Malarone. This is Kira's strongest translational finding.

**SmTGR: 0 compounds with dual-species data.** Nobody has tested any SmTGR inhibitor against human thioredoxin reductase 1. Selectivity is unknown, not confirmed. The biology is favorable — SmTGR's unique fusion architecture has no direct human equivalent — but the data doesn't exist.

This selectivity finding is the core of the Kira publication. It was not previously documented as a systematic analysis across the entire ChEMBL chemical series for these targets.

---

# PART TWO: THE DATABASES AND DATA SOURCES

Kira doesn't generate data from scratch. It integrates, filters, cross-references, and adjudicates data from multiple public databases. Understanding what each database contains and what it misses is essential to understanding what Kira can and cannot do.

---

## LEVEL 6: PrimeKG — THE KNOWLEDGE GRAPH

PrimeKG, built by Harvard's Zitnik Lab, is a biomedical knowledge graph containing approximately 8 million relationships between biomedical entities. Each entry says: "Entity A has relationship R with Entity B."

Types of entities: diseases, genes/proteins, drugs, biological processes, pathways, anatomical structures. Types of relationships: disease-protein association, drug-disease indication, protein-protein interaction, disease-disease similarity, drug-gene interaction.

When Kira queries PrimeKG for schistosomiasis, it finds 52 edges: 6 human genes (NOS2, SOD1, TLR9, PRDX1, TIMP3, RASSF1), 3 known drugs (praziquantel, oxamniquine, stibophen), and a taxonomy of disease subtypes (intestinal, urinary, neuro).

**The critical lesson from PrimeKG:** It maps human biology, not parasite biology. The 6 genes are human immune response genes, not worm proteins. NOS2 is your nitric oxide synthase that kills parasites. SOD1 is your antioxidant defense during infection. TLR9 is your immune receptor that detects pathogen DNA. PrimeKG says nothing about SmTGR, SmHDAC8, or any protein inside the worm.

This is the first evidence pathway limitation: a human-centric knowledge graph cannot provide parasite-specific drug target information. That sent us to ChEMBL.


## LEVEL 7: ChEMBL — THE EXPERIMENTAL DATABASE

ChEMBL, maintained by the European Bioinformatics Institute, is fundamentally different from PrimeKG. It stores experimental measurements: "Compound X was tested against protein target Y in assay Z and showed IC50 of W nanomolar."

The data comes from published scientific papers. When a researcher synthesizes a compound, tests it against a purified protein, measures the IC50, and publishes, that data gets curated into ChEMBL. This is primary experimental evidence, not inferred relationships.

Kira queries ChEMBL in three ways:

**Target-based activity (Script 02).** Query: "What protein targets from Schistosoma mansoni exist in ChEMBL?" Result: 10 parasite proteins. For 5 of them, there are dose-response activity measurements — 226 quality-filtered data points across 176 compounds.

**Whole-organism activity (Script 06).** Query: "What compounds have been tested against live Schistosoma mansoni worms?" Result: 494 raw records. After filtering for dose-response measurements in consistent units, 27 records across 6 compounds survive. Praziquantel has 74 raw records (by far the most) with a best activity of 77 nM. This is the evidence that rescued praziquantel in the ranking — it works against the whole worm even though its molecular target isn't characterized.

**Human orthologue activity (Script 10).** Query: "For each compound tested against a parasite target, has it also been tested against the human version of the same protein?" Result: 517 records. 73 compounds have data against both SmHDAC8 and human HDAC8. 18 have data against both SmDHODH and human DHODH. Zero have data against both SmTGR and human TrxR1.

**What ChEMBL misses:** Drugs with phenotypic evidence but no target-based IC50 (like praziquantel — it kills worms but nobody knows exactly which protein it hits). Compounds tested in assays that don't produce clean dose-response curves. Negative results (compounds that were tested and failed — these rarely get published or deposited). This publication bias means ChEMBL over-represents positive results and under-represents failures, which affects evaluation set construction.


## LEVEL 8: PubMed — THE LITERATURE

PubMed, maintained by the US National Library of Medicine, indexes over 36 million biomedical publications. Kira uses it for novelty screening: for each top-ranked compound, query "[compound name] AND schistosomiasis" and count publications.

Praziquantel: 4,030 papers. Obviously well-known. Vorinostat: 4 papers. Atovaquone: 4 papers. Quisinostat: 0 papers. Most unnamed research compounds: 0 papers.

**What PubMed tells you and what it doesn't:** Zero publications for a compound means the exact string query didn't match. It does not definitively mean nobody has ever considered the compound for schistosomiasis — they might have used a different name, discussed the compound class without naming individual molecules, or published in a journal not indexed by PubMed. The novelty filter is a useful first screen, not definitive proof of novelty.


## LEVEL 9: WHO ESSENTIAL MEDICINES LIST — THE DEPLOYMENT REALITY

The WHO Model List of Essential Medicines (EML) is a curated list of medications considered most important for a basic health system. It is updated every two years. If a drug is on the EML, it is more likely to be available in government health systems across the developing world, including sub-Saharan Africa.

This is domain knowledge that cannot be queried from any API. It was manually curated for Kira from the WHO EML 2023 edition:

- Praziquantel: on the EML (Section 6.1, Antihelminthics). Donated by Merck through WHO. Widely available across sub-Saharan Africa. Oral tablet, 600mg. Cost: essentially free for endemic country programs.
- Atovaquone: on the EML (Section 6.5.3, Antimalarials). Available as Malarone in many African countries. Oral. Cost: moderate (more expensive than artemisinin combinations). Limited availability compared to praziquantel.
- Vorinostat, panobinostat, quisinostat, idebenone: NOT on the EML. Not available in African health systems. Cancer or specialty drugs with costs that exclude deployment in NTD settings.

This layer is what makes Kira contextual rather than generic. A drug repurposing pipeline built in Boston would tell you "vorinostat has a 590 nM IC50 against SmHDAC8." Kira tells you "vorinostat has a 590 nM IC50 against SmHDAC8, is counter-selective against human HDAC8, is not on the WHO EML, is not available in African health systems, and costs over $100 per treatment. It is not a viable schistosomiasis drug candidate." That additional context is the deployment reality layer, and it requires knowing what a clinic in Kigali actually looks like — knowledge that comes from growing up in Rwanda, not from querying a database.

---

# PART THREE: THE COMPUTATIONAL METHODS

Kira uses four computational methods, each providing a different kind of evidence about drug candidates.

---

## LEVEL 10: MOLECULAR FINGERPRINTING AND STRUCTURAL SIMILARITY

Every molecule has a structure — a specific arrangement of atoms connected by chemical bonds. This structure is encoded as a SMILES string (Simplified Molecular Input Line Entry System), a text representation of the molecule's connectivity. Praziquantel's SMILES:

    O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1

This encodes: a carbonyl group (O=C) attached to a cyclohexane ring (C1CCCCC1), connected through nitrogen to a bicyclic system containing a benzene ring (c3ccccc3) and a lactam.

RDKit, the open-source cheminformatics toolkit, converts SMILES into molecular fingerprints. A Morgan fingerprint (also called ECFP) works by examining each atom and recording the chemical environment around it out to a specified radius. At radius 2, it captures each atom plus everything within 2 bonds. These local environments are hashed into a bit vector — in our case, 2048 bits. Each bit represents the presence or absence of a particular molecular substructure.

Two molecules with similar fingerprints share similar substructures. The similarity is measured by the Tanimoto coefficient:

    Tanimoto = |bits set in both A and B| / |bits set in either A or B|

Range: 0 (completely different) to 1 (identical). In drug discovery, Tanimoto > 0.4 is considered structurally similar. In Kira, each compound's maximum Tanimoto similarity to any known potent active (IC50 < 1000 nM) is computed as the structural similarity signal.

This is what rescued praziquantel from rank 193. It had no target-based IC50, so the potency signal was zero. But its fingerprint has Tanimoto 0.524 to the nearest known active — it shares enough substructural features to be recognized as relevant. Structural similarity is Kira's second evidence pathway.


## LEVEL 11: ADMET AND DRUG-LIKENESS

ADMET stands for Absorption, Distribution, Metabolism, Excretion, Toxicity — the five pharmacokinetic hurdles every drug must clear.

**Absorption.** Can the drug get from the gut into the bloodstream? If swallowed as a pill, it must dissolve in intestinal fluid, cross the gut wall (a lipid membrane), and enter the blood. Lipinski's Rule of Five predicts oral absorption from four molecular properties:
- Molecular weight < 500 daltons (larger molecules don't cross the gut wall efficiently)
- LogP < 5 (LogP measures lipophilicity — how "greasy" a molecule is; too greasy and it won't dissolve in blood)
- Hydrogen bond donors < 5 (too many H-bond donors prevent membrane crossing)
- Hydrogen bond acceptors < 10 (same logic)

A compound violating more than one of these rules is unlikely to work as an oral drug. Oral administration is the only realistic option for mass drug administration in rural East Africa — you cannot set up IV infusions in village health posts.

**Distribution.** Does the drug reach the site of infection? Schistosome worms live in portal blood vessels. The drug must circulate through the hepatic portal system.

**Metabolism.** The liver actively destroys foreign chemicals using cytochrome P450 enzymes. A drug metabolized too quickly never reaches therapeutic levels.

**Excretion.** The kidneys clear waste. A drug that accumulates causes toxicity.

**Toxicity.** Does the drug damage human cells? This is where selectivity becomes a pharmacokinetic reality, not just a biochemical ratio.

RDKit computes all of these properties from the SMILES string: molecular weight, LogP, hydrogen bond donors and acceptors, topological polar surface area (TPSA — if above 140 Å², oral absorption is poor), rotatable bonds (more than 10 means the molecule is too floppy), and QED (Quantitative Estimate of Drug-likeness, a composite score from 0 to 1 where higher means more drug-like).

The key finding: 181 of 194 compounds in our evaluation set pass Lipinski with zero violations. Most compounds in ChEMBL are designed to be drug-like. The ADMET filter did not dramatically reshape the ranking because the chemical space was already drug-like. The selectivity filter was far more impactful.


## LEVEL 12: COMPOSITE SCORING — HOW KIRA RANKS COMPOUNDS

Kira's v3 ranking integrates seven normalized (0-1) signals into a single composite score via weighted summation:

**Signal 1: Potency (weight 0.25).** IC50 converted to pIC50, scaled to 0-1. 10 nM → 0.80. 100 nM → 0.60. 1,000 nM → 0.40.

**Signal 2: Target essentiality (weight 0.15).** Domain-knowledge score per target. SmTGR = 1.0 (single point of failure). SmHDAC8 = 0.85. SmDHODH = 0.80. SmVKR2 = 0.60. These are human judgments encoding what a parasitologist would tell you about each target's importance. They are not learned from data.

**Signal 3: Data confidence (weight 0.05).** More independent measurements → more confidence. Log-scaled, saturating. 1 measurement = 0.20. 10+ measurements = 0.80.

**Signal 4: Drug development stage (weight 0.10).** Approved = 1.0. Phase III = 0.7. Phase II = 0.5. Phase I = 0.3. Research compound = 0.1. Approved drugs are immediately actionable.

**Signal 5: Multi-target activity (weight 0.05).** Compounds hitting 2+ parasite targets get a bonus. Multi-mechanism attack is harder for the parasite to resist.

**Signal 6: Structural similarity (weight 0.15).** Maximum Tanimoto to any known potent active. Rescues compounds without IC50 data.

**Signal 7: Whole-organism activity (weight 0.25).** Whole-worm IC50/EC50 converted to same pIC50 scale as target-based potency. Captures drugs like praziquantel.

The composite score is then multiplied by:
- **Selectivity multiplier:** counter-selective (0.10), poor (0.50), moderate (0.85), selective (1.10), unknown (0.70-0.80)
- **ADMET penalty:** based on Lipinski violations and TPSA

The final product — the deployment score — also incorporates supply chain data for approved drugs: WHO EML status, African availability, cost tier, and route of administration.


## LEVEL 13: EVALUATION — HOW WE KNOW IF THE PIPELINE WORKS

An evaluation set is a labeled dataset where you know the correct answer. You hide the answers from the pipeline, run it, and check whether the ranking matches reality.

Kira's evaluation set v2 contains 228 compounds:
- 110 actives (IC50 < 1000 nM against a parasite target)
- 68 weak (1000-10000 nM)
- 50 inactives (5 experimental negatives from whole-worm screens + 45 assumed negatives from unrelated drug classes)

The set is split 80/20 into development (182 compounds) and held-out test (50 compounds) with a fixed random seed. Curated drugs (praziquantel, etc.) appear in both for sanity checking.

**AUROC** (Area Under the Receiver Operating Characteristic curve) measures whether the algorithm ranks actives above inactives. 0.5 = random guessing. 1.0 = perfect separation. Kira achieves 1.0 with bootstrap 95% CI [1.0, 1.0].

**Why 1.0 is honest but insufficient:** The discrimination task is easy. Active compounds have ChEMBL activity data and structural similarity to other actives. Inactive compounds are drugs for depression and diabetes with completely different structures. Distinguishing them is like telling chefs from plumbers by kitchen-related features. The hard test — distinguishing active compounds from structurally similar but inactive compounds — requires more experimental negatives. Kira identified 27 hard negatives (counter-selective compounds with Tanimoto > 0.3 to a selective compound) but these were discovered during the analysis, not used for evaluation.

**The honest statement:** "The pipeline perfectly separates known anti-schistosomal compounds from unrelated approved drugs. Discrimination against structurally similar but non-selective compounds remains undertested."

---

# PART FOUR: THE 13-SCRIPT PIPELINE

Each script built on the previous one. The git history records every step. Here is what each script does, what it discovered, and what lesson it taught.

---

## LEVEL 14: THE PIPELINE, SCRIPT BY SCRIPT

**Script 01 — PrimeKG Exploration.** Downloaded 8 million biomedical knowledge graph edges. Queried schistosomiasis. Found 52 edges: 6 human genes, 3 known drugs, disease taxonomy. **Lesson learned:** PrimeKG maps human biology, not parasite biology. The 6 genes are human immune response genes. The pipeline needs a different data source for parasite-specific targets.

**Script 02 — ChEMBL Parasite Targets.** Queried ChEMBL for Schistosoma mansoni protein targets. Found 10 parasite proteins. Retrieved 226 quality-filtered activity measurements across 176 compounds against 5 data-rich targets. Identified 4 approved drugs with measured parasite activity: atovaquone (430 nM vs SmDHODH), panobinostat (450 nM vs SmHDAC8), vorinostat (590 nM vs SmHDAC8), idebenone (1900 nM vs sirtuin). **Lesson learned:** The parasite targets exist and have real experimental data. But coverage is uneven — SmHDAC8 has 100 compounds while SmTPx has zero.

**Script 03 — Evaluation Set v1.** Built a 194-compound ground-truth dataset. 110 actives, 68 weak, 16 inactives (1 experimental + 15 assumed). Added praziquantel, oxamniquine, mefloquine as curated positives based on clinical evidence despite lacking target-based IC50. **Lesson learned:** Class imbalance (110 actives vs 16 inactives) limits evaluation reliability. Negative controls need expansion.

**Script 04 — Ranking Algorithm v1.** Built a five-signal composite scoring function. AUROC 0.985. **Critical discovery:** Praziquantel ranked 193 out of 194 — dead last, indistinguishable from metformin and sertraline. The algorithm was measuring "has ChEMBL data" vs "doesn't have ChEMBL data," not actual drug quality. **Lesson learned:** A high AUROC can hide a fundamental flaw. If the standard of care ranks last, the pipeline is broken regardless of the aggregate metric.

**Script 05 — Structural Similarity.** Added Morgan fingerprints and Tanimoto similarity as a second evidence pathway. Praziquantel rose from rank 193 to 132 (Tanimoto 0.524 to nearest active). AUROC improved to 0.999. **Lesson learned:** Fingerprint similarity helps when the unknown drug resembles known actives. It cannot fully rescue drugs with unique scaffolds and uncharacterized mechanisms.

**Script 06 — Whole-Organism Activity.** Added phenotypic worm-killing data as a third evidence pathway. Found 74 records for praziquantel (best activity 77 nM). Praziquantel rose from rank 132 to 44. **Lesson learned:** Multiple evidence pathways cover each other's blind spots. Target-based IC50, structural similarity, and whole-organism activity provide three independent lines of evidence.

**Script 07 — ADMET Filtering.** Added Lipinski Rule of Five, TPSA, QED. 181/194 compounds pass with zero violations. Atovaquone dropped from rank 5 to 61 due to LogP at the borderline. **Lesson learned:** ADMET catches real pharmacokinetic limitations (atovaquone genuinely has bioavailability challenges), but the penalty was too aggressive for a drug already proven to work in humans.

**Script 08 — Benchmark Hardening.** Expanded inactives from 16 to 50. Created train/test split. Computed bootstrap confidence intervals. Discovered that new negatives had been assigned zero scores (not run through the pipeline), artificially inflating AUROC. Fixed by scoring all compounds through the same pipeline. **Lesson learned:** Evaluation infrastructure is as important as the ranking algorithm. A perfect AUROC built on a flawed benchmark is worse than a modest AUROC built on a credible benchmark, because the former creates false confidence.

**Script 09 — Novelty Filter.** Queried PubMed for each top-50 compound. Found quisinostat as apparently novel (0 schistosomiasis publications, 106 general publications, Phase II drug). 45 of 50 "novel" but almost all are unnamed research compounds with zero general publications. **Lesson learned:** "Zero PubMed hits" means different things for named drugs (genuinely unexplored) vs unnamed research molecules (nobody writes papers about individual CHEMBL IDs). Novelty must be interpreted with context.

**Script 10 — Selectivity Analysis.** Queried ChEMBL for human orthologue activity data. This was the most important script. **Central finding:** 90.4% of SmHDAC8 chemical matter is non-selective. SmDHODH has a favorable selectivity profile. SmTGR selectivity is unknown. Quisinostat — the "novel" finding from Script 09 — has a selectivity ratio of 0.02 (50x more potent against human HDAC8 than parasite HDAC8). It was killed. **Lesson learned:** Selectivity is the translational filter that separates real candidates from shiny traps. Novelty without selectivity is dangerous. This script changed the entire project.

**Script 11 — Selectivity-Adjusted Re-Ranking.** Rebuilt the shortlist with counter-selective compounds eliminated (0.10x multiplier). 40 compounds removed. The top 4 are all genuinely selective. 137 compounds survive all filters. **Lesson learned:** When you add a hard translational constraint, the ranking reorganizes around viability rather than raw potency. Compounds with 30x selectivity and 23 nM IC50 now outrank compounds with 28 nM IC50 but 0.02x selectivity. That's correct behavior.

**Script 12 — Publication Analysis.** Resolved all compound identities (fixing "nan" names). Per-target selectivity statistics. Hard negatives identified (27 counter-selective compounds structurally similar to selective ones). Sanity anchor analysis. Discovery vs translational tier separation. Publication-ready tables. Working title and abstract. **Lesson learned:** The finding needs packaging. A pipeline is not a paper. Specific numbers, honest limitations, and clear scope statements convert a pipeline into a contribution.

**Script 13 — Supply Chain + Definitive Shortlist.** Fixed max_phase permanently by querying ChEMBL directly. Added WHO EML status, African availability, cost tier. Produced deployment-adjusted scores. Praziquantel ranks #1 among translational candidates (supply chain score 1.000). Atovaquone ranks #2 (score 0.690). **Lesson learned:** The supply chain layer is what makes Kira contextual. Anyone can query ChEMBL. Almost nobody encodes East African health system reality into the scoring function.

### Phase 2: Molecular Docking

**Script 14 — SmTGR Molecular Docking.** Downloaded crystal structures from the Protein Data Bank: SmTGR (PDB 2X99) and human thioredoxin reductase 1 (PDB 2ZZC). Docked all 43 SmTGR compounds against both targets using AutoDock Vina. **First run failed:** human receptor coordinates were estimated at (65, 45, 35) — 85 Angstroms from the actual FAD active site at (-20.0, 19.4, -43.3). Every human docking score was 0.0 kcal/mol. All 43 compounds appeared "parasite-selective" with ratios of 500,000x. This was an artifact, not biology. **Diagnosis:** inspected HETATM records in both PDB files to locate co-crystallized ligands. SmTGR had GSH (glutathione) at (34.4, 24.9, 6.5) and FAD at (19.9, 4.7, 15.7). Human TrxR1 had FAD at (-20.0, 19.4, -43.3). **Fix:** recentered docking boxes on actual ligand coordinates. **Second run produced real results:** SmTGR energies -5.0 to -8.5 kcal/mol, human TrxR1 energies -5.2 to -10.0 kcal/mol. **Central finding:** 38/43 compounds (88%) showed equivalent or preferential binding to the human enzyme. Only 5 showed marginally favorable parasite binding (delta +0.0 to +0.2 kcal/mol — essentially noise). The most potent compound in the dataset (CHEMBL3322287, 10 nM) had delta -1.3 kcal/mol — it computationally prefers the human target. **Validation:** Spearman rho = 0.367 (p = 0.016) between SmTGR docking energy and experimental IC50 — the docking captures real binding signal, confirming the non-selectivity is genuine, not a setup artifact. **Lesson learned:** When a result looks too good to be true (43/43 selective, ratios of 500,000x), it is. Diagnosing the failure (wrong coordinates), fixing it, and accepting the less exciting but honest result (88% non-selective) is how real science works. The failed first run taught more than the successful second run.

---

# PART FIVE: THE FINDINGS

## LEVEL 15: WHAT DID KIRA ACTUALLY DISCOVER?

Strip away all the infrastructure. What does the pipeline tell us about the world that we didn't know, or didn't know in this systematic form, before?

**Finding 1: SmHDAC8 is a selectivity trap.** 90.4% of the chemical matter tested against this target is non-selective or counter-selective against the human orthologue. The median selectivity ratio is 1.0x — the typical compound hits human and parasite HDAC8 equally. All three clinically advanced HDAC inhibitors (vorinostat, panobinostat, quisinostat) are severely counter-selective. This target, despite being the most extensively studied in schistosomiasis drug discovery, has a fundamental translational barrier that had not been systematically documented across the full ChEMBL chemical series.

**Finding 2: SmDHODH is the most tractable repurposing target.** 61.1% of compounds with dual-species data show at least 3x selectivity. The median ratio is 4.7x. Three compounds exceed 10x. The best — CHEMBL155771 — combines 23 nM potency, 30.8x selectivity, and QED 0.89 with zero Lipinski violations. This compound has no prior schistosomiasis literature. It is a novel, potent, selective, drug-like SmDHODH inhibitor.

**Finding 3: Atovaquone is the strongest translational candidate.** An approved antimalarial already on the WHO Essential Medicines List, available in sub-Saharan Africa as Malarone, with 6x selectivity over human DHODH and 430 nM parasite potency. It has been proposed in a few publications but not systematically prioritized through a selectivity-first pipeline.

**Finding 4: SmTGR is computationally non-selective.** The most biologically compelling target (single point of failure, unique fusion enzyme) with the most potent hits (10 nM) was assessed computationally through molecular docking against both SmTGR (PDB 2X99) and human TrxR1 (PDB 2ZZC). 88% of compounds showed equivalent or preferential binding to the human enzyme. The most potent compound (CHEMBL3322287, 10 nM) computationally prefers the human target (delta -1.3 kcal/mol). Docking scores correlated significantly with experimental IC50 (Spearman rho = 0.37, p = 0.016), validating the setup. The remaining selectivity opportunity lies in SmTGR's unique TrxR-GR fusion interface — a structural feature absent from human TrxR1 — which was not targeted in this analysis and would require molecular dynamics or co-crystallography to evaluate.

**Finding 5: Multi-pathway evidence integration is necessary.** A pipeline using only target-based IC50 ranked praziquantel last (193/194). Adding structural similarity raised it to 132. Adding whole-organism data raised it to 44. Adding supply chain data raised it to #1 among translational candidates. No single evidence pathway captures the full picture. The pipeline's architecture — not any individual signal — is what produces correct rankings.

**Finding 6: 27 hard negatives reveal structural selectivity determinants.** Counter-selective compounds with high Tanimoto similarity to selective compounds (up to 0.755 similarity) point to specific structural features responsible for selectivity. These compound pairs are starting points for rational design of selective inhibitors.


## LEVEL 16: WHAT IS THIS NOT?

Kira is not a clinically validated repurposing engine. No compound from its shortlist has been tested in a living worm based on Kira's recommendation. The pipeline ranks possibilities; it does not adjudicate translational reality.

Kira is not evidence of novel therapeutic biology. The selectivity finding is a new systematic analysis, but the underlying data is public. The contribution is systematization and integration, not generation of new experimental data.

Kira's benchmark is structured but probably circular. Active compounds have ChEMBL data. Inactive compounds do not. The AUROC measures this association, not general predictive power.

Kira's supply chain data is manually curated and limited. A comprehensive African formulary database does not exist in queryable form.

Every limitation is documented in the pipeline's own output. That honesty is the most important feature. A system that reports 0.999 AUROC without noting that the discrimination task is easy is misleading. A system that reports its limitations alongside its findings is trustworthy.

---

# PART SIX: THE TECHNICAL STACK

## LEVEL 17: WHAT RUNS ON THE MACHINE

**Hardware:** Apple M4 Pro MacBook, 48 GB unified memory. All computation runs on CPU. No GPU was used. No cloud compute was required.

**Python 3.11** in the bio-builder conda environment.

**pandas** — Tabular data manipulation. Every CSV is loaded, filtered, merged, and analyzed with DataFrames.

**RDKit** — Cheminformatics. Morgan fingerprints, Tanimoto similarity, SMILES parsing, Lipinski properties, TPSA, QED, molecular formula computation.

**chembl_webresource_client** — Python API for ChEMBL queries. HTTP requests to ChEMBL servers returning structured JSON.

**scikit-learn** — Evaluation metrics only (AUROC, AUPRC). No models were trained.

**requests** — HTTP library for PubMed E-utilities API.

**NumPy** — Numerical computation. Bootstrap confidence intervals, logarithmic transformations.

**git** — Version control. Every major step committed with descriptive message.

**AutoDock Vina** — Molecular docking engine. Computes binding energies (kcal/mol) for compound-protein interactions by placing a ligand into a receptor's active site and scoring the fit.

**meeko** — Ligand preparation for Vina. Converts RDKit molecules to PDBQT format with AutoDock atom types and torsion trees.

**gemmi** — Macromolecular structure library. Required by meeko for chemical template generation.

**Total pipeline runtime:** ~40 minutes including docking (mostly network queries, cached after first run).

No deep learning. No JAX. No GPU. No cloud. The power is integrative, not computational.


## LEVEL 18: THE DATA ARCHITECTURE

```
~/kira/
├── data/
│   ├── raw/primekg.csv                           ~250 MB knowledge graph
│   ├── processed/
│   │   ├── schistosomiasis_subgraph.csv           52 edges, human context
│   │   ├── schisto_parasite_targets.csv           10 parasite proteins
│   │   ├── schisto_filtered_activities.csv        226 measurements
│   │   ├── schisto_repurposing_candidates.csv     + drug approval
│   │   ├── schisto_whole_organism_raw.csv         494 phenotypic records
│   │   ├── human_orthologue_activities.csv        517 human target records
│   │   ├── kira_selectivity_analysis.csv          91 selectivity ratios
│   │   ├── kira_ranked_candidates_v1-v3.csv       Progressive rankings
│   │   ├── kira_shortlist_v2.csv                  Selectivity-adjusted
│   │   ├── kira_definitive_shortlist_v1.csv       Final with supply chain
│   │   ├── smtgr_docking_results.csv              43 docking selectivity results
│   │   └── kira_novelty_analysis.csv              PubMed checks
│   ├── eval/
│   │   ├── evaluation_set_v1.csv                  194 compounds
│   │   ├── evaluation_set_v2.csv                  228 compounds
│   │   ├── eval_v2_dev.csv / eval_v2_test.csv     80/20 split
│   │   └── eval_v2_inclusion_criteria.txt         Frozen documentation
│   ├── publication/
│   │   ├── table1_target_selectivity.csv          Per-target statistics
│   │   ├── table2_selective_compounds.csv         Top selective hits
│   │   ├── table3_hard_negatives.csv              Structural pairs
│   │   ├── table4_smtgr_docking_selectivity.csv   Docking selectivity
│   │   ├── translational_candidates.csv           Approved drugs
│   │   ├── discovery_candidates.csv               Research compounds
│   │   └── publication_analysis.txt               Full analysis report
│   └── reports/                                   Human-readable reports
├── docs/
│   └── kira-ground-to-god.md                      This document
├── data/docking/
│   ├── receptors/                                  PDB files, cleaned PDBs, PDBQT receptors
│   └── ligands/                                    43 SmTGR compound PDBQT files
├── mission.md / plan.md / agents.md / architecture.md
└── Scripts 01-14
```

---

# PART SEVEN: THE META-LESSONS

## LEVEL 19: WHAT THE PIPELINE TAUGHT THAT TRANSCENDS THE PIPELINE

**Lesson 1: Aggregate metrics lie when case-level behavior fails.** AUROC 0.985 in Script 04 looked excellent. Praziquantel ranked last. The metric and the reality contradicted each other. Always check load-bearing examples, not just aggregate numbers.

**Lesson 2: Every data source has a blind spot.** PrimeKG misses parasite biology. ChEMBL misses phenotypic evidence. PubMed misses compound synonyms. WHO EML misses actual in-country availability. No single source tells the full story. The pipeline's value is in integration, not in any individual query.

**Lesson 3: The hardest test is not "does it separate good from bad?" but "does it rank the known correct answer correctly?"** Praziquantel — the actual standard of care — was the pipeline's hardest case. Its journey from rank 193 to rank 44 to translational #1 is the story of the entire project: each improvement came from understanding why a known-good drug was being ranked incorrectly and adding the right evidence pathway to fix it.

**Lesson 4: Selectivity is not a bonus feature — it is the gatekeeping filter.** Before Script 10, the pipeline rewarded potency. Quisinostat (28 nM, 0 schistosomiasis publications) looked like a breakthrough. After Script 10, quisinostat was revealed as 50x more potent against the human enzyme. Selectivity turned the most exciting compound into the most dangerous one. This applies to all of drug discovery: potency without selectivity is toxicity.

**Lesson 5: Deployment context is a competitive advantage, not a cosmetic addition.** Someone at a Boston lab could replicate Kira's ChEMBL queries. They could not replicate the supply chain scoring because they don't know that praziquantel is donated by Merck, that Malarone is available in East Africa but expensive, that cancer HDAC inhibitors don't exist in African formularies at any price. That contextual knowledge is irreplaceable.

**Lesson 6: Documenting limitations is a feature, not a weakness.** Every Kira report explicitly states what the pipeline can and cannot do. This is not self-deprecation — it is the difference between a trustworthy system and a demo. A supervisor, reviewer, or hiring manager who reads "AUROC 1.0 but we note this reflects an easy discrimination task" thinks: this person understands evaluation. A person who reads "AUROC 1.0" with no qualification thinks: this person doesn't.

**Lesson 7: The iterative process is the product.** The pipeline didn't spring into existence as a 14-script system. It grew one script at a time, each responding to a specific failure or gap in the previous version. Script 04 revealed the praziquantel problem. Script 05 partially fixed it. Script 06 finished fixing it. Script 10 revealed the selectivity problem. Script 11 fixed it. Script 13 added the deployment layer. Script 14 produced a false result (43/43 selective due to wrong coordinates), was diagnosed, fixed, and produced the honest result (88% non-selective). That iterative, error-driven development — visible in the git history — is more impressive than any final output. It demonstrates the ability to diagnose, prioritize, and fix problems systematically.

**Lesson 8: When a result looks too good, it is.** Script 14's first run showed every SmTGR compound as massively parasite-selective (ratios of 500,000x). Instead of celebrating, we diagnosed: the human receptor docking box was pointed at empty space 85 Angstroms from the active site. The corrected run showed 88% non-selectivity — a less exciting but honest and scientifically useful result. The ability to be skeptical of your own best-looking findings is the single most important scientific skill.


## LEVEL 20: WHAT HAPPENS NEXT

The pipeline is a scaffold. The findings are real but preliminary. Four "not yet" statements define the path forward.

**Not yet a moat → Moat.** Build a proprietary NTD ontology that compounds over time: target essentiality scores refined by wet-lab results, selectivity metadata accumulated across multiple parasites, African supply chain constraints encoded as a queryable database, evaluation datasets that grow with each project. The moat is not the code — it is the accumulated, validated, context-specific knowledge that the code operates on.

**Not yet a platform → Platform.** Generalize the selectivity-first triage methodology beyond schistosomiasis to a second NTD (Trypanosoma or Leishmania). If the same approach reveals a similar pattern on a different parasite — most chemical matter non-selective, one target with a favorable selectivity landscape — the methodology is validated as a reusable engine, not a one-off analysis.

**Not yet novel structural biology → Novel structural biology.** The SmTGR docking showed that current chemical matter is non-selective at the GSH/FAD binding sites. The remaining hypothesis is that SmTGR's unique TrxR-GR fusion interface — a structural feature absent from human TrxR1 — could be a selectivity pocket. Testing this requires molecular dynamics simulation or AlphaFold 3 complex prediction to capture the conformational dynamics of the fusion enzyme.

**Not yet publishable → Published.** The preprint is drafted with all 14 scripts' findings including the docking results. What remains: submit to bioRxiv, get the DOI, and send it to the labs that matter.

**Not yet clinically meaningful → Clinically meaningful.** Test CHEMBL155771 and CHEMBL4452960 in whole-worm assays. Verify atovaquone against adult S. mansoni. Even one experimental validation converts the pipeline from a computational analysis into a discovery platform. Reach out to parasitology groups at universities with S. mansoni assay capability.

The preprint is written. The README is ready. The GitHub repository is live at https://github.com/Danny2045/Kira. The selectivity finding across all three major targets is complete — two experimental, one computational.

The pipeline built the road. The selectivity analysis drove somewhere new. The docking confirmed where the road doesn't go (SmTGR GSH site) and pointed toward where it might (the fusion interface). The publication makes it visible. And the person who built it — alone, on a MacBook, in Nashville — is the one who walks through the door it opens.

That is Kira, from ground level to god level.
