# The Emerging Stack of Programmable Medicine
## A Technical Manual from Physics to Clinic

*Daniel Ngabonziza — March 2026*

---

## One Synthesis

The deepest unifying idea is this: biology is not just chemistry, and not just genetics. It is a hierarchy of information-processing and control systems spanning atoms, cells, tissues, and bodies.

Michael Levin's work supplies the conceptual language for that claim. The Xenobot memory preprint gives a mesoscopic experimental probe of it. AlphaFold 3 gives an atomic-scale inference engine for biomolecular structure. IsoDDE tries to push that atomic engine toward true drug-design usefulness in novel chemical space. Read together, these are not disconnected curiosities. They are pieces of one emerging science: programmable living matter.

A useful way to orient everything is to separate three layers that are often confused. First is the philosophical and conceptual layer: what counts as agency, memory, goal-directedness, and intelligence in nonstandard substrates. Second is the mechanistic physiological layer: how cell collectives actually sense, coordinate, store state, and produce behavior without neurons. Third is the molecular and design layer: how proteins, ligands, antibodies, nucleic acids, and pockets interact at atomic resolution, and how one predicts or engineers those interactions. Levin's work emphasizes the first two layers. The Xenobot paper experimentally probes the second. AF3 and IsoDDE dominate the third. The real frontier is coupling all three.

The right mental model is not "AI helps medicine." It is medicine becoming a cyber-physical production system. The old world was serial and artisanal: hypothesis, experiment, formulation, trial, paperwork, reimbursement. The new world is increasingly parallel and machine-mediated: foundation models generate hypotheses, cloud labs execute experiments, formulation systems tune biodistribution, virtual-cell systems reduce blind search, clinical copilots compress information friction, and real-world deployment feeds back into the next design cycle.

What follows is the full stack, from physics to clinic.

---

# I. Layer 0 — The Invariants: Why Medicine Has Been Hard

The base truth is that medicine is not one optimization problem. It is a cascade of inverse problems under uncertainty. You are trying to infer causal mechanisms from noisy biological systems, then design interventions, then deliver them into the right spatiotemporal compartments, then prove benefit and acceptable risk in humans, then scale and pay for them.

Eroom's law gives the historical macro-view: the number of new drugs approved per billion dollars of R&D spend has halved roughly every nine years since 1950. A precise way to think about this:

```
T_total = T_target + T_design + T_preclinical + T_delivery + T_clinical + T_manufacturing/regulatory
```

and

```
P_launch = ∏(p_i) for i = 1 to n
```

where each p_i is the conditional probability of surviving one stage. AI can reduce several T_i and improve some p_i, but if one stage remains hard — delivery, recruitment, long-term safety, CMC, reimbursement — the whole product still stalls. This is why naive "AI cures everything in 10 years" narratives fail: they assume the front-end search problem dominates the entire pipeline. In medicine, it often does not.

Rare adverse-event statistics show the irreducibility of some clinical burdens. If an adverse event has true incidence p, then the probability of seeing at least one event in n patients is:

```
P(≥1 event) = 1 − (1 − p)^n
```

For p = 10⁻³ and 95% detection probability, you need roughly n ≈ 3,000. For p = 10⁻⁴, you need roughly n ≈ 30,000. That is before subgroup heterogeneity, adherence failures, delayed toxicity, and endpoint ambiguity. FDA guidance for gene therapies exists precisely because some risks are delayed and durable, with long-term follow-up required for as long as 15 years in some cases.

So the first principle of this whole stack is: faster design is necessary but not sufficient. The stack only bends the curve when every downstream bottleneck is also being attacked.


# II. Layer 1 — Biological Information Theory: Sequence, Structure, Function, Programmability

The foundational axiom: biological information is encoded in sequence, structure emerges from physical interactions, and function emerges from structure. Synthetic biology adds standardization, modularity, and abstraction so that living systems become engineerable information-processing substrates.

At the molecular level, the key chain is:

```
Sequence → Conformational ensemble → Binding / catalysis / assembly → Cellular phenotype
```

But the modern update is that this chain is not deterministic in the trivial sense. The true mapping is from sequence to a distribution over conformational and interaction states, modulated by solvent, post-translational modifications, cofactors, membranes, crowding, and partner molecules.

Structure AI therefore matters because it approximates a previously intractable map from symbolic biological information to spatially grounded biochemical possibility. This is also why synthetic biology and structural biology are converging. Structural biology gives you constraints and mechanistic priors; synthetic biology gives you controllable substrates and design freedom.


# III. Layer 2 — Design Intelligence: AlphaFold 3 and Interaction-Space Inference

## What AF3 actually is

AlphaFold 3 is not "AlphaFold 2 but larger." The Nature paper describes a substantially updated diffusion-based architecture that predicts the joint structure of complexes including proteins, nucleic acids, ligands, ions, and modified residues. It explicitly de-emphasizes heavy MSA processing, replaces AF2's evoformer dominance with a pairformer trunk, and replaces the AF2 structure module with a diffusion module that predicts raw atom coordinates directly.

The key architectural elements: input embedder, template module, 4-block MSA module, 48-block pairformer, diffusion module (3 + 24 + 3 blocks), recycling, and a confidence module. The MSA representation is no longer the central carrier of state; information passes predominantly through pair and single representations.

AF3 directly predicts raw atom coordinates rather than operating on amino-acid-specific frames and side-chain torsions. The diffusion process lets the network improve local structure at low noise and global structure at higher noise. That means AF3 is a general conditional generative model over heterogeneous biomolecular complexes, not a protein-only fold engine. The input can be "protein + nucleic acid + ligand + ion + modified residue + bonding constraints," and the output is a sampled coordinate realization plus confidence estimates. That is a category change.

## Why diffusion matters

Three advantages. First, a common output language: atomic coordinates, removing special-case machinery for mixed molecular systems. Second, denoising across noise scales means the model learns both local stereochemical regularities and large-scale interface geometry. Third, prediction becomes distributional, so multiple samples and seeds are meaningful ways to explore uncertainty.

Diffusion does not magically solve physics. AF3 is a learned structural prior, not a full molecular dynamics engine. It produces plausible structures, not thermodynamic ensembles. The paper is candid about limitations in stereochemistry, disordered regions, and state selection.

## Clinical consequence

AF3 compresses the front-end of therapeutic reasoning: ligand pose hypotheses, antibody-antigen interfaces, nucleic-acid complexes, covalent modifications, and mixed biomolecular assemblies become accessible to machine-guided exploration. AF3 is the first strong expression of a new design-layer principle: therapeutic search begins with an interaction-space prior, not with blind combinatorics.


# IV. Layer 3 — Beyond AF3: IsoDDE and Drug-Design Engines

The Isomorphic Labs technical report names AF3's residual problem directly: accurate structure prediction is not yet real-world drug design. Despite AF3-class advances, limitations persist in generalizing to unexplored molecular space, estimating binding affinity, and detecting binding sites on previously uncharacterized protein surfaces.

## What IsoDDE claims

Four capabilities: protein-ligand structure prediction, antibody-antigen structure prediction, binding affinity prediction, and pocket identification. The report says IsoDDE more than doubles AF3 accuracy on a hard protein-ligand generalization benchmark, improves antibody-antigen interface prediction, exceeds gold-standard physics-based methods on certain affinity benchmarks, and outperforms P2Rank for pocket identification.

This is the move from structure engine to decision engine. AF3: here is a plausible 3D arrangement. IsoDDE: here is a more experimentally useful simulation layer for discovery, especially in out-of-distribution chemistry.

## Why out-of-distribution is the real game

Medicinal chemistry rarely fails on easy, in-distribution cases. It fails on novel pockets, induced-fit systems, unseen chemotypes, cryptic binding sites, and dynamic interfaces. The report defines "globally cryptic" pockets via apo-holo RMSD criteria and benchmarks against P2Rank — much more discovery-relevant than simple pose refinement in a known site.

## Honest caution

AF3 gives architecture, training logic, and limitations in a peer-reviewed paper. IsoDDE gives benchmark claims in a company report without the same mechanistic disclosure. The design layer is clearly moving toward integrated drug-design engines, but the most ambitious post-AF3 claims deserve independent stress testing.


# V. Layer 4 — Autonomous Experimentation: Closed-Loop Lab Control

## What the Ginkgo/OpenAI system did

GPT-5-driven autonomous laboratory optimizing cell-free protein synthesis in Ginkgo's cloud lab. Across six optimization steps: 480 384-well plate designs, 29,527 unique reaction compositions, six months. Reduced cost from $698/g to $422/g and increased titer by 27%. Only ~1% of plate designs were fundamentally flawed.

GPT-5 produced designs under schema-validated interface, executed physically in a real laboratory, with measurements feeding the next round.

## The control-theoretic picture

```
x_{t+1} = Update(x_t, a_t, y_t)
a_t = argmax_{a ∈ A} E[U(a; x_t)]
```

The system turns the reasoning model into a bounded experimental policy generator. Autonomous driving became real not when a model could describe roads, but when perception, planning, validation, and actuation were tightly coupled. The same is now happening in biology.

## Why this changes medicine

Structure models emit thousands of candidate designs. Without autonomous validation, that creates a downstream bottleneck. The cloud-lab layer absorbs that overflow. The rate limit moves from human experimental labor to assay kinetics, reagent logistics, and model quality.


# VI. Layer 5 — Delivery Physics: Why Most Payloads Die Here

The chokepoint that medicine people respect and software people underestimate.

## The problem

Delivery must solve circulation, biodistribution, tissue penetration, cell entry, endosomal escape, intracellular localization, dose, immunogenicity, and manufacturability simultaneously.

## Why LNPs love the liver

Standard LNPs exhibit strong liver tropism because serum ApoE adsorbs onto the particle and mediates hepatocyte uptake via LDL receptors. The problem of extrahepatic delivery: how do you break or re-route that default biodistribution?

## Delivery as optimization

```
P_effective = P_circulation · P_tissue_access · P_cell_binding · P_internalization · P_endosomal_escape · P_intracellular_action

TI = desired_target_exposure / (off_target_exposure + toxicity)
```

Every design parameter — peptide, lipid, polymer, PEG density, particle size, charge, route, dosing — perturbs these terms. The design space is combinatorial, coupled, and experimentally measured. This is exactly why delivery is now an AI + autonomous-lab problem.


# VII. Layer 6 — Payloads: What the Stack Delivers

## CRISPR: from break to write

CRISPR-Cas9 made programmable DNA targeting practical. Base editors allowed single-base changes without double-strand breaks. Prime editors expanded to search-and-replace. The most dramatic proof of concept: a personalized in vivo CRISPR therapy for an infant developed in six months.

## Platformizing rare-disease cures

The Center for Pediatric CRISPR Cures: design, preclinical safety, manufacturing, clinical, and regulatory infrastructure under one coordinated program, initially aiming to treat eight patients.

## Epigenetic reprogramming as state reset

ER-100: controlled expression of OCT4, SOX2, KLF4. FDA clearance in January 2026 for human trials in optic neuropathies. First cellular rejuvenation therapy using epigenetic reprogramming to reach human testing. Classical pharmacology: "modulate a target." Reprogramming: "change cell state." The payload layer becomes state-transition operators.

## mRNA and transient expression

Transient-expression payloads are ideal for delivering reprogramming factors, genome editors, immune instructions, or proteins without permanent integration. Nonviral delivery and LNPs are now central to the gene-therapy future.


# VIII. Layer 7 — Virtual Cells and Perturbation Atlases

Autonomous labs need better priors than brute-force search.

Arc's Virtual Cell Atlas: observational and perturbational datasets from over 600 million cells. First-generation model State predicts how stem cells, cancer cells, and immune cells respond to perturbations, trained on ~170 million observational and ~100 million perturbational cells across 70 cell lines.

Virtual-cell models approximate the transition kernel:

```
P(x_{t+1} | x_t, u_t)
```

where x_t is cell state and u_t is intervention. Better approximations mean better experimental triage and fewer wasted assays.

The future preclinical sequence:

```
design model → virtual-cell screening → autonomous wet-lab validation → delivery optimization → animal/human evidence
```

Not replacement for wet lab. A more data-efficient acquisition strategy for wet lab.


# IX. Layer 8 — Clinical Intelligence

The real deployment problem is clinical information orchestration under uncertainty.

## Evaluation

HealthBench: 262 physicians, 60 countries, 5,000 conversations, 48,562 rubric criteria. Tests uncertainty handling, context gathering, communication, and aligned reasoning — not multiple-choice medical trivia.

## Clinical function

A model becomes clinically valuable when it retrieves evidence, aligns to institutional pathways, and grounds outputs in patient-specific context simultaneously. That is a workflow coupling problem.

## The deep role

Clinical AI reduces the human information tax between scientific possibility and actual care. Even a perfect therapy underperforms if the patient is routed late, the chart is fragmented, documentation burden is high, prior auth is delayed, or follow-up is missed. Clinical intelligence is the deployment layer of the medical stack.


# X. Layer 9 — Neurotechnology and Read/Write Interfaces

Medicine shifting from passive observation to closed-loop interfaces.

PRIMA: subretinal photovoltaic implant restoring central form vision in geographic atrophy patients. Medicine becoming a sensor-decoder-stimulator system.

Intracranial BCIs with implantable arrays, wireless telemetry, and decoders: existence proof that some parts of medicine are moving from pharmacology toward real-time read/write control of physiological systems.


# XI. Layer 10 — Levin, Xenobots, and the Frontier Beyond Molecular Medicine

The deepest frontier. What comes after molecule-centric medicine.

## The Xenobot contribution

March 2026 preprint: basal Xenobots (no scaffolds, no synthetic circuits, no genetic editing) retain distinct, stimulus-specific, long-term signatures in behavior, transcription, and calcium physiology after brief chemical exposure. A synthetic non-neural cellular collective storing at least two distinct persistent internal state changes.

This expands the architecture of medicine. Most of the stack assumes the target is a molecule, receptor, pathway, or cell state. Levin asks whether tissues themselves are distributed information-processing systems whose bioelectric, calcium, and morphogenetic states can be measured, modeled, and steered.

## Medical implications

If non-neural collectives can store state, discriminate stimuli, and coordinate behavior, the endpoint of medicine may not be "deliver the right molecule" but "communicate with and reprogram multicellular control systems." That is morphophysiological control.

The current stack:

```
molecular design → delivery → cellular perturbation → clinical deployment
```

The future extension:

```
molecular design → cell-state control → tissue-state control → anatomical / regenerative control
```

Beyond today's clinical reality. Not science fiction. A programmatic expansion of what counts as a controllable medical variable.


# XII. Worked End-to-End Example

Extrahepatic in vivo editing for a monogenic pediatric disease:

**Step 1 — Target/mechanism.** Genetics + transcriptomics + perturbation atlases + structure models define the optimal edit class. AF3/IsoDDE map constraints. Virtual-cell systems predict state transitions.

**Step 2 — Payload.** Cas nuclease, base editor, prime editor, CRISPRa/i, reprogramming cassette, siRNA, or mRNA. Depends on edit class, durability, risk tolerance, manufacturability.

**Step 3 — Formulation/targeting.** Optimize ionizable lipid, helper lipid, cholesterol, PEG-lipid, N:P ratio, peptide ligand, route, dosing. Autonomous labs indispensable — design space too large for artisanal search.

**Step 4 — Autonomous preclinical loop.** Robotic screening: on-target editing, viability, cytokine induction, delivery efficiency, endosomal escape, manufacturability. Design → assay → update belief → redesign at industrial cadence.

**Step 5 — Animal/translational validation.** Biodistribution, PK, immunogenicity, germline risk, tumor risk, persistence, histopathology. No AI stack abolishes this.

**Step 6 — Clinical protocol.** Trial eligibility, phenotyping, consent, AE workflows, monitoring, documentation. Clinical intelligence compresses information friction.

**Step 7 — Long-term safety flywheel.** Post-treatment evidence flows back into model retraining, delivery redesign, patient-selection logic. The stack becomes self-improving.


# XIII. What Remains Hard

**Causality under heterogeneity.** Mixed etiologies, comorbid states, history-dependent physiology. Better priors do not eliminate perturbation-grounded evidence.

**Delivery external validity.** Protein coronas change across species and patients. Receptor expression changes with inflammation, fibrosis, tumor context, age.

**Long-term safety and economics.** Scientifically elegant, commercially fragile. Faster science does not solve reimbursement.

**Governance of autonomous systems.** Lab: structured control and validation. Clinic: auditability, escalation, liability, calibration.

**Multicellular control theory.** The stack focuses on molecules, cells, workflows. The next frontier: regenerative setpoints, bioelectric state, non-neural memory, tissue-level decision-making. No mature medical stack for that layer yet. But we are starting to see it.


# XIV. Final Synthesis

Medicine is being rebuilt as a multi-layer control stack.

Structural and synthetic biology provide the information-theoretic and physical substrate. AF3-class models turn biomolecular interaction space into a tractable search problem. IsoDDE-class systems push toward drug-design inference in novel chemical space. Autonomous cloud labs convert model proposals into real data at machine speed. Delivery engineering determines whether payloads reach the right tissue. Modern payloads — CRISPR, base/prime editing, mRNA, reprogramming constructs — make root-cause and state-reset interventions plausible. Virtual-cell systems reduce preclinical blind search. Clinical intelligence compresses the information tax that dilutes scientific progress at the point of care.

And Levin/Xenobot-style work hints that the eventual endpoint is not only precise molecular intervention, but direct control of multicellular physiological and morphogenetic state.

The stack is not complete. The delivery layer is brutally empirical. Clinical evidence requires time, patients, and money no model can simulate away. Multicellular control theory is in its infancy. But the trajectory is unmistakable: medicine is becoming a cyber-physical production system, and the people who understand all the layers — from atomic physics to tissue control to clinical deployment — are the ones who will build what comes next.

---

*This document is a companion to the Kira drug repurposing pipeline (github.com/Danny2045/Kira), which operates at Layers 2-3 of this stack: molecular design intelligence applied to neglected tropical disease pharmacology, constrained by cross-species selectivity and East African deployment reality. The selectivity-first methodology Kira demonstrates is one instantiation of the broader principle: therapeutic search must be constrained by translational physics at every layer, not just target potency.*
