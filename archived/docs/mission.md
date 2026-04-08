# Kira — Mission

## The Problem

Schistosomiasis infects approximately 200 million people, overwhelmingly in
sub-Saharan Africa. The global treatment strategy depends almost entirely on a
single drug: praziquantel. The WHO has flagged this single-drug dependency as a
resistance vulnerability. If praziquantel resistance emerges at scale, there is
no adequate backup.

New drug development takes 10–15 years and over $1 billion. Pharmaceutical
companies do not invest in diseases that primarily affect populations who cannot
pay. This is the market failure.

## What Kira Does

Kira is a computational drug repurposing engine. It systematically searches the
space of existing approved drugs to find candidates that may also be effective
against schistosomiasis.

The logic: approved drugs have already passed safety testing. Their
pharmacokinetic profiles are known. If one of them also inhibits a critical
schistosome protein, the path to clinical use is dramatically shorter and
cheaper than de novo drug design.

## What "Good" Means

A good output from Kira is a ranked shortlist of repurposing candidates where:

- Each candidate has a mechanistic rationale (which target, what binding mode)
- Each candidate has an ADMET profile compatible with the target population
- Each candidate is assessed for real-world availability in East African supply chains
- Confidence levels and provenance are attached to every claim
- The ranking is validated against known outcomes (retrospective evaluation)

## What Is Out of Scope

- Kira does not replace wet-lab validation or clinical trials
- Kira does not generate novel molecular structures (that is de novo design)
- Kira does not make clinical recommendations
- Kira is a triage and prioritization engine, not an oracle

## Target Disease (Phase 1)

Schistosomiasis (Schistosoma mansoni as primary species)

## Key Targets Under Investigation

- SmTGR (thioredoxin glutathione reductase) — essential for parasite redox defense
- SmHDAC8 (histone deacetylase 8) — involved in parasite gene regulation
- SmCA (carbonic anhydrase) — pH regulation
- Additional targets to be identified through knowledge graph exploration

## Deployment Context

- Resource-constrained clinical environments in East Africa
- Intermittent connectivity
- Integration with existing health information infrastructure
- Outputs must be interpretable by researchers and clinicians, not just ML engineers
