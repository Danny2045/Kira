# Kira — Phase 1 Plan

## Goal

Produce a validated, reproducible drug repurposing triage result for
schistosomiasis that can be published or presented.

## Pipeline Architecture

```
[Knowledge Graph] → [Target Identification] → [Candidate Retrieval]
        ↓                                            ↓
[Disease-Target Links]                    [Molecular Docking]
                                                     ↓
                                          [ADMET Filtering]
                                                     ↓
                                          [Supply Chain Check]
                                                     ↓
                                          [Composite Ranking]
                                                     ↓
                                          [Retrospective Eval]
```

## Milestones

### M1: Data Foundation (Week 1-2)
- Download and load PrimeKG knowledge graph
- Query schistosomiasis subgraph: disease → targets → drugs → pathways
- Identify characterized drug targets with available 3D structures
- Build ground-truth evaluation set from literature

### M2: Candidate Generation (Week 3-4)
- Retrieve approved drugs connected to identified targets
- Expand candidates via pathway and mechanism similarity
- Implement basic molecular docking against priority targets
- Score and rank initial candidate list

### M3: Filtering Pipeline (Week 5-6)
- ADMET prediction for top candidates
- Supply chain feasibility check (WHO Essential Medicines List, availability)
- Composite scoring function combining all signals

### M4: Evaluation and Validation (Week 7-8)
- Test pipeline against known outcomes:
  - Does praziquantel rank highly? (positive control)
  - Does oxamniquine appear? (known active, older drug)
  - Do known-inactive compounds rank low? (negative controls)
- Quantify precision, recall, ranking quality
- Document methodology for publication

## Acceptance Criteria

- Pipeline runs end-to-end from PrimeKG query to ranked candidate list
- Retrospective evaluation shows pipeline correctly separates known actives
  from known inactives above chance level
- Every candidate in output carries: target rationale, docking score,
  ADMET summary, availability status, confidence level
- All results are reproducible from a single command

## Risk Boundaries

- Do not claim clinical validity — this is computational triage
- Do not skip evaluation — a ranked list without validation is a demo, not science
- Do not optimize for pretty outputs before the pipeline is correct
