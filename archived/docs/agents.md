# Kira — Coding Rules (agents.md)

## Repository Map

```
kira/
├── mission.md              # Why this project exists
├── plan.md                 # Phase 1 milestones and acceptance criteria
├── agents.md               # THIS FILE — coding rules and repo conventions
├── architecture.md         # Technical design decisions
├── data/
│   ├── raw/                # Downloaded datasets (NOT in git)
│   ├── processed/          # Pipeline outputs (NOT in git)
│   └── eval/               # Ground-truth evaluation sets
├── src/kira/               # Main Python package
│   ├── data/               # Data loading and processing
│   ├── graph/              # Knowledge graph queries
│   ├── dock/               # Molecular docking
│   ├── filter/             # ADMET and supply chain filtering
│   └── eval/               # Evaluation and scoring
├── tests/                  # All tests
├── notebooks/              # Jupyter exploration (not production code)
├── harness/                # Linters, validators, constraints
└── docs/                   # Ontology definitions, decision records
```

## Rules

1. Every function has a docstring explaining WHAT it does and WHY.
2. No magic numbers — all thresholds and parameters are named constants.
3. Data files in data/raw/ are never modified. Processing creates new files
   in data/processed/.
4. Every pipeline step logs: what input it received, what it produced, how
   long it took.
5. Type hints on all function signatures.
6. Tests exist for every module in src/kira/.

## Commands

- Run all tests: `pytest tests/ -v`
- Run the exploration script: `python 01_explore_primekg.py`
- Start Jupyter: `jupyter lab` (from the kira/ directory)

## What Not To Do

- Do not put analysis code in src/. Exploration goes in notebooks/.
  Only stable, tested code goes in src/kira/.
- Do not commit data/raw/ files to git (they are too large).
- Do not skip writing tests because "it's just exploration."
  If it influences the pipeline, it needs a test.
