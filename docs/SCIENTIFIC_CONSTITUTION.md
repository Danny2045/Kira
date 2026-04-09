# Kira Scientific Constitution

## Formal Claim Boundary

### What Kira Is Allowed To Claim Right Now

Kira is allowed to claim that it is:

- a computational repository for retrospective selectivity analysis across neglected tropical disease targets
- a repository of curated parasite-human comparator mappings, target-pair pocket definitions, and translational ranking heuristics
- a proof-of-concept benchmark showing that curated pocket-divergence descriptors can outperform an ESM-2 baseline on a small leave-one-disease-out target-pair feature benchmark
- a mechanistic hypothesis generator that highlights residue-level pocket differences and residue-energy heuristics associated with observed selectivity
- an approximate structure-checking toolkit that reports steric clashes, Lennard-Jones summaries, and backbone-dihedral counts on protein structures
- an auditable scientific record of historical analyses, including failed approaches, revised benchmarks, and archived docking studies

### What Kira Is Not Allowed To Claim Right Now

Kira is not allowed to claim that it is:

- a causal discovery system in the scientific sense
- a validated drug discovery engine
- a validated force field or validated physics engine
- a compound-conditioned selectivity predictor
- a binding-affinity predictor
- a protein-ligand binding explanation engine
- a full docking engine in the active package
- an eight-check active validation engine
- a system with active mmCIF parsing support
- a system that explains why a specific compound is selective through explicit ligand-conditioned energetic modeling
- a source of experimentally validated causal or mechanistic evidence

### Canonical Present-Tense Claim

Kira is a computational repository for retrospective translational analysis, pair-level selectivity benchmarking, curated biological comparison, and residue-level mechanistic hypothesis generation.

## Evidence-Tier Policy

### Tier Definitions

- `Tier 1`: descriptive / curation / representation
- `Tier 2`: retrospective translational analysis
- `Tier 3`: benchmark correlation / predictive proof-of-concept
- `Tier 4`: mechanistic hypothesis generation
- `Tier 5`: validated causal or experimental evidence

### Repository Rule

No current module in this repository qualifies for `Tier 5`.

Every public-facing:

- README claim
- technical or narrative documentation claim
- docstring
- CLI command help text
- experiment description

must declare its evidence tier in plain language.

### Minimum Declaration Template

- `Evidence tier: Tier X`
- `Scientific status: descriptive / retrospective / benchmark proof-of-concept / mechanistic hypothesis / validated experimental`
- `Not a claim of: compound-level prediction / binding affinity / causal proof / experimental validation`

### Constitutional Constraint

No public-facing text may describe a module above its assigned evidence tier.

## Controlled Scientific Vocabulary

### Forbidden Overclaim Terms

The following terms are forbidden unless directly backed by code and evidence at the necessary tier:

- causal discovery
- causality module
- explains why
- explanation of selectivity
- validated engine
- validated force field
- validated physics engine
- Trust Report
- eight checks
- composite trust score
- binding energy differential
- binding explanation
- full pipeline
- end-to-end pipeline
- full causality pipeline
- docking engine
- selectivity predictor
- compound-level prediction
- predicts selectivity
- validated on real proteins
- mmCIF support
- mechanistic proof
- genuine non-selectivity
- captures the mechanism
- identifies what drives binding

### Approved Replacement Terms

- causal discovery -> `retrospective association analysis` or `mechanistic hypothesis generation`
- causality module -> `mechanistic hypothesis module`
- explains why -> `proposes a residue-level rationale for`
- explanation of selectivity -> `heuristic selectivity attribution`
- validated engine -> `benchmark-tested proof-of-concept`
- validated force field -> `approximate Lennard-Jones model`
- validated physics engine -> `approximate structure-checking toolkit`
- Trust Report -> `heuristic structure-quality report`
- eight checks -> `documented check families` or `planned check families`
- composite trust score -> `heuristic score`
- binding energy differential -> `residue-energy differential` or `docking-score difference`
- binding explanation -> `binding-site hypothesis`
- full pipeline -> `curated workflow` or `partial analytical workflow`
- end-to-end pipeline -> `case-study workflow`
- full causality pipeline -> `curated hypothesis workflow`
- docking engine -> `archived docking workflow`
- selectivity predictor -> `selectivity benchmark model`
- compound-level prediction -> `pair-level feature benchmark`
- predicts selectivity -> `shows benchmark association with selectivity labels`
- validated on real proteins -> `tested for plausibility on selected real proteins`
- mmCIF support -> `PDB-format parsing`
- mechanistic proof -> `mechanistic hypothesis`
- genuine non-selectivity -> `suggestive evidence of non-selectivity`
- captures the mechanism -> `captures one proposed pocket-level signal`
- identifies what drives binding -> `highlights residues associated with the heuristic score`

## Allowed Claims / Forbidden Claims

### Allowed Claims

- Kira performs retrospective selectivity analysis on public and curated data.
- Kira benchmarks target-pair-level pocket descriptors against historical selectivity labels.
- Kira generates residue-level mechanistic hypotheses from curated pockets and approximate residue-energy heuristics.
- Kira performs approximate structural quality checks using currently implemented geometry and steric routines.
- Kira preserves an auditable historical record of archived translational and docking analyses.

### Forbidden Claims

- Kira discovers causal mechanisms.
- Kira validates drug candidates experimentally.
- Kira predicts compound-conditioned selectivity from explicit protein-ligand physics.
- Kira computes true binding affinities or binding free energies.
- Kira actively implements a complete docking engine in the live package.
- Kira currently performs a fully implemented eight-check validation workflow.
- Kira parses mmCIF files in the active parser.

## Closing Principle

Kira must always separate:

- what was curated
- what was benchmarked
- what was hypothesized
- what was experimentally shown

The repository is scientifically trustworthy only if every claim stays at its real evidence tier.
