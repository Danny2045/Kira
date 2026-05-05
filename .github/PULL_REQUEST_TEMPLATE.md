# Summary

Describe the branch in one or two paragraphs. Include the PR title and the
scientific purpose in plain language.

## Changed files

List the exact files changed and why each changed.

- `path/to/file.py`:
- `docs/path.md`:
- `tests/path_test.py`:

## Scientific purpose

State the scientific or process problem this PR addresses.

Examples:

- data substrate repair
- benchmark repair
- contrast definition
- evidence audit
- experiment or data ticket generation
- returned-data audit
- claim-boundary documentation

## Validation

List exact commands run and their outcomes. Include focused tests when the PR
touches a specific module or workflow.

Required full-project checks unless explicitly not applicable:

```bash
ruff check .
pytest -q
```

Focused checks:

```bash
# command:
# result:
```

## Non-claims

This PR does not claim that:

- Kira discovered drugs
- Kira solved AMR
- Kira solved regeneration
- Kira gives clinical recommendations
- Kira gives prescribing advice
- Kira wet-lab validated tickets
- Kira proved a new model-performance result unless explicitly supported by
  results

Add any PR-specific non-claims here:

- 

## Scientific risk

Describe the main scientific, data, benchmark, measurement, or claim-boundary
risk that remains after this PR.

## Reproducibility / Carmack check

Confirm:

- exact commands are listed
- generated artifacts are either committed intentionally or excluded
- local-only state is not required to reproduce the result
- diffs are limited to the stated purpose
- validation does not depend on untracked private files

## Claim-boundary check

Confirm:

- conclusions are bounded by the actual data and tests in this PR
- no unsupported clinical, therapeutic, public-health, benchmark, or wet-lab
  claims were added
- institution, scientist, or partner names are not used to imply endorsement,
  partnership, or equivalence

## Next step

State the single most useful next branch, measurement, ticket, audit, or
benchmark repair action.
