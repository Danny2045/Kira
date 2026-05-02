# Regeneration / CRISPR contrast scout

This branch adds a narrow scout adapter for regenerative biology, CRISPR perturbation, partial reprogramming, and Michael Levin-style bioelectric pattern-control examples.

It is not a broad regenerative-medicine platform claim. It only turns a small set of frontier systems into Kira's contrast-core shape:

```text
intervention
→ desired biological context
→ control / safety / failure context
→ measurable readout
→ evidence status
→ missing measurement
→ experiment ticket
```

## Scope

The scout covers six seed contrasts:

1. Bioelectric membrane-voltage / proton-pump modulation in a Xenopus tail-regeneration setting.
2. Planarian target-morphology / gap-junction / bioelectric-gradient perturbation.
3. Human forebrain assembloid CRISPR perturbation mapped to interneuron generation or migration phenotypes.
4. CRISPR knockout / CRISPRi / CRISPRa / single-cell perturbation in primary human 3D organoids.
5. Yamanaka-factor partial reprogramming as a bounded safety-window contrast.
6. Morphogen or pathway-control organoid morphology phenotyping as a proposed future assay shape.

## Non-claims

The scout does not claim that Kira has discovered drugs, validated a wet-lab assay, solved regeneration, solved CRISPR, produced a clinical therapeutic, or proven a new model-performance result.

## Evidence-status discipline

The first branch intentionally includes multiple statuses:

| Status | Meaning in this scout |
|---|---|
| `paired` | Desired and control/failure sides are represented well enough to form a first contrast. |
| `bounded` | The contrast is useful only when dose, timing, or morphology boundaries are explicit. |
| `single-side-only` | A useful perturbation side exists, but Kira should demand a missing safety/control side. |
| `conflicting` | The domain has promising and failure/safety evidence that must be paired before stronger claims. |
| `proposed` | The record is only an executable experiment-ticket shape. |

## Code surface

New package:

```text
src/kira/regeneration/__init__.py
src/kira/regeneration/scout.py
tests/test_regeneration_scout.py
docs/REGENERATION_CONTRAST_SCOUT.md
```

Public helpers:

```python
from kira.regeneration import (
    build_contrast_specs,
    build_evidence_records,
    build_experiment_tickets,
    build_regeneration_scout,
    regeneration_scout_seeds,
    summarize_evidence_statuses,
    tickets_as_dicts,
)
```

## Expected test commands

```bash
pytest tests/test_regeneration_scout.py
pytest tests/test_contrast_schemas.py tests/test_contrast_tickets.py tests/test_regeneration_scout.py
ruff check src/kira/regeneration tests/test_regeneration_scout.py
```

## Source anchors

These source anchors are used only to motivate seed shapes, not to assert therapeutic validation:

- Michael Levin, *Bioelectric signaling: Reprogrammable circuits underlying embryogenesis, regeneration, and cancer*, Cell, 2021. DOI: `10.1016/j.cell.2021.02.034`.
- Adams, Masi, and Levin, H+ pump-dependent membrane-voltage changes in Xenopus tail regeneration, Development, 2007. DOI: `10.1242/dev.02812`.
- Emmons-Bell et al., gap-junction blockade and planarian head-anatomy outcomes, International Journal of Molecular Sciences, 2015. DOI: `10.3390/ijms161126061`.
- Meng et al., assembloid CRISPR screens in human neurodevelopment, Nature, 2023. DOI: `10.1038/s41586-023-06564-w`.
- Lo et al., large-scale CRISPR screening in primary human 3D gastric organoids, Nature Communications, 2025. DOI: `10.1038/s41467-025-62818-3`.
- Yücel and colleagues, partial reprogramming-induced rejuvenation review, Nature Communications, 2024. DOI: `10.1038/s41467-024-46020-5`.
