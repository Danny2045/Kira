"""Rwanda AMR contrast scout adapter."""

from kira.amr.scout import (
    AMR_SCOUT_VERSION,
    AmrScoutSeed,
    amr_scout_seeds,
    build_amr_scout,
    build_contrast_specs,
    build_evidence_records,
    build_experiment_tickets,
    summarize_evidence_statuses,
    tickets_as_dicts,
    validate_amr_seed,
)

__all__ = [
    "AMR_SCOUT_VERSION",
    "AmrScoutSeed",
    "build_amr_scout",
    "build_contrast_specs",
    "build_evidence_records",
    "build_experiment_tickets",
    "amr_scout_seeds",
    "summarize_evidence_statuses",
    "tickets_as_dicts",
    "validate_amr_seed",
]
