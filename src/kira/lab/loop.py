"""Closed-loop optimization for the design-build-test-learn cycle.

Implements the iterative experimental cycle:
  1. DESIGN — generate candidate compounds/formulations
  2. BUILD — specify synthesis/procurement
  3. TEST — design and execute experiments
  4. LEARN — analyze results, update priors, generate next round

This is the computational backbone of autonomous experimentation.
The module does not execute experiments directly — it generates
structured experimental plans and processes results to inform
the next iteration.

References:
    - Ginkgo/OpenAI cloud lab (2024) — autonomous experiment design
    - Recursion OS LOWE — workflow orchestration
    - Generate Bio — generate-build-measure-learn loop
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum


class CyclePhase(Enum):
    """Current phase in the DBTL cycle."""

    DESIGN = "design"
    BUILD = "build"
    TEST = "test"
    LEARN = "learn"
    COMPLETE = "complete"


class CompoundStatus(Enum):
    """Status of a compound in the pipeline."""

    DESIGNED = "designed"
    PROCURED = "procured"
    TESTED = "tested"
    HIT = "hit"
    LEAD = "lead"
    ELIMINATED = "eliminated"


@dataclass
class ExperimentResult:
    """Result from a single experiment.

    Attributes
    ----------
    compound_id : str
        Compound identifier.
    target_name : str
        Target tested against.
    assay_type : str
        Type of assay.
    ic50_nm : float | None
        IC50 in nM if determined.
    max_inhibition_pct : float
        Maximum inhibition observed (%).
    hill_coefficient : float
        Hill coefficient from curve fit.
    z_prime : float
        Z' factor for plate quality.
    is_active : bool
        Whether compound meets activity threshold.
    raw_data_path : str
        Path to raw data file.
    notes : str
        Analysis notes.
    """

    compound_id: str
    target_name: str
    assay_type: str = ""
    ic50_nm: float | None = None
    max_inhibition_pct: float = 0.0
    hill_coefficient: float = 1.0
    z_prime: float = 0.0
    is_active: bool = False
    raw_data_path: str = ""
    notes: str = ""

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "compound_id": self.compound_id,
            "target_name": self.target_name,
            "assay_type": self.assay_type,
            "ic50_nm": self.ic50_nm,
            "max_inhibition_pct": self.max_inhibition_pct,
            "hill_coefficient": self.hill_coefficient,
            "z_prime": self.z_prime,
            "is_active": self.is_active,
            "notes": self.notes,
        }


@dataclass
class SelectivityResult:
    """Selectivity result for a compound across target pairs.

    Attributes
    ----------
    compound_id : str
        Compound identifier.
    parasite_target : str
        Parasite target name.
    human_target : str
        Human orthologue name.
    parasite_ic50_nm : float | None
        IC50 against parasite target.
    human_ic50_nm : float | None
        IC50 against human target.
    selectivity_ratio : float | None
        human_ic50 / parasite_ic50.
    is_selective : bool
        Whether selectivity ratio > threshold.
    selectivity_threshold : float
        Threshold used for classification.
    """

    compound_id: str
    parasite_target: str
    human_target: str
    parasite_ic50_nm: float | None = None
    human_ic50_nm: float | None = None
    selectivity_ratio: float | None = None
    is_selective: bool = False
    selectivity_threshold: float = 10.0

    def __post_init__(self) -> None:
        if (
            self.parasite_ic50_nm is not None
            and self.human_ic50_nm is not None
            and self.parasite_ic50_nm > 0
        ):
            self.selectivity_ratio = self.human_ic50_nm / self.parasite_ic50_nm
            self.is_selective = self.selectivity_ratio >= self.selectivity_threshold

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "compound_id": self.compound_id,
            "parasite_target": self.parasite_target,
            "human_target": self.human_target,
            "parasite_ic50_nm": self.parasite_ic50_nm,
            "human_ic50_nm": self.human_ic50_nm,
            "selectivity_ratio": self.selectivity_ratio,
            "is_selective": self.is_selective,
            "selectivity_threshold": self.selectivity_threshold,
        }


@dataclass
class CycleIteration:
    """A single iteration of the design-build-test-learn cycle.

    Attributes
    ----------
    iteration : int
        Iteration number (1-indexed).
    phase : CyclePhase
        Current phase.
    compounds_tested : list[str]
        Compound IDs tested in this iteration.
    results : list[ExperimentResult]
        Experiment results.
    selectivity_results : list[SelectivityResult]
        Selectivity results.
    hits : list[str]
        Compound IDs classified as hits.
    eliminated : list[str]
        Compound IDs eliminated.
    key_findings : list[str]
        Key findings from this iteration.
    next_actions : list[str]
        Recommended actions for next iteration.
    started_at : str
        Timestamp.
    completed_at : str
        Timestamp.
    """

    iteration: int
    phase: CyclePhase = CyclePhase.DESIGN
    compounds_tested: list[str] = field(default_factory=list)
    results: list[ExperimentResult] = field(default_factory=list)
    selectivity_results: list[SelectivityResult] = field(default_factory=list)
    hits: list[str] = field(default_factory=list)
    eliminated: list[str] = field(default_factory=list)
    key_findings: list[str] = field(default_factory=list)
    next_actions: list[str] = field(default_factory=list)
    started_at: str = ""
    completed_at: str = ""

    def __post_init__(self) -> None:
        if not self.started_at:
            self.started_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

    @property
    def n_tested(self) -> int:
        return len(self.compounds_tested)

    @property
    def n_hits(self) -> int:
        return len(self.hits)

    @property
    def hit_rate(self) -> float:
        if self.n_tested == 0:
            return 0.0
        return self.n_hits / self.n_tested

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "iteration": self.iteration,
            "phase": self.phase.value,
            "n_tested": self.n_tested,
            "n_hits": self.n_hits,
            "hit_rate": round(self.hit_rate, 3),
            "hits": self.hits,
            "eliminated": self.eliminated,
            "key_findings": self.key_findings,
            "next_actions": self.next_actions,
            "results": [r.to_dict() for r in self.results],
            "selectivity_results": [s.to_dict() for s in self.selectivity_results],
            "started_at": self.started_at,
            "completed_at": self.completed_at,
        }


@dataclass
class DiscoveryCampaign:
    """A complete discovery campaign tracking the DBTL cycle.

    Manages the state of an iterative drug discovery campaign,
    tracking compounds through the pipeline and maintaining
    the history of experiments and decisions.

    Attributes
    ----------
    campaign_id : str
        Unique campaign identifier.
    name : str
        Campaign name.
    target_pair : tuple[str, str]
        (parasite_target, human_orthologue).
    disease : str
        Target disease.
    iterations : list[CycleIteration]
        History of DBTL iterations.
    active_compounds : dict[str, CompoundStatus]
        Current status of all compounds.
    selectivity_threshold : float
        Minimum selectivity ratio to qualify as selective.
    activity_threshold_nm : float
        Maximum IC50 (nM) to qualify as active.
    """

    campaign_id: str
    name: str
    target_pair: tuple[str, str]
    disease: str = ""
    iterations: list[CycleIteration] = field(default_factory=list)
    active_compounds: dict[str, CompoundStatus] = field(default_factory=dict)
    selectivity_threshold: float = 10.0
    activity_threshold_nm: float = 1000.0

    @property
    def current_iteration(self) -> int:
        return len(self.iterations)

    @property
    def total_compounds_tested(self) -> int:
        return len(self.active_compounds)

    @property
    def total_hits(self) -> int:
        return sum(
            1
            for s in self.active_compounds.values()
            if s in {CompoundStatus.HIT, CompoundStatus.LEAD}
        )

    @property
    def total_eliminated(self) -> int:
        return sum(
            1
            for s in self.active_compounds.values()
            if s == CompoundStatus.ELIMINATED
        )

    def start_iteration(self) -> CycleIteration:
        """Start a new DBTL iteration."""
        iteration = CycleIteration(
            iteration=self.current_iteration + 1,
            phase=CyclePhase.DESIGN,
        )
        self.iterations.append(iteration)
        return iteration

    def record_results(
        self,
        results: list[ExperimentResult],
        selectivity_results: list[SelectivityResult] | None = None,
    ) -> CycleIteration:
        """Record experiment results for the current iteration.

        Parameters
        ----------
        results : list[ExperimentResult]
            Experiment results to record.
        selectivity_results : list[SelectivityResult] or None
            Selectivity results if available.

        Returns
        -------
        CycleIteration
            Updated current iteration with analysis.
        """
        if not self.iterations:
            self.start_iteration()

        current = self.iterations[-1]
        current.results.extend(results)
        current.phase = CyclePhase.LEARN

        if selectivity_results:
            current.selectivity_results.extend(selectivity_results)

        # Classify compounds
        for result in results:
            cid = result.compound_id
            current.compounds_tested.append(cid)

            if result.is_active and result.ic50_nm is not None:
                if result.ic50_nm <= self.activity_threshold_nm:
                    self.active_compounds[cid] = CompoundStatus.HIT
                    current.hits.append(cid)
                else:
                    self.active_compounds[cid] = CompoundStatus.TESTED
            else:
                self.active_compounds[cid] = CompoundStatus.ELIMINATED
                current.eliminated.append(cid)

        # Process selectivity results
        if selectivity_results:
            for sel in selectivity_results:
                if sel.is_selective:
                    # Upgrade from HIT to LEAD if selective
                    if self.active_compounds.get(sel.compound_id) == CompoundStatus.HIT:
                        self.active_compounds[sel.compound_id] = CompoundStatus.LEAD
                        current.key_findings.append(
                            f"{sel.compound_id}: selective "
                            f"({sel.selectivity_ratio:.1f}x) — promoted to LEAD"
                        )
                elif sel.selectivity_ratio is not None:
                    current.key_findings.append(
                        f"{sel.compound_id}: non-selective "
                        f"({sel.selectivity_ratio:.1f}x) — flagged for off-target risk"
                    )

        # Generate recommendations
        current.next_actions = self._recommend_next_actions(current)
        current.completed_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

        return current

    def _recommend_next_actions(
        self, iteration: CycleIteration
    ) -> list[str]:
        """Generate recommendations for the next iteration."""
        actions: list[str] = []

        if iteration.n_hits == 0:
            actions.append(
                "No hits found. Consider expanding chemical space or "
                "lowering activity threshold."
            )
        else:
            actions.append(
                f"{iteration.n_hits} hits identified. "
                "Proceed to selectivity counter-screen if not already done."
            )

        leads = [
            cid
            for cid, status in self.active_compounds.items()
            if status == CompoundStatus.LEAD
        ]
        if leads:
            actions.append(
                f"{len(leads)} leads with confirmed selectivity. "
                "Recommend ADMET profiling and cytotoxicity counter-screen."
            )

        if iteration.hit_rate < 0.05:
            actions.append(
                f"Low hit rate ({iteration.hit_rate:.1%}). "
                "Consider SAR analysis of tested compounds to guide "
                "next round of selections."
            )

        return actions

    def summary(self) -> str:
        """Campaign summary."""
        return (
            f"Campaign: {self.name} ({self.campaign_id})\n"
            f"  Target: {self.target_pair[0]} vs {self.target_pair[1]}\n"
            f"  Iterations: {self.current_iteration}\n"
            f"  Compounds: {self.total_compounds_tested} tested, "
            f"{self.total_hits} hits, {self.total_eliminated} eliminated\n"
            f"  Selectivity threshold: {self.selectivity_threshold}x\n"
            f"  Activity threshold: {self.activity_threshold_nm} nM"
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "campaign_id": self.campaign_id,
            "name": self.name,
            "target_pair": list(self.target_pair),
            "disease": self.disease,
            "current_iteration": self.current_iteration,
            "total_compounds_tested": self.total_compounds_tested,
            "total_hits": self.total_hits,
            "total_eliminated": self.total_eliminated,
            "selectivity_threshold": self.selectivity_threshold,
            "activity_threshold_nm": self.activity_threshold_nm,
            "iterations": [it.to_dict() for it in self.iterations],
            "compound_statuses": {
                cid: status.value
                for cid, status in self.active_compounds.items()
            },
        }
