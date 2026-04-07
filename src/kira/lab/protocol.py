"""Structured experiment protocol design.

Generates experimental protocols for selectivity assays, dose-response
curves, and compound screening — the building blocks of the
design-build-test-learn cycle.

Protocols are structured as machine-readable specifications that can be:
  1. Reviewed by a scientist before execution
  2. Translated into cloud lab API calls (Strateos, ECL, etc.)
  3. Compared against previous protocols for reproducibility
  4. Used as training data for closed-loop optimization

References:
    - Ginkgo/OpenAI cloud lab collaboration (2024)
    - Strateos API specification
    - Emerald Cloud Lab protocol format
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum


class AssayType(Enum):
    """Standard assay types for drug discovery."""

    DOSE_RESPONSE = "dose_response"
    SELECTIVITY_PANEL = "selectivity_panel"
    ADMET = "admet"
    CYTOTOXICITY = "cytotoxicity"
    WHOLE_ORGANISM = "whole_organism"
    BINDING_AFFINITY = "binding_affinity"
    ENZYMATIC_ACTIVITY = "enzymatic_activity"


class ReadoutType(Enum):
    """Measurement readout types."""

    FLUORESCENCE = "fluorescence"
    LUMINESCENCE = "luminescence"
    ABSORBANCE = "absorbance"
    MASS_SPEC = "mass_spectrometry"
    CELL_VIABILITY = "cell_viability"
    IMAGING = "imaging"
    RADIOACTIVE = "radioactive"


class PlateFormat(Enum):
    """Microplate formats."""

    PLATE_96 = 96
    PLATE_384 = 384
    PLATE_1536 = 1536


@dataclass
class ConcentrationSeries:
    """A series of concentrations for dose-response experiments.

    Attributes
    ----------
    top_concentration_um : float
        Highest concentration in micromolar.
    dilution_factor : float
        Fold dilution between each step (typically 3 or 10).
    n_points : int
        Number of concentration points.
    n_replicates : int
        Number of replicates per concentration.
    """

    top_concentration_um: float = 100.0
    dilution_factor: float = 3.0
    n_points: int = 8
    n_replicates: int = 3

    @property
    def concentrations_um(self) -> list[float]:
        """Generate the concentration series in micromolar."""
        return [
            self.top_concentration_um / (self.dilution_factor ** i)
            for i in range(self.n_points)
        ]

    @property
    def concentrations_nm(self) -> list[float]:
        """Generate the concentration series in nanomolar."""
        return [c * 1000 for c in self.concentrations_um]

    @property
    def total_wells(self) -> int:
        """Total wells needed for this series."""
        return self.n_points * self.n_replicates

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "top_concentration_um": self.top_concentration_um,
            "dilution_factor": self.dilution_factor,
            "n_points": self.n_points,
            "n_replicates": self.n_replicates,
            "concentrations_um": [round(c, 4) for c in self.concentrations_um],
            "total_wells": self.total_wells,
        }


@dataclass
class CompoundSpec:
    """Specification for a compound in an experiment.

    Attributes
    ----------
    compound_id : str
        Identifier (ChEMBL ID or internal).
    name : str
        Human-readable name.
    smiles : str
        SMILES string for the compound.
    stock_concentration_mm : float
        Stock solution concentration in mM.
    volume_ul : float
        Required volume in microliters.
    """

    compound_id: str
    name: str
    smiles: str = ""
    stock_concentration_mm: float = 10.0
    volume_ul: float = 50.0


@dataclass
class TargetSpec:
    """Specification for a target in an experiment.

    Attributes
    ----------
    target_name : str
        Target identifier (e.g., 'SmDHODH', 'HsDHODH').
    organism : str
        Source organism.
    uniprot_id : str
        UniProt accession.
    protein_concentration_nm : float
        Protein concentration in nM for enzymatic assays.
    substrate : str
        Substrate name if enzymatic assay.
    substrate_concentration_um : float
        Substrate concentration in micromolar.
    """

    target_name: str
    organism: str = ""
    uniprot_id: str = ""
    protein_concentration_nm: float = 50.0
    substrate: str = ""
    substrate_concentration_um: float = 100.0


@dataclass
class ControlSpec:
    """Control specifications for an experiment.

    Attributes
    ----------
    positive_control : str
        Positive control compound or condition.
    negative_control : str
        Negative control (DMSO or vehicle).
    reference_compound : str
        Reference inhibitor with known IC50.
    reference_ic50_nm : float
        Known IC50 of reference compound.
    n_control_wells : int
        Number of wells per control.
    """

    positive_control: str = "known_inhibitor"
    negative_control: str = "DMSO"
    reference_compound: str = ""
    reference_ic50_nm: float = 0.0
    n_control_wells: int = 8


@dataclass
class ExperimentProtocol:
    """Complete experiment protocol specification.

    A machine-readable protocol that specifies everything needed
    to execute an experiment: compounds, targets, concentrations,
    controls, readout, and quality criteria.

    Attributes
    ----------
    protocol_id : str
        Unique protocol identifier.
    name : str
        Protocol name.
    assay_type : AssayType
        Type of assay.
    objective : str
        Scientific objective of the experiment.
    compounds : list[CompoundSpec]
        Compounds to test.
    targets : list[TargetSpec]
        Targets to test against.
    concentration_series : ConcentrationSeries
        Dose-response parameters.
    controls : ControlSpec
        Control specifications.
    readout_type : ReadoutType
        Measurement type.
    plate_format : PlateFormat
        Microplate format.
    incubation_hours : float
        Incubation time.
    temperature_c : float
        Incubation temperature.
    created_at : str
        Timestamp.
    notes : str
        Additional notes.
    quality_criteria : dict[str, float]
        Quality thresholds (e.g., Z' factor > 0.5).
    """

    protocol_id: str = ""
    name: str = ""
    assay_type: AssayType = AssayType.DOSE_RESPONSE
    objective: str = ""
    compounds: list[CompoundSpec] = field(default_factory=list)
    targets: list[TargetSpec] = field(default_factory=list)
    concentration_series: ConcentrationSeries = field(
        default_factory=ConcentrationSeries
    )
    controls: ControlSpec = field(default_factory=ControlSpec)
    readout_type: ReadoutType = ReadoutType.FLUORESCENCE
    plate_format: PlateFormat = PlateFormat.PLATE_384
    incubation_hours: float = 1.0
    temperature_c: float = 37.0
    created_at: str = ""
    notes: str = ""
    quality_criteria: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.protocol_id:
            self.protocol_id = f"KIRA-{uuid.uuid4().hex[:8].upper()}"
        if not self.created_at:
            self.created_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        if not self.quality_criteria:
            self.quality_criteria = {
                "z_prime_min": 0.5,
                "cv_max_pct": 20.0,
                "signal_to_background_min": 3.0,
            }

    @property
    def n_compounds(self) -> int:
        return len(self.compounds)

    @property
    def n_targets(self) -> int:
        return len(self.targets)

    @property
    def total_data_wells(self) -> int:
        """Total data wells (excluding controls)."""
        return (
            self.n_compounds
            * self.n_targets
            * self.concentration_series.total_wells
        )

    @property
    def total_control_wells(self) -> int:
        """Total control wells."""
        return self.controls.n_control_wells * 2 * self.n_targets

    @property
    def total_wells(self) -> int:
        """Total wells needed."""
        return self.total_data_wells + self.total_control_wells

    @property
    def n_plates(self) -> int:
        """Number of plates needed."""
        wells_per_plate = self.plate_format.value
        return max(1, -(-self.total_wells // wells_per_plate))  # ceiling division

    @property
    def estimated_cost_usd(self) -> float:
        """Rough cost estimate based on plate count and assay type.

        Very approximate — for planning purposes only.
        """
        cost_per_plate = {
            AssayType.DOSE_RESPONSE: 150.0,
            AssayType.SELECTIVITY_PANEL: 200.0,
            AssayType.ADMET: 500.0,
            AssayType.CYTOTOXICITY: 100.0,
            AssayType.WHOLE_ORGANISM: 300.0,
            AssayType.BINDING_AFFINITY: 250.0,
            AssayType.ENZYMATIC_ACTIVITY: 120.0,
        }
        per_plate = cost_per_plate.get(self.assay_type, 200.0)
        return self.n_plates * per_plate

    def summary(self) -> str:
        """Human-readable summary."""
        return (
            f"Protocol {self.protocol_id}: {self.name}\n"
            f"  Type: {self.assay_type.value}\n"
            f"  Objective: {self.objective}\n"
            f"  Compounds: {self.n_compounds}, Targets: {self.n_targets}\n"
            f"  Concentrations: {self.concentration_series.n_points} points, "
            f"{self.concentration_series.n_replicates} replicates\n"
            f"  Plates: {self.n_plates} x {self.plate_format.value}-well\n"
            f"  Total wells: {self.total_wells}\n"
            f"  Estimated cost: ${self.estimated_cost_usd:,.0f}"
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "protocol_id": self.protocol_id,
            "name": self.name,
            "assay_type": self.assay_type.value,
            "objective": self.objective,
            "compounds": [
                {
                    "compound_id": c.compound_id,
                    "name": c.name,
                    "smiles": c.smiles,
                    "stock_concentration_mm": c.stock_concentration_mm,
                }
                for c in self.compounds
            ],
            "targets": [
                {
                    "target_name": t.target_name,
                    "organism": t.organism,
                    "uniprot_id": t.uniprot_id,
                    "protein_concentration_nm": t.protein_concentration_nm,
                }
                for t in self.targets
            ],
            "concentration_series": self.concentration_series.to_dict(),
            "controls": {
                "positive_control": self.controls.positive_control,
                "negative_control": self.controls.negative_control,
                "reference_compound": self.controls.reference_compound,
                "reference_ic50_nm": self.controls.reference_ic50_nm,
                "n_control_wells": self.controls.n_control_wells,
            },
            "readout_type": self.readout_type.value,
            "plate_format": self.plate_format.value,
            "incubation_hours": self.incubation_hours,
            "temperature_c": self.temperature_c,
            "n_plates": self.n_plates,
            "total_wells": self.total_wells,
            "estimated_cost_usd": self.estimated_cost_usd,
            "quality_criteria": self.quality_criteria,
            "created_at": self.created_at,
            "notes": self.notes,
        }


# ---------------------------------------------------------------------------
# Protocol generators for common NTD assay types
# ---------------------------------------------------------------------------

def design_selectivity_assay(
    compounds: list[CompoundSpec],
    parasite_target: TargetSpec,
    human_target: TargetSpec,
    top_concentration_um: float = 100.0,
    n_dose_points: int = 10,
    n_replicates: int = 3,
) -> ExperimentProtocol:
    """Design a selectivity assay comparing parasite vs human target.

    The core experiment of Kira's selectivity analysis: measure
    dose-response curves against both the parasite target and its
    human orthologue to determine the selectivity ratio.

    Parameters
    ----------
    compounds : list[CompoundSpec]
        Compounds to test.
    parasite_target : TargetSpec
        Parasite target specification.
    human_target : TargetSpec
        Human orthologue specification.
    top_concentration_um : float
        Top concentration for dose-response.
    n_dose_points : int
        Number of concentration points.
    n_replicates : int
        Replicates per condition.

    Returns
    -------
    ExperimentProtocol
        Complete selectivity assay protocol.
    """
    conc_series = ConcentrationSeries(
        top_concentration_um=top_concentration_um,
        dilution_factor=3.0,
        n_points=n_dose_points,
        n_replicates=n_replicates,
    )

    controls = ControlSpec(
        positive_control=f"known_{parasite_target.target_name}_inhibitor",
        negative_control="DMSO (0.1% v/v)",
        n_control_wells=8,
    )

    protocol = ExperimentProtocol(
        name=f"Selectivity: {parasite_target.target_name} vs {human_target.target_name}",
        assay_type=AssayType.SELECTIVITY_PANEL,
        objective=(
            f"Determine selectivity ratios for {len(compounds)} compounds "
            f"between {parasite_target.target_name} ({parasite_target.organism}) "
            f"and {human_target.target_name} ({human_target.organism})"
        ),
        compounds=compounds,
        targets=[parasite_target, human_target],
        concentration_series=conc_series,
        controls=controls,
        readout_type=ReadoutType.FLUORESCENCE,
        plate_format=PlateFormat.PLATE_384,
        incubation_hours=1.0,
        temperature_c=37.0,
        notes=(
            "Run parasite and human target in parallel on the same plate "
            "to minimize inter-plate variability. Include reference inhibitor "
            "on every plate for quality control."
        ),
    )

    return protocol


def design_dose_response(
    compounds: list[CompoundSpec],
    target: TargetSpec,
    top_concentration_um: float = 100.0,
    n_dose_points: int = 8,
    n_replicates: int = 3,
) -> ExperimentProtocol:
    """Design a standard dose-response assay.

    Parameters
    ----------
    compounds : list[CompoundSpec]
        Compounds to test.
    target : TargetSpec
        Target specification.
    top_concentration_um : float
        Top concentration.
    n_dose_points : int
        Number of concentration points.
    n_replicates : int
        Replicates per condition.

    Returns
    -------
    ExperimentProtocol
        Dose-response protocol.
    """
    conc_series = ConcentrationSeries(
        top_concentration_um=top_concentration_um,
        dilution_factor=3.0,
        n_points=n_dose_points,
        n_replicates=n_replicates,
    )

    return ExperimentProtocol(
        name=f"Dose-response: {target.target_name}",
        assay_type=AssayType.DOSE_RESPONSE,
        objective=(
            f"Determine IC50 values for {len(compounds)} compounds "
            f"against {target.target_name}"
        ),
        compounds=compounds,
        targets=[target],
        concentration_series=conc_series,
        readout_type=ReadoutType.FLUORESCENCE,
        plate_format=PlateFormat.PLATE_384,
    )


def design_cytotoxicity_counter_screen(
    compounds: list[CompoundSpec],
    cell_type: str = "HEK293",
    top_concentration_um: float = 100.0,
    n_dose_points: int = 8,
    n_replicates: int = 3,
    incubation_hours: float = 48.0,
) -> ExperimentProtocol:
    """Design a cytotoxicity counter-screen.

    Essential for distinguishing selective target inhibition from
    general cytotoxicity. Compounds that kill cells non-specifically
    will look active in any assay.

    Parameters
    ----------
    compounds : list[CompoundSpec]
        Compounds to test.
    cell_type : str
        Cell line for counter-screen.
    top_concentration_um : float
        Top concentration.
    n_dose_points : int
        Number of concentration points.
    n_replicates : int
        Replicates per condition.
    incubation_hours : float
        Incubation time (longer for cytotoxicity).

    Returns
    -------
    ExperimentProtocol
        Cytotoxicity counter-screen protocol.
    """
    conc_series = ConcentrationSeries(
        top_concentration_um=top_concentration_um,
        dilution_factor=3.0,
        n_points=n_dose_points,
        n_replicates=n_replicates,
    )

    target = TargetSpec(
        target_name=f"{cell_type}_viability",
        organism="Homo sapiens",
    )

    controls = ControlSpec(
        positive_control="staurosporine (1 uM)",
        negative_control="DMSO (0.1% v/v)",
        n_control_wells=8,
    )

    return ExperimentProtocol(
        name=f"Cytotoxicity: {cell_type}",
        assay_type=AssayType.CYTOTOXICITY,
        objective=(
            f"Determine cytotoxicity (CC50) for {len(compounds)} compounds "
            f"in {cell_type} cells to calculate selectivity index"
        ),
        compounds=compounds,
        targets=[target],
        concentration_series=conc_series,
        controls=controls,
        readout_type=ReadoutType.CELL_VIABILITY,
        plate_format=PlateFormat.PLATE_384,
        incubation_hours=incubation_hours,
        temperature_c=37.0,
        notes=(
            f"CellTiter-Glo readout. 48h incubation in {cell_type} cells. "
            "Include staurosporine as positive control for cell death."
        ),
    )
