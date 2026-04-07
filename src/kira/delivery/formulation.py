"""LNP formulation design space.

Models the combinatorial design space of lipid nanoparticle formulations:
ionizable lipid, helper lipid, cholesterol analog, PEG-lipid, N:P ratio,
particle size, and targeting peptides.

Each formulation parameter affects multiple delivery outcomes simultaneously:
circulation half-life, tissue tropism, cell binding, endosomal escape,
intracellular action, immunogenicity, and manufacturability.

References:
    - Hou et al. (2021) Nature Reviews Materials — LNP design principles
    - Patel et al. (2023) Nature Nanotechnology — organ-selective LNPs
    - Cheng et al. (2020) Nature Nanotechnology — selective organ targeting (SORT)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class IonizableLipidClass(Enum):
    """Major ionizable lipid families for mRNA delivery."""

    DLIN_MC3 = "DLin-MC3-DMA"  # FDA-approved (Onpattro), liver-tropic
    SM102 = "SM-102"  # Moderna COVID vaccine
    ALC0315 = "ALC-0315"  # BioNTech COVID vaccine
    LIPID5 = "Lipid 5"  # Optimized for low-dose mRNA delivery
    C12_200 = "C12-200"  # Potent, hepatocyte-selective
    CUSTOM = "custom"


class HelperLipidClass(Enum):
    """Helper lipids that modulate membrane properties."""

    DSPC = "DSPC"  # Standard saturated phospholipid
    DOPE = "DOPE"  # Fusogenic, enhances endosomal escape
    DPPC = "DPPC"  # Alternative saturated phospholipid
    DOPC = "DOPC"  # Unsaturated, more fluid membranes


class PEGLipidClass(Enum):
    """PEG-lipid variants controlling circulation and immunogenicity."""

    DMG_PEG2000 = "DMG-PEG2000"  # Short acyl chain, faster shedding
    DSPE_PEG2000 = "DSPE-PEG2000"  # Long acyl chain, stable coating
    C14_PEG2000 = "C14-PEG2000"  # Intermediate shedding rate
    C18_PEG2000 = "C18-PEG2000"  # Slow shedding, extended circulation


class AdministrationRoute(Enum):
    """Route of administration affecting biodistribution."""

    IV = "intravenous"  # Systemic, strong liver tropism
    IM = "intramuscular"  # Local + draining lymph nodes
    SC = "subcutaneous"  # Local depot, slow release
    IN = "intranasal"  # Mucosal, lung-accessible
    IT = "intrathecal"  # CNS-directed
    INTRATUMORAL = "intratumoral"  # Direct tumor injection
    NEBULIZED = "nebulized"  # Lung epithelium


@dataclass
class TargetingPeptide:
    """A targeting peptide conjugated to the LNP surface.

    Peptides can redirect LNP tropism away from the liver
    toward specific tissues or cell types.

    Attributes
    ----------
    name : str
        Peptide identifier or name.
    sequence : str
        Amino acid sequence.
    target_receptor : str
        Known or putative receptor target.
    tissue_affinity : dict[str, float]
        Tissue name -> relative affinity score (0-1).
    kd_nm : float | None
        Binding affinity in nM if known.
    source : str
        How the peptide was identified (phage display, computational, etc.).
    """

    name: str
    sequence: str
    target_receptor: str
    tissue_affinity: dict[str, float] = field(default_factory=dict)
    kd_nm: float | None = None
    source: str = "unknown"


@dataclass
class LNPFormulation:
    """Complete LNP formulation specification.

    Encodes all tunable parameters of a lipid nanoparticle.
    The design space is combinatorial: each parameter affects
    multiple delivery outcomes simultaneously.

    Attributes
    ----------
    name : str
        Formulation identifier.
    ionizable_lipid : IonizableLipidClass
        Primary ionizable lipid component.
    ionizable_lipid_mol_pct : float
        Molar percentage of ionizable lipid (typically 30-50%).
    helper_lipid : HelperLipidClass
        Structural/fusogenic helper lipid.
    helper_lipid_mol_pct : float
        Molar percentage of helper lipid (typically 10-20%).
    cholesterol_mol_pct : float
        Molar percentage of cholesterol (typically 30-45%).
    peg_lipid : PEGLipidClass
        PEG-lipid for steric stabilization.
    peg_lipid_mol_pct : float
        Molar percentage of PEG-lipid (typically 1-5%).
    np_ratio : float
        Nitrogen-to-phosphate ratio (typically 3-10).
    target_size_nm : float
        Target particle diameter in nanometers (typically 60-150).
    route : AdministrationRoute
        Intended route of administration.
    targeting_peptides : list[TargetingPeptide]
        Optional targeting peptides conjugated to the surface.
    payload_type : str
        Type of payload (mRNA, siRNA, CRISPR-Cas9 RNP, etc.).
    payload_size_nt : int
        Payload size in nucleotides (affects encapsulation).
    """

    name: str
    ionizable_lipid: IonizableLipidClass = IonizableLipidClass.SM102
    ionizable_lipid_mol_pct: float = 50.0
    helper_lipid: HelperLipidClass = HelperLipidClass.DSPC
    helper_lipid_mol_pct: float = 10.0
    cholesterol_mol_pct: float = 38.5
    peg_lipid: PEGLipidClass = PEGLipidClass.DMG_PEG2000
    peg_lipid_mol_pct: float = 1.5
    np_ratio: float = 6.0
    target_size_nm: float = 80.0
    route: AdministrationRoute = AdministrationRoute.IV
    targeting_peptides: list[TargetingPeptide] = field(default_factory=list)
    payload_type: str = "mRNA"
    payload_size_nt: int = 1000

    @property
    def total_lipid_mol_pct(self) -> float:
        """Total molar percentage of all lipid components."""
        return (
            self.ionizable_lipid_mol_pct
            + self.helper_lipid_mol_pct
            + self.cholesterol_mol_pct
            + self.peg_lipid_mol_pct
        )

    @property
    def is_balanced(self) -> bool:
        """Check if lipid molar percentages sum to ~100%."""
        return abs(self.total_lipid_mol_pct - 100.0) < 1.0

    @property
    def has_targeting(self) -> bool:
        """Whether the formulation includes targeting peptides."""
        return len(self.targeting_peptides) > 0

    def validate(self) -> list[str]:
        """Return a list of validation warnings for this formulation."""
        warnings: list[str] = []

        if not self.is_balanced:
            warnings.append(
                f"Lipid molar percentages sum to {self.total_lipid_mol_pct:.1f}%, "
                f"expected ~100%"
            )

        if self.ionizable_lipid_mol_pct < 20 or self.ionizable_lipid_mol_pct > 60:
            warnings.append(
                f"Ionizable lipid at {self.ionizable_lipid_mol_pct:.1f}% is outside "
                f"typical range (20-60%)"
            )

        if self.np_ratio < 2 or self.np_ratio > 15:
            warnings.append(
                f"N:P ratio {self.np_ratio:.1f} is outside typical range (2-15)"
            )

        if self.target_size_nm < 40 or self.target_size_nm > 200:
            warnings.append(
                f"Target size {self.target_size_nm:.0f} nm is outside typical range "
                f"(40-200 nm)"
            )

        if self.peg_lipid_mol_pct > 5:
            warnings.append(
                f"PEG-lipid at {self.peg_lipid_mol_pct:.1f}% may reduce cellular uptake"
            )

        if self.payload_size_nt > 5000 and self.target_size_nm < 80:
            warnings.append(
                f"Large payload ({self.payload_size_nt} nt) may require larger particle "
                f"(current target: {self.target_size_nm:.0f} nm)"
            )

        return warnings

    def to_dict(self) -> dict:
        """Serialize formulation to a dictionary."""
        return {
            "name": self.name,
            "ionizable_lipid": self.ionizable_lipid.value,
            "ionizable_lipid_mol_pct": self.ionizable_lipid_mol_pct,
            "helper_lipid": self.helper_lipid.value,
            "helper_lipid_mol_pct": self.helper_lipid_mol_pct,
            "cholesterol_mol_pct": self.cholesterol_mol_pct,
            "peg_lipid": self.peg_lipid.value,
            "peg_lipid_mol_pct": self.peg_lipid_mol_pct,
            "np_ratio": self.np_ratio,
            "target_size_nm": self.target_size_nm,
            "route": self.route.value,
            "targeting_peptides": [
                {
                    "name": p.name,
                    "sequence": p.sequence,
                    "target_receptor": p.target_receptor,
                    "tissue_affinity": p.tissue_affinity,
                    "kd_nm": p.kd_nm,
                    "source": p.source,
                }
                for p in self.targeting_peptides
            ],
            "payload_type": self.payload_type,
            "payload_size_nt": self.payload_size_nt,
            "total_lipid_mol_pct": self.total_lipid_mol_pct,
            "is_balanced": self.is_balanced,
            "has_targeting": self.has_targeting,
        }


# ---------------------------------------------------------------------------
# Reference formulations from literature
# ---------------------------------------------------------------------------

ONPATTRO_FORMULATION = LNPFormulation(
    name="Onpattro (patisiran)",
    ionizable_lipid=IonizableLipidClass.DLIN_MC3,
    ionizable_lipid_mol_pct=50.0,
    helper_lipid=HelperLipidClass.DSPC,
    helper_lipid_mol_pct=10.0,
    cholesterol_mol_pct=38.5,
    peg_lipid=PEGLipidClass.DMG_PEG2000,
    peg_lipid_mol_pct=1.5,
    np_ratio=3.0,
    target_size_nm=80.0,
    route=AdministrationRoute.IV,
    payload_type="siRNA",
    payload_size_nt=21,
)

MODERNA_COVID_FORMULATION = LNPFormulation(
    name="Spikevax (mRNA-1273)",
    ionizable_lipid=IonizableLipidClass.SM102,
    ionizable_lipid_mol_pct=50.0,
    helper_lipid=HelperLipidClass.DSPC,
    helper_lipid_mol_pct=10.0,
    cholesterol_mol_pct=38.5,
    peg_lipid=PEGLipidClass.DMG_PEG2000,
    peg_lipid_mol_pct=1.5,
    np_ratio=6.0,
    target_size_nm=100.0,
    route=AdministrationRoute.IM,
    payload_type="mRNA",
    payload_size_nt=4284,
)

BIONTECH_COVID_FORMULATION = LNPFormulation(
    name="Comirnaty (BNT162b2)",
    ionizable_lipid=IonizableLipidClass.ALC0315,
    ionizable_lipid_mol_pct=46.3,
    helper_lipid=HelperLipidClass.DSPC,
    helper_lipid_mol_pct=9.4,
    cholesterol_mol_pct=42.7,
    peg_lipid=PEGLipidClass.DMG_PEG2000,
    peg_lipid_mol_pct=1.6,
    np_ratio=6.0,
    target_size_nm=80.0,
    route=AdministrationRoute.IM,
    payload_type="mRNA",
    payload_size_nt=4284,
)
