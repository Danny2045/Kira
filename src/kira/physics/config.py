"""Default configuration for Physics Auditor.

All thresholds, weights, and parameters are defined here.
Override via YAML config file passed to CLI.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ClashConfig:
    """Steric clash detection parameters."""

    vdw_tolerance: float = 0.4  # Angstroms subtracted from sum of vdW radii
    ignore_bonded_neighbors: int = 3  # Skip 1-2 and 1-3 bonded pairs
    severe_clash_threshold: float = 0.2  # Below this overlap = severe


@dataclass(frozen=True)
class BondGeometryConfig:
    """Bond length and angle validation parameters."""

    zscore_warning: float = 3.0  # Flag bonds/angles beyond this z-score
    zscore_severe: float = 5.0  # Severe violation threshold


@dataclass(frozen=True)
class RamachandranConfig:
    """Backbone dihedral validation parameters."""

    # Fraction thresholds for scoring
    favored_weight: float = 1.0
    allowed_weight: float = 0.5
    generous_weight: float = 0.1
    outlier_weight: float = 0.0


@dataclass(frozen=True)
class PeptideConfig:
    """Peptide bond planarity parameters."""

    omega_trans_ideal: float = 180.0  # degrees
    omega_cis_ideal: float = 0.0  # degrees (rare, mostly Pro)
    strain_threshold: float = 20.0  # degrees deviation = strained
    severe_threshold: float = 30.0  # degrees deviation = severe


@dataclass(frozen=True)
class ChiralityConfig:
    """Chirality validation parameters."""

    # Improper dihedral at Cα: positive for L-amino acids
    l_amino_acid_sign: float = 1.0  # Expected sign of improper dihedral
    tolerance_degrees: float = 30.0  # Tolerance around ideal


@dataclass(frozen=True)
class RotamerConfig:
    """Sidechain rotamer validation parameters."""

    outlier_threshold: float = 1.0  # % probability below which = outlier


@dataclass(frozen=True)
class LennardJonesConfig:
    """Lennard-Jones energy computation parameters."""

    cutoff_angstrom: float = 12.0  # Distance cutoff for LJ computation
    epsilon_default: float = 0.05  # kcal/mol, default well depth
    sigma_default: float = 3.4  # Angstroms, default collision diameter
    energy_cap: float = 1000.0  # Cap per-pair energy to avoid infinity


@dataclass(frozen=True)
class DisulfideConfig:
    """Disulfide bond validation parameters."""

    ss_bond_length_ideal: float = 2.05  # Angstroms
    ss_bond_length_tolerance: float = 0.15  # Angstroms
    ss_detection_cutoff: float = 2.5  # Angstroms, max distance to consider
    dihedral_ideal: float = 90.0  # |χ3| should be ~90°
    dihedral_tolerance: float = 30.0  # degrees


@dataclass(frozen=True)
class CompositeConfig:
    """Weights for composite trust score and recommendation thresholds."""

    # Check weights (must sum to 1.0)
    weight_clashes: float = 0.20
    weight_bond_geometry: float = 0.15
    weight_ramachandran: float = 0.15
    weight_lj_energy: float = 0.15
    weight_peptide_planarity: float = 0.10
    weight_chirality: float = 0.10
    weight_rotamers: float = 0.10
    weight_disulfides: float = 0.05

    # Recommendation thresholds
    accept_threshold: float = 0.85
    short_md_threshold: float = 0.60
    # Below short_md_threshold = discard


@dataclass(frozen=True)
class AuditorConfig:
    """Top-level configuration aggregating all sub-configs."""

    clashes: ClashConfig = field(default_factory=ClashConfig)
    bond_geometry: BondGeometryConfig = field(default_factory=BondGeometryConfig)
    ramachandran: RamachandranConfig = field(default_factory=RamachandranConfig)
    peptide: PeptideConfig = field(default_factory=PeptideConfig)
    chirality: ChiralityConfig = field(default_factory=ChiralityConfig)
    rotamers: RotamerConfig = field(default_factory=RotamerConfig)
    lennard_jones: LennardJonesConfig = field(default_factory=LennardJonesConfig)
    disulfides: DisulfideConfig = field(default_factory=DisulfideConfig)
    composite: CompositeConfig = field(default_factory=CompositeConfig)


def load_config(yaml_path: str | None = None) -> AuditorConfig:
    """Load configuration, optionally overriding defaults from a YAML file.

    Parameters
    ----------
    yaml_path : str or None
        Path to YAML config file. If None, returns defaults.

    Returns
    -------
    AuditorConfig
        Complete configuration with all parameters.
    """
    if yaml_path is None:
        return AuditorConfig()

    import yaml

    with open(yaml_path) as f:
        overrides = yaml.safe_load(f)

    if overrides is None:
        return AuditorConfig()

    # Build sub-configs from overrides
    kwargs = {}
    config_map = {
        "clashes": ClashConfig,
        "bond_geometry": BondGeometryConfig,
        "ramachandran": RamachandranConfig,
        "peptide": PeptideConfig,
        "chirality": ChiralityConfig,
        "rotamers": RotamerConfig,
        "lennard_jones": LennardJonesConfig,
        "disulfides": DisulfideConfig,
        "composite": CompositeConfig,
    }

    for key, cls in config_map.items():
        if key in overrides:
            kwargs[key] = cls(**overrides[key])

    return AuditorConfig(**kwargs)
