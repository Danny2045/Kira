"""Cell-type expression context for drug targets.

Maps drug targets to their expression profiles across cell types
and tissues. This is critical for understanding:
  1. Where a target is expressed (on-target tissues)
  2. Where the human orthologue is expressed (off-target risk)
  3. Whether delivery to the right tissue will reach the right cells

Integrates with public expression data (Human Protein Atlas,
GTEx, CZ CELLxGENE) to provide expression context for Kira's
selectivity analysis.

References:
    - Human Protein Atlas (proteinatlas.org)
    - GTEx Consortium (2020) Science — tissue expression
    - CZ CELLxGENE Census — single-cell expression
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class ExpressionLevel(Enum):
    """Qualitative expression level categories."""

    NOT_DETECTED = "not_detected"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"

    @classmethod
    def from_tpm(cls, tpm: float) -> ExpressionLevel:
        """Classify TPM value into expression level.

        Thresholds based on Human Protein Atlas conventions:
        - Not detected: < 1 TPM
        - Low: 1-10 TPM
        - Medium: 10-100 TPM
        - High: > 100 TPM
        """
        if tpm < 1.0:
            return cls.NOT_DETECTED
        if tpm < 10.0:
            return cls.LOW
        if tpm < 100.0:
            return cls.MEDIUM
        return cls.HIGH


@dataclass
class TissueExpression:
    """Expression of a gene in a specific tissue.

    Attributes
    ----------
    tissue : str
        Tissue name (e.g., 'liver', 'brain', 'lung').
    tpm : float
        Transcripts per million (TPM) expression value.
    level : ExpressionLevel
        Qualitative expression level.
    cell_types : dict[str, float]
        Cell-type breakdown within the tissue (cell_type -> TPM).
        Empty if single-cell data not available.
    source : str
        Data source (e.g., 'HPA', 'GTEx', 'CELLxGENE').
    """

    tissue: str
    tpm: float
    level: ExpressionLevel
    cell_types: dict[str, float] = field(default_factory=dict)
    source: str = "curated"


@dataclass
class GeneExpressionProfile:
    """Complete expression profile for a gene across tissues.

    Attributes
    ----------
    gene_name : str
        Gene symbol (e.g., 'DHODH', 'HDAC8', 'TXNRD1').
    organism : str
        Source organism (e.g., 'Homo sapiens', 'Schistosoma mansoni').
    uniprot_id : str
        UniProt accession if available.
    tissue_expression : list[TissueExpression]
        Expression across tissues.
    primary_tissues : list[str]
        Tissues where expression is HIGH.
    ubiquitous : bool
        Whether the gene is expressed in most tissues.
    """

    gene_name: str
    organism: str
    uniprot_id: str = ""
    tissue_expression: list[TissueExpression] = field(default_factory=list)
    primary_tissues: list[str] = field(default_factory=list)
    ubiquitous: bool = False

    @property
    def n_tissues_expressed(self) -> int:
        """Number of tissues with detectable expression."""
        return sum(
            1 for te in self.tissue_expression
            if te.level != ExpressionLevel.NOT_DETECTED
        )

    @property
    def max_tpm(self) -> float:
        """Maximum TPM across all tissues."""
        if not self.tissue_expression:
            return 0.0
        return max(te.tpm for te in self.tissue_expression)

    @property
    def tissue_specificity(self) -> float:
        """Tissue specificity index (0 = ubiquitous, 1 = tissue-specific).

        Computed as 1 - (entropy / max_entropy) over tissue TPM values.
        """
        import math

        tpms = [te.tpm for te in self.tissue_expression if te.tpm > 0]
        if len(tpms) <= 1:
            return 1.0

        total = sum(tpms)
        probs = [t / total for t in tpms]
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        max_entropy = math.log2(len(tpms))

        if max_entropy == 0:
            return 0.0

        return 1.0 - (entropy / max_entropy)

    def get_expression_in(self, tissue: str) -> TissueExpression | None:
        """Get expression data for a specific tissue."""
        for te in self.tissue_expression:
            if te.tissue.lower() == tissue.lower():
                return te
        return None

    def summary(self) -> str:
        """Human-readable summary."""
        top_tissues = sorted(
            self.tissue_expression, key=lambda t: t.tpm, reverse=True
        )[:3]
        top_str = ", ".join(
            f"{t.tissue}: {t.tpm:.1f} TPM" for t in top_tissues
        )
        return (
            f"{self.gene_name} ({self.organism}): "
            f"expressed in {self.n_tissues_expressed} tissues, "
            f"specificity={self.tissue_specificity:.2f}. "
            f"Top: {top_str}"
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "gene_name": self.gene_name,
            "organism": self.organism,
            "uniprot_id": self.uniprot_id,
            "n_tissues_expressed": self.n_tissues_expressed,
            "max_tpm": self.max_tpm,
            "tissue_specificity": round(self.tissue_specificity, 3),
            "ubiquitous": self.ubiquitous,
            "primary_tissues": self.primary_tissues,
            "tissue_expression": [
                {
                    "tissue": te.tissue,
                    "tpm": te.tpm,
                    "level": te.level.value,
                    "source": te.source,
                }
                for te in self.tissue_expression
            ],
        }


# ---------------------------------------------------------------------------
# Curated expression data for Kira's NTD targets
# ---------------------------------------------------------------------------

# Human orthologue expression profiles (from Human Protein Atlas / GTEx)
# These are the targets whose expression determines off-target risk
# when delivering anti-parasitic compounds.

def _build_hsdhodh_profile() -> GeneExpressionProfile:
    """Human DHODH expression profile.

    DHODH (dihydroorotate dehydrogenase) is a mitochondrial enzyme
    in the de novo pyrimidine biosynthesis pathway. It is ubiquitously
    expressed but with highest levels in rapidly dividing cells.

    Source: Human Protein Atlas + GTEx consensus.
    """
    tissues = [
        TissueExpression("liver", 45.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("bone_marrow", 120.0, ExpressionLevel.HIGH, source="HPA"),
        TissueExpression("lymph_node", 85.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("small_intestine", 65.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("spleen", 55.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("thymus", 95.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("skin", 30.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("lung", 25.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("kidney", 20.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("heart", 15.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("brain", 12.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("muscle", 8.0, ExpressionLevel.LOW, source="HPA"),
        TissueExpression("adipose", 5.0, ExpressionLevel.LOW, source="HPA"),
    ]
    return GeneExpressionProfile(
        gene_name="DHODH",
        organism="Homo sapiens",
        uniprot_id="Q02127",
        tissue_expression=tissues,
        primary_tissues=["bone_marrow", "thymus"],
        ubiquitous=True,
    )


def _build_hshdac8_profile() -> GeneExpressionProfile:
    """Human HDAC8 expression profile.

    HDAC8 (histone deacetylase 8) is a class I HDAC involved in
    chromatin remodeling. Broadly expressed with enrichment in
    smooth muscle and brain.

    Source: Human Protein Atlas + GTEx consensus.
    """
    tissues = [
        TissueExpression("brain", 35.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("smooth_muscle", 50.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("lung", 28.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("liver", 22.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("kidney", 20.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("heart", 18.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("spleen", 15.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("bone_marrow", 12.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("small_intestine", 10.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("skin", 8.0, ExpressionLevel.LOW, source="HPA"),
        TissueExpression("muscle", 6.0, ExpressionLevel.LOW, source="HPA"),
        TissueExpression("adipose", 4.0, ExpressionLevel.LOW, source="HPA"),
    ]
    return GeneExpressionProfile(
        gene_name="HDAC8",
        organism="Homo sapiens",
        uniprot_id="Q9BY41",
        tissue_expression=tissues,
        primary_tissues=["smooth_muscle", "brain"],
        ubiquitous=True,
    )


def _build_hstxnrd1_profile() -> GeneExpressionProfile:
    """Human TXNRD1 (thioredoxin reductase 1) expression profile.

    TXNRD1 is the human orthologue relevant for SmTGR selectivity.
    Ubiquitously expressed — essential antioxidant enzyme.

    Source: Human Protein Atlas + GTEx consensus.
    """
    tissues = [
        TissueExpression("liver", 150.0, ExpressionLevel.HIGH, source="HPA"),
        TissueExpression("kidney", 95.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("lung", 60.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("heart", 45.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("brain", 40.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("spleen", 35.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("small_intestine", 55.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("muscle", 30.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("skin", 25.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("bone_marrow", 40.0, ExpressionLevel.MEDIUM, source="HPA"),
        TissueExpression("adipose", 20.0, ExpressionLevel.MEDIUM, source="HPA"),
    ]
    return GeneExpressionProfile(
        gene_name="TXNRD1",
        organism="Homo sapiens",
        uniprot_id="Q16881",
        tissue_expression=tissues,
        primary_tissues=["liver"],
        ubiquitous=True,
    )


# Registry of curated expression profiles
EXPRESSION_PROFILES: dict[str, GeneExpressionProfile] = {
    "DHODH": _build_hsdhodh_profile(),
    "HDAC8": _build_hshdac8_profile(),
    "TXNRD1": _build_hstxnrd1_profile(),
}


def get_expression_profile(gene_name: str) -> GeneExpressionProfile | None:
    """Look up a curated expression profile by gene name.

    Parameters
    ----------
    gene_name : str
        Gene symbol (case-insensitive).

    Returns
    -------
    GeneExpressionProfile or None
        Curated profile if available, None otherwise.
    """
    return EXPRESSION_PROFILES.get(gene_name.upper())


# ---------------------------------------------------------------------------
# Selectivity-expression integration
# ---------------------------------------------------------------------------

@dataclass
class SelectivityExpressionContext:
    """Integrates selectivity data with expression context.

    For a parasite target and its human orthologue, this answers:
    "If we deliver a selective compound to tissue X, how much
    off-target inhibition will occur based on orthologue expression?"

    Attributes
    ----------
    parasite_target : str
        Parasite target name (e.g., 'SmDHODH').
    human_orthologue : str
        Human orthologue gene name (e.g., 'DHODH').
    selectivity_ratio : float
        IC50(human) / IC50(parasite). Higher = more selective.
    human_profile : GeneExpressionProfile
        Expression profile of the human orthologue.
    tissue_risk : dict[str, float]
        Tissue -> off-target risk score (0-1).
        Risk = expression_level * (1/selectivity_ratio).
    safest_delivery_tissues : list[str]
        Tissues where delivery is safest (lowest off-target expression).
    riskiest_tissues : list[str]
        Tissues where off-target inhibition risk is highest.
    """

    parasite_target: str
    human_orthologue: str
    selectivity_ratio: float
    human_profile: GeneExpressionProfile
    tissue_risk: dict[str, float] = field(default_factory=dict)
    safest_delivery_tissues: list[str] = field(default_factory=list)
    riskiest_tissues: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Human-readable summary."""
        safe_str = ", ".join(self.safest_delivery_tissues[:3])
        risky_str = ", ".join(self.riskiest_tissues[:3])
        return (
            f"{self.parasite_target} ({self.selectivity_ratio:.1f}x selective): "
            f"human orthologue {self.human_orthologue} is "
            f"{'ubiquitous' if self.human_profile.ubiquitous else 'tissue-specific'}. "
            f"Safest delivery: {safe_str}. "
            f"Highest risk: {risky_str}."
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "parasite_target": self.parasite_target,
            "human_orthologue": self.human_orthologue,
            "selectivity_ratio": self.selectivity_ratio,
            "tissue_risk": self.tissue_risk,
            "safest_delivery_tissues": self.safest_delivery_tissues,
            "riskiest_tissues": self.riskiest_tissues,
            "human_profile": self.human_profile.to_dict(),
        }


def compute_selectivity_expression_context(
    parasite_target: str,
    human_orthologue: str,
    selectivity_ratio: float,
) -> SelectivityExpressionContext | None:
    """Compute selectivity-expression context for a target pair.

    Integrates Kira's selectivity ratios with human orthologue
    expression data to identify which tissues are safest/riskiest
    for compound delivery.

    Parameters
    ----------
    parasite_target : str
        Parasite target name.
    human_orthologue : str
        Human orthologue gene symbol.
    selectivity_ratio : float
        IC50(human) / IC50(parasite).

    Returns
    -------
    SelectivityExpressionContext or None
        Context if expression data available, None otherwise.
    """
    profile = get_expression_profile(human_orthologue)
    if profile is None:
        return None

    # Compute per-tissue risk: normalized expression * (1/selectivity)
    # Higher expression + lower selectivity = higher risk
    max_tpm = profile.max_tpm
    if max_tpm == 0:
        return None

    tissue_risk: dict[str, float] = {}
    for te in profile.tissue_expression:
        # Normalized expression (0-1)
        norm_expr = te.tpm / max_tpm
        # Risk is proportional to expression and inversely to selectivity
        risk = norm_expr / max(selectivity_ratio, 0.1)
        tissue_risk[te.tissue] = round(min(1.0, risk), 4)

    # Sort tissues by risk
    sorted_tissues = sorted(tissue_risk.items(), key=lambda x: x[1])
    safest = [t for t, r in sorted_tissues[:3]]
    riskiest = [t for t, r in sorted_tissues[-3:]]
    riskiest.reverse()

    return SelectivityExpressionContext(
        parasite_target=parasite_target,
        human_orthologue=human_orthologue,
        selectivity_ratio=selectivity_ratio,
        human_profile=profile,
        tissue_risk=tissue_risk,
        safest_delivery_tissues=safest,
        riskiest_tissues=riskiest,
    )
