"""Target pair definitions for selectivity analysis.

Each target pair has:
- Parasite and human protein names + ChEMBL IDs
- UniProt IDs for sequence retrieval
- Known active site / binding pocket residue positions (0-indexed into aligned sequence)
- Disease context

These pocket definitions come from crystal structure analysis and
literature-reported active sites. They define the LOCAL binding-site
positions where divergence matters for selectivity — the positions
that ESM-2 global embeddings cannot see.

Pocket positions should be updated as better structural data becomes available.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TargetPair:
    """A parasite-target / human-ortholog pair for selectivity analysis."""

    name: str  # Short identifier
    disease: str
    parasite_name: str
    human_name: str
    parasite_chembl: str
    human_chembl: str
    parasite_uniprot: str
    human_uniprot: str
    # Pocket residue positions (0-indexed into the ALIGNED sequences)
    # These are the positions where binding-site divergence is computed.
    # Source of pocket definitions noted in comments.
    parasite_pocket_sequence: str  # One-letter codes of pocket residues
    human_pocket_sequence: str  # Aligned one-letter codes of pocket residues
    notes: str = ""


# ---------------------------------------------------------------------------
# All target pairs across the three diseases
# ---------------------------------------------------------------------------

TARGET_PAIRS: dict[str, TargetPair] = {
    # === SCHISTOSOMIASIS ===

    "SmDHODH": TargetPair(
        name="SmDHODH",
        disease="Schistosomiasis",
        parasite_name="Dihydroorotate dehydrogenase (quinone), mitochondrial",
        human_name="Dihydroorotate dehydrogenase (human)",
        parasite_chembl="CHEMBL4523950",
        human_chembl="CHEMBL1966",
        parasite_uniprot="G4VFD7",
        human_uniprot="Q02127",
        # Ubiquinone binding site residues from crystal structures
        # SmDHODH: modeled from human DHODH (PDB 1D3H) + homology
        # Key residues lining the ubiquinone tunnel
        # From Mori et al. 2021 (FEBS J), PDB 6UY4 (SmDHODH) / PDB 1D3H (HsDHODH)
        # 9 aligned ubiquinone binding site positions:
        # G46/M43, A49/A55, H50/H56, S53/L59, F92/F98, I128/V134, R130/R136, F357/F362, V358/P364
        # Key differences: S53→L59 (hydrophobicity flip), V358→P364 (flexibility), G46→M43 (volume)
        parasite_pocket_sequence="GAHSFIRFV",
        human_pocket_sequence="MAHLFVRFP",
        notes="ESM-2 cosine=0.9897 yet 30.8x selectivity exists. "
              "Pocket identity 55.6%. Key: Ser53→Leu59 flips pocket polarity. "
              "Central example of global-local divergence gap.",
    ),

    "SmHDAC8": TargetPair(
        name="SmHDAC8",
        disease="Schistosomiasis",
        parasite_name="Histone deacetylase 8",
        human_name="Histone deacetylase 8 (human)",
        parasite_chembl="CHEMBL3797017",
        human_chembl="CHEMBL3192",
        parasite_uniprot="A0A3Q0KTZ8",
        human_uniprot="Q9BY41",
        # HDAC8 catalytic channel residues
        # From PDB 6HSF (SmHDAC8) and PDB 1T69 (HsHDAC8)
        # Active site: Zn-binding residues + catalytic channel
        # From PDB 6HSF (SmHDAC8), PDB 1T69 (HsHDAC8), PMC9486128
        # Active site Zn-binding residues + catalytic channel + Trp198 allosteric region
        # Catalytic site is almost identical; Trp198(Sm)→Ser(Hs) is the key difference
        parasite_pocket_sequence="DDHHDYCFWL",
        human_pocket_sequence="DDHHDYCFSL",  # Only Trp→Ser at allosteric site
        notes="90.4% non-selective. Catalytic site conserved (90% pocket identity). "
              "Selectivity potential only at Trp198 allosteric site, which differs "
              "from human (PMC9486128).",
    ),

    "SmTGR": TargetPair(
        name="SmTGR",
        disease="Schistosomiasis",
        parasite_name="Thioredoxin glutathione reductase",
        human_name="Thioredoxin reductase 1 (human)",
        parasite_chembl="CHEMBL6110",
        human_chembl="CHEMBL3952",
        parasite_uniprot="Q86LC0",
        human_uniprot="Q16881",
        # GSH binding site from PDB 2X99 (SmTGR) / PDB 2ZZC (HsTrxR1)
        # Used in Kira Script 14 docking
        parasite_pocket_sequence="CVNVGCVIPAKVLRN",
        human_pocket_sequence="CVNVGCTIPAKILRN",  # Fusion enzyme = unique architecture
        notes="88% non-selective at GSH site by docking. SmTGR is a fusion "
              "enzyme (TrxR+GR) unique to parasite — selectivity potential "
              "is in the fusion junction, not the conserved active site.",
    ),

    # === TRYPANOSOMIASIS ===

    "TbCathB": TargetPair(
        name="TbCathB",
        disease="Trypanosomiasis",
        parasite_name="Cathepsin B-like cysteine protease",
        human_name="Cathepsin B (human)",
        parasite_chembl="CHEMBL612545",
        human_chembl="CHEMBL3837",
        parasite_uniprot="Q95PM0",
        human_uniprot="P07858",
        # Cysteine protease active site (Cys-His-Asn catalytic triad)
        # + S1/S2/S3 substrate binding subsites
        parasite_pocket_sequence="CHNWGLGAVVLAQEL",
        human_pocket_sequence="CHNWGLGGVVLAQEL",
        notes="94.4% non-selective. Cysteine protease active sites are "
              "extremely conserved across species.",
    ),

    "TbPDEB1": TargetPair(
        name="TbPDEB1",
        disease="Trypanosomiasis",
        parasite_name="Class 1 phosphodiesterase PDEB1",
        human_name="Phosphodiesterase (human)",
        parasite_chembl="CHEMBL613698",
        human_chembl="CHEMBL1952",  # Human PDE4
        parasite_uniprot="Q38F42",
        human_uniprot="Q07343",
        # PDE catalytic pocket — metal-binding residues + substrate pocket
        # PDE catalytic pocket — metal-binding + substrate site
        # Parasite PDEB1 has a unique P-pocket extension not in human PDE4
        # Catalytic residues conserved; P-pocket residues differ
        parasite_pocket_sequence="HDHDYTHFIMQTV",
        human_pocket_sequence="HDHDYTHFIMNAL",  # Q→N, T→A, V→L at P-pocket
        notes="10.3% selective. Parasite PDE has a unique P-pocket extension "
              "not found in human PDEs. 3 residue differences at P-pocket.",
    ),

    # === LEISHMANIASIS ===

    "LmPTR1": TargetPair(
        name="LmPTR1",
        disease="Leishmaniasis",
        parasite_name="Pteridine reductase 1",
        human_name="Dihydrofolate reductase (human)",
        parasite_chembl="CHEMBL1855",
        human_chembl="CHEMBL202",
        parasite_uniprot="Q01782",
        human_uniprot="P00374",
        # PTR1 active site: NADPH binding + pteridine substrate site
        # From PDB 1E7W (LmPTR1)
        # PTR1 is NOT a direct ortholog of DHFR — different fold entirely
        # This is why selectivity is so high: different enzyme, different pocket
        parasite_pocket_sequence="DSRFYNAKELVPKMF",
        human_pocket_sequence="QNLIVNAKELVPRTW",  # Very different — different enzyme family
        notes="75.6% selective, median 68.2x. BEST target across all 3 diseases. "
              "PTR1 and DHFR are different enzyme families that happen to "
              "reduce the same substrate. Different fold = different pocket = "
              "high selectivity. This is the gold standard.",
    ),

    "LmDHFR": TargetPair(
        name="LmDHFR",
        disease="Leishmaniasis",
        parasite_name="Bifunctional dihydrofolate reductase-thymidylate synthase",
        human_name="Dihydrofolate reductase (human)",
        parasite_chembl="CHEMBL2413",
        human_chembl="CHEMBL202",
        parasite_uniprot="P07382",
        human_uniprot="P00374",
        # DHFR folate binding site
        # From PDB structures of Leishmania DHFR-TS
        # DHFR folate binding site from PDB structures of Leishmania DHFR-TS
        # Same enzyme family as human DHFR — mostly conserved
        # But bifunctional TS domain creates interface differences
        parasite_pocket_sequence="DLGQNLIVNDKRFWP",
        human_pocket_sequence="DLGQNLIVNDRRFWP",  # K→R at one position (conservative)
        notes="42.0% selective. Same enzyme family as human DHFR (unlike PTR1). "
              "Selectivity comes from subtle active site differences + "
              "the bifunctional TS domain interaction.",
    ),
}


def get_target_pair(name: str) -> TargetPair:
    """Look up a target pair by short name."""
    if name not in TARGET_PAIRS:
        available = ", ".join(TARGET_PAIRS.keys())
        raise KeyError(f"Unknown target pair '{name}'. Available: {available}")
    return TARGET_PAIRS[name]


def get_pairs_for_disease(disease: str) -> list[TargetPair]:
    """Get all target pairs for a given disease."""
    return [p for p in TARGET_PAIRS.values() if p.disease == disease]


def list_all_pairs() -> list[str]:
    """List all target pair names."""
    return list(TARGET_PAIRS.keys())
