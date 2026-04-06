"""Drug data — curated SMILES, WHO Essential Medicines, stage mappings.

Contains reference data for known drugs relevant to NTD repurposing,
including manually curated SMILES strings for compounds that may not
have structural data in ChEMBL activity exports.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Curated SMILES for known drugs
# ---------------------------------------------------------------------------
# Standard SMILES representations from ChEMBL/PubChem for drugs that
# appear in the pipeline but may lack canonical_smiles in activity data.

CURATED_SMILES: dict[str, str] = {
    # Anti-schistosomal / anti-parasitic drugs
    "PRAZIQUANTEL": "O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1",
    "OXAMNIQUINE": "CC(CO)Nc1ccc2c(c1)[C@@H](C)C[C@H](C)N2O",
    "MEFLOQUINE": "OC(c1cc(C(F)(F)F)nc2c(C(F)(F)F)cccc12)C1CCCCN1",
    "ATOVAQUONE": "O=C1C(=O)c2ccccc2C(=O)C1C1CCC(c2ccc(Cl)cc2)CC1",

    # Repurposing candidates (from WHO Essential Medicines List)
    "METFORMIN": "CN(C)C(=N)NC(=N)N",
    "LISINOPRIL": (
        "NCCCC[C@@H](N[C@@H](CCc1ccccc1)C(=O)O)"
        "C(=O)N1CCCC1C(=O)O"
    ),
    "SERTRALINE": "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21",
    "OMEPRAZOLE": "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    "ATORVASTATIN": (
        "CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(=O)O)"
        "c(-c2ccccc2)c(-c2ccc(F)cc2)c1C(=O)Nc1ccccc1"
    ),
    "AMLODIPINE": (
        "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl"
    ),
    "LEVOTHYROXINE": (
        "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O"
    ),
    "MONTELUKAST": (
        "CC(C)(O)c1ccccc1CC[C@@H](SCC1(CC(=O)O)CC1)"
        "c1cccc(-c2cccc3ccc(Cl)cc23)c1"
    ),
    "CLOPIDOGREL": "COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1",
    "TAMSULOSIN": (
        "CCOc1ccc(CC(C)NCC[C@@H](O)c2ccc(OC)"
        "c(S(N)(=O)=O)c2)cc1"
    ),
    "ESCITALOPRAM": "N#Cc1ccc2c(c1)C(CCCCN1CCC1)(OC2)c1ccc(F)cc1",
    "GABAPENTIN": "NCC1(CC(=O)O)CCCCC1",
    "PANTOPRAZOLE": (
        "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC"
    ),
    "ROSUVASTATIN": (
        "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)"
        "c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O"
    ),
    "TRAMADOL": "COc1cccc(C2(O)CCCCC2CN(C)C)c1",
}


def get_smiles(drug_name: str) -> str | None:
    """Look up a curated SMILES string by drug name.

    Parameters
    ----------
    drug_name : str
        Drug name (case-sensitive, uppercase expected).

    Returns
    -------
    str or None
        SMILES string if found, None otherwise.
    """
    return CURATED_SMILES.get(drug_name)
