"""Tests for kira.drugs — curated SMILES and drug data."""

from kira.drugs import CURATED_SMILES, get_smiles


class TestCuratedSmiles:
    """Verify the CURATED_SMILES dictionary is well-formed."""

    def test_praziquantel_present(self):
        assert "PRAZIQUANTEL" in CURATED_SMILES

    def test_oxamniquine_present(self):
        assert "OXAMNIQUINE" in CURATED_SMILES

    def test_mefloquine_present(self):
        assert "MEFLOQUINE" in CURATED_SMILES

    def test_all_smiles_are_nonempty_strings(self):
        for name, smi in CURATED_SMILES.items():
            assert isinstance(smi, str), f"{name} SMILES is not a string"
            assert len(smi) > 0, f"{name} SMILES is empty"

    def test_at_least_four_drugs(self):
        assert len(CURATED_SMILES) >= 4

    def test_no_whitespace_only_smiles(self):
        for name, smi in CURATED_SMILES.items():
            assert smi.strip() == smi, f"{name} SMILES has leading/trailing whitespace"


class TestGetSmiles:
    def test_known_drug(self):
        smi = get_smiles("PRAZIQUANTEL")
        assert smi is not None
        assert len(smi) > 0

    def test_unknown_drug(self):
        assert get_smiles("NONEXISTENT_DRUG") is None

    def test_case_sensitive(self):
        # CURATED_SMILES keys are uppercase
        assert get_smiles("praziquantel") is None
