from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

PRIMARY_V4_PAIRS = {
    "SmDHODH",
    "SmHDAC8",
    "TbCathB",
    "TbPDEB1",
    "LmPTR1",
    "LmDHFR",
}

TARGET_NAME_MAP_V4: dict[str, str] = {
    # Trypanosomiasis
    "Cathepsin B-like cysteine protease": "TbCathB",
    "Class 1 phosphodiesterase PDEB1": "TbPDEB1",
    "Class I phosphodiesterase PDEB1": "TbPDEB1",
    # Leishmaniasis
    "Pteridine reductase 1": "LmPTR1",
    "Bifunctional dihydrofolate reductase-thymidylate synthase": "LmDHFR",
    # Schistosomiasis
    "Histone deacetylase 8": "SmHDAC8",
    "Dihydroorotate dehydrogenase (quinone), mitochondrial": "SmDHODH",
    "Dihydroorotate dehydrogenase (quino": "SmDHODH",
    "Thioredoxin glutathione reductase": "SmTGR",
}

SCHISTO_RELEVANT_TARGETS = {
    "Histone deacetylase 8",
    "Dihydroorotate dehydrogenase (quinone), mitochondrial",
    "Dihydroorotate dehydrogenase (quino",
    "Thioredoxin glutathione reductase",
}

STANDARD_COLUMNS = [
    "molecule_chembl_id",
    "compound_name",
    "parasite_ic50",
    "human_ic50",
    "selectivity_ratio",
    "parasite_target",
    "human_target",
    "human_target_id",
    "disease",
    "source_file",
    "context_source",
    "source_mw",
    "source_logp",
    "source_tpsa",
]


FINAL_COLUMNS = STANDARD_COLUMNS + [
    "target_pair_id",
    "is_selective",
    "label_quality",
    "structure_status",
    "smiles_source",
    "canonical_smiles",
    "inchikey",
    "murcko_scaffold",
    "row_status",
    "drop_reason",
    "compound_pair_key",
    "n_obs_for_compound_pair",
    "is_duplicate_compound_pair",
]

def find_data_dir() -> Path:
    """Find the Kira data directory."""
    env_data = os.environ.get("KIRA_DATA")
    if env_data:
        p = Path(env_data).expanduser().resolve()
        if p.is_dir():
            return p

    candidates = [
        Path(__file__).resolve().parents[3] / "data",
        Path.home() / "Kira" / "data",
        Path.home() / "kira" / "data",
        Path("data"),
    ]
    for p in candidates:
        if p.is_dir():
            return p

    raise FileNotFoundError(
        "Could not find Kira data directory. Run from repo root or set KIRA_DATA."
    )


def _safe_col(df: pd.DataFrame, name: str, default: object = pd.NA) -> pd.Series:
    """Return a column if present; otherwise a same-length default Series."""
    if name in df.columns:
        return df[name]
    return pd.Series([default] * len(df), index=df.index)


def _standardize_frame(
    df: pd.DataFrame,
    *,
    source_file: str,
    context_source: str,
    disease_default: str | None = None,
) -> pd.DataFrame:
    """Project a source-specific frame into the common v4 raw schema."""
    out = pd.DataFrame(
        {
            "molecule_chembl_id": _safe_col(df, "molecule_chembl_id"),
            "compound_name": _safe_col(df, "pref_name"),
            "parasite_ic50": _safe_col(df, "parasite_ic50"),
            "human_ic50": _safe_col(df, "human_ic50"),
            "selectivity_ratio": _safe_col(df, "selectivity_ratio"),
            "parasite_target": _safe_col(df, "parasite_target"),
            "human_target": _safe_col(df, "human_target"),
            "human_target_id": _safe_col(df, "human_target_id"),
            "disease": _safe_col(
                df, "disease", disease_default if disease_default else pd.NA
            ),
            "source_file": source_file,
            "context_source": context_source,
            "source_mw": _safe_col(df, "MW"),
            "source_logp": _safe_col(df, "LogP"),
            "source_tpsa": _safe_col(df, "TPSA"),
        }
    )
    return out


def load_leish_rows(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "leishmania" / "leish_selectivity.csv"
    if not path.exists():
        return pd.DataFrame(columns=STANDARD_COLUMNS)

    df = pd.read_csv(path)
    return _standardize_frame(
        df,
        source_file=path.name,
        context_source="leish_selectivity",
        disease_default="Leishmaniasis",
    )


def load_tryp_rows(data_dir: Path) -> pd.DataFrame:
    expanded = data_dir / "trypanosoma" / "tryp_selectivity_expanded.csv"
    fallback = data_dir / "trypanosoma" / "tryp_selectivity.csv"

    path = expanded if expanded.exists() else fallback
    if not path.exists():
        return pd.DataFrame(columns=STANDARD_COLUMNS)

    df = pd.read_csv(path)
    return _standardize_frame(
        df,
        source_file=path.name,
        context_source=path.stem,
        disease_default="Trypanosomiasis",
    )


def load_schisto_rows(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "publication" / "discovery_candidates.csv"
    if not path.exists():
        return pd.DataFrame(columns=STANDARD_COLUMNS)

    df = pd.read_csv(path).copy()

    if "best_target" in df.columns:
        df = df[df["best_target"].isin(SCHISTO_RELEVANT_TARGETS)].copy()

    df["parasite_target"] = _safe_col(df, "best_target")
    df["selectivity_ratio"] = _safe_col(df, "best_selectivity_ratio")
    df["compound_name"] = _safe_col(df, "pref_name")
    df["disease"] = "Schistosomiasis"

    out = _standardize_frame(
        df,
        source_file=path.name,
        context_source="discovery_candidates",
        disease_default="Schistosomiasis",
    )
    return out


def build_v4_raw_table(data_dir: Path | None = None) -> pd.DataFrame:
    """Load and concatenate all v4-relevant raw rows."""
    data_dir = data_dir or find_data_dir()

    frames = [
        load_leish_rows(data_dir),
        load_tryp_rows(data_dir),
        load_schisto_rows(data_dir),
    ]
    frames = [f for f in frames if len(f) > 0]

    if not frames:
        raise FileNotFoundError(f"No selectivity source files found under {data_dir}")

    raw = pd.concat(frames, ignore_index=True)

    for col in ["compound_name", "parasite_target", "human_target", "human_target_id"]:
        raw[col] = raw[col].replace("", pd.NA)

    return raw



def annotate_v4_rows(raw: pd.DataFrame) -> pd.DataFrame:
    """Add v4 benchmark bookkeeping fields without resolving structures yet."""
    df = raw.copy()

    for col in [
        "parasite_ic50",
        "human_ic50",
        "selectivity_ratio",
        "source_mw",
        "source_logp",
        "source_tpsa",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["target_pair_id"] = df["parasite_target"].map(TARGET_NAME_MAP_V4)

    df["is_selective"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
    has_ratio = df["selectivity_ratio"].notna()
    df.loc[has_ratio, "is_selective"] = (
        df.loc[has_ratio, "selectivity_ratio"] >= 10.0
    ).astype("Int64")

    df["label_quality"] = "missing_ratio"
    df.loc[has_ratio, "label_quality"] = "numeric_exact"

    df["structure_status"] = "unresolved"
    df["smiles_source"] = pd.NA
    df["canonical_smiles"] = pd.NA
    df["inchikey"] = pd.NA
    df["murcko_scaffold"] = pd.NA

    df["row_status"] = "pending_structure"
    df["drop_reason"] = pd.NA

    unmapped = df["target_pair_id"].isna()
    missing_ratio = df["selectivity_ratio"].isna()

    df.loc[unmapped, "row_status"] = "unmapped_pair"
    df.loc[unmapped, "drop_reason"] = "unmapped_pair"

    only_missing_ratio = (~unmapped) & missing_ratio
    df.loc[only_missing_ratio, "row_status"] = "missing_ratio"
    df.loc[only_missing_ratio, "drop_reason"] = "missing_ratio"

    df["compound_pair_key"] = (
        df["molecule_chembl_id"].fillna("NA").astype(str)
        + "|"
        + df["target_pair_id"].fillna("NA").astype(str)
    )

    df["n_obs_for_compound_pair"] = (
        df.groupby("compound_pair_key")["compound_pair_key"].transform("size").astype(int)
    )
    df["is_duplicate_compound_pair"] = df["n_obs_for_compound_pair"] > 1

    return df[FINAL_COLUMNS].copy()

def prepare_primary_v4_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Primary v4 candidate rows:
    - mapped to one of the six current benchmark pairs
    - have a numeric selectivity ratio
    - still pending structure resolution in this first pass
    """
    keep = df["target_pair_id"].isin(PRIMARY_V4_PAIRS) & df["selectivity_ratio"].notna()
    out = df.loc[keep].copy()
    return out
def _processed_dir(data_dir: Path) -> Path:
    processed = data_dir / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    return processed


def _load_chembl_cache(cache_path: Path) -> dict[str, dict[str, str | None]]:
    if not cache_path.exists():
        return {}
    try:
        return json.loads(cache_path.read_text())
    except Exception:
        return {}


def _save_chembl_cache(cache: dict[str, dict[str, str | None]], cache_path: Path) -> None:
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))


def _canonicalize_smiles(smiles: str | None) -> tuple[str | None, str | None, str | None]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None, None, None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None

    canonical = Chem.MolToSmiles(mol, canonical=True)
    inchikey = Chem.MolToInchiKey(mol)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or None
    return canonical, inchikey, scaffold


def _fetch_chembl_record(
    molecule_chembl_id: str,
    cache: dict[str, dict[str, str | None]],
) -> dict[str, str | None]:
    cid = str(molecule_chembl_id).strip()
    if not cid or cid.lower() == "nan":
        return {"canonical_smiles": None, "pref_name": None}

    if cid in cache:
        return cache[cid]

    try:
        rec = new_client.molecule.get(cid)
    except Exception:
        rec = None

    structures = (rec or {}).get("molecule_structures") or {}
    out = {
        "canonical_smiles": structures.get("canonical_smiles"),
        "pref_name": (rec or {}).get("pref_name"),
    }
    cache[cid] = out
    return out


def resolve_structures(df: pd.DataFrame, data_dir: Path | None = None) -> pd.DataFrame:
    """
    Resolve canonical structures primarily from ChEMBL ID, then canonicalize with RDKit.
    """
    data_dir = data_dir or find_data_dir()
    cache_path = _processed_dir(data_dir) / "chembl_smiles_cache.json"
    cache = _load_chembl_cache(cache_path)

    out = df.copy()

    for idx, row in out.iterrows():
        chembl_id = row["molecule_chembl_id"]
        compound_name = row["compound_name"]

        smiles = None
        smiles_source = None

        # Primary route: ChEMBL ID
        if pd.notna(chembl_id):
            rec = _fetch_chembl_record(str(chembl_id), cache)
            smiles = rec.get("canonical_smiles")
            if pd.isna(compound_name) and rec.get("pref_name"):
                out.at[idx, "compound_name"] = rec["pref_name"]
            if smiles:
                smiles_source = "chembl_id"

        # Fallback route: curated local name map
        if not smiles and pd.notna(out.at[idx, "compound_name"]):
            try:
                from kira.drugs import get_smiles

                fallback = get_smiles(str(out.at[idx, "compound_name"]).upper())
                if fallback:
                    smiles = fallback
                    smiles_source = "curated_drugs"
            except Exception:
                pass

        canonical, inchikey, scaffold = _canonicalize_smiles(smiles)

        if canonical:
            out.at[idx, "structure_status"] = "resolved"
            out.at[idx, "smiles_source"] = smiles_source
            out.at[idx, "canonical_smiles"] = canonical
            out.at[idx, "inchikey"] = inchikey
            out.at[idx, "murcko_scaffold"] = scaffold

            if out.at[idx, "row_status"] == "pending_structure":
                out.at[idx, "row_status"] = "resolved_candidate"

        elif smiles:
            out.at[idx, "structure_status"] = "bad_smiles"
            out.at[idx, "smiles_source"] = smiles_source
            if out.at[idx, "row_status"] == "pending_structure":
                out.at[idx, "row_status"] = "bad_smiles"

        else:
            if out.at[idx, "row_status"] == "pending_structure":
                out.at[idx, "row_status"] = "missing_smiles"

    _save_chembl_cache(cache, cache_path)
    return out


def prepare_primary_v4_trainable(df: pd.DataFrame) -> pd.DataFrame:
    """
    First-pass trainable set:
    - six primary benchmark pairs
    - numeric selectivity ratio
    - resolved structure
    - excludes repeated compound_pair_key rows
    """
    keep = (
        df["target_pair_id"].isin(PRIMARY_V4_PAIRS)
        & df["selectivity_ratio"].notna()
        & (df["structure_status"] == "resolved")
        & (~df["is_duplicate_compound_pair"])
    )
    return df.loc[keep].copy()

def _counts(series: pd.Series) -> dict[str, int]:
    vc = series.value_counts(dropna=False)
    return {str(k): int(v) for k, v in vc.items()}



def make_summary(
    df_all: pd.DataFrame,
    df_primary: pd.DataFrame,
    df_trainable: pd.DataFrame,
) -> dict:
    summary = {
        "rows_all": int(len(df_all)),
        "rows_primary_candidate": int(len(df_primary)),
        "rows_primary_trainable": int(len(df_trainable)),
        "mapped_rows_all": int(df_all["target_pair_id"].notna().sum()),
        "labeled_rows_all": int(df_all["selectivity_ratio"].notna().sum()),
        "resolved_rows_all": int((df_all["structure_status"] == "resolved").sum()),
        "resolved_rows_primary_candidate": int(
            (df_primary["structure_status"] == "resolved").sum()
        ),
        "rows_by_status": _counts(df_all["row_status"]),
        "rows_by_structure_status": _counts(df_all["structure_status"]),
        "rows_by_disease_all": _counts(df_all["disease"]),
        "rows_by_pair_all": _counts(df_all["target_pair_id"].fillna("NA")),
        "rows_by_pair_primary": _counts(df_primary["target_pair_id"]),
        "rows_by_pair_trainable": _counts(df_trainable["target_pair_id"]),
        "duplicate_compound_pair_rows_all": int(df_all["is_duplicate_compound_pair"].sum()),
        "duplicate_compound_pair_rows_primary": int(df_primary["is_duplicate_compound_pair"].sum()),
    }
    return summary


def main() -> None:
    data_dir = find_data_dir()
    processed_dir = _processed_dir(data_dir)

    raw = build_v4_raw_table(data_dir)
    df_all = annotate_v4_rows(raw)
    df_all = resolve_structures(df_all, data_dir)

    df_primary = prepare_primary_v4_candidates(df_all)
    df_trainable = prepare_primary_v4_trainable(df_all)

    summary = make_summary(df_all, df_primary, df_trainable)

    all_path = processed_dir / "selectivity_v4_rows_all.csv"
    primary_path = processed_dir / "selectivity_v4_rows_primary_candidate.csv"
    trainable_path = processed_dir / "selectivity_v4_rows_primary_trainable.csv"
    summary_path = processed_dir / "selectivity_v4_summary.json"

    df_all.to_csv(all_path, index=False)
    df_primary.to_csv(primary_path, index=False)
    df_trainable.to_csv(trainable_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2))

    print("=" * 80)
    print("selectivity_v4_data: summary")
    print("=" * 80)
    print(json.dumps(summary, indent=2))
    print("\nWrote:")
    print(f"  - {all_path}")
    print(f"  - {primary_path}")
    print(f"  - {trainable_path}")
    print(f"  - {summary_path}")

if __name__ == "__main__":
    main()
