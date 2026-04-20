from kira.experiments.selectivity_v4_data import (
    PRIMARY_V4_PAIRS,
    annotate_v4_rows,
    build_v4_raw_table,
    find_data_dir,
    prepare_primary_v4_candidates,
    prepare_primary_v4_trainable,
    resolve_structures,
)


def _tables():
    data_dir = find_data_dir()
    raw = build_v4_raw_table(data_dir)
    df_all = annotate_v4_rows(raw)
    df_all = resolve_structures(df_all, data_dir)
    df_primary = prepare_primary_v4_candidates(df_all)
    df_trainable = prepare_primary_v4_trainable(df_all)
    return df_all, df_primary, df_trainable


def test_v4_row_counts_are_stable():
    df_all, df_primary, df_trainable = _tables()

    assert len(df_all) == 310
    assert len(df_primary) == 237
    assert len(df_trainable) == 235
    assert int(df_all["target_pair_id"].notna().sum()) == 310
    assert int(df_all["selectivity_ratio"].notna().sum()) == 237
    assert int((df_all["structure_status"] == "resolved").sum()) == 310


def test_v4_status_counts_are_stable():
    df_all, _, _ = _tables()

    assert df_all["row_status"].value_counts().to_dict() == {
        "resolved_candidate": 237,
        "missing_ratio": 73,
    }
    assert df_all["structure_status"].value_counts().to_dict() == {
        "resolved": 310,
    }


def test_primary_pair_counts_are_stable():
    _, df_primary, _ = _tables()

    assert df_primary["target_pair_id"].value_counts().to_dict() == {
        "TbPDEB1": 70,
        "LmDHFR": 69,
        "LmPTR1": 45,
        "TbCathB": 36,
        "SmDHODH": 10,
        "SmHDAC8": 7,
    }
    assert set(df_primary["target_pair_id"].unique()) == PRIMARY_V4_PAIRS


def test_trainable_pair_counts_are_stable():
    _, _, df_trainable = _tables()

    assert df_trainable["target_pair_id"].value_counts().to_dict() == {
        "LmDHFR": 69,
        "TbPDEB1": 68,
        "LmPTR1": 45,
        "TbCathB": 36,
        "SmDHODH": 10,
        "SmHDAC8": 7,
    }


def test_sm_tgr_is_excluded_from_primary_and_trainable():
    _, df_primary, df_trainable = _tables()

    assert "SmTGR" not in set(df_primary["target_pair_id"])
    assert "SmTGR" not in set(df_trainable["target_pair_id"])


def test_primary_rows_are_labeled_and_resolved():
    _, df_primary, _ = _tables()

    assert df_primary["selectivity_ratio"].notna().all()
    assert (df_primary["label_quality"] == "numeric_exact").all()
    assert (df_primary["structure_status"] == "resolved").all()
    assert (df_primary["row_status"] == "resolved_candidate").all()


def test_binary_label_matches_threshold():
    _, df_primary, _ = _tables()

    expected = (df_primary["selectivity_ratio"] >= 10.0).astype(int).tolist()
    actual = df_primary["is_selective"].astype(int).tolist()
    assert actual == expected


def test_duplicate_compound_pair_rows_are_explicit_and_excluded_from_trainable():
    df_all, df_primary, df_trainable = _tables()

    assert int(df_all["is_duplicate_compound_pair"].sum()) == 2
    assert int(df_primary["is_duplicate_compound_pair"].sum()) == 2
    assert int(df_trainable["is_duplicate_compound_pair"].sum()) == 0
