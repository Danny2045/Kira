import numpy as np
import pandas as pd

from kira.experiments.selectivity_v4_features import (
    block_slices,
    build_feature_names,
    build_feature_vector,
    build_v4_blocks,
    featurize_dataframe,
    select_feature_blocks,
)


def test_different_compounds_same_pair_change_full_vector():
    v1 = build_feature_vector("CCO", "LmDHFR", n_bits=64)
    v2 = build_feature_vector("c1ccccc1", "LmDHFR", n_bits=64)
    assert not np.array_equal(v1, v2)


def test_same_compound_different_pairs_change_full_vector():
    v1 = build_feature_vector("CCO", "LmDHFR", n_bits=64)
    v2 = build_feature_vector("CCO", "TbCathB", n_bits=64)
    assert not np.array_equal(v1, v2)


def test_pair_delta_block_is_15d():
    blocks = build_v4_blocks("CCO", "LmDHFR", n_bits=64)
    assert blocks["pair_delta"].shape == (15,)


def test_compatibility_block_changes_across_pairs():
    b1 = build_v4_blocks("CCO", "LmDHFR", n_bits=64)
    b2 = build_v4_blocks("CCO", "TbCathB", n_bits=64)
    assert not np.array_equal(b1["compat_diff"], b2["compat_diff"])


def test_feature_name_count_matches_vector_width():
    names = build_feature_names(n_bits=64)
    vec = build_feature_vector("CCO", "LmDHFR", n_bits=64)
    assert len(names) == len(vec)


def test_block_slices_cover_full_vector():
    slices = block_slices(n_bits=64)
    vec = build_feature_vector("CCO", "LmDHFR", n_bits=64)

    total = sum(s.stop - s.start for s in slices.values())
    assert total == len(vec)
    assert slices["pair_delta"].stop - slices["pair_delta"].start == 15


def test_same_pair_different_compounds_share_pair_block_but_not_compound_block():
    b1 = build_v4_blocks("CCO", "LmDHFR", n_bits=64)
    b2 = build_v4_blocks("c1ccccc1", "LmDHFR", n_bits=64)

    assert np.array_equal(b1["pair_delta"], b2["pair_delta"])
    assert not np.array_equal(b1["compound_fp"], b2["compound_fp"])


def test_featurize_dataframe_shape():
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CCO", "c1ccccc1"],
            "target_pair_id": ["LmDHFR", "TbCathB"],
        }
    )
    X, names = featurize_dataframe(df, n_bits=64)
    assert X.shape == (2, len(names))


def test_select_feature_blocks_returns_expected_width():
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CCO", "c1ccccc1"],
            "target_pair_id": ["LmDHFR", "TbCathB"],
        }
    )
    X, _ = featurize_dataframe(df, n_bits=64)

    X_pair = select_feature_blocks(X, ["pair_delta"], n_bits=64)
    X_compound = select_feature_blocks(X, ["compound_fp", "compound_desc"], n_bits=64)

    assert X_pair.shape == (2, 15)
    assert X_compound.shape[1] == 64 + 9
