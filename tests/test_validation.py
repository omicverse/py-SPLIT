from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import splitst as split


def test_non_unit_weights_are_normalized_with_warning():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    scaled = weights * 2

    with pytest.warns(RuntimeWarning, match="do not sum to 1"):
        split.purify(adata, scaled, reference, primary_cell_type=adata.obs["primary_cell_type"])

    assert "split_purified" in adata.layers


@pytest.mark.parametrize(
    "mutator,error",
    [
        (lambda adata, weights, reference: weights.__setitem__(weights.columns[0], -1.0), "negative"),
        (lambda adata, weights, reference: weights.__setitem__(weights.columns[0], np.nan), "finite"),
        (lambda adata, weights, reference: reference.__setitem__(reference.columns[0], -1.0), "negative"),
        (lambda adata, weights, reference: reference.__setitem__(reference.columns[0], np.inf), "finite"),
    ],
)
def test_purify_rejects_invalid_weights_and_reference(mutator, error):
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    weights = weights.copy()
    reference = reference.copy()
    mutator(adata, weights, reference)

    with pytest.raises(ValueError, match=error):
        split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])


def test_purify_rejects_negative_counts():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    adata.layers["counts"][0, 0] = -1

    with pytest.raises(ValueError, match="negative counts"):
        split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])


def test_purify_rejects_unknown_primary_cell_type():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    primary = adata.obs["primary_cell_type"].copy()
    primary.iloc[0] = "missing_type"

    with pytest.raises(ValueError, match="absent from weights"):
        split.purify(adata, weights, reference, primary_cell_type=primary)


def test_purify_rejects_unaligned_cells_to_purify():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)

    with pytest.raises(ValueError, match="cells absent"):
        split.purify(
            adata,
            weights,
            reference,
            primary_cell_type=adata.obs["primary_cell_type"],
            cells_to_purify=["not_a_cell"],
        )


def test_purify_rejects_bad_chunk_size():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)

    with pytest.raises(ValueError, match="chunk_size"):
        split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"], chunk_size=0)


def test_spatial_score_rejects_bad_spatial_inputs():
    adata, weights, _ = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)

    with pytest.raises(ValueError, match="k"):
        split.spatial_score(adata, weights, adata.obs["primary_cell_type"], k=0)

    with pytest.raises(ValueError, match="radius"):
        split.spatial_score(adata, weights, adata.obs["primary_cell_type"], radius=-1)

    adata_bad = adata.copy()
    adata_bad.obsm["spatial"][0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        split.spatial_score(adata_bad, weights, adata_bad.obs["primary_cell_type"])


def test_spatial_score_rejects_unknown_secondary_cell_type():
    adata, weights, _ = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    secondary = adata.obs["secondary_cell_type"].copy()
    secondary.iloc[0] = "missing_type"

    with pytest.raises(ValueError, match="secondary_cell_type"):
        split.spatial_score(adata, weights, adata.obs["primary_cell_type"], secondary_cell_type=secondary)


def test_balance_rejects_missing_or_bad_score():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])

    with pytest.raises(KeyError, match="neighborhood_weights_second_type"):
        split.balance(adata)

    adata.obs["neighborhood_weights_second_type"] = 0.0
    adata.obs.iloc[0, adata.obs.columns.get_loc("neighborhood_weights_second_type")] = np.nan
    with pytest.raises(ValueError, match="finite"):
        split.balance(adata)

    with pytest.raises(ValueError, match="threshold"):
        split.balance(adata, threshold=np.inf)


def test_reassign_residuals_rejects_bad_arguments_even_without_senders():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])
    adata.obs["neighborhood_weights_second_type"] = 0.0
    split.balance(adata, threshold=0.15)

    with pytest.raises(ValueError, match="mode"):
        split.reassign_residuals(adata, mode="bad")

    with pytest.raises(ValueError, match="self_keep"):
        split.reassign_residuals(adata, self_keep=1.5)


def test_dataframe_alignment_is_by_labels_not_order():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    shuffled_weights = weights.sample(frac=1.0, random_state=3)
    shuffled_reference = reference.sample(frac=1.0, axis=1, random_state=4)
    shuffled_reference = shuffled_reference.sample(frac=1.0, random_state=5)

    split.purify(adata, shuffled_weights, shuffled_reference, primary_cell_type=adata.obs["primary_cell_type"])
    expected = split.purify(
        adata.copy(),
        weights,
        reference,
        primary_cell_type=adata.obs["primary_cell_type"],
        copy=True,
    )

    assert np.allclose(adata.layers["split_purified"], expected.layers["split_purified"])


def test_numpy_inputs_are_accepted_by_position():
    adata, weights, reference = split.create_split_fixture(n_cells=24, n_genes=12, n_cell_types=3)
    split.purify(
        adata,
        weights.to_numpy(),
        reference.to_numpy(),
        primary_cell_type=None,
    )

    assert "split_purified" in adata.layers
