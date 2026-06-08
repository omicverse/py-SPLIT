from __future__ import annotations

import numpy as np
from scipy import sparse

import splitst as split


def _manual_split(counts, weights, reference, primary):
    counts_arr = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    w = weights.to_numpy(dtype=float)
    ref = reference.to_numpy(dtype=float)
    pidx = reference.index.get_indexer(primary.to_numpy())
    primary_weight = w[np.arange(w.shape[0]), pidx]
    n_types = np.maximum((w > 0).sum(axis=1), 1)
    denom = w @ ref + 1e-10
    numer = primary_weight[:, None] * ref[pidx] + 1e-10 / n_types[:, None]
    return counts_arr * numer / denom


def test_purify_matches_reference_formula_dense():
    adata, weights, reference = split.create_split_fixture(n_cells=36, n_genes=18, n_cell_types=3)

    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])

    expected = _manual_split(adata.layers["counts"], weights, reference, adata.obs["primary_cell_type"])
    assert np.allclose(adata.layers["split_purified"], expected, rtol=1e-10, atol=1e-10)
    assert set(adata.obs["purification_status"]).issubset({"raw", "purified"})


def test_purify_accepts_sparse_and_infers_primary():
    dense, weights, reference = split.create_split_fixture(n_cells=36, n_genes=18, n_cell_types=3)
    sparse_adata = dense.copy()
    sparse_adata.X = sparse.csr_matrix(sparse_adata.X)
    sparse_adata.layers["counts"] = sparse.csr_matrix(sparse_adata.layers["counts"])

    split.purify(dense, weights, reference, primary_cell_type=dense.obs["primary_cell_type"])
    split.purify(sparse_adata, weights, reference)

    assert sparse.issparse(sparse_adata.layers["split_purified"])
    assert np.allclose(sparse_adata.layers["split_purified"].toarray(), dense.layers["split_purified"])


def test_cells_to_purify_keeps_unselected_raw():
    adata, weights, reference = split.create_split_fixture(n_cells=36, n_genes=18, n_cell_types=3)
    selected = adata.obs_names[:10]

    split.purify(
        adata,
        weights,
        reference.T,
        primary_cell_type=adata.obs["primary_cell_type"],
        cells_to_purify=selected,
    )

    assert np.allclose(adata.layers["split_purified"][10:], adata.layers["counts"][10:])
    assert set(adata.obs.iloc[10:]["purification_status"]) == {"raw"}


def test_spatial_score_writes_diffusion_metrics():
    adata, weights, _ = split.create_split_fixture(n_cells=48, n_genes=20, n_cell_types=4)

    split.spatial_score(
        adata,
        weights,
        primary_cell_type=adata.obs["primary_cell_type"],
        secondary_cell_type=adata.obs["secondary_cell_type"],
        k=5,
        radius=8.0,
    )

    scores = adata.obs["neighborhood_weights_second_type"].to_numpy()
    assert "split_spatial_neighbors" in adata.uns
    assert np.isfinite(scores).all()
    assert ((scores >= 0) & (scores <= 1)).all()
    assert adata.uns["split_spatial_neighbors"]["indices"].shape[0] == adata.n_obs


def test_balance_selects_purified_profiles_by_score():
    adata, weights, reference = split.create_split_fixture(n_cells=36, n_genes=18, n_cell_types=3)

    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])
    adata.obs["neighborhood_weights_second_type"] = 0.0
    adata.obs.iloc[:5, adata.obs.columns.get_loc("neighborhood_weights_second_type")] = 0.5
    split.balance(adata, threshold=0.15)

    balanced = adata.layers["split_balanced"].toarray()
    raw = np.asarray(adata.layers["counts"])
    purified = np.asarray(adata.layers["split_purified"])
    assert np.allclose(balanced[:5], purified[:5])
    assert np.allclose(balanced[5:], raw[5:])
    assert list(adata.obs["split_balance_status"][:5]) == ["purified"] * 5
    assert set(adata.obs["split_balance_status"][5:]) == {"raw"}


def test_split_shift_marks_label_swaps():
    adata, weights, reference = split.create_split_fixture(n_cells=36, n_genes=18, n_cell_types=3)

    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])
    adata.obs["neighborhood_weights_second_type"] = 1.0
    adata.obs["first_type_neighborhood"] = adata.obs["second_type"].astype(str)
    split.balance(adata, threshold=0.15, swap_labels=True)

    assert adata.obs["split_shift_swap"].any()


def test_reassign_residuals_preserves_nonnegative_counts_and_operator():
    adata, weights, reference = split.create_split_fixture(n_cells=48, n_genes=20, n_cell_types=4)
    adata.X = sparse.csr_matrix(adata.X)
    adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])

    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])
    adata.obs["neighborhood_weights_second_type"] = 1.0
    split.balance(adata, threshold=0.15)
    split.reassign_residuals(adata, k=5, mode="uniform")

    reassigned = adata.layers["split_reassigned"]
    assert sparse.issparse(reassigned)
    assert reassigned.data.min() >= 0
    op = adata.uns["split_reassignment_operator"]
    row_sums = np.asarray(op.sum(axis=1)).ravel()
    assert np.all(row_sums[row_sums > 0] <= 1.0 + 1e-12)


def test_purification_improves_distance_to_clean_target():
    adata, weights, reference = split.create_split_fixture(n_cells=48, n_genes=20, n_cell_types=4)

    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])

    contaminated = np.asarray(adata.layers["counts"])
    purified = np.asarray(adata.layers["split_purified"])
    clean = np.asarray(adata.layers["clean_counts"])
    assert np.linalg.norm(purified - clean) < np.linalg.norm(contaminated - clean)

