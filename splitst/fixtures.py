"""Deterministic fixtures for SPLIT tests, examples, and benchmarks."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def create_split_fixture(
    n_cells: int = 240,
    n_genes: int = 120,
    n_cell_types: int = 4,
    contamination_rate: float = 0.35,
    random_state: int = 7,
):
    """Create a spatial-like AnnData fixture with known SPLIT inputs.

    Returns
    -------
    adata
        AnnData with `counts`, `clean_counts`, `spatial`, primary/secondary labels,
        and simulated contamination flags.
    weights
        Cell-by-cell-type deconvolution weights.
    reference
        Cell-type-by-gene reference expression profiles.
    """
    rng = np.random.default_rng(random_state)
    cell_types = pd.Index([f"CellType_{i + 1}" for i in range(n_cell_types)])
    genes = pd.Index([f"Gene_{i + 1:04d}" for i in range(n_genes)])
    cells = pd.Index([f"Cell_{i + 1:04d}" for i in range(n_cells)])

    reference = rng.negative_binomial(n=8, p=0.35, size=(n_cell_types, n_genes)).astype(float)
    for ct in range(n_cell_types):
        start = (ct * max(1, n_genes // n_cell_types)) % n_genes
        marker_idx = np.arange(start, min(start + max(5, n_genes // 8), n_genes))
        reference[ct, marker_idx] += rng.uniform(10, 25, size=len(marker_idx))
    reference = pd.DataFrame(reference, index=cell_types, columns=genes)

    primary_codes = np.arange(n_cells) % n_cell_types
    rng.shuffle(primary_codes)
    primary = pd.Series(cell_types[primary_codes].to_numpy(), index=cells, name="primary_cell_type")
    secondary = primary.map(lambda ct: cell_types[(cell_types.get_loc(ct) + 1) % len(cell_types)])
    suspicious = rng.random(n_cells) < contamination_rate

    weights = pd.DataFrame(0.0, index=cells, columns=cell_types)
    for cell in cells:
        weights.loc[cell, primary.loc[cell]] = 1.0
    for cell in cells[suspicious]:
        weights.loc[cell, primary.loc[cell]] = 0.8
        weights.loc[cell, secondary.loc[cell]] = 0.2

    raw_counts = weights.to_numpy() @ reference.to_numpy()
    primary_idx = reference.index.get_indexer(primary.to_numpy())
    primary_weight = weights.to_numpy()[np.arange(n_cells), primary_idx]
    clean_counts = primary_weight[:, None] * reference.to_numpy()[primary_idx]

    angles = np.linspace(0, 2 * np.pi, n_cell_types, endpoint=False)
    centers = {ct: np.array([np.cos(a), np.sin(a)]) * 10 for ct, a in zip(cell_types, angles)}
    coords = np.vstack([centers[ct] + rng.normal(scale=1.1, size=2) for ct in primary])
    for i, cell in enumerate(cells):
        if suspicious[i]:
            coords[i] = 0.45 * centers[primary.loc[cell]] + 0.55 * centers[secondary.loc[cell]] + rng.normal(scale=0.7, size=2)

    adata = ad.AnnData(
        X=raw_counts.copy(),
        obs=pd.DataFrame(
            {
                "primary_cell_type": primary,
                "secondary_cell_type": secondary,
                "simulated_contamination": suspicious,
            },
            index=cells,
        ),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts"] = raw_counts.copy()
    adata.layers["clean_counts"] = clean_counts
    adata.obsm["spatial"] = coords
    return adata, weights, reference

