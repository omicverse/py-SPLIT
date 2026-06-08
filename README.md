# SPLIT-ST

Pure-Python SPLIT for spatial transcriptomics purification.

`py-SPLIT` is a standalone, AnnData-native reimplementation of the SPLIT workflow. It accepts generic deconvolution inputs instead of RCTD/Seurat-specific objects:

- spatial expression counts in `AnnData`
- cell-by-cell-type deconvolution weights
- cell-type-by-gene reference profiles
- primary cell-type labels, or inferred primary labels from the max weight

The initial workflow mirrors the SPLIT tutorial logic:

1. default SPLIT purification
2. spatially aware second-type scoring
3. selective raw/purified balancing
4. SPLIT-shift metadata swaps
5. residual transcript reassignment

## Install

From a local checkout:

```bash
pip install -e ".[dev]"
```

## Quick Start

```python
import splitst as split

adata, weights, reference = split.create_split_fixture()

split.purify(
    adata,
    deconvolution_weights=weights,
    reference=reference,
    primary_cell_type=adata.obs["primary_cell_type"],
)

split.spatial_score(
    adata,
    deconvolution_weights=weights,
    primary_cell_type=adata.obs["primary_cell_type"],
    secondary_cell_type=adata.obs["secondary_cell_type"],
)

split.balance(adata, threshold=0.15)
split.reassign_residuals(adata)
```

Results are written back to AnnData:

- `adata.layers["split_purified"]`
- `adata.obs["neighborhood_weights_second_type"]`
- `adata.layers["split_balanced"]`
- `adata.layers["split_reassigned"]`
- `adata.uns["split_reassignment_operator"]`

## API

```python
split.purify(...)
split.spatial_score(...)
split.balance(...)
split.reassign_residuals(...)
```

Compatibility aliases are also exported:

```python
split.split_purify(...)
split.split_spatial_score(...)
split.split_balance(...)
split.split_reassign_residuals(...)
```

## Benchmarks

Run the Python workflow benchmark:

```bash
python benchmarks/benchmark_split.py
```

Run tests:

```bash
pytest
```

R parity scripts are kept under `scripts/` for environments with R and the original SPLIT package installed.

## Relationship To OmicVerse

This repository develops SPLIT as a standalone Python package first. After the API and parity benchmarks settle, `omicverse.space.split_*` can wrap this package in the main `omicverse` repository.

