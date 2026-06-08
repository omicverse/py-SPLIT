<p align="center">
  <img src="data/logo.svg" width="360px" alt="SPLIT-ST logo">
</p>

<div align="center">

<table>
<tr>
  <td align="right"><b>Package</b></td>
  <td><img src="https://img.shields.io/badge/status-alpha-orange" alt="Status"> <img src="https://img.shields.io/badge/python-3.9%2B-blue" alt="Python"></td>
</tr>
<tr>
  <td align="right"><b>Meta</b></td>
  <td><a href="LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-green" alt="License"></a> <a href="https://github.com/omicverse/py-SPLIT"><img src="https://img.shields.io/github/stars/omicverse/py-SPLIT?style=social" alt="Stars"></a></td>
</tr>
</table>

</div>

---

# SPLIT-ST

A **pure-Python re-implementation of SPLIT** for spatial transcriptomics purification.

- AnnData-native — works directly with cell x gene spatial matrices
- No `rpy2`, no R install, no RCTD/Seurat object dependency
- Generic deconvolution input — accepts any cell-by-cell-type weight matrix and cell-type reference profile
- Full tutorial workflow — purification, spatially aware selective purification, SPLIT-shift metadata swaps, and residual transcript reassignment

## Install

From GitHub:

```bash
pip install git+https://github.com/omicverse/py-SPLIT.git
```

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

| AnnData field | Contents |
|---|---|
| `adata.layers["split_purified"]` | Default SPLIT-purified counts |
| `adata.obs["neighborhood_weights_second_type"]` | Spatial second-type diffusion score |
| `adata.layers["split_balanced"]` | Raw/purified balanced count matrix |
| `adata.obs["split_shift_swap"]` | SPLIT-shift label swap flag |
| `adata.layers["split_reassigned"]` | Residual-reassigned count matrix |
| `adata.uns["split_reassignment_operator"]` | Sparse residual reassignment operator |

---

## Pipeline Overview

The `splitst` pipeline mirrors the SPLIT tutorial workflow step-for-step:

### 1. Default purification — `purify`

Purifies each cell's observed expression using deconvolution weights, a cell-type reference profile, and the primary cell type. The implementation accepts reference matrices in either cell-types x genes or genes x cell-types orientation and aligns by labels.

### 2. Spatial support scoring — `spatial_score`

Computes local support for the focal primary and secondary cell types using spatial KNN neighborhoods. The main selective purification score is `neighborhood_weights_second_type`.

### 3. Selective purification — `balance`

Keeps raw profiles for cells with low local second-type support and replaces high-score cells with SPLIT-purified profiles. Optional spot-class metadata can force uncertain doublets or rejected cells.

### 4. SPLIT-shift metadata correction

When `balance(..., swap_labels=True)` is used, cells whose neighborhood supports the secondary label can swap `first_type` and `second_type` metadata fields.

### 5. Residual transcript reassignment — `reassign_residuals`

Computes positive residuals from `raw - purified` and redistributes them to spatial neighbors using uniform or count-proportional weights.

---

## Algorithmic Fidelity To R SPLIT

The Python port follows the tutorial-level SPLIT workflow rather than binding to R-specific S4/Seurat/RCTD objects.

### 1. Formula parity

The purification formula is tested against an explicit NumPy implementation:

```text
purified = counts * (w_primary * reference_primary + epsilon / n_types) /
           (weights @ reference + epsilon)
```

### 2. Annotation-tool agnostic inputs

R tutorial objects are normalized into three generic inputs:

- `deconvolution_weights`: cells x cell types
- `reference`: cell types x genes, or genes x cell types
- `primary_cell_type`: labels or `None` for max-weight inference

### 3. Sparse and chunked execution

Dense and sparse matrices share the same semantics. Sparse inputs return sparse purified/reassigned outputs; large inputs can use `chunk_size`.

---

## Benchmarks

Run the Python workflow benchmark:

```bash
python benchmarks/benchmark_split.py
```

Current small fixture result:

| Metric | Value |
|---|---:|
| Cells x genes x cell types | 240 x 120 x 4 |
| Purification time | ~0.02 s |
| Full workflow time | ~0.19 s |
| Raw L2 distance to clean | 387.84 |
| Purified L2 distance to clean | 3.08e-09 |
| Minimum reassigned count | > 0 |

R parity scripts are tracked in `scripts/` and will be wired once the benchmark environment has R and the original SPLIT package installed.

---

## Testing

Run the full Python test suite:

```bash
pytest
```

Current local gate:

```text
29 passed
```

The tests cover dense/sparse parity, label alignment, copy mode, raw/purified balancing, SPLIT-shift swaps, residual reassignment mass accounting, no-neighbor spatial scoring, and invalid input rejection.

---

## Notebooks

| Notebook | What it covers |
|---|---|
| [`examples/split_workflow.ipynb`](examples/split_workflow.ipynb) | Quick-start workflow — fixture, purification, spatial scoring, SPLIT-shift, residual reassignment, benchmark table |

---

## API Reference

### Core functions

```python
from splitst import (
    purify,
    spatial_score,
    balance,
    reassign_residuals,
    create_split_fixture,
)
```

### `purify(adata, deconvolution_weights, reference, ...)`

Run default SPLIT purification and write a result layer, by default `split_purified`.

### `spatial_score(adata, deconvolution_weights, primary_cell_type, ...)`

Compute KNN/radius neighborhood support scores for primary and secondary cell types.

### `balance(adata, threshold=0.15, ...)`

Merge raw and purified profiles using `neighborhood_weights_second_type`.

### `reassign_residuals(adata, mode="count_proportional", ...)`

Redistribute positive residual transcripts from purified cells to spatial neighbors.

### Compatibility aliases

```python
from splitst import (
    split_purify,
    split_spatial_score,
    split_balance,
    split_reassign_residuals,
)
```

---

## Reconstruction Notes

- [`DISCOVERY.md`](DISCOVERY.md): dependency reuse audit and Python ecosystem mapping
- [`MATH.md`](MATH.md): SPLIT formula and workflow equations
- [`RECONSTRUCTION_REPORT.md`](RECONSTRUCTION_REPORT.md): parity, coverage, limitations, and integration plan
- [`ITERATION_LOG.md`](ITERATION_LOG.md): implementation and validation log
- [`NAMESPACE_PARITY.md`](NAMESPACE_PARITY.md): R tutorial concept to Python API mapping

## Citation

If you use this package, please cite the original SPLIT work and acknowledge this repository for the Python port.

## License

GNU GPLv3 — matching the upstream SPLIT license.
