# DISCOVERY.md - Dependency Reuse Audit

## Target Package

- **R Package**: SPLIT
- **Description**: Spatial transcript purification for contamination-aware spatial transcriptomics workflows
- **License**: GPL-3
- **Repository**: https://github.com/bdsc-tds/SPLIT

## Existing py- Mirror Check

- **py-SPLIT exists in omicverse?**: No - this is a new port.
- **Checked via**: Manual inspection of the `omicverse` organization's py-* reconstruction repositories.

## R Dependency Audit

| R Dependency / Object | Python Equivalent | Reuse Decision | Notes |
|---|---|---|---|
| Matrix | scipy.sparse | **hard dep** | Sparse matrices and residual operators |
| base matrix algebra | numpy | **hard dep** | Core purification formula |
| data.frame / tibble | pandas | **hard dep** | Weights, reference profiles, metadata |
| RCTD output objects | pandas / AnnData inputs | **normalize** | Core package accepts generic weights and reference |
| Seurat objects | AnnData | **skip direct support** | Convert upstream before calling `splitst` |
| spatial KNN utilities | scipy.spatial.cKDTree | **hard dep** | Fast KNN and radius pruning |
| ggplot2 plotting | matplotlib | **optional example dep** | Used only in notebooks |

## Ecosystem Savings

- AnnData replaces R object-specific S4/Seurat handling.
- scipy.sparse and cKDTree remove the need for custom sparse and spatial neighbor code.
- pandas alignment prevents order bugs in cell, gene, and cell-type inputs.

## Port Boundary

The first Python port targets the SPLIT tutorial workflow:

1. generic deconvolution input
2. default purification
3. spatially aware selective purification
4. SPLIT-shift metadata swaps
5. residual transcript reassignment

RCTD-specific and Seurat-specific post-processing are intentionally outside the first core package.

