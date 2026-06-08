# RECONSTRUCTION_REPORT.md - splitst

## 1. Identity

| Field | Value |
|---|---|
| Port name | splitst |
| Repository | https://github.com/omicverse/py-SPLIT |
| Upstream R package | SPLIT |
| Upstream repository | https://github.com/bdsc-tds/SPLIT |
| Algorithm class | Deterministic numerical workflow |
| Parity metric | Formula parity and workflow sanity on deterministic AnnData fixtures |
| Current parity | Python formula parity green; R parity pending local R SPLIT environment |
| Audit classification | A for formula implementation; B pending R object workflow parity |
| License | GPL-3 |
| Python version | 0.1.0 |

## 2. R Tutorial Coverage Audit

| SPLIT tutorial concept | Python equivalent | Status | Notes |
|---|---|---|---|
| Generic purification | `splitst.purify()` | Ported | Main entry point |
| Primary label inference | `primary_cell_type=None` | Ported | Uses max deconvolution weight |
| Spatial second-type score | `splitst.spatial_score()` | Ported | KNN/radius support |
| Spatially aware selective purification | `splitst.balance()` | Ported | Raw/purified merge |
| SPLIT-shift | `splitst.balance(..., swap_labels=True)` | Ported | Metadata swap logic |
| Residual transcript reassignment | `splitst.reassign_residuals()` | Ported | Sparse operator |
| RCTD S4 object handling | Not implemented | Deferred | Use generic weights/reference |
| Seurat object handling | Not implemented | Deferred | Convert to AnnData first |

## 3. Parity Evidence

Python tests cover:

- dense formula parity against explicit NumPy implementation
- sparse input/output equivalence
- gene/cell-type orientation alignment
- primary label inference
- selective `cells_to_purify`
- spatial KNN/radius scoring
- balance threshold behavior
- SPLIT-shift label swaps
- residual reassignment non-negativity and operator row sums
- biological sanity: purified profiles move closer to clean profiles

Current test gate:

```bash
pytest
# 29 passed
```

## 4. Benchmark Evidence

Small fixture:

| Metric | Value |
|---|---:|
| Cells x genes x cell types | 240 x 120 x 4 |
| Purification time | ~0.02 s |
| Full workflow time | ~0.19 s |
| Raw L2 distance to clean | 387.84 |
| Purified L2 distance to clean | 3.08e-09 |
| Minimum reassigned count | > 0 |

Run:

```bash
python benchmarks/benchmark_split.py
```

## 5. Code Quality Audit

| Check | Status |
|---|---|
| `pip install -e ".[dev]"` succeeds | Pass |
| `pytest` green | Pass (29 tests) |
| Quickstart notebook pre-executed | Pass |
| License matches upstream class | Pass |
| CI workflow added | Pass |

## 6. Known Limitations

1. R parity is not yet wired because the current development machine has no `Rscript`.
2. RCTD and Seurat object adapters are outside the core package; users should pass generic weights and reference profiles.
3. Plotting is kept in notebooks; core package has no plotting API yet.
4. Residual reassignment implements practical neighbor-based redistribution, pending exact comparison against the upstream tutorial outputs.

## 7. Integration Into OmicVerse

- **Package location**: `splitst/` as a standalone package
- **Public API**: `purify`, `spatial_score`, `balance`, `reassign_residuals`
- **OmicVerse wrapper plan**: `omicverse.space.split_*` can call `splitst` after the standalone package stabilizes
- **Tutorial slot**: spatial deconvolution / purification tutorials

## 8. Sign-off

| Field | Value |
|---|---|
| Author | splitst port (agent-assisted) |
| Date | 2026-06-08 |
| Audit classification | A/B pending R parity |
