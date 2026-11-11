# RECONSTRUCTION_REPORT.md - splitst

## Native R decision validation (2026-11-11)

Eight native-R balance conditions now cover four thresholds and both swap settings. On the shared graph, maximum score error is 5.56e-16 and retained expression error is 5.69e-13, with zero swap/status disagreements. This exposed and repaired residual-profile substitution during label exchange, missing-class swaps and single-type status propagation. R removes rejected observations; Python intentionally retains zeroed rows. The newer evidence supersedes the earlier unverified-score/balance/shift statement below, but still does not validate independently constructed graphs or biological decontamination accuracy. Real R outputs and controlled metadata are preserved under tests/data, exercised by test_native_r_decisions.py.


## Current validation update (2026-11-10)

This section supersedes historical R-pending statements below. Direct R evidence now covers purification and residual redistribution, not every R/S4 workflow. examples/real_dlpfc_R_parity.ipynb is executed with spatial figures and numerical comparisons.

The real input is 240 spots x 512 genes from retained DLPFC 151673 counts. Reference profiles are within-slice cortical-layer means; NNLS weights are illustrative inputs, not independent cell-type truth or RCTD outputs. Agreement tests numerical kernels, not biological purification accuracy.

Unmodified R functions at bdsc-tds/SPLIT@f226987f0fecf2847646bda50a64886e0ea432a9 give maximum absolute errors of 5.12e-13 for purification and 5.12e-13 for reassignment. Reassignment shares the neighbor graph and balance decisions to isolate its kernel. Independent spatial scoring, automatic balance decisions, and label-swap parity remain unverified.

The default recipient rule was repaired: receiver primary type must match sender secondary type; rejected recipients are excluded when status is available. This follows R's second_type_neighbors_no_reject rule. receiver_policy="all" explicitly selects the former Python-only behavior. Missing labels produce an actionable error. The embedded OmicVerse implementation has the same repair.

validation/ contains CSV inputs, the R runner, hashes and R versions. tests/test_real_reference.py compares actual stored R outputs; test_receivers.py covers rejected and mismatched neighbors.

```sh
Rscript validation/run_r_reference.R /path/to/SPLIT validation
python -m pytest -q
```

The historical sections below document the initial synthetic-only delivery. Their test counts, timings, installation state, and pending-R statements are not current claims.


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

1. R parity is not yet wired because the current development machine has no `Rscript`; the repository intentionally ships no R runtime code.
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
