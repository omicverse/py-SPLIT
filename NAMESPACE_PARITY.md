# SPLIT Namespace Parity

This file tracks parity between the original R SPLIT workflow and `splitst`.

| R SPLIT tutorial concept | `splitst` API | Status |
| --- | --- | --- |
| Generic purification from deconvolution weights and reference | `splitst.purify()` | Implemented |
| Neighborhood second-type evidence | `splitst.spatial_score()` | Implemented |
| Spatially aware selective purification | `splitst.balance()` | Implemented |
| SPLIT-shift label swap metadata | `splitst.balance(..., swap_labels=True)` | Implemented |
| Residual transcript reassignment | `splitst.reassign_residuals()` | Implemented |
| R object-specific RCTD/Seurat post-processing | Not planned for core package | Out of scope |
| R numeric parity fixture | `scripts/run_r_parity.R` | Pending R environment |

Acceptance target for R parity:

- `max_abs_err <= 1e-8` for dense fixture parity
- sparse float path: `rtol=1e-7`, `atol=1e-10`
- identical cell/gene order
- documented sparse nonzero pattern differences, if any

