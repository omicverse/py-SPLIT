# ITERATION_LOG.md

## 2026-06-08 - Initial Scaffold

- Created standalone `py-SPLIT` repository under `omicverse`.
- Implemented AnnData-native SPLIT workflow in `splitst`.
- Added deterministic fixture generator.
- Added tests for purification, spatial scoring, balancing, SPLIT-shift, and residual reassignment.
- Added benchmark script and pre-executed example notebook.

Validation:

```bash
pytest
# 8 passed
```

## 2026-06-08 - py-Augur Style Alignment

- Switched from `src/` layout to top-level package layout (`splitst/`) to mirror `py-Augur`.
- Added reconstruction docs: `DISCOVERY.md`, `MATH.md`, `RECONSTRUCTION_REPORT.md`, `ITERATION_LOG.md`.
- Updated README to include badges, pipeline overview, algorithmic fidelity, benchmark table, notebooks, API reference, and citation/license sections.
- Added simple `data/logo.svg` for README branding.

## 2026-06-08 - Test Hardening

- Added stricter validation for non-finite/negative weights, reference profiles, counts, spatial coordinates, thresholds, and radius.
- Added explicit secondary-label validation before neighborhood scoring.
- Added residual reassignment mass accounting in `adata.uns["split_residual_stats"]`.
- Expanded tests from 8 to 29 cases, covering copy mode, X-layer input, no-neighbor spatial scoring, spot-class balancing, self-keep residual accounting, reproducible fixtures, input alignment, NumPy inputs, and invalid argument paths.

Validation:

```bash
pytest
# 29 passed
```
