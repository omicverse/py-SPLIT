from __future__ import annotations

import argparse
import time

import numpy as np

import splitst as split


def run_case(n_cells: int, n_genes: int, n_cell_types: int):
    adata, weights, reference = split.create_split_fixture(
        n_cells=n_cells,
        n_genes=n_genes,
        n_cell_types=n_cell_types,
        random_state=11,
    )
    start = time.perf_counter()
    split.purify(adata, weights, reference, primary_cell_type=adata.obs["primary_cell_type"])
    purify_s = time.perf_counter() - start

    start = time.perf_counter()
    split.spatial_score(adata, weights, adata.obs["primary_cell_type"], adata.obs["secondary_cell_type"], k=20)
    split.balance(adata, threshold=0.15)
    split.reassign_residuals(adata, k=20)
    workflow_s = time.perf_counter() - start

    clean = np.asarray(adata.layers["clean_counts"])
    raw = np.asarray(adata.layers["counts"])
    purified = np.asarray(adata.layers["split_purified"])
    return {
        "cells": n_cells,
        "genes": n_genes,
        "cell_types": n_cell_types,
        "purify_seconds": purify_s,
        "workflow_seconds": workflow_s,
        "cells_per_second_purify": n_cells / purify_s if purify_s else float("inf"),
        "raw_l2_to_clean": float(np.linalg.norm(raw - clean)),
        "purified_l2_to_clean": float(np.linalg.norm(purified - clean)),
        "min_reassigned": float(adata.layers["split_reassigned"].min()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=240)
    parser.add_argument("--genes", type=int, default=120)
    parser.add_argument("--cell-types", type=int, default=4)
    args = parser.parse_args()
    result = run_case(args.cells, args.genes, args.cell_types)
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()

