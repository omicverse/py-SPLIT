# MATH.md - SPLIT-ST Formula Notes

## Input Convention

`splitst` uses AnnData cell x gene orientation.

| Symbol | Shape | Meaning |
|---|---:|---|
| `X` | cells x genes | raw contaminated counts |
| `W` | cells x cell types | deconvolution weights |
| `R` | cell types x genes | reference expression profiles |
| `p_i` | scalar | primary cell type for cell `i` |
| `w_i,p` | scalar | primary weight for cell `i` |

Rows of `W` are normalized to sum to one when needed.

## Default Purification

For cell `i` and gene `g`:

```text
D_i,g = sum_c W_i,c R_c,g + epsilon
N_i,g = W_i,p_i R_p_i,g + epsilon / n_i
Y_i,g = X_i,g * N_i,g / D_i,g
```

where `n_i` is the number of positive cell-type weights for cell `i`.

`Y` is stored in `adata.layers["split_purified"]` by default.

## Spatial Score

For each focal cell, neighbors are found with `scipy.spatial.cKDTree`.

The local composition vector accumulates weighted primary/secondary support from neighboring cells and normalizes to sum to one. The key score is:

```text
neighborhood_weights_second_type_i = local_composition_i[secondary_i]
```

## Selective Balance

Cells with second-type score above a threshold use purified profiles:

```text
balanced_i = purified_i if score_i > threshold else raw_i
```

Optional spot classes can force `doublet_uncertain` cells to purified and `reject` cells to zero.

## SPLIT-shift

When local neighborhood evidence supports the secondary label rather than the primary label, `balance(..., swap_labels=True)` records `split_shift_swap=True` and swaps first/second type metadata.

## Residual Reassignment

Positive residuals are:

```text
residual = max(raw - balanced, 0)
```

A sparse operator `A` redistributes sender residuals to spatial neighbors by uniform or count-proportional weights:

```text
reassigned = balanced + A.T @ residual
```

The operator is stored in `adata.uns["split_reassignment_operator"]`.

