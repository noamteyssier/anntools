import sys
from typing import Literal

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix


def _downsample_binomial(matrix: csr_matrix, frac: float) -> csr_matrix:
    """
    Downsample a sparse count matrix using binomial sampling.

    Args:
        matrix: sparse count matrix (cells x genes)
        downsample_frac: fraction to downsample to (0 < frac <= 1)

    Returns:
        Downsampled sparse matrix
    """
    # Work with the sparse matrix data directly
    downsampled_matrix = matrix.copy()

    # For each non-zero entry, sample from binomial
    downsampled_matrix.data = np.random.binomial(
        matrix.data.astype(int),
        frac,
    )

    # Remove zeros that resulted from downsampling
    downsampled_matrix.eliminate_zeros()

    return downsampled_matrix


def _downsample_multinomial(matrix: csr_matrix, frac: float) -> csr_matrix:
    """
    Downsample a sparse count matrix using multinomial sampling.

    For each cell (row), creates a probability vector from the guide counts,
    then resamples using multinomial with n = rowsum * frac draws.

    Args:
        matrix: sparse count matrix (cells x guides)
        frac: fraction to downsample to (0 < frac <= 1)

    Returns:
        Downsampled sparse matrix
    """
    n_cells, _ = matrix.shape  # type: ignore

    # Convert to LIL format for efficient row-wise operations
    downsampled_matrix = matrix.tolil()

    for i in range(n_cells):
        # Get the row as a dense array
        row = matrix.getrow(i).toarray().flatten()

        # Skip empty rows
        row_sum = row.sum()
        if row_sum == 0:
            continue

        # Calculate probability vector (normalize counts)
        probs = row / row_sum

        # Calculate number of draws for this cell
        n_draws = int(row_sum * frac)

        if n_draws == 0:
            # Set row to zeros
            downsampled_matrix[i, :] = 0
            continue

        # Sample from multinomial
        new_counts = np.random.multinomial(n_draws, probs)

        # Update the row
        downsampled_matrix[i, :] = new_counts

    # Convert back to CSR format and eliminate zeros
    downsampled_matrix = downsampled_matrix.tocsr()
    downsampled_matrix.eliminate_zeros()

    return csr_matrix(downsampled_matrix)


def _downsample_cells(matrix: csr_matrix, fraction: float) -> np.ndarray:
    mask = (
        np.random.random(
            matrix.shape[0]  # type: ignore
        )
        < fraction
    )
    return mask


def downsample_anndata(
    adata: ad.AnnData,
    fraction: float,
    method: Literal["binomial", "multinomial"],
    which: Literal["umis", "cells"],
    seed: int | None = None,
) -> ad.AnnData:
    if seed:
        np.random.seed(seed)

    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)

    pre_shape = adata.X.shape
    total_elements = np.prod(
        pre_shape  # type: ignore
    )
    pre_nonzero_fraction = adata.X.nnz / total_elements
    pre_nonzero_mean = adata.X.data.mean()

    if which == "umis":
        if method == "binomial":
            matrix = _downsample_binomial(adata.X, fraction)
        elif method == "multinomial":
            matrix = _downsample_multinomial(adata.X, fraction)
        else:
            raise ValueError(f"Unknown method {method}")
        adata.X = matrix
    else:
        mask = _downsample_cells(adata.X, fraction)
        adata = adata[mask, :]
        matrix: csr_matrix = adata.X  # type: ignore
        total_elements = np.prod(
            matrix.shape  # type: ignore
        )

    post_shape = adata.shape
    post_nonzero_fraction = matrix.nnz / total_elements
    post_nonzero_mean = matrix.data.mean()

    print(f"Pre-filter shape: {pre_shape}", file=sys.stderr)
    print(f"Pre-filter nonzero sparsity: {pre_nonzero_fraction}", file=sys.stderr)
    print(f"Pre-filter nonzero mean: {pre_nonzero_mean}", file=sys.stderr)
    print(f"Post-filter shape: {post_shape}", file=sys.stderr)
    print(f"Post-filter nonzero sparsity: {post_nonzero_fraction}", file=sys.stderr)
    print(f"Post-filter nonzero mean: {post_nonzero_mean}", file=sys.stderr)
    return adata
