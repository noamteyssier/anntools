import sys
from pathlib import Path

import anndata as ad
import numpy as np
from anndata.abc import CSRDataset
from scipy.sparse import issparse


def _format_size(size_bytes: int | float) -> str:
    """Format bytes into human-readable size."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


def _get_dtype_info(matrix) -> str:
    """Get dtype information from matrix."""
    if issparse(matrix):
        return str(matrix.dtype)
    elif hasattr(matrix, "dtype"):
        return str(matrix.dtype)
    return "unknown"


def _get_sparsity_info(matrix) -> tuple[str, float]:
    """Get sparsity information from matrix."""
    if isinstance(matrix, CSRDataset):
        matrix = matrix.to_memory()

    if issparse(matrix):
        total_elements = matrix.shape[0] * matrix.shape[1]
        nonzero = matrix.nnz
        sparsity = 1 - (nonzero / total_elements)
        return "sparse", sparsity
    return "dense", 0.0


def _summarize_column(series, col_name: str, max_display: int = 5) -> str:
    """Summarize a pandas Series column."""
    dtype = str(series.dtype)

    if series.dtype == "object" or series.dtype.name == "category":
        n_unique = series.nunique()
        if n_unique <= max_display:
            values = sorted(series.unique())
            return f"{col_name} ({dtype}): {n_unique} unique - {values}"
        else:
            sample_values = sorted(series.unique())[:max_display]
            return f"{col_name} ({dtype}): {n_unique} unique - {sample_values}..."
    elif np.issubdtype(series.dtype, np.number):
        return f"{col_name} ({dtype}): min={series.min():.2f}, max={series.max():.2f}, mean={series.mean():.2f}, median={series.median():.2f}"
    elif np.issubdtype(series.dtype, np.bool_):
        n_true = series.sum()
        return f"{col_name} ({dtype}): {n_true} True, {len(series) - n_true} False"
    else:
        return f"{col_name} ({dtype})"


def display_info(h5ad_path: str, verbose: bool = False):
    """Display comprehensive information about an h5ad file."""

    # Get file size
    file_path = Path(h5ad_path)
    if not file_path.exists():
        print(f"Error: File not found: {h5ad_path}", file=sys.stderr)
        sys.exit(1)

    file_size = file_path.stat().st_size

    print("=" * 80)
    print(f"AnnData File Info: {file_path.name}")
    print("=" * 80)

    # Load in backed mode for efficiency
    adata = ad.read_h5ad(h5ad_path, backed="r")

    # Basic information
    print(f"\nFile size: {_format_size(file_size)}")
    print(f"Shape: {adata.shape[0]:,} observations × {adata.shape[1]:,} variables")

    # Matrix information
    if adata.X is not None:
        matrix_type, sparsity = _get_sparsity_info(adata.X)
        dtype = _get_dtype_info(adata.X)
        print("\nMain matrix (X):")
        print(f"  Type: {matrix_type}")
        print(f"  Dtype: {dtype}")
        if matrix_type == "sparse":
            print(f"  Sparsity: {sparsity:.2%}")
    else:
        print("\nMain matrix (X): None")

    # Layers
    if adata.layers and len(adata.layers) > 0:
        print(f"\nLayers ({len(adata.layers)}):")
        for layer_name in adata.layers.keys():
            layer = adata.layers[layer_name]
            matrix_type, sparsity = _get_sparsity_info(layer)
            dtype = _get_dtype_info(layer)
            print(f"  - {layer_name}: {matrix_type}, {dtype}", end="")
            if matrix_type == "sparse":
                print(f", sparsity={sparsity:.2%}")
            else:
                print()

    # Observations metadata
    if adata.obs is not None and len(adata.obs.columns) > 0:
        print(f"\nObservations metadata (obs) - {len(adata.obs.columns)} columns:")
        if verbose:
            for col in adata.obs.columns:
                summary = _summarize_column(adata.obs[col], col)
                print(f"  - {summary}")
        else:
            cols = list(adata.obs.columns)
            if len(cols) <= 10:
                print(f"  {', '.join(cols)}")
            else:
                print(f"  {', '.join(cols[:10])}... (+{len(cols) - 10} more)")

    # Variables metadata
    if adata.var is not None and len(adata.var.columns) > 0:
        print(f"\nVariables metadata (var) - {len(adata.var.columns)} columns:")
        if verbose:
            for col in adata.var.columns:
                summary = _summarize_column(adata.var[col], col)
                print(f"  - {summary}")
        else:
            cols = list(adata.var.columns)
            if len(cols) <= 10:
                print(f"  {', '.join(cols)}")
            else:
                print(f"  {', '.join(cols[:10])}... (+{len(cols) - 10} more)")

    # Unstructured annotations
    if adata.uns and len(adata.uns) > 0:
        print(f"\nUnstructured annotations (uns) - {len(adata.uns)} keys:")
        keys = list(adata.uns.keys())
        if len(keys) <= 15:
            for key in keys:
                value = adata.uns[key]
                if isinstance(value, dict):
                    print(f"  - {key}: dict with {len(value)} keys")
                elif isinstance(value, (list, tuple)):
                    print(f"  - {key}: {type(value).__name__} of length {len(value)}")
                elif isinstance(value, np.ndarray):
                    print(f"  - {key}: array of shape {value.shape}")
                else:
                    print(f"  - {key}: {type(value).__name__}")
        else:
            print(f"  {', '.join(keys[:15])}... (+{len(keys) - 15} more)")

    # obsm (multi-dimensional observations)
    if adata.obsm and len(adata.obsm) > 0:
        print(f"\nMulti-dimensional observations (obsm) - {len(adata.obsm)} keys:")
        for key in adata.obsm.keys():
            value = adata.obsm[key]
            if hasattr(value, "shape"):
                print(f"  - {key}: shape {value.shape}")
            else:
                print(f"  - {key}: {type(value).__name__}")

    # varm (multi-dimensional variables)
    if adata.varm and len(adata.varm) > 0:
        print(f"\nMulti-dimensional variables (varm) - {len(adata.varm)} keys:")
        for key in adata.varm.keys():
            value = adata.varm[key]
            if hasattr(value, "shape"):
                print(f"  - {key}: shape {value.shape}")
            else:
                print(f"  - {key}: {type(value).__name__}")

    # obsp (pairwise observations)
    if adata.obsp and len(adata.obsp) > 0:
        print(f"\nPairwise observations (obsp) - {len(adata.obsp)} keys:")
        for key in adata.obsp.keys():
            value = adata.obsp[key]
            matrix_type, sparsity = _get_sparsity_info(value)
            print(f"  - {key}: {matrix_type}, shape {value.shape}", end="")
            if matrix_type == "sparse":
                print(f", sparsity={sparsity:.2%}")
            else:
                print()

    # varp (pairwise variables)
    if adata.varp and len(adata.varp) > 0:
        print(f"\nPairwise variables (varp) - {len(adata.varp)} keys:")
        for key in adata.varp.keys():
            value = adata.varp[key]
            matrix_type, sparsity = _get_sparsity_info(value)
            print(f"  - {key}: {matrix_type}, shape {value.shape}", end="")
            if matrix_type == "sparse":
                print(f", sparsity={sparsity:.2%}")
            else:
                print()

    print("=" * 80)

    if not verbose:
        print("\nUse --verbose (-v) flag for detailed metadata summaries")
