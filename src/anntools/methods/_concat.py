import sys
from pathlib import Path
from typing import Literal

import anndata as ad


def concat_anndata(
    h5ad_files: list[str],
    output: str,
    join: Literal["inner", "outer"] = "inner",
    batch_key: str | None = None,
    batch_categories: str | None = None,
):
    """
    Concatenate multiple h5ad files along the observation axis.

    Args:
        h5ad_files: List of paths to h5ad files
        output: Output path for concatenated h5ad file
        join: Join type - 'inner' keeps only common variables, 'outer' keeps all
        batch_key: Optional key to store batch labels in obs
        batch_categories: Optional comma-separated batch labels (default: file basenames)
    """
    if len(h5ad_files) < 2:
        print("Error: Need at least 2 files to concatenate", file=sys.stderr)
        sys.exit(1)

    # Validate join type
    if join not in ["inner", "outer"]:
        print(f"Error: join must be 'inner' or 'outer', got '{join}'", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {len(h5ad_files)} files...", file=sys.stderr)

    # Load all files
    adatas = []
    for i, h5ad_file in enumerate(h5ad_files):
        file_path = Path(h5ad_file)
        if not file_path.exists():
            print(f"Error: File not found: {h5ad_file}", file=sys.stderr)
            sys.exit(1)

        print(
            f"  [{i + 1}/{len(h5ad_files)}] Loading {file_path.name}...",
            file=sys.stderr,
        )
        adata = ad.read_h5ad(h5ad_file)
        print(
            f"      Shape: {adata.shape[0]:,} obs × {adata.shape[1]:,} vars",
            file=sys.stderr,
        )
        adatas.append(adata)

    # Parse batch categories if provided
    if batch_categories:
        batch_labels = [label.strip() for label in batch_categories.split(",")]
        if len(batch_labels) != len(h5ad_files):
            print(
                f"Error: Number of batch categories ({len(batch_labels)}) must match "
                f"number of files ({len(h5ad_files)})",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        # Use file basenames without extension as default batch labels
        batch_labels = [Path(f).stem for f in h5ad_files]

    print(f"\nConcatenating with join='{join}'...", file=sys.stderr)

    # Concatenate
    if batch_key:
        print(f"Adding batch labels to obs['{batch_key}']", file=sys.stderr)
        for i, label in enumerate(batch_labels):
            print(f"  - {label}: {adatas[i].shape[0]:,} cells", file=sys.stderr)

        result = ad.concat(
            adatas,
            join=join,
            label=batch_key,
            keys=batch_labels,
            index_unique="_",
        )
    else:
        result = ad.concat(
            adatas,
            join=join,
            index_unique="_",
        )

    print(
        f"\nConcatenated shape: {result.shape[0]:,} obs × {result.shape[1]:,} vars",
        file=sys.stderr,
    )

    # Report on variable handling
    if join == "inner":
        original_var_counts = [adata.shape[1] for adata in adatas]
        max_vars = max(original_var_counts)
        min_vars = min(original_var_counts)
        final_vars = result.shape[1]

        if max_vars != min_vars:
            print("\nVariable statistics:", file=sys.stderr)
            print(
                f"  Original variable counts: {min_vars:,} to {max_vars:,}",
                file=sys.stderr,
            )
            print(f"  Final (common) variables: {final_vars:,}", file=sys.stderr)
            print(f"  Variables dropped: {max_vars - final_vars:,}", file=sys.stderr)
    elif join == "outer":
        total_unique_vars = result.shape[1]
        original_var_counts = [adata.shape[1] for adata in adatas]
        total_original_vars = sum(original_var_counts)

        if total_unique_vars < total_original_vars:
            print("\nVariable statistics:", file=sys.stderr)
            print(
                f"  Total variables across files: {total_original_vars:,}",
                file=sys.stderr,
            )
            print(f"  Unique variables (union): {total_unique_vars:,}", file=sys.stderr)
            overlaps = total_original_vars - total_unique_vars
            print(f"  Overlapping variables: {overlaps:,}", file=sys.stderr)

    # Save result
    print(f"\nWriting to {output}...", file=sys.stderr)
    result.write_h5ad(output)
    print("Done!", file=sys.stderr)
