import typer

app = typer.Typer()


def main():
    app()


@app.command()
def downsample(
    h5ad: str,
    fraction: float,
    output: str | None = typer.Option(None, help="Output file path"),
    method: str = typer.Option(
        "binomial", help="Downsampling method used for UMIs [binomial, multinomial]"
    ),
    which: str = typer.Option("umis", help="Downsampling method [umis, cells]"),
    seed: int | None = typer.Option(None, help="Random seed for reproducibility"),
):
    """Downsample UMIs or Cells in a given h5ad to a specified fraction."""
    import anndata as ad

    from anntools.methods._downsample import downsample_anndata

    adata = ad.read_h5ad(h5ad)
    if adata.X is None:
        raise ValueError("Input file does not contain data")
    adata = downsample_anndata(
        adata,
        fraction=fraction,
        method=method,  # type: ignore
        which=which,  # type: ignore
        seed=seed,
    )
    output_path = output or h5ad.replace(".h5ad", f".ds_{which}_{fraction:.2f}.h5ad")
    adata.write_h5ad(output_path)


@app.command()
def sparse(
    h5ad: str,
    output: str | None = typer.Option(None, help="Output file path"),
    replace: bool = typer.Option(False, help="Replace existing file"),
):
    """Convert data to CSR sparse format."""
    import sys

    import anndata as ad
    from scipy.sparse import csr_matrix

    adata = ad.read_h5ad(h5ad)
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)
    else:
        print("Data is already in CSR sparse format - doing nothing", file=sys.stderr)
        return

    if replace:
        if output is not None:
            raise ValueError("Cannot specify output path when replacing existing file")
        output_path = h5ad  # set to overwrite existing file
    else:
        output_path = output or h5ad.replace(".h5ad", "_sparse.h5ad")
    adata.write_h5ad(output_path)


@app.command()
def qc(
    h5ad: str,
    output: str | None = typer.Option(None, help="Output file path"),
    replace: bool = typer.Option(False, help="Replace existing file"),
):
    """Calculate quality control metrics for an h5ad file."""
    import anndata as ad
    import scanpy as sc

    adata = ad.read_h5ad(h5ad)
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    if replace:
        if output is not None:
            raise ValueError("Cannot specify output path when replacing existing file")
        output_path = h5ad  # set to overwrite existing file
    else:
        output_path = output or h5ad.replace(".h5ad", "_qc.h5ad")
    adata.write_h5ad(output_path)


@app.command()
def view_obs(h5ad: str):
    """View the .obs dataframe of an h5ad file."""
    import sys

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(h5ad, backed="r")
    assert adata.obs is not None, "Input file does not contain observation metadata"
    if not isinstance(adata.obs, pd.DataFrame):
        obs_df = adata.obs.to_memory()
    else:
        obs_df = adata.obs
    obs_df.reset_index().to_csv(sys.stdout, sep="\t", index=False)


@app.command()
def view_var(h5ad: str):
    """View the .var dataframe of an h5ad file."""
    import sys

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(h5ad, backed="r")
    assert adata.var is not None, "Input file does not contain variable metadata"
    if not isinstance(adata.var, pd.DataFrame):
        var_df = adata.var.to_memory()
    else:
        var_df = adata.var
    var_df.reset_index().to_csv(sys.stdout, sep="\t", index=False)


@app.command()
def info(
    h5ad: str,
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Show detailed information"
    ),
):
    """Display comprehensive information about an h5ad file."""
    from anntools.methods._info import display_info

    display_info(h5ad, verbose=verbose)


@app.command()
def concat(
    output: str = typer.Argument(..., help="Output file path for concatenated h5ad"),
    h5ad_files: list[str] = typer.Argument(..., help="Input h5ad files to concatenate"),
    join: str = typer.Option("inner", help="Join type: 'inner' or 'outer'"),
    batch_key: str | None = typer.Option(None, help="Key to store batch labels in obs"),
    batch_categories: str | None = typer.Option(
        None, help="Comma-separated batch labels (default: file basenames)"
    ),
):
    """Concatenate multiple h5ad files along the observation axis."""
    from anntools.methods._concat import concat_anndata

    concat_anndata(
        h5ad_files=h5ad_files,
        output=output,
        join=join,  # type: ignore
        batch_key=batch_key,
        batch_categories=batch_categories,
    )


@app.command()
def pseudobulk(
    h5ad: str = typer.Argument(..., help="Input h5ad file to pseudobulk"),
    groupby: list[str] = typer.Argument(..., help="Keys to group by"),
    output: str | None = typer.Option(
        None, help="Output file path for pseudobulked h5ad"
    ),
    layer: str | None = typer.Option(None, help="Layer to pseudobulk"),
    method: str = typer.Option(
        "mean", help="Aggregation method to use [mean, median, sum]"
    ),
):
    """Pseudobulk multiple h5ad files along the observation axis."""
    import anndata as ad

    from anntools.methods._pseudobulk import pseudobulk

    adata = ad.read_h5ad(h5ad)
    bulked = pseudobulk(
        adata,
        groupby=groupby,
        layer=layer,
        method=method,  # type: ignore
    )
    output_path = output or h5ad.replace(".h5ad", ".pseudobulk.h5ad")
    bulked.write_h5ad(output_path)
