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
        "binomial", help="Downsampling method [binomial, multinomial]"
    ),
):
    import anndata as ad

    from adtk.methods._downsample import downsample_anndata

    adata = ad.read_h5ad(h5ad)
    if adata.X is None:
        raise ValueError("Input file does not contain data")
    adata.X = downsample_anndata(
        adata,
        fraction=fraction,
        method=method,  # type: ignore
    )
    output_path = output or h5ad.replace(".h5ad", f"_dsf{fraction:.2f}.h5ad")
    adata.write_h5ad(output_path)


@app.command()
def sparse(h5ad: str, output: str | None = typer.Option(None, help="Output file path")):
    import sys

    import anndata as ad
    from scipy.sparse import csr_matrix

    adata = ad.read_h5ad(h5ad)
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)
    else:
        print("Data is already in CSR sparse format - doing nothing", file=sys.stderr)
        return
    output_path = output or h5ad.replace(".h5ad", "_sparse.h5ad")
    adata.write_h5ad(output_path)


@app.command()
def view_obs(h5ad: str):
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
