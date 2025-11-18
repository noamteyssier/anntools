import sys
from typing import Literal

import anndata as ad
import pandas as pd
from adpbulk import ADPBulk
from scipy.sparse import csr_matrix


def pseudobulk(
    adata: ad.AnnData,
    groupby: str | list[str],
    method: Literal["mean", "sum", "median"] = "mean",
    layer: str | None = None,
) -> ad.AnnData:
    adpb = ADPBulk(adata, groupby=groupby, method=method, layer=layer)
    bulked = adpb.fit_transform()
    meta = adpb.get_meta().rename(columns={"SampleName": "bulk-name"})
    print(f"Input shape: {adata.shape}", file=sys.stderr)  # type: ignore
    print(f"Pseudobulk shape: {bulked.shape}", file=sys.stderr)  # type: ignore
    return ad.AnnData(
        X=csr_matrix(bulked.values),
        obs=meta.set_index("bulk-name"),
        var=pd.DataFrame(index=bulked.columns),
    )
