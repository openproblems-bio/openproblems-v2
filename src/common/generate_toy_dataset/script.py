# https://github.com/scverse/anndata-tutorials/blob/master/getting-started.ipynb
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


## VIASH START
par = {
    "output": "foo_data.h5ad"
}
## VIASH END

counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency

adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))

adata.uns["random"] = [1, 2, 3]

adata.layers["log_transformed"] = np.log1p(adata.X)

print(">> Generating random toy data")
adata.write(par['output'], compression="gzip")
