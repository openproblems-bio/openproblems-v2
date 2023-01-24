import anndata
import numpy as np
from scipy import sparse

## VIASH START
par = {
    "input_mod1": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.mod1.h5ad",
    "output": "output/output_prediction.h5ad",
    "n_dims": 100,
}
## VIASH END

print("Load and prepare data")
adata_mod1 = anndata.read_h5ad(par["input_mod1"])

X = np.random.randn(adata_mod1.shape[0], par["n_dims"])
print("Saving output")
adata_out = anndata.AnnData(
    X=X,
    obsm= {"X_emb": sparse.csr_matrix(X)},
    obs=adata_mod1.obs[["batch"]],
    uns={"dataset_id": adata_mod1.uns["dataset_id"], "method_id": "random_embed"},
)
del adata_out.X

adata_out.write_h5ad(par["output"], compression="gzip")
