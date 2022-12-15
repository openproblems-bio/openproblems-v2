import anndata
import numpy as np

## VIASH START
par = {
    "input_mod1": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.mod1.h5ad",
    "input_mod2": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.mod2.h5ad",
    "output": "tmp/output_prediction.h5ad",
    "n_dims": 1,
}
## VIASH END

print("Load and prepare data")
adata_mod1 = anndata.read_h5ad(par["input_mod1"])

X = np.zeros((adata_mod1.shape[0], par["n_dims"]))
print("Saving output")
adata_out = anndata.AnnData(
    X=X,
    obs=adata_mod1.obs,
    uns={"dataset_id": adata_mod1.uns["dataset_id"], "method_id": "dummy_zeros"},
)
adata_out.write_h5ad(par["output"], compression="gzip")
