import scanpy as sc
import numpy as np

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output': "output.h5ad",
    'layer_output': "sqrt_cpm",
    'obs_size_factors': "size_factors_sqrt_cpm",
    'n_cp': 1e6,
}
meta = {
    "functionality_name": "normalize_sqrt_cpm"
}
## VIASH END

print(">> Load data")
adata = sc.read_h5ad(par['input'])

print(">> Normalize data")
norm = sc.pp.normalize_total(
    adata, 
    target_sum=par['n_cp'], 
    layer="counts", 
    inplace=False
)
lognorm = np.sqrt(norm['X'])

print(">> Store output in adata")
adata.layers[par["layer_output"]] = lognorm
adata.obs[par["obs_size_factors"]] = norm["norm_factor"]
adata.uns["normalization_id"] = par["norm_id"]

print(">> Write data")
adata.write_h5ad(par['output'], compression="gzip")
