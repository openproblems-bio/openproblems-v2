import scanpy as sc
import numpy as np

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output': "output.h5ad",
    'layer_output': "sqrt_cpm",
    'obs_size_factors': "size_factors_sqrt_cpm"
}
meta = {
    "functionality_name": "normalize_sqrt_cpm"
}
## VIASH END

print(">> Load data", flush=True)
adata = sc.read_h5ad(par['input'])

print(">> Normalize data", flush=True)
norm = sc.pp.normalize_total(
    adata, 
    target_sum=1e6, 
    layer="counts", 
    inplace=False
)
lognorm = np.sqrt(norm["X"])

print(">> Store output in adata", flush=True)
adata.layers[par["layer_output"]] = lognorm
adata.obs[par["obs_size_factors"]] = norm["norm_factor"]
adata.uns["normalization_id"] = meta["functionality_name"]

print(">> Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
