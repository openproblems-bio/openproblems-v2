import anndata as ad
import sklearn.decomposition


## VIASH START
par = {
    "input": "resources_test/common/multimodal/normalized_mod1.h5ad",
    "output": "output.h5ad",
    "layer_input": "normalized",
    "obsm_embedding": "X_svd",
    "num_components": 100,
}
## VIASH END

print(">> Load data", flush=True)
adata = ad.read(par["input"])

print(">> check parameters", flush=True)
n_svd = min([par["num_components"], min(adata.layers[par["layer_input"]].shape) - 1])

print(">> Run SVD", flush=True)
svd = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.layers[par["layer_input"]])

print(">> Storing output", flush=True)
adata.obsm[par["obsm_embedding"]] = svd

print(">> Writing data", flush=True)
adata.write_h5ad(par["output"])
