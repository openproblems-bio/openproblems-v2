import anndata as ad
from muon import atac as ac

## VIASH START
par = {
    'input': "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod2.h5ad",
    'output': "output_norm.h5ad"
}
meta = {
    'functionality_name': "tfidf"
}
## VIASH END

print("Load data", flush=True)
adata = ad.read_h5ad(par['input'])

print("Normalize data", flush=True)
ac.pp.tfidf(
    adata, 
    inplace=True,
    from_layer="counts",
    to_layer=par["layer_output"]
)

print("Store output in adata", flush=True)
adata.uns["normalization_id"] = par["normalization_id"] or meta['functionality_name']

print("Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
