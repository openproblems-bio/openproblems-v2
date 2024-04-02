import anndata as ad
from scvi.model import SCVI

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

ad_out = adata.copy()

print('Run scvi', flush=True)
SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

# Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
n_latent = 30
n_hidden = 128
n_layers = 2

vae = SCVI(
    adata,
    gene_likelihood="nb",
    n_layers=n_layers,
    n_latent=n_latent,
    n_hidden=n_hidden,
)
train_kwargs = {"train_size": 1.0}
vae.train(**train_kwargs)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": vae.get_latent_representation(),
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
