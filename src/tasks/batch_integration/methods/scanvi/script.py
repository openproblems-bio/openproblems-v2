import anndata as ad
from scvi.model import SCVI, SCANVI

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
    'n_latent': 30,
    'n_hidden': 128,
    'n_layers': 2,
    'max_epochs_scvi': 400,
    'max_epochs_scanvi': 400
}
meta = {
    'functionality_name' : 'scanvi',
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    adata = adata[:, idx].copy()

# based on scib
# -> https://github.com/theislab/scib/blob/main/scib/integration.py#L290-L297
if not par["max_epochs_scvi"]:
    par["max_epochs_scvi"] = min(int(round(20000 / adata.n_obs) * 400), 400)
    print(f"Setting max_epochs_scvi to {par['max_epochs_scvi']}", flush=True)
if not par["max_epochs_scanvi"]:
    par["max_epochs_scanvi"] = min(max(2, int(round(par["max_epochs_scvi"] / 3.0))), 10)
    print(f"Setting max_epochs_scanvi to {par['max_epochs_scanvi']}", flush=True)

print("Processing data", flush=True)
SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

print("Run scVI", flush=True)
model_kwargs = {
    key: par[key]
    for key in ["n_latent", "n_hidden", "n_layers"]
    if par[key] is not None
}

vae = SCVI(adata, **model_kwargs)

vae.train(max_epochs=par["max_epochs_scvi"], train_size=1.0)

print('Run SCANVI', flush=True)
scanvae = SCANVI.from_scvi_model(
    scvi_model=vae,
    labels_key="label",
    unlabeled_category="UnknownUnknown", # pick anything definitely not in a dataset
)
scanvae.train(max_epochs=par["max_epochs_scanvi"], train_size=1.0)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": scanvae.get_latent_representation(),
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
