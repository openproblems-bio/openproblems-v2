import anndata as ad
import numpy as np
from scvi.model import SCVI, SCANVI

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END



print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])


# Set the max epochs for training
n_epochs_scVI = int(np.min([round((20000 / adata.n_obs) * 400), 400]))  # 400
n_epochs_scANVI = int(np.min([10, np.max([2, round(n_epochs_scVI / 3.0)])]))

print('Run SCVI', flush=True)
SCVI.setup_anndata(adata, layer="normalized", batch_key="batch")

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
train_kwargs["max_epochs"] = n_epochs_scVI

vae.train(**train_kwargs)

print('Run SCANVI', flush=True)
scanvae = SCANVI.from_scvi_model(
    scvi_model=vae,
    labels_key="label",
    unlabeled_category="UnknownUnknown",  # pick anything definitely not in a dataset
)
scanvae.train(max_epochs=n_epochs_scANVI, train_size=1.0)
adata.obsm["X_emb"] = scanvae.get_latent_representation()
del adata.X

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
