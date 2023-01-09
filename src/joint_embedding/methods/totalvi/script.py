import anndata
import scanpy as sc
from scipy import sparse
from scvi.model import TOTALVI

## VIASH START
par = {
    "input_mod1": "output/public_datasets/joint_embedding/totalvi_spleen_lymph_111/totalvi_spleen_lymph_111.censor_dataset.output_mod1.h5ad",
    "input_mod2": "output/public_datasets/joint_embedding/totalvi_spleen_lymph_111/totalvi_spleen_lymph_111.censor_dataset.output_mod2.h5ad",
    "output": "tmp/output_prediction.h5ad",
    "hvg_number": 4000,
    "max_epochs": 20
}

meta = {
    'funcionality_name': "foo"
}
## VIASH END

print("Load and prepare data", flush=True)
adata_mod1 = anndata.read_h5ad(par['input_mod1'])
adata_mod2 = anndata.read_h5ad(par['input_mod2'])
adata_mod1.obsm['protein_expression'] = adata_mod2.layers["counts"].toarray()

print('Select highly variable genes', flush=True)
sc.pp.highly_variable_genes(
    adata_mod1,
    n_top_genes=par['hvg_number'],
    flavor="seurat_v3",
    batch_key="batch",
    subset=True
)

print("Set up model", flush=True)
TOTALVI.setup_anndata(
    adata_mod1,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)

print('Train totalVI with', par['max_epochs'], 'epochs', flush=True)
vae = TOTALVI(adata_mod1, latent_distribution="normal")
vae.train(max_epochs = par['max_epochs'])

print("Postprocessing and saving output", flush=True)
adata_out = anndata.AnnData(
    X=vae.get_latent_representation(),
    obs=adata_mod1.obs[['batch']],
    uns={
        "dataset_id": adata_mod1.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    },
    obsm = {"X_emb": sparse.csr_matrix(vae.get_latent_representation())}
)

del adata_out.X

adata_out.write_h5ad(par['output'], compression = "gzip")
