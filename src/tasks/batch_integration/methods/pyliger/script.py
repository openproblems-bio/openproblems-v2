import anndata as ad
import pyliger

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])
print(adata)
adata.X = adata.layers['counts']
adata.obs.index.name = 'cell_barcode'
adata.var.index.name = 'gene_id'

adata_per_batch = []
for batch in adata.obs['batch'].unique():
  ad = adata[adata.obs['batch'] == batch]
  ad.uns['sample_name'] = batch
  adata_per_batch.append(ad)

print('Run pyliger', flush=True)
lobj = pyliger.create_liger(
  adata_per_batch,
  remove_missing=False,
)

# normalise
pyliger.normalize(lobj, remove_missing=False)

# do not select genes
lobj.var_genes = adata.var_names

# scale data
pyliger.scale_not_center(lobj, remove_missing=False)

# run model
pyliger.optimize_ALS(
  lobj,
  k=20,
  thresh=3,
  nrep=5e-5,
)
pyliger.quantile_norm(
  lobj,
  resolution=0.4,
  small_clust_thresh=20,
)
print(lobj, flush=True)

adata.obsm['X_emb'] = lobj.H_norm

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
