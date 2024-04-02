import anndata as ad
import mnnpy

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

ad_out = adata.copy()

print('Run mnn', flush=True)
adata.X = adata.layers['normalized']
split = []
batch_categories = adata.obs['batch'].cat.categories
for i in batch_categories:
    split.append(adata[adata.obs['batch'] == i].copy())
corrected, _, _ = mnnpy.mnn_correct(
        *split,
        batch_key='batch',
        batch_categories=batch_categories,
        index_unique=None
    )

ad_out.layers['corrected_counts'] = corrected.X


print("Store outputs", flush=True)
ad_out.uns['method_id'] = meta['functionality_name']
ad_out.write_h5ad(par['output'], compression='gzip')
