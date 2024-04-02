import anndata as ad
import bbknn

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

print('Run BBKNN', flush=True)
kwargs = dict(batch_key='batch', copy=True)
if adata.n_obs >= 1e5:
    kwargs['neighbors_within_batch'] = 25

ad_bbknn = bbknn.bbknn(adata, **kwargs)

ad_out = adata.copy()
ad_out.obsp['connectivities'] = ad_bbknn.obsp['connectivities']
ad_out.obsp['distances'] = ad_bbknn.obsp['distances']
ad_out.uns['neighbors'] = ad_bbknn.uns['neighbors']
ad_out.uns['method_id'] = meta['functionality_name']

print("Store outputs", flush=True)
ad_out.write_h5ad(par['output'], compression='gzip')
