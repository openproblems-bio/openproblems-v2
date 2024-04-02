import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True
}

meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}

## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

ad_out = adata.copy()

print('Run Combat', flush=True)
adata.X = adata.layers['normalized']
adata.X = sc.pp.combat(adata, key='batch', inplace=False)

ad_out.layers['corrected_counts'] = csr_matrix(adata.X)

print("Store outputs", flush=True)
ad_out.uns['method_id'] = meta['functionality_name']
ad_out.write_h5ad(par['output'], compression='gzip')
