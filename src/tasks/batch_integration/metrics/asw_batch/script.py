import anndata as ad
from scib.metrics import silhouette_batch

from read_anndata_partial import read_anndata

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print('Read input', flush=True)
input_solution = read_anndata(par['input_integrated'], obs='obsm', uns='uns')
input_solution.obs = read_anndata(par['input_solution'], obs='obs').obs

print('compute score', flush=True)
score = silhouette_batch(
    input_solution,
    batch_key='batch',
    label_key='label',
    embed='X_emb',
)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': input_solution.uns['dataset_id'],
        'normalization_id': input_solution.uns['normalization_id'],
        'method_id': input_solution.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
