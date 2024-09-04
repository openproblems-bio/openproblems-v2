import sys
import anndata as ad
import mnnpy

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)

if par['n_hvg']:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:par['n_hvg']]
    adata = adata[:, idx].copy()

print('Run mnn', flush=True)
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

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    },
    layers={
        'corrected_counts': corrected.X,
    }
)


print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
