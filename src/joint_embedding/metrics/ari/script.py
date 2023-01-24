import pprint
import scanpy as sc
import anndata as ad
import scib

## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    input_solution="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
    output="openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END

if par['debug']:
    pprint.pprint(par)

print("Read prediction anndata")
adata = ad.read_h5ad(par['input_prediction'])
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = ad.read_h5ad(par['input_solution'])

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]

print('Preprocessing')
sc.pp.neighbors(adata, use_rep='X_emb')

print('Clustering')
scib.cl.opt_louvain(
    adata,
    label_key='cell_type',
    cluster_key='cluster',
    plot=False,
    inplace=True,
    force=True
)

print('Compute score')
score = scib.me.ari(adata, group1='cluster', group2='cell_type')

# store adata with metrics
print("Create output object")
out = ad.AnnData(
    uns=dict(
        dataset_id=adata.uns['dataset_id'],
        method_id=adata.uns['method_id'],
        metric_ids=["ari"],
        metric_values=[score]
    )
)

print("Write output to h5ad file")
out.write(par['output'], compression='gzip')
