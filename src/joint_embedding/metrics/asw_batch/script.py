import pprint
import scanpy as sc
import anndata as ad
import scib

## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    input_solution="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
    output="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.asw_batch.tsv",
    debug=True
)

## VIASH END

if par['debug']:
    pprint.pprint(par)

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = ad.read(input_prediction)
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = ad.read(input_solution)

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]

print('Compute score')
score = scib.me.silhouette_batch(
    adata,
    batch_key='batch',
    group_key='cell_type',
    embed='X_emb',
    verbose=False
)

# store adata with metrics
print("Create output object")
out = ad.AnnData(
    uns=dict(
        dataset_id=adata.uns['dataset_id'],
        method_id=adata.uns['method_id'],
        metric_ids=['asw_batch'],
        metric_values=[score]
    )
)

print("Write output to h5ad file")
out.write(output, compression='gzip')
