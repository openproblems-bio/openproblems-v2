

print('Importing libraries')
import pprint
import numpy as np
import scanpy as sc
import anndata as ad
import scib

## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    input_solution="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
    output="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.ti_cons.h5ad",
    debug=True
)
## VIASH END

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'ti_cons'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = ad.read_h5ad(input_prediction)
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = ad.read_h5ad(input_solution)

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]
adt_atac_trajectory = 'pseudotime_order_ATAC' if 'pseudotime_order_ATAC' in adata_solution.obs else 'pseudotime_order_ADT'

print('Preprocessing')
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute scores')
obs_keys = adata_solution.obs_keys()

if 'pseudotime_order_GEX' in obs_keys:
    score_rna = scib.me.trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key='pseudotime_order_GEX'
    )
else:
    score_rna = np.nan

if adt_atac_trajectory in obs_keys:
    score_adt_atac = scib.me.trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key=adt_atac_trajectory
    )
else:
    score_adt_atac = np.nan

score_mean = (score_rna + score_adt_atac) / 2

# store adata with metrics
print("Create output object")
out = ad.AnnData(
    uns=dict(
        dataset_id=adata.uns['dataset_id'],
        method_id=adata.uns['method_id'],
        metric_ids=['ti_cons_RNA', 'ti_cons_ADT_ATAC', 'ti_cons_mean'],
        metric_values=[score_rna, score_adt_atac, score_mean]
    )
)

print("Write output to h5ad file")
out.write(output, compression='gzip')