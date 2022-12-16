import pprint
import anndata as ad
import scib

## VIASH START
par = dict(
    input_prediction="resources_test/common/joint_embedding/cite_random_prediction.h5ad",
    input_solution="resources_test/common/joint_embedding/cite_solution.h5ad",
    output="resources_test/common/joint_embedding/score_cc_cons.h5ad",
    debug=False
)
## VIASH END


if par['debug']:
    pprint.pprint(par)

print("Read prediction anndata")
adata = ad.read_h5ad(par['input_prediction'])
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = ad.read_h5ad(par['input_solution'])
organism = adata_solution.uns['organism']

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]
recompute_cc = 'S_score' not in adata_solution.obs_keys() or \
               'G2M_score' not in adata_solution.obs_keys()

print('Compute score')
score = scib.me.cell_cycle(
    adata_pre=adata_solution,
    adata_post=adata,
    batch_key='batch',
    embed='X_emb',
    recompute_cc=recompute_cc,
    organism=organism
)

# store adata with metrics
print("Create output object")
out = ad.AnnData(
    uns= {
        "dataset_id":adata.uns['dataset_id'],
        "method_id":adata.uns['method_id'],
        "metric_ids":['cc_cons'],
        "metric_values":[score],
    }
)

print("Write output to h5ad file")
out.write(par['output'], compression='gzip')
