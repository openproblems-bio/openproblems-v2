import anndata as ad
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import adjusted_rand_score

## VIASH START
par = {
  'input_embedding': 'resources_test/dimensionality_reduction/pancreas/embedding.h5ad',
  'input_solution': 'resources_test/dimensionality_reduction/pancreas/solution.h5ad',
  'output': 'output.h5ad',
  'nmi_avg_method': 'arithmetic'
}
meta = {
  'functionality_name': 'clustering_perfrmance'
}
## VIASH END

print('Reading input files', flush=True)
input_embedding = ad.read_h5ad(par['input_embedding'])
input_solution = ad.read_h5ad(par['input_solution'])

print('Compute metrics', flush=True)

# Perform Leiden clustering on dimensionlity reduction embedding
n = 20
resolutions = [2 * x / n for x in range(1, n + 1)]
score_max = 0
res_max = resolutions[0]
clustering = None
score_all = []

try:
  input_embedding.uns["neighbors"]
except KeyError:
  print("Computing KNN")
  sc.pp.neighbors(input_embedding, use_rep="X_emb")

for res in resolutions:
  sc.tl.leiden(input_embedding, resolution=res)
  score = normalized_mutual_info_score(input_solution.obs["cell_type"], input_embedding.obs["leiden"], average_method = par['nmi_avg_method'])
  score_all.append(score)

  if score_max < score:
    score_max = score
    res_max = res
    clustering = input_embedding.obs["leiden"]

input_embedding.obs["leiden"] = clustering

# Compute NMI scores
nmi = normalized_mutual_info_score(input_solution.obs["cell_type"], input_embedding.obs["leiden"], average_method = par['nmi_avg_method'])

# Compute ARI scores
ari = adjusted_rand_score(input_solution.obs["cell_type"], input_embedding.obs["leiden"])

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'dataset_id': input_embedding.uns['dataset_id'],
    'normalization_id': input_embedding.uns['normalization_id'],
    'method_id': input_embedding.uns['method_id'],
    'metric_ids': [ 'normalized_mutual_information', 'adjusted_rand_index' ],
    'metric_values': [ nmi, ari ]
  }
)
output.write_h5ad(par['output'], compression='gzip')
