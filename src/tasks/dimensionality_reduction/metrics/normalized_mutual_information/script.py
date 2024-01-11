import anndata as ad
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score

## VIASH START
par = {
  'input_embedding': 'resources_test/dimensionality_reduction/pancreas/embedding.h5ad',
  'input_solution': 'resources_test/dimensionality_reduction/pancreas/solution.h5ad',
  'output': 'output.h5ad',
  'n_clusters': 10,
  'seed': 42, 
  'nmi_avg_method': 'arithmetic'
}
meta = {
  'functionality_name': 'normalized_mutual_information'
}
## VIASH END

print('Reading input files', flush=True)
input_embedding = ad.read_h5ad(par['input_embedding'])
input_solution = ad.read_h5ad(par['input_solution'])

high_dim = input_solution.layers["normalized"]
X_emb = input_embedding.obsm["X_emb"]

print('Compute metrics', flush=True)

# k-means clustering for high-dimensional data and reduced data
clusters_high = KMeans(n_clusters = par['n_clusters'], random_state = par['seed'], n_init='auto').fit_predict(high_dim)
clusters_reduced = KMeans(n_clusters = par['n_clusters'], random_state = par['seed'], n_init='auto').fit_predict(X_emb)

# Compute NMI scores
nmi = normalized_mutual_info_score(clusters_high, clusters_reduced, average_method = par['nmi_avg_method'])

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'dataset_id': input_embedding.uns['dataset_id'],
    'normalization_id': input_embedding.uns['normalization_id'],
    'method_id': input_embedding.uns['method_id'],
    'metric_ids': [ 'normalized_mutual_information' ],
    'metric_values': [ nmi ]
  }
)
output.write_h5ad(par['output'], compression='gzip')
