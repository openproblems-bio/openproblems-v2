import anndata as ad
import sklearn.metrics
import pandas as pd

## VIASH START
par = {
  'input_method': 'resources_test/spatially_variable_genes/10x_Visium_mouse_brain/output.csv',
  'input_solution': 'resources_test/spatially_variable_genes/10x_Visium_mouse_brain/solution.csv',
  'output': 'score.h5ad'
}
meta = {
  'functionality_name': 'correlation'
}
## VIASH END

print('Reading input files', flush=True)
input_method = pd.read_csv(par['input_method'])
input_solution = pd.read_csv(par['input_solution'])


print('Compute metrics', flush=True)
df = pd.merge(input_method, input_solution, how='left')
groupby = df.groupby('gene_name', observed=True)
corr = groupby.apply(lambda x: x['pred_spatial_var_score'].corr(x['true_spatial_var_score'], method='kendall'))

uns_metric_ids = [ 'correlation' ]
uns_metric_values = [ corr.mean() ]

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'dataset_id': input_method.uns['dataset_id'],
    'method_id': input_method.uns['method_id'],
    'metric_ids': uns_metric_ids,
    'metric_values': uns_metric_values
  }
)
output.write_h5ad(par['output'], compression='gzip')

