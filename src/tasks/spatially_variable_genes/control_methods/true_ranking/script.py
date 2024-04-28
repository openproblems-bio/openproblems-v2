import anndata as ad
import pandas as pd

## VIASH START
par = {
  'input_data': 'resources_test/spatially_variable_genes/10x_Visium_mouse_brain/input.h5ad',
  'input_solution': 'resources_test/spatially_variable_genes/10x_Visium_mouse_brain/solution.csv',
  'output': 'output.csv'
}
meta = {
  'functionality_name': 'true_ranking'
}
## VIASH END

print('Reading input files', flush=True)
input_data = ad.read_h5ad(par['input_data'])

print('Generate predictions', flush=True)
output = pd.read_csv(par['input_solution'])

print("Write output to file", flush=True)
output.columns = ['feature_name', 'gene_name', 'pred_spatial_var_score']

output.to_csv(par['output'], index=False)