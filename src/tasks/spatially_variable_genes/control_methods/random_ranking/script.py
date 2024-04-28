import anndata as ad
import numpy as np

## VIASH START
par = {
  'input_data': 'resources_test/spatially_variable_genes/10x_Visium_mouse_brain/input.h5ad',
  'output': 'output.csv'
}
meta = {
  'functionality_name': 'random_ranking'
}
## VIASH END

print('Reading input files', flush=True)
input_data = ad.read_h5ad(par['input_data'])

print('Generate predictions', flush=True)
input_data.var['ranking'] = np.random.permutation(input_data.n_vars)

print("Write output AnnData to file", flush=True)
output = input_data.var

output.to_csv(par['output'])