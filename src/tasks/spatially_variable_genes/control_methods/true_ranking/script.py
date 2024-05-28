import anndata as ad

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/10x_visium_mouse_brain/dataset.h5ad',
    'input_solution': 'resources_test/spatially_variable_genes/10x_visium_mouse_brain/solution.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'true_ranking'
}
# VIASH END

print('Reading input files', flush=True)
input_data = ad.read_h5ad(par['input_data'])

print('Generate predictions', flush=True)
input_solution = ad.read_h5ad(par['input_solution'])

df = input_solution.var
df.columns = ['feature_name', 'gene_name', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': '10x_visium_mouse_brain',
                         'method_id': 'true_ranking'})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
