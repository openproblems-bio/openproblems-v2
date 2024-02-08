import anndata as ad

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/pancreas/single_cell.h5ad',
  'input_spatial': 'resources_test/spatial_decomposition/pancreas/spatial.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'functionality_name': 'cell2location'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial'])

print('Preprocess data', flush=True)
# ... preprocessing ...

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  obsm={
    'coordinates': obsm_coordinates,
    'proportions_pred': obsm_proportions_pred
  },
  layers={
    'counts': layers_counts
  },
  uns={
    'cell_type_names': uns_cell_type_names,
    'dataset_id': input_single_cell.uns['dataset_id'],
    'method_id': meta['functionality_name']
  }
)
output.write_h5ad(par['output'], compression='gzip')
