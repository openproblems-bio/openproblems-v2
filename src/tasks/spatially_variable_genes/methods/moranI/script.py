import anndata as ad
import squidpy as sq

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/10x_visium_mouse_brain/input_data.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'moranI'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run moranI', flush=True)
sq.gr.spatial_neighbors(adata,
                        coord_type="generic",
                        delaunay=True)

sq.gr.spatial_autocorr(adata,
                       mode="moran",
                       layer='normalized',
                       n_perms=100,
                       genes=adata.var_names)

# save results
df = adata.uns["moranI"]
df = df.loc[adata.var_names][['I']]
df = df.reset_index()
df.columns = ['feature_name', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
