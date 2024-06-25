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

print('Generate predictions', flush=True)
adata = ad.read_h5ad(par['input_data'])

sq.gr.spatial_neighbors(adata,
                        coord_type="generic",
                        delaunay=True)

sq.gr.spatial_autocorr(adata,
                       mode="moran",
                       layer="normalized_data",
                       n_perms=100,
                       n_jobs=10,
                       genes=adata.var_names)

# save results
df = adata.uns["moranI"]
df = df.loc[adata.var_names][['I']]
df = df.reset_index()
df.columns = ['feature_name', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': '10x_visium_mouse_brain',
                         'method_id': 'moranI'})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
