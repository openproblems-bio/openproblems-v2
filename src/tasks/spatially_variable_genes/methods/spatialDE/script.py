import scanpy as sc
import anndata as ad
import NaiveDE
import SpatialDE

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

# run spatialDE
sc.pp.calculate_qc_metrics(adata, 
                           layer='counts', 
                           inplace=True, 
                           percent_top=[10])
    
counts = sc.get.obs_df(adata, 
                       keys=list(adata.var_names), 
                       use_raw=False, layer='counts')

total_counts = sc.get.obs_df(adata, keys=["total_counts"])
norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(total_counts, 
                                 norm_expr.T, 
                                 "np.log(total_counts)").T
    
df = SpatialDE.run(adata.obsm["spatial"], resid_expr)

# save results
df.set_index("g", inplace=True)
df = df.loc[adata.var_names][['FSV']]
df = df.reset_index()
df.columns = ['feature_name', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': '10x_visium_mouse_brain',
                         'method_id': 'spatialDE'})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
