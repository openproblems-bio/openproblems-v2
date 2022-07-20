## VIASH START
par = {
    "input": "input_data.h5ad",
    "output": "output_data.h5ad",
    "max_iter": 100,
    # "not_required": "foo"
}
## VIASH END
import anndata
import scanpy as sc

# Read data based on variables passed to the component
adata = sc.read(par['input'])

# Do something with the other variables
adata.uns["iteractions"] = par['max_iter']

# Do something if the required parameter was passed
if par.get('not_required'):
    adata.uns["not_required"] = par['not_required']

# Write resulted data in the path passed to the component
adata.write(par['output'], compression="gzip")
