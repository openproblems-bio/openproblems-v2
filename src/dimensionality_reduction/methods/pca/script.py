## VIASH START
par = {
    'input': "../../../../resources_test/dimensionality_reduction/tenx_5k_pbmc/toy_data.h5ad",
    'output': "/tmp/pca.h5ad"
}
## VIASH END

import sys
sys.path.append(meta['resources_dir'])
import scanpy as sc
from utils import check_version




print(">> Load data")
adata = sc.read(par['input'])

print(f">> Running {meta['functionality_name']} method")
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
adata.obsm["X_emb"] = adata.obsm["X_pca"][:, :2]
adata.uns["method_id"] = meta["functionality_name"]
adata.uns["method_code_version"] = check_version("scikit-learn")

print(">> Write data")
adata.write(par['output'], compression="gzip")
