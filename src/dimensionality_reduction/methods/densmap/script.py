## VIASH START
par = {
    'input': "../../../../resources_test/dimensionality_reduction/tenx_5k_pbmc/toy_data.h5ad",
    'pca': True
    'output': "/tmp/pca.h5ad"
}
## VIASH END

import sys
sys.path.append(meta['resources_dir'])
import scanpy as sc
from utils import check_version
from umap import UMAP

print(">> Load data")
adata = sc.read(par['input'])

print(f">> Running {meta['functionality_name']} method")
if par['pca']:
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    X = adata.obsm["X_pca"]
else:
    X = adata.X

adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(X)
adata.uns["method_id"] = meta["functionality_name"]
adata.uns["method_code_version"] = check_version("umap-learn")

print(">> Write data")
adata.write(par['output'], compression="gzip")
