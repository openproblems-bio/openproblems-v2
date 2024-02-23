import sys
import anndata as ad
import scanpy as sc

from utils import generate_synthetic_dataset
from utils import filter_genes_cells

def _pancreas_synthetic(
    adata: ad.AnnData,
    alpha: float,
    n_obs: int = 100 
) -> ad.AnnData:

    adata.X = adata.layers["counts"]
    sc.pp.filter_genes(adata, min_counts=10)
    adata_merged = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=alpha)
    adata.uns["spatial_data_summary"] = f"Dirichlet alpha={alpha}"
    filter_genes_cells(adata_merged)
    adata_merged.X = None
    adata_merged.obs['is_primary_data'] = adata_merged.obs['is_primary_data'].fillna(False)

    return adata_merged

# def pancreas_alpha_1(adata, n_obs=100):
#     return _pancreas_synthetic(adata, n_obs=n_obs, alpha=1)

# def pancreas_alpha_5(adata, n_obs=100):
#     return _pancreas_synthetic(adata, n_obs=n_obs, alpha=5)

# def pancreas_alpha_0_5(adata, n_obs=100):
#     return _pancreas_synthetic(adata, n_obs=n_obs, alpha=0.5)

if __name__ == "__main__":
    adata = ad.read_h5ad("resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad")
    pancreas = _pancreas_synthetic(adata, n_obs=100, alpha=1)
