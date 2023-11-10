import anndata as ad

## VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    "output_spatial": "dataset.h5ad",
    "output_single_cell": "single_cell.h5ad",
    "output_solution": "solution.h5ad",
}
## VIASH END

print(">> Load common dataset", flush=True)
adata = ad.read_h5ad(par["input"])

def _pancreas_synthetic(
    alpha: float,
    test: bool = False,
    n_obs: int = 100,
    keep_techs: Optional[List[str]] = None,
):
    import scanpy as sc

    adata = load_pancreas(test=test, keep_techs=keep_techs or ["inDrop3"])
    sc.pp.filter_genes(adata, min_counts=10)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(
        adata, n_obs=n_obs, alpha=alpha, test=test
    )
    filter_genes_cells(merged_adata)
    return merged_adata


def pancreas_alpha_1(test=False, n_obs=100, keep_techs: Optional[List[str]] = None):
    return _pancreas_synthetic(test=test, n_obs=n_obs, alpha=1, keep_techs=keep_techs)

pancreas = pancreas_alpha_1()

print(">> Create dataset for methods", flush=True)
output_dataset = ad.AnnData(
  # <- add the data you wish to store in the output_dataset object here
)

print(">> Create solution object for metrics", flush=True)
output_solution = ad.AnnData(
  # <- add the data you wish to store in the output_solution object here
)

print(">> Write to disk", flush=True)
output_dataset.write_h5ad(par["output_dataset"])
output_solution.write_h5ad(par["output_solution"])