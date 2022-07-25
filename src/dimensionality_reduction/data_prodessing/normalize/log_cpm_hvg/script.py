##VIASH START
par = {
    'input': "../../../../resources_test/label_projection/pancreas/toy_preprocessed_data.h5ad",
    'output': "output.h5ad"
}
##VIASH END

import scanpy as sc

print(">> Load data")
adata = sc.read(par['input'])

print(">> Normalize data")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")
sc.pp.log1p(adata)

if adata.n_vars < par['n_genes']:
    log.warning(
        f"Less than {n_genes} genes, setting 'n_genes' to {int(adata.n_vars * 0.5)}"
    )
    n_genes = int(adata.n_vars * 0.5)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]].copy()

adata.uns["normalization_method"] = "log_cpm"

print(">> Write data")
adata.write(par['output'], compression="gzip")
