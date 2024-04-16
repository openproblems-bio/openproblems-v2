import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

## VIASH START
par = {
  "input_mod1": "cite_rna_merged.h5ad",
  "input_mod2": "cite_prot_merged.h5ad",
  "mod1": "GEX",
  "mod2": "ATAC",
  "dataset_id": "openproblems/neurips2022_pbmc",
  "dataset_name": "Kaggle22 PBMC (CITE-seq)",
  "dataset_url_mod1": "s3://openproblems-nextflow/datasets_private/neurips2022/cite_rna_merged.h5ad",
  "dataset_url_mod2": "s3://openproblems-nextflow/datasets_private/neurips2022/cite_prot_merged.h5ad",
  "dataset_reference": "Neurips22",
  "dataset_summary": "Neurips22 competition dataset",
  "dataset_description": "value",
  "dataset_organism": "homo_sapiens",
  "output_mod1": "output/mod1.h5ad",
  "output_mod2": "output/mod2.h5ad"
}
meta = {
  "functionality_name": "openproblems_neurips2022_pbmc",
  "resources_dir": "/tmp/viash_inject_openproblems_neurips2021_bmmc14365472827677740971", # to be adjusted?
}
## VIASH END


def convert_matrix(adata):
  for key in adata:
      if isinstance(adata[key], sparse.csr_matrix):
        adata[key] = sparse.csc_matrix(adata[key])
      

print("load dataset modality 1 file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])

print("load dataset modality 2 file", flush=True)
adata_mod2 = ad.read_h5ad(par["input_mod2"])

# Convert to sparse csc_matrix
convert_matrix(adata_mod1.layers)
convert_matrix(adata_mod1.obsm)
convert_matrix(adata_mod2.layers)
convert_matrix(adata_mod2.obsm)


print("Add metadata to uns", flush=True)
metadata_fields = [
    "dataset_id", "dataset_name", "dataset_url", "dataset_reference",
    "dataset_summary", "dataset_description", "dataset_organism"
]
for key in metadata_fields:
    if key in par:
        print(f"  Setting .uns['{key}']", flush=True)
        adata_mod1.uns[key] = par[key]
        adata_mod2.uns[key] = par[key]


print("Writing adata to file", flush=True)
adata_mod1.write_h5ad(par["output_mod1"], compression="gzip")
adata_mod2.write_h5ad(par["output_mod2"], compression="gzip")




