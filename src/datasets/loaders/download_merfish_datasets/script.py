import subprocess
import tempfile
import scanpy as sc

# VIASH START
par = {
    "input_data": "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.250.expand.rep1_data.h5ad?download=1",
    "dataset_id": "spatial_merfish/human_cortex_1",
    "dataset_name": "Human Cortex",
    "dataset_url": "https://www.science.org/doi/10.1126/science.abm1741",
    "dataset_summary": "Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH)",
    "dataset_organism": "Homo sapiens",
    "dataset": "dataset.h5ad",
    "spot_filter_min_genes": 10,
    "gene_filter_min_spots": 100,
    "remove_mitochondrial": True
}
meta = {
    "functionality_name": "download_merfish_datasets"
}
# VIASH END

print(f"Downloading data", flush=True)
with tempfile.TemporaryDirectory() as tempdir:
    input_data = "input_data.h5ad"
    epx_data = subprocess.run(["wget", "-O", f"{tempdir}/{input_data}", par['input_data']], stderr=subprocess.STDOUT)
    adata = sc.read_h5ad(filename=f"{tempdir}/{input_data}")

# Make variable names unique
adata.var_names_make_unique()

sc.pp.calculate_qc_metrics(adata, inplace=True)

if par["remove_mitochondrial"]:
    print("Removing mitochondrial genes")
    non_mito_genes_list = [name for name in adata.var_names if not (
        name.startswith('MT-') or name.startswith('mt-'))]
    adata = adata[:, non_mito_genes_list]

print("Filtering spots or genes")
t0 = adata.shape
# remove cells with few counts
if par["spot_filter_min_counts"]:
    sc.pp.filter_cells(
        adata, min_counts=par["spot_filter_min_counts"], inplace=True)

# remove cells with few genes
if par["spot_filter_min_genes"]:
    sc.pp.filter_cells(
        adata, min_genes=par["spot_filter_min_genes"], inplace=True)

# remove genes that have few counts
if par["gene_filter_min_counts"]:
    sc.pp.filter_genes(
        adata, min_counts=par["gene_filter_min_counts"], inplace=True)

# remove genes that are found in few cells
if par["gene_filter_min_spots"]:
    sc.pp.filter_genes(
        adata, min_cells=par["gene_filter_min_spots"], inplace=True)

t1 = adata.shape
print(f"Removed {t0[0] - t1[0]} cells and {(t0[1] - t1[1])} genes.")

# Rename .var columns
adata.var['feature_name'] = adata.var_names

# Move counts to .layers
print("Add metadata to uns", flush=True)
adata.layers["counts"] = adata.X
adata.X = None

# Add metadata
print("Add metadata to uns", flush=True)
metadata_fields = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for key in metadata_fields:
    if key in par:
        print(f"Setting .uns['{key}']", flush=True)
        adata.uns[key] = par[key]

print("Writing adata to file", flush=True)
adata.write_h5ad(par["dataset"], compression="gzip")
