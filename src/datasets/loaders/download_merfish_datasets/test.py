import os
import subprocess
import anndata as ad

input_data ="https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.250.expand.rep1_data.h5ad?download=1"
dataset_id = "spatial_merfish/human_cortex_1"
dataset_name = "human_cortex_1"
dataset_url = "https://www.science.org/doi/10.1126/science.abm1741"
dataset_summary = "Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH)"
dataset_description = "Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH)"
dataset_organism = "Homo sapiens"
dataset = "dataset.h5ad"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta['executable'],
        "--input_data",  input_data, 
        "--dataset_id", dataset_id, 
        "--dataset_name", dataset_name, 
        "--dataset_url", dataset_url, 
        "--dataset_summary", dataset_summary, 
        "--dataset_description", dataset_description, 
        "--dataset_organism", dataset_organism, 
        "--dataset", dataset
    ],
    stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout, flush=True)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.", flush=True)
    exit(out.returncode)

print(">> Checking whether output file exists", flush=True)
assert os.path.exists(dataset), "Output does not exist"

print(">> Read output anndata", flush=True)
adata = ad.read_h5ad(dataset)

print(adata)

print(">> Check that output fits expected API", flush=True)
assert adata.X is None, "adata.X should be None/empty"
assert "counts" in adata.layers, "Counts layer not found in .layers"
assert adata.uns["dataset_id"] == dataset_id, f"Expected {dataset_id} as value"
assert adata.uns["dataset_name"] == dataset_name, f"Expected {dataset_name} as value"
assert adata.uns["dataset_url"] == dataset_url, f"Expected {dataset_url} as value"
assert adata.uns["dataset_summary"] == dataset_summary, f"Expected {dataset_summary} as value"
assert adata.uns["dataset_organism"] == dataset_organism, f"Expected {dataset_organism} as value"
assert 'spatial' in adata.obsm, "Spatial spot coordinates not found in .obsm"

print(">> All tests passed successfully", flush=True)
