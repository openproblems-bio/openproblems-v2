import os
import subprocess
import anndata as ad

input_data ="https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_5.6.h5ad?download=1"
dataset_id = "spatial_stereo_seq/drosophila_embryo_e5_6"
dataset_name = "drosophila_embryo_e5_6"
dataset_url = "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
dataset_summary = "Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution"
dataset_description = "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity"
dataset_organism = "Drosophila"
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
