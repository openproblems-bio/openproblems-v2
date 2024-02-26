import anndata as ad
import sys 

## VIASH START
par = {
    "input": "resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad",
    # "input": "resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/dataset_simulated.h5ad",
    "output_spatial_masked": "spatial_masked.h5ad",
    "output_single_cell": "single_cell_ref.h5ad",
    "output_solution": "solution.h5ad",
    "generate_dataset": True, 
    "alpha": 1,
    "simulated_data": "dataset_simulated.h5ad"
}
meta = {
    "functionality_name": "process_dataset",
    "resources_dir": "src/tasks/spatial_decomposition/process_dataset",
    "config": "target/nextflow/spatial_decomposition/process_dataset/.config.vsh.yaml"
}
## VIASH END

sys.path.append(meta['resources_dir'])
from subset_anndata import read_config_slots_info, subset_anndata

print(">> Load dataset", flush=True)

if par['generate_dataset']:
    from sample_datasets import _pancreas_synthetic
    raw_input = ad.read_h5ad(par["input"])
    adata = _pancreas_synthetic(raw_input, n_obs=100, alpha=par['alpha'])
    adata.write_h5ad(par['simulated_data'], compression='gzip')
else:
    adata = ad.read_h5ad(par["input"])

print(">> Figuring out which data needs to be copied to which output file", flush=True)
slot_info = read_config_slots_info(meta["config"])

print(">> Split dataset by modality", flush=True)
is_sp = adata.obs["modality"] == "sp"
adata_sp = adata[is_sp, :].copy()
adata_sc = adata[~is_sp, :].copy()

print(">> Create dataset for methods", flush=True)
output_spatial_masked = subset_anndata(adata_sp, slot_info['output_spatial_masked'])
output_single_cell = subset_anndata(adata_sc, slot_info['output_single_cell'])

print(">> Create solution object for metrics", flush=True)
output_solution = subset_anndata(adata_sp, slot_info['output_solution'])

print(">> Write to disk", flush=True)
output_spatial_masked.write_h5ad(par["output_spatial_masked"])
output_single_cell.write_h5ad(par["output_single_cell"])
output_solution.write_h5ad(par["output_solution"])
