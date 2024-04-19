import sys
import anndata as ad
# from openproblems.tasks._cell_cell_communication._common.utils import ligand_receptor_resource
import openproblems as op
import pandas as pd

## VIASH START
par = {
    "input": "resources_test/cell_cell_communication_source_target/allen_brain_atlas",
    'output_train': 'train.h5ad',
    'output_test': 'test.h5ad',
    'output_solution': 'solution.h5ad'
}
meta = {
    "functionality_name": "split_data",
    "config": "src/tasks/cell_cell_communication_source_target/process_dataset/.config.vsh.yaml"
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_anndata import read_config_slots_info, subset_anndata

print(">> Load Data", flush=True)
adata = ad.read_h5ad(par["input"])

if 'ccc_target' in adata.uns:
    if 'source' in adata.uns['ccc_target'].columns and 'target' in adata.uns['ccc_target'].columns:
        adata.uns["ligand_receptor_resource"] = op.tasks._cell_cell_communication._common.utils.ligand_receptor_resource(
                adata.uns["target_organism"]
            )
        adata.uns['ccc_source'] = adata.uns['ligand_receptor_resource']['source']
        adata.uns['ccc_target'] = adata.uns['ligand_receptor_resource']['target']
        adata.uns['ccc_ligand'] = adata.uns['ligand_receptor_resource']['ligand_genesymbol']
        adata.uns['ccc_receptor'] = adata.uns['ligand_receptor_resource']['receptor_genesymbol']
        del adata.uns['ligand_receptor_resource']
        for key, value in adata.uns.items():
            if isinstance(value, pd.Series):
                arr = adata.uns[key].to_numpy()
                adata.uns[key] = arr

print(">> Figuring out which data needs to be copied to which output file", flush=True)
slot_info = read_config_slots_info(meta["config"])

print(">> Creating train data", flush=True)
output_train = subset_anndata(adata, slot_info["output_train"])

print(">> Creating test data", flush=True)
output_test = subset_anndata(adata, slot_info["output_test"])

print(">> Creating solution data", flush=True)
output_solution = subset_anndata(adata, slot_info["output_solution"])

if 'ccc_source' in adata.uns:
    print("ccc_source exists\n")

print(">> Writing", flush=True)
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
output_solution.write_h5ad(par["output_solution"])
