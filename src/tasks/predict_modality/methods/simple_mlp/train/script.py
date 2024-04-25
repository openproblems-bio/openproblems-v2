import os
import logging
import anndata as ad

logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    # 'input_train_mod1': 'resources/predict_modality/datasets/openproblems_neurips2021/bmmc_multiome/normal/log_cp10k/train_mod1.h5ad',
    # 'input_train_mod2': 'resources/predict_modality/datasets/openproblems_neurips2021/bmmc_multiome/normal/log_cp10k/train_mod2.h5ad',
    # 'input_test_mod1': 'resources/predict_modality/datasets/openproblems_neurips2021/bmmc_multiome/normal/log_cp10k/test_mod1.h5ad',
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad',
    'output': 'output/model'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/simple_mlp',
    'cpus': 10
}
## VIASH END

resources_dir = f"{meta['resources_dir']}/resources"

import sys
sys.path.append(resources_dir)
from train import train

input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])

mod_1 = input_train_mod1.uns["modality"]
mod_2 = input_train_mod2.uns["modality"]

task = f'{mod_1}2{mod_2}'
yaml_path = f'{resources_dir}/yaml/mlp_{task}.yaml'

os.makedirs(par['output'], exist_ok=True)

if not os.path.exists(yaml_path):
    logging.warning(f"No configuration file found for task '{task}'. Skipping training.")
    sys.exit(0)

train(
    task=task,
    cp=resources_dir,
    wp=par['output'],
    tr1=input_train_mod1,
    tr2=input_train_mod2,
    num_workers=meta["cpus"]
)
