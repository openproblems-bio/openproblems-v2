#!/bin/bash

# fail on error
set -e

viash run src/tasks/predict_modality/process_dataset/config.vsh.yaml -- \
  --input_rna "resources/datasets/openproblems_neurips2021/bmmc_cite/log_cp10k/dataset_rna.had" \
  --input_other_mod "resources/datasets/openproblems_neurips2021/bmmc_cite/log_cp10k/dataset_other_mod.h5ad" \
  --output_train_mod1 "output/predict_modality/neurips2021_bmmc_cite/train_mod1.h5ad" \
  --output_train_mod2 "output/predict_modality/neurips2021_bmmc_cite/train_mod2.h5ad" \
  --output_test_mod1 "output/predict_modality/neurips2021_bmmc_cite/test_mod1.h5ad" \
  --output_test_mod2 "output/predict_modality/neurips2021_bmmc_cite/test_mod2.h5ad" \
  --swap FALSE \
  --seed 1 \