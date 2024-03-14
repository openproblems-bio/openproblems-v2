#!/bin/bash

viash run src/tasks/predict_modality/methods/novel/config.vsh.yaml -- \
    --input_train_mod1 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod1.h5ad' \
    --input_train_mod2 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod2.h5ad' \
    --input_test_mod1 'resources_test/predict_modality/neurips2021_bmmc_cite/test_mod1.h5ad' \
    --input_test_mod2 'resources_test/predict_modality/neurips2021_bmmc_cite/test_mod2.h5ad' \
    --pretrain 'model.pt' \
    --output 'output/novel_test.h5ad'