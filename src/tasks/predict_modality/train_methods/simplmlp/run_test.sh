#!/bin/bash

DATASET_DIR=resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome

viash run src/tasks/predict_modality/train_methods/simplmlp/config.vsh.yaml -- \
    --input_train_mod1 $DATASET_DIR/train_mod1.h5ad \
    --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --output 'output/models'