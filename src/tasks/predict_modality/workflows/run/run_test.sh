#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/predict_modality/bmmc_cite_starter

# run benchmark
export NXF_VER=23.04.2

nextflow \
  run . \
  -main-script src/tasks/predict_modality/workflows/run/main.nf \
  -profile docker \
  -resume \
  --id pancreas \
  --input_train_mod1 $DATASET_DIR/train_mod1.h5ad \
  --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
  --input_test_mod1 $DATASET_DIR/test_mod1.h5ad \
  --input_test_mod2 $DATASET_DIR/test_mod2.h5ad \
  --output scores.tsv \
  --publish_dir output/label_projection/