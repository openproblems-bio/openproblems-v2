#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=src/label_projection/starter_kits/resources/

mkdir -p $DATASET_DIR

target/docker/common/generate_toy_dataset/generate_toy_dataset\
    --output $DATASET_DIR/toy_data.h5ad
