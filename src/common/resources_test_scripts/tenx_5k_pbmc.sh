#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/dimensionality_reduction/tenx_5k_pbmc

mkdir -p $DATASET_DIR

target/docker/common/dataset_loader/download/download\
    --url "https://ndownloader.figshare.com/files/25555739"\
    --name "tenx_5k_pbmc"\
    --output $DATASET_DIR/raw_data.h5ad

target/docker/common/data_processing/subsample/subsample\
    --input $DATASET_DIR/raw_data.h5ad\
    --output $DATASET_DIR/toy_data.h5ad

target/docker/dimensionality_reduction/data_processing/normalize/log_cpm_hvg/log_cpm_hvg \
    --input $DATASET_DIR/toy_data.h5ad \
    --output $DATASET_DIR/toy_log_cpm_hvg.h5ad
