#!/bin/bash

# make sure the following command has been executed
# viash_build -q 'spatial_decomposition|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/spatial_decomposition/pancreas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# generate synthetic spatial data
python3 src/spatial_decomposition/datasets/sample_datasets.py $RAW_DATA > $DATASET_DIR/dataset_synthetic.h5ad

# split dataset
viash run src/tasks/label_projection/process_dataset/config.vsh.yaml -- \
    --input $DATASET_DIR/dataset_synthetic.h5ad \
    --output_spatial $DATASET_DIR/spatial.h5ad \
    --output_single_cell $DATASET_DIR/single_cell_ref.h5ad \
    --output_solution $DATASET_DIR/solution.h5ad \
