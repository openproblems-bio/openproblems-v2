#!/bin/bash

# make sure the following command has been executed
# viash ns build -q 'spatial_decomposition|common'

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
SYNTHETIC_DATA=$DATASET_DIR/dataset_synthetic.h5ad
python3 src/tasks/spatial_decomposition/datasets/sample_datasets.py $RAW_DATA $SYNTHETIC_DATA

# split dataset
viash run src/tasks/label_projection/process_dataset/config.vsh.yaml -- \
    --input $SYNTHETIC_DATA \
    --output_spatial $DATASET_DIR/spatial.h5ad \
    --output_single_cell $DATASET_DIR/single_cell_ref.h5ad \
    --output_solution $DATASET_DIR/solution.h5ad \

