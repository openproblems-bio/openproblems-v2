#!/bin/bash

# make sure the following command has been executed
# viash ns build -q 'spatial_decomposition|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad
DATASET_DIR=resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# generate synthetic spatial data
SIMULATED_DATA=$DATASET_DIR/dataset_simulated.h5ad
python3 src/tasks/spatial_decomposition/datasets/sample_datasets.py $RAW_DATA $SIMULATED_DATA

# process dataset
viash run src/tasks/spatial_decomposition/process_dataset/config.vsh.yml -- \
    --input $SIMULATED_DATA \
    --output_spatial_masked $DATASET_DIR/spatial_masked.h5ad \
    --output_single_cell $DATASET_DIR/single_cell_ref.h5ad \
    --output_solution $DATASET_DIR/solution.h5ad

# process dataset
# echo Running process_dataset
# nextflow run . \
#   -main-script target/nextflow/spatial_decomposition/workflows/process_datasets/main.nf \
#   -profile docker \
#   -entry auto \
#   --input_states "$RAW_DATA/**/state.yaml" \
#   --rename_keys 'input:output_dataset' \
#   --settings '{"output_spatial_masked": "$id/spatial_masked.h5ad", "output_single_cell": "$id/single_cell_ref.h5ad", "output_solution": "$id/solution.h5ad"}' \
#   --publish_dir "$DATASET_DIR" \
#   --output_state '$id/state.yaml'