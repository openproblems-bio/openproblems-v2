#!/bin/bash

# make sure the following command has been executed
# viash ns build -q 'spatial_decomposition|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# RAW_DATA=resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad
# DATASET_DIR=resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/spatial_decomposition

mkdir -p $DATASET_DIR

# process dataset
# viash run src/tasks/spatial_decomposition/process_dataset/config.vsh.yml -- \
#     --input $RAW_DATA \
#     --output_spatial_masked $DATASET_DIR/spatial_masked.h5ad \
#     --output_single_cell $DATASET_DIR/single_cell_ref.h5ad \
#     --output_solution $DATASET_DIR/solution.h5ad \
#     --generate_dataset true \
#     --alpha 1 \
#     --simulated_data $DATASET_DIR/dataset_simulated.h5ad 

echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/spatial_decomposition/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_spatial_masked": "$id/spatial_masked.h5ad", "output_single_cell": "$id/single_cell_ref.h5ad", "output_solution": "$id/solution.h5ad", "generate_dataset": true, "alpha": 1, "simulated_data": "$id/dataset_simulated.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'

# run one method
viash run src/tasks/spatial_decomposition/methods/rctd/config.vsh.yaml -- \
    --input_single_cell $DATASET_DIR/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad \
    --input_spatial $DATASET_DIR/cxg_mouse_pancreas_atlas/spatial_masked.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/output.h5ad

# run one metric
viash run src/tasks/spatial_decomposition/metrics/r2/config.vsh.yaml -- \
    --input_method $DATASET_DIR/cxg_mouse_pancreas_atlas/output.h5ad \
    --input_solution $DATASET_DIR/cxg_mouse_pancreas_atlas/solution.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/score.h5ad
