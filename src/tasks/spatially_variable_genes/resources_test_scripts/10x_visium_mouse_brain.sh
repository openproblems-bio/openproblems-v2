#!/bin/bash

# make sure the following command has been executed
# viash ns build -q 'spatially_variable_genes|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/spatially_variable_genes

mkdir -p $DATASET_DIR

echo "Running process_dataset"
nextflow run . \
  -main-script target/nextflow/spatially_variable_genes/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input resources_test/common/10x_visium_mouse_brain/common_dataset.h5ad \
  --output_dataset 10x_visium_mouse_brain/dataset.h5ad \
  --output_solution 10x_visium_mouse_brain/solution.h5ad \
  --publish_dir resources_test/spatially_variable_genes

# run control method
viash run src/tasks/spatially_variable_genes/control_methods/true_ranking/config.vsh.yaml -- \
    --input_dataset $DATASET_DIR/10x_visium_mouse_brain/dataset.h5ad \
    --input_solution $DATASET_DIR/10x_visium_mouse_brain/solution.h5ad \
    --output $DATASET_DIR/10x_visium_mouse_brain/output.h5ad

# run one metric
viash run src/tasks/spatially_variable_genes/metrics/correlation/config.vsh.yaml -- \
    --input_method $DATASET_DIR/10x_visium_mouse_brain/output.h5ad \
    --input_solution $DATASET_DIR/10x_visium_mouse_brain/solution.h5ad \
    --output $DATASET_DIR/10x_visium_mouse_brain/score.h5ad
