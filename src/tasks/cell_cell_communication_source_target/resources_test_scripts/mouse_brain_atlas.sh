#!/bin/bash
#make sure the following command has been executed
#viash ns build -q 'cell_cell_communication_source_target|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/cell_cell_communication_source_target

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/cell_cell_communication_source_target/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_train": "$id/train.h5ad", "output_test": "$id/test.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved


# # run one method
# viash run src/tasks/cell_cell_communication_source_target/methods/method_to_be_added/config.vsh.yaml -- \
#     --input_train $DATASET_DIR/mouse_brain_atlas/train.h5ad \
#     --input_test $DATASET_DIR/mouse_brain_atlas/test.h5ad \
#     --output $DATASET_DIR/mouse_brain_atlas/prediction.h5ad

# # run one metric
# viash run src/tasks/cell_cell_communication_source_target/metrics/metric_to_be_added/config.vsh.yaml -- \
#     --input_prediction $DATASET_DIR/mouse_brain_atlas/prediction.h5ad \
#     --input_solution $DATASET_DIR/mouse_brain_atlas/solution.h5ad \
#     --output $DATASET_DIR/mouse_brain_atlas/score.h5ad

# # run benchmark
# export NXF_VER=22.04.5

# # after having added a split dataset component
# nextflow \
#   run . \
#   -main-script src/tasks/dimensionality_reduction/workflows/run/main.nf \
#   -profile docker \
#   --id pancreas \
#   --input_dataset $DATASET_DIR/dataset.h5ad \
#   --input_solution $DATASET_DIR/solution.h5ad \
#   --output scores.tsv \
#   --publish_dir $DATASET_DIR/