#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/match_modalities

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/match_modalities/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input_mod1:output_dataset_mod1,input_mod2:output_dataset_mod2' \
  --settings '{"output_mod1": "$id/output_mod1.h5ad", "output_mod2": "$id/output_mod2.h5ad", "output_solution_mod1": "$id/output_solution_mod1.h5ad", "output_solution_mod2": "$id/output_solution_mod2.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved

# # run one method
# viash run src/tasks/match_modalities/methods/knn/config.vsh.yaml -- \
#     --input_train $DATASET_DIR/train.h5ad \
#     --input_test $DATASET_DIR/test.h5ad \
#     --output $DATASET_DIR/knn.h5ad

# # run one metric
# viash run src/tasks/match_modalities/metrics/accuracy/config.vsh.yaml -- \
#     --input_prediction $DATASET_DIR/knn.h5ad \
#     --input_solution $DATASET_DIR/solution.h5ad \
#     --output $DATASET_DIR/knn_accuracy.h5ad

# # # run benchmark
# # export NXF_VER=22.04.5

# # nextflow \
# #   run . \
# #   -main-script src/tasks/match_modalities/workflows/run/main.nf \
# #   -profile docker \
# #   -resume \
# #   --id pancreas \
# #   --input_train $DATASET_DIR/train.h5ad \
# #   --input_test $DATASET_DIR/test.h5ad \
# #   --input_solution $DATASET_DIR/solution.h5ad \
# #   --output scores.tsv \
# #   --publish_dir $DATASET_DIR/