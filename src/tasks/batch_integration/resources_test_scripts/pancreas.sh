#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/batch_integration

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
export NXF_VER=22.04.5
nextflow run . \
  -main-script target/nextflow/batch_integration/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c src/wf_utils/labels_ci.config \
  -entry auto \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input_dataset:output_dataset,input_solution:output_solution' \
  --settings '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state "state.yaml"
# output_state should be moved to settings once workaround is solved

# echo Running BBKNN
# viash run src/tasks/batch_integration/methods/bbknn/config.vsh.yaml -- \
#   --input $DATASET_DIR/pancreas/dataset.h5ad \
#   --output $DATASET_DIR/pancreas/integrated_graph.h5ad

# echo Running SCVI
# viash run src/tasks/batch_integration/methods/scvi/config.vsh.yaml -- \
#   --input $DATASET_DIR/pancreas/dataset.h5ad \
#   --output $DATASET_DIR/pancreas/integrated_embedding.h5ad

# echo Running combat
# viash run src/tasks/batch_integration/methods/combat/config.vsh.yaml -- \
#   --input $DATASET_DIR/pancreas/dataset.h5ad \
#   --output $DATASET_DIR/pancreas/integrated_feature.h5ad