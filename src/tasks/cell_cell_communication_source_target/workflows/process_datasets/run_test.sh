#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'batch_integration'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/cell_cell_communication_source_target/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --id resources_test \
  --input_states "resources_test/common/**/state.yaml" \
  --rename_keys 'input_train:output_train,input_test:output_test,input_solution:output_solution' \
  --settings '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml", "output_task_info": "task_info.yaml"}' \
  --publish_dir "resources_test/cell_cell_communication_source_target"