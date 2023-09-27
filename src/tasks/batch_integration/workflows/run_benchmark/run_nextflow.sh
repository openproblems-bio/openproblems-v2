#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'batch_integration'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -xe

DATASET_DIR=resources_test/batch_integration

# run benchmark
export NXF_VER=22.04.5

  # -profile docker \
nextflow run . \
  -main-script src/tasks/batch_integration/workflows/run_benchmark/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  -resume \
  --id foo \
  --input_states "$DATASET_DIR/**/state.yaml" \
  --publish_dir "output"