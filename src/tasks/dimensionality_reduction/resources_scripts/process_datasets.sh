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
  -main-script target/nextflow/dimensionality_reduction/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --id run_test \
  --input_states "resources/common/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "dataset.h5ad", "output_solution": "solution.h5ad"}' \
  --publish_dir "resources/dimensionality_reduction"