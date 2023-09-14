#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'batch_integration'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -xe

# run benchmark
export NXF_VER=22.04.5

  # -profile docker \
nextflow run . \
  -main-script src/tasks/batch_integration/workflows/process_datasets/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  -resume \
  --id pancreas \
  --input resources_test/common/pancreas/dataset.h5ad \
  --schema src/tasks/batch_integration/api/file_dataset.yaml \
  --publish_dir resources_test/batch_integration/pancreas_foo/ \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  $@