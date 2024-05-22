#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'spatially_variable_genes'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

nextflow run . \
  -main-script target/nextflow/spatially_variable_genes/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input resources_test/common/10x_visium_mouse_brain/common_dataset.h5ad \
  --output_dataset 10x_visium_mouse_brain/dataset.h5ad \
  --output_solution 10x_visium_mouse_brain/solution.h5ad \
  --publish_dir resources_test/spatially_variable_genes