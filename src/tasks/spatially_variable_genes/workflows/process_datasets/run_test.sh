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
  -c src/wf_utils/labels_ci.config \
  --id 10x_visium_mouse_brain \
  --input "resources_test/common/10x_visium_mouse_brain/dataset.h5ad" \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  --publish_dir "resources_test/spatially_variable_genes/10x_visium_mouse_brain" \
  --output_state "state.yaml"