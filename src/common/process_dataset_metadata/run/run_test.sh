#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
INPUT_DIR="output/meta.yaml"
OUTPUT_DIR="output/temp/metadata"

# # temp sync
# aws s3 sync $INPUT_DIR output/temp

# start the run
NXF_VER=23.10.0 nextflow run . \
  -main-script target/nextflow/common/process_dataset_metadata/run/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  --id "process" \
  --input "$INPUT_DIR" \
  --publish_dir "$OUTPUT_DIR"

# cause quarto rerender to index page when in preview mode
# touch ../website/results/$TASK/index.qmd
