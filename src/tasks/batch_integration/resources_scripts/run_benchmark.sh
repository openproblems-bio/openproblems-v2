#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/batch_integration/datasets/openproblems_v1"
OUTPUT_DIR="resources/batch_integration/benchmarks/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script src/tasks/batch_integration/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  --id resources \
  --input_dir "$DATASETS_DIR" \
  --rename_keys 'input_dataset:output_dataset,input_solution:output_solution' \
  --settings '{"output": "scores.tsv"}' \
  --publish_dir "$OUTPUT_DIR"