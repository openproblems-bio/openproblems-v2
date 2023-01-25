#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/joint_embedding/datasets/openproblems_v1"
OUTPUT_DIR="resources/joint_embedding/benchmarks/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  python << HERE
import yaml
import os

dataset_dir = "$DATASETS_DIR"
output_dir = "$OUTPUT_DIR"

# read split datasets yaml
with open(dataset_dir + "/params.yaml", "r") as file:
  split_list = yaml.safe_load(file)
datasets = split_list['param_list']

# figure out where train/test files were stored
param_list = []

for dataset in datasets:
  id = dataset["id"]
  input_train = dataset_dir + "/" + id + ".train.h5ad"
  input_test = dataset_dir + "/" + id + ".test.h5ad"
  
  if os.path.exists(input_test):
    obj = {
      'id': id, 
    'id': id, 
      'id': id, 
      'dataset_id': dataset["dataset_id"],
      'input_train': input_train,
      'input_test': input_test
    }
    param_list.append(obj)

# write as output file
output = {
  "param_list": param_list,
}

with open(output_dir + "/params.yaml", "w") as file:
  yaml.dump(output, file)
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script src/joint_embedding/workflows/run/main.nf \
  -profile docker \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR" \
  -with-tower

bin/tools/docker/nextflow/process_log/process_log \
  --output "$OUTPUT_DIR/nextflow_log.tsv"