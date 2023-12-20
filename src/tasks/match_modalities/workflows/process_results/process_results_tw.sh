#!/bin/bash

DATASETS_DIR="s3://openproblems-data/resources/match_modalities/results/"

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: match_modalities
input_scores: "$DATASETS_DIR/scores.yaml"
input_dataset_info: "$DATASETS_DIR/dataset_info.yaml"
input_method_configs: "$DATASETS_DIR/method_configs.yaml"
input_metric_configs: "$DATASETS_DIR/metric_configs.yaml"
input_execution: "$DATASETS_DIR/trace.txt"
input_task_info: "$DATASETS_DIR/task_info.yaml"
task_id: "match_modalities" 
output_scores: "results.json"
output_method_info: "method_info.json"
output_metric_info: "metric_info.json"
output_dataset_info: "dataset_info.json"
output_task_info: "task_info.json"
publish_dir: $DATASETS_DIR
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}


HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/common/workflows/transform_meta/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config