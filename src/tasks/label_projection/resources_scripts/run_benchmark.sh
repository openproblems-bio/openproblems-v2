#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/label_projection/results/${RUN_ID}"

cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/label_projection/datasets/**/state.yaml
rename_keys: 'input_train:output_train,input_test:output_test,input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

cat src/wf_utils/labels_tw.config > /tmp/nextflow.config

cat >> /tmp/nextflow.config << HERE
trace {
    enabled = true
    overwrite = true
    file = "$publish_dir/trace.txt"
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/label_projection/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  --labels label_projection,full