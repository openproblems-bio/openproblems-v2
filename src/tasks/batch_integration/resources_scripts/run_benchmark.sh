#!/bin/bash


# try running on nf tower
cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/batch_integration/datasets/**/state.yaml
rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
settings: '{"output": "scores.tsv"}'
output_state: "state.yaml"
publish_dir: s3://openproblems-data/resources/batch_integration/results
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/batch_integration/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config