#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/openproblems_v1_multimodal"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  cat > "$params_file" << 'HERE'
param_list:
  - id: citeseq_cbmc
    dataset_id: citeseq_cbmc
    layer_counts: counts

  - id: scicar_cell_lines
    dataset_id: scicar_cell_lines
    obs_celltype: cell_name
    layer_counts: counts

  - id: scicar_mouse_kidney
    dataset_id: scicar_mouse_kidney
    obs_celltype: cell_name
    obs_batch: replicate
    layer_counts: counts

output_dataset_mod1: '$id/dataset_mod1.h5ad'
output_dataset_mod1: '$id/dataset_mod2.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod1.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod2.h5ad'
output_state: '$id/state.yaml'
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1_multimodal/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
