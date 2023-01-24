#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash run src/datasets/loaders/openproblems_v1_multimodal/config.vsh.yaml -- \
  --id "scicar_mouse_kidney" \
  --obs_celltype "cell_name" \
  --obs_batch "replicate" \
  --layer_counts "counts" \
  --output_mod1 "foo_mod1.h5ad" \
  --output_mod2 "foo_mod2.h5ad"

# export NXF_VER=22.04.5

# nextflow \
#   run . \
#   -main-script target/nextflow/datasets/loaders/openproblems_v1/main.nf \
#   -resume \
#   -profile docker \
#   --param_list src/datasets/loaders/openproblems_v1/datasets.csv \
#   --publish_dir output/datasets