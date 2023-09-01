#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'denoising|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/denoising/pancreas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# # split dataset
# viash run src/tasks/denoising/process_dataset/config.vsh.yaml -- \
#     --input $RAW_DATA \
#     --output_train $DATASET_DIR/train.h5ad \
#     --output_test $DATASET_DIR/test.h5ad \
#     --seed 123

# # run one method
# viash run src/tasks/denoising/methods/magic/config.vsh.yaml -- \
#     --input_train $DATASET_DIR/train.h5ad \
#     --output $DATASET_DIR/denoised.h5ad

# # run one metric
# viash run src/tasks/denoising/metrics/poisson/config.vsh.yaml -- \
#     --input_denoised $DATASET_DIR/denoised.h5ad \
#     --input_test $DATASET_DIR/test.h5ad \
#     --output $DATASET_DIR/score.h5ad

# run benchmark
export NXF_VER=22.04.5

nextflow \
  run . \
  -main-script src/tasks/denoising/workflows/run/main.nf \
  -profile docker \
  -resume \
  --id pancreas \
  --input_train $DATASET_DIR/train.h5ad \
  --input_test $DATASET_DIR/test.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/