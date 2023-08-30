#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'denoising|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

MOD_1_DATA=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad
MOD_2_DATA=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_mod2.h5ad
DATASET_DIR=resources_test/joint_embedding/bmmc_cite

if [ ! -f $MOD_1_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# split dataset
bin/viash run src/joint_embedding/mask_dataset/config.vsh.yaml -- \
    --input_mod1 $MOD_1_DATA \
    --input_mod2 $MOD_2_DATA \
    --output_mod1 $DATASET_DIR/cite_mod1.h5ad \
    --output_mod2 $DATASET_DIR/cite_mod2.h5ad \
    --output_solution $DATASET_DIR/cite_solution.h5ad

# run one method
bin/viash run src/joint_embedding/methods/pca/config.vsh.yaml -- \
    --input_mod1 $DATASET_DIR/cite_mod1.h5ad \
    --input_mod2 $DATASET_DIR/cite_mod2.h5ad \
    --output $DATASET_DIR/pca.h5ad

# run one metric
bin/viash run src/joint_embedding/metrics/ari/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/pca.h5ad \
    --input_solution $DATASET_DIR/cite_solution.h5ad \
    --output $DATASET_DIR/ari.h5ad

# run benchmark
export NXF_VER=22.04.5

bin/nextflow \
  run . \
  -main-script src/joint_embedding/workflows/run/main.nf \
  -profile docker \
  -resume \
  --id bmmc_cite \
  --dataset_id bmmc_site \
  --input_mod1 $DATASET_DIR/cite_mod1.h5ad \
  --input_mod2 $DATASET_DIR/cite_mod2.h5ad \
  --input_solution $DATASET_DIR/cite_solution.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/