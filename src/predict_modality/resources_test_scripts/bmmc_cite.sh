#!/bin/bash
#
#make sure the following command has been executed
#viash ns build -q 'match_modality|common' --parallel --setup cb

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

MOD_1_DATA=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad
MOD_2_DATA=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_mod2.h5ad
DATASET_DIR=resources_test/predict_modality/bmmc_cite

if [ ! -f $MOD_1_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# maskdataset
viash run src/predict_modality/mask_dataset/config.vsh.yaml -- \
    --input_mod1 $MOD_1_DATA \
    --input_mod2 $MOD_2_DATA \
    --output_train_mod1 $DATASET_DIR/cite_train_mod1.h5ad \
    --output_train_mod2 $DATASET_DIR/cite_train_mod2.h5ad \
    --output_test_mod1 $DATASET_DIR/cite_test_mod1.h5ad \
    --output_solution $DATASET_DIR/cite_solution.h5ad

# run one method
viash run src/predict_modality/methods/knnr_py/config.vsh.yaml -- \
    --input_train_mod1 $DATASET_DIR/cite_train_mod1.h5ad \
    --input_train_mod2 $DATASET_DIR/cite_train_mod2.h5ad \
    --input_test_mod1 $DATASET_DIR/cite_test_mod1.h5ad \
    --output $DATASET_DIR/knnr_py.h5ad

# run one metric
viash run src/predict_modality/metrics/mse/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/knnr_py.h5ad \
    --output_solution $DATASET_DIR/cite_solution.h5ad
    --output $DATASET_DIR/mse.h5ad

# run benchmark
export NXF_VER=22.04.5

nextflow \
  run . \
  -main-script src/predict_modality/workflows/run/main.nf \
  -profile docker \
  --id bmmc_cite \
  --dataset_id bmmc_site \
  --input_train_mod1 $DATASET_DIR/cite_train_mod1.h5ad \
  --input_train_mod2 $DATASET_DIR/cite_train_mod2.h5ad \
  --input_train_sol $DATASET_DIR/cite_train_sol.h5ad \
  --input_test_mod1 $DATASET_DIR/cite_test_mod1.h5ad \
  --input_test_mod2 $DATASET_DIR/cite_test_mod2.h5ad \
  --output_solution $DATASET_DIR/cite_solution.h5ad
  --output scores.tsv \
  --publish_dir $DATASET_DIR/