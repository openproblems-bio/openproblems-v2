#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/label_projection/pancreas

mkdir -p $DATASET_DIR

target/docker/common/dataset_loader/download/download\
    --url "https://ndownloader.figshare.com/files/24539828"\
    --obs_celltype "celltype"\
    --obs_batch "tech"\
    --name "pancreas"\
    --output $DATASET_DIR/raw_data.h5ad

target/docker/label_projection/data_processing/subsample/subsample\
    --input $DATASET_DIR/raw_data.h5ad\
    --celltype_categories "acinar:beta"\
    --tech_categories "celseq:inDrop4:smarter"\
    --output $DATASET_DIR/toy_data.h5ad

target/docker/label_projection/data_processing/randomize/randomize\
    --input $DATASET_DIR/toy_data.h5ad\
    --output $DATASET_DIR/toy_preprocessed_data.h5ad

target/docker/label_projection/data_processing/normalize/log_cpm/log_cpm\
    --input $DATASET_DIR/toy_preprocessed_data.h5ad\
    --output $DATASET_DIR/toy_normalized_log_cpm_data.h5ad

target/docker/label_projection/data_processing/normalize/scran/log_scran_pooling/log_scran_pooling\
    --input $DATASET_DIR/toy_preprocessed_data.h5ad\
    --output $DATASET_DIR/toy_normalized_log_scran_pooling_data.h5ad

target/docker/label_projection/methods/baseline/majority_vote/majority_vote\
    --input $DATASET_DIR/toy_normalized_log_cpm_data.h5ad\
    --output $DATASET_DIR/toy_baseline_pred_data.h5ad
