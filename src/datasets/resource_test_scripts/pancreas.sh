#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/pancreas

set -e

mkdir -p $DATASET_DIR

# download dataset
viash run src/datasets/loaders/openproblems_v1/config.vsh.yaml -- \
    --obs_celltype "celltype" \
    --obs_batch "tech" \
    --layer_counts "counts" \
    --dataset_id pancreas \
    --dataset_name "Human pancreas" \
    --data_url "https://theislab.github.io/scib-reproducibility/dataset_pancreas.html" \
    --data_reference "luecken2022benchmarking" \
    --dataset_summary "Human pancreas cells dataset from the scIB benchmarks" \
    --dataset_description "Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq)." \
    --dataset_organism "homo_sapiens" \
    --output $DATASET_DIR/temp_dataset_full.h5ad

wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh_hm.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt
wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh_hm.txt -O $DATASET_DIR/temp_s_genes_tirosh_hm.txt
KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt $DATASET_DIR/temp_s_genes_tirosh_hm.txt | paste -sd ":" -`

# subsample
viash run src/datasets/processors/subsample/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset_full.h5ad \
    --keep_celltype_categories "acinar:beta" \
    --keep_batch_categories "celseq:inDrop4:smarter" \
    --keep_features "$KEEP_FEATURES" \
    --output $DATASET_DIR/raw.h5ad \
    --seed 123

# run log cpm normalisation
viash run src/datasets/normalization/log_cp/config.vsh.yaml -- \
    --input $DATASET_DIR/raw.h5ad \
    --output $DATASET_DIR/normalized.h5ad

# run pca
viash run src/datasets/processors/pca/config.vsh.yaml -- \
    --input $DATASET_DIR/normalized.h5ad \
    --output $DATASET_DIR/pca.h5ad

# run hvg
viash run src/datasets/processors/hvg/config.vsh.yaml -- \
    --input $DATASET_DIR/pca.h5ad \
    --output $DATASET_DIR/hvg.h5ad

# run knn
viash run src/datasets/processors/knn/config.vsh.yaml -- \
    --input $DATASET_DIR/hvg.h5ad \
    --output $DATASET_DIR/dataset.h5ad

# run log cp10k normalisation
viash run src/datasets/normalization/log_cp/config.vsh.yaml -- \
    --input $DATASET_DIR/raw.h5ad \
    --n_cp 10000 \
    --output $DATASET_DIR/cp10k_normalized.h5ad

# run pca
viash run src/datasets/processors/pca/config.vsh.yaml -- \
    --input $DATASET_DIR/cp10k_normalized.h5ad \
    --output $DATASET_DIR/cp10k_pca.h5ad

# run hvg
viash run src/datasets/processors/hvg/config.vsh.yaml -- \
    --input $DATASET_DIR/cp10k_pca.h5ad \
    --output $DATASET_DIR/cp10k_hvg.h5ad

# run knn
viash run src/datasets/processors/knn/config.vsh.yaml -- \
    --input $DATASET_DIR/cp10k_hvg.h5ad \
    --output $DATASET_DIR/cp10k_dataset.h5ad

rm -r $DATASET_DIR/temp_*