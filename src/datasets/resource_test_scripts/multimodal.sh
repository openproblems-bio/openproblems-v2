#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/multimodal

set -e

mkdir -p $DATASET_DIR

# download dataset
viash run src/datasets/loaders/openproblems_v1_multimodal/config.vsh.yaml -- \
    --obs_tissue "source" \
    --layer_counts "counts" \
    --dataset_id scicar_cell_lines \
    --dataset_name "sci-CAR cell lines" \
    --data_url "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089" \
    --data_reference "cao2018joint" \
    --dataset_summary "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells" \
    --dataset_description "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells. Here, we use two sciCAR datasets that were obtained from the same study. The first dataset contains 4,825 cells from three cell lines (HEK293T cells, NIH/3T3 cells, and A549 cells) at multiple timepoints (0, 1 hour, 3 hours) after dexamethasone treatment. The second dataset contains 11,233 cells from wild-type adult mouse kidney." \
    --dataset_organism "[homo_sapiens, mus_musculus]" \
    --output_mod1 $DATASET_DIR/temp_mod1_full.h5ad \
    --output_mod2 $DATASET_DIR/temp_mod2_full.h5ad

# wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh_hm.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt
# wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh_hm.txt -O $DATASET_DIR/temp_s_genes_tirosh_hm.txt
# KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt $DATASET_DIR/temp_s_genes_tirosh_hm.txt | paste -sd ":" -`

# subsample
viash run src/datasets/processors/subsample/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_mod1_full.h5ad \
    --input_mod2 $DATASET_DIR/temp_mod2_full.h5ad \
    --n_obs 600 \
    --n_vars 1500 \
    --output $DATASET_DIR/raw_mod1.h5ad \
    --output_mod2 $DATASET_DIR/raw_mod2.h5ad
    --seed 123

# # run log cpm normalisation
# viash run src/datasets/normalization/log_cpm/config.vsh.yaml -- \
#     --input $DATASET_DIR/raw.h5ad \
#     --output $DATASET_DIR/normalized.h5ad

# # run pca
# viash run src/datasets/processors/pca/config.vsh.yaml -- \
#     --input $DATASET_DIR/normalized.h5ad \
#     --output $DATASET_DIR/pca.h5ad

# # run hvg
# viash run src/datasets/processors/hvg/config.vsh.yaml -- \
#     --input $DATASET_DIR/pca.h5ad \
#     --output $DATASET_DIR/hvg.h5ad

# # run knn
# viash run src/datasets/processors/knn/config.vsh.yaml -- \
#     --input $DATASET_DIR/hvg.h5ad \
#     --output $DATASET_DIR/dataset.h5ad

rm -r $DATASET_DIR/temp_*