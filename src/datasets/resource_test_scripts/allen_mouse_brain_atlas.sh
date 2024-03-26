#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common

set -e

mkdir -p $DATASET_DIR

# wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh_hm.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt
# wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh_hm.txt -O $DATASET_DIR/temp_s_genes_tirosh_hm.txt
# KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt $DATASET_DIR/temp_s_genes_tirosh_hm.txt | paste -sd ":" -`

# download dataset
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  -resume \
  --id allen_brain_atlas \
  --input_id allen_brain_atlas \
  --obs_cell_type "celltype" \
  --obs_batch "tech" \
  --var_feature_name "index" \
  --layer_counts "counts" \
  --dataset_name "Mouse Brain Atlas" \
  --dataset_url "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585" \
  --dataset_reference "tasic2016adult" \
  --dataset_summary "Adult mouse primary visual cortex" \
  --dataset_description "A murine brain atlas with adjacent cell types as assumed benchmark truth, inferred from deconvolution proportion correlations using matching 10x Visium slides (see Dimitrov et al., 2022)." \
  --dataset_organism "mus_musculus" \
  --keep_cell_type_categories "acinar:beta" \
  --keep_batch_categories "celseq:inDrop4:smarter" \
  --seed 123 \
  --normalization_methods log_cp10k \
  --do_subsample true \
  --output_raw '$id/raw.h5ad' \
  --output_normalized '$id/normalized.h5ad' \
  --output_hvg '$id/hvg.h5ad' \
  --output_pca '$id/pca.h5ad' \
  --output_knn '$id/knn.h5ad' \
  --output_dataset '$id/dataset.h5ad' \
  --output_meta '$id/dataset_meta.yaml' \
  --output_state '$id/state.yaml' \
  --publish_dir "$DATASET_DIR"

# rm -r $DATASET_DIR/temp_*

# # run task process dataset components
# src/tasks/cell_cell_communication_source_target/resources_test_scripts/mouse_brain_atlas.sh