#!/bin/bash

NEURIPS2021_URL="https://github.com/openproblems-bio/neurips2021_multimodal_viash/raw/main/resources_test/common"
DATASET_DIR="resources_test/common"

SUBDIR="$DATASET_DIR/bmmc_cite_starter"
mkdir -p "$SUBDIR"
wget "$NEURIPS2021_URL/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad" \
  -O "$SUBDIR/dataset_rna.h5ad"
wget "$NEURIPS2021_URL/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_mod2.h5ad" \
  -O "$SUBDIR/dataset_adt.h5ad"

cat > "$SUBDIR/state.yaml" << HERE
id: bmmc_cite_starter
output_dataset_rna: !file dataset_rna.h5ad
output_dataset_other_mod: !file dataset_adt.h5ad
HERE

SUBDIR="$DATASET_DIR/bmmc_multiome_starter"
mkdir -p "$SUBDIR"
wget "$NEURIPS2021_URL/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad" \
  -O "$SUBDIR/dataset_rna.h5ad"
wget "$NEURIPS2021_URL/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_mod2.h5ad" \
  -O "$SUBDIR/dataset_atac.h5ad"

cat > "$SUBDIR/state.yaml" << HERE
id: bmmc_multiome_starter
output_dataset_rna: !file dataset_rna.h5ad
output_dataset_other_mod: !file dataset_atac.h5ad
HERE

# run task process dataset components
src/tasks/predict_modality/resources_test_scripts/bmmc_x_starter.sh