#!/bin/bash

# fail on error
set -e

viash run src/datasets/loaders/openproblems_neurips2021_bmmc/config.vsh.yaml -- \
  --input resources_test/common/openproblems_neurips2021/neurips2021_bmmc_cite.h5ad \
  --mod1 GEX \
  --mod2 ADT \
  --dataset_name neurips2021_bmmc_cite \
  --dataset_url https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122 \
  --dataset_reference Neurips \
  --dataset_summary value \
  --dataset_description value \
  --dataset_organism homo_sapiens \
  --output_mod1 output/mod1_gex.h5ad \
  --output_mod2 output/mod2_adt.h5ad