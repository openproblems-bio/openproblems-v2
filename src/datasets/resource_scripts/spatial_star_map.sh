#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_star_map/mouse_brain_2d_zstep10_0
    input_data: "https://zenodo.org/records/12785822/files/STARmap_Wang2018three_data_2D_zstep10_0_data.h5ad?download=1"
    dataset_name: STARmap - Mouse Brain 1
    dataset_url: "https://www.science.org/doi/10.1126/science.aat5691"
    dataset_summary: Three-dimensional intact-tissue sequencing of single-cell transcriptional states.
    dataset_description: "3D architecture of cell types in visual cortex volumes."
    dataset_organism: Mus musculus
    dataset_reference: wang2018three
    spot_filter_min_genes: 1
    gene_filter_min_spots: 1
    num_reference_genes: 25
    select_top_variable_genes: 25
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_star_map/mouse_brain_2d_zstep15_0
    input_data: "https://zenodo.org/records/12785822/files/STARmap_Wang2018three_data_2D_zstep15_0_data.h5ad?download=1"
    dataset_name: STARmap - Mouse Brain 2
    dataset_url: "https://www.science.org/doi/10.1126/science.aat5691"
    dataset_summary: Three-dimensional intact-tissue sequencing of single-cell transcriptional states.
    dataset_description: "3D architecture of cell types in visual cortex volumes."
    dataset_organism: Mus musculus
    dataset_reference: wang2018three
    spot_filter_min_genes: 1
    gene_filter_min_spots: 1
    num_reference_genes: 25
    select_top_variable_genes: 25
    remove_mitochondrial: true
    coord_type: generic

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: s3://openproblems-data/resources/datasets
remove_mitochondrial: true
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_spatial_from_zenodo/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
