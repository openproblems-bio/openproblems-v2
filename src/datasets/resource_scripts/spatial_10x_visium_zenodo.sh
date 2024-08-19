#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_10x_visium/human_heart_myocardial_infarction_1
    input_data: "https://zenodo.org/records/13328275/files/10X0018.h5ad?download=1"
    dataset_name: 10X Visium - Human Heart MI 1
    dataset_url: "https://www.nature.com/articles/s41586-022-05060-x"
    dataset_summary: Gene expression library of human heart using 10x Visium.
    dataset_description: "Frozen heart samples were embedded in OCT (Tissue-Tek) and cryosectioned (Thermo Cryostar). The 10-µm section was placed on the pre-chilled Optimization slides (Visium, 10X Genomics, PN-1000193) and the optimal lysis time was determined. The tissues were treated as recommended by 10X Genomics and the optimization procedure showed an optimal permeabilization time of 12 or 18 min of digestion and release of RNA from the tissue slide. Spatial gene expression slides (Visium, 10X Genomics, PN-1000187) were used for spatial transcriptomics following the Visium User Guides"
    dataset_reference: kuppe2022spatial
    dataset_organism: Homo sapiens
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type_proc: grid
    coord_type_moran_i: generic
    coord_type_sepal: grid
    max_neighs_speal: 6
    n_cp: -1

  - id: spatial_10x_visium/human_heart_myocardial_infarction_2
    input_data: "https://zenodo.org/records/13328275/files/10X009.h5ad?download=1"
    dataset_name: 10X Visium - Human Heart MI 2
    dataset_url: "https://www.nature.com/articles/s41586-022-05060-x"
    dataset_summary: Gene expression library of human heart using 10x Visium.
    dataset_description: "Frozen heart samples were embedded in OCT (Tissue-Tek) and cryosectioned (Thermo Cryostar). The 10-µm section was placed on the pre-chilled Optimization slides (Visium, 10X Genomics, PN-1000193) and the optimal lysis time was determined. The tissues were treated as recommended by 10X Genomics and the optimization procedure showed an optimal permeabilization time of 12 or 18 min of digestion and release of RNA from the tissue slide. Spatial gene expression slides (Visium, 10X Genomics, PN-1000187) were used for spatial transcriptomics following the Visium User Guides"
    dataset_reference: kuppe2022spatial
    dataset_organism: Homo sapiens
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type_proc: grid
    coord_type_moran_i: generic
    coord_type_sepal: grid
    max_neighs_speal: 6
    n_cp: -1

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
