#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_dbit_seq/mouse_e10_brain
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_brain_gene_25um_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Brain (E10)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E10 whole mouse embryo tissue (brain in early-stage organogenesis) profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_dbit_seq/mouse_e10_eye
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_eye_and_nearby_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Eye (E10)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E10 whole mouse embryo tissue (eye in early-stage organogenesis) profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_dbit_seq/mouse_e10_whole_body
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_whole_gene_best_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Whole Body (E10)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E10 whole mouse embryo tissue profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_dbit_seq/mouse_e11_lower_body
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E11_lower_body_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Lower Body (E11)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E11 whole mouse embryo tissue (lower body in early-stage organogenesis) profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_dbit_seq/mouse_e11_1
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_GSM4364244_E11-FL-1L_gene_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Whole Body 1 (E11)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E11 whole mouse embryo tissue profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_dbit_seq/mouse_e11_2
    input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_GSM4364245_E11-FL-2L_gene_data.h5ad?download=1"
    dataset_name: DBiT-seq - Mouse Whole Body 2 (E11)
    dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
    dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
    dataset_description: "Gene expression library of an E11 whole mouse embryo tissue profiled using DBiT-seq."
    dataset_organism: Mus musculus
    dataset_reference: liu2020high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 200
    select_top_variable_genes: 50
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
