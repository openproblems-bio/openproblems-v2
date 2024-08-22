#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_seqfish/mouse_organogenesis
    input_data: "https://zenodo.org/records/12785822/files/seqfish.h5ad?download=1"
    dataset_name: Seqfish - Mouse Organogenesis
    dataset_url: "https://www.nature.com/articles/s41587-021-01006-2"
    dataset_summary: Single-cell spatial expression of mouse organogenesis.
    dataset_description: "Sagittal sections from mouse embryo at the 8-12 ss was profiled by seqFISH."
    dataset_organism: Mus musculus
    dataset_reference: lohoff2021integration
    spot_filter_min_genes: 10
    gene_filter_min_spots: 10
    remove_mitochondrial: true

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: resources/datasets
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
