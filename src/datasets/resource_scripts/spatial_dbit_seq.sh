#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_slide_tags/mouse_embryo
    input_data: "Slide-tags enables single-nucleus barcoding for multimodal spatial genomics"
    dataset_name: mouse_embryo
    dataset_url: "https://www.nature.com/articles/s41586-023-06837-4"
    dataset_summary: 
    dataset_description: ""
    dataset_organism: Mus musculus
    spot_filter_min_genes: 100
    gene_filter_min_spots: 500

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
  --main-script target/nextflow/datasets/workflows/process_10x_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
