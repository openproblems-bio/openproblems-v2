#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_slideseq_v2/mouse_olfactory_bulb_puck
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_SlideSeqV2_Mouse_Olfactory_bulb_Puck_200127_15_data_whole.h5ad?download=1"
    dataset_name: Slide-seqV2 - Mouse Olfactory Bulb Puck
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
    dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
    dataset_description: "Gene expression library of mouse olfactory bulk puck profiled using Slide-seq V2."
    dataset_reference: stickels2020highly
    dataset_organism: Mus musculus
    spot_filter_min_genes: 100
    gene_filter_min_spots: 500

  - id: spatial_slideseq_v2/mouse_cortex
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_palla2021squidpy_Slide-seqV2_Mouse_Cortex_data_whole.h5ad?download=1"
    dataset_name: Slide-seqV2 - Mouse Cortex
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
    dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
    dataset_description: "Gene expression library of Mouse cortex profiled using Slide-seq V2."
    dataset_reference: stickels2020highly
    dataset_organism: Mus musculus
    spot_filter_min_genes: 10
    gene_filter_min_spots: 500

  - id: spatial_slideseq_v2/mouse_cerebellum
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_Cerebellum_SCP948_data_whole.h5ad?download=1"
    dataset_name: Slide-seqV2 - Mouse Cerebellum
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
    dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
    dataset_description: "Gene expression library of mouse cerebellum profiled using Slide-seq V2."
    dataset_reference: stickels2020highly
    dataset_organism: Mus musculus
    spot_filter_min_genes: 100
    gene_filter_min_spots: 500

  - id: spatial_slideseq_v2/mouse_hippocampus_puck
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_Hippocampus_Puck_200115_08_data_whole.h5ad?download=1"
    dataset_name: Slide-seqV2 - Mouse Hippocampus Puck
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
    dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
    dataset_description: "Gene expression library of mouse hippocampus puck profiled using Slide-seq V2."
    dataset_reference: stickels2020highly
    dataset_organism: Mus musculus
    spot_filter_min_genes: 200
    gene_filter_min_spots: 500

  - id: spatial_slideseq_v2/mouse_somatosensory_cortex_puck
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_SomatosensoryCortex_Puck_200306_03_data_whole.h5ad?download=1"
    dataset_name: Slide-seqV2 - Mouse Somatosensory Cortex Puck
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
    dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
    dataset_description: "Gene expression library of mouse somatosensory cortex puck profiled using Slide-seq V2."
    dataset_reference: stickels2020highly
    dataset_organism: Mus musculus
    spot_filter_min_genes: 200
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
  --main-script target/nextflow/datasets/workflows/process_spatial_from_zenodo/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
