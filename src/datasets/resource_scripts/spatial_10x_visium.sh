#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_10x_visium/mouse_brain_coronal_section1
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_spatial.tar.gz"
    dataset_name: Mouse Brain Coronal Section 1 (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard"
    dataset_summary: Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set
    dataset_description: "FFPE Mouse Brain tissue blocks sectioned as described in Visium CytAssist Spatial Gene Expression for FFPE - Tissue Preparation Guide Demonstrated Protocol. The H&E stained glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression slide. The probe extension and library construction steps follow the standard Visium for FFPE workflow outside of the instrument. The H&E image was acquired using Olympus VS200 Slide Scanning Microscope. Sequencing depth was 53,497 reads per spot. Sequencing configuration: 28bp read 1 (16bp Visium spatial barcode, 12bp UMI), 90bp read 2 (transcript), 10bp i7 sample barcode and 10bp i5 sample barcode. Key metrics include: 2,310 spots detected under tissue; 6,736 median genes per spot; 24,862 median UMI counts per spot."
    dataset_organism: Mus musculus

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_hvg: force_null
publish_dir: s3://openproblems-data/resources/datasets
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
  --main-script target/nextflow/datasets/workflows/process_10x_visium/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
