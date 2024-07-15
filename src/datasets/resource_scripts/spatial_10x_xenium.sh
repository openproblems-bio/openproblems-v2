#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_10x_xenium/human_colon_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_spatial.tar.gz"
    dataset_name: CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1
    dataset_url: "https://www.10xgenomics.com/datasets/visium-cytassist-gene-expression-libraries-of-post-xenium-human-colon-cancer-ffpe-using-the-human-whole-transcriptome-probe-set-2-standard"
    dataset_summary: Gene expression library of Post Xenium Human Colon Cancer (CytAssist FFPE) using the Human Whole Transcriptome Probe Set - Replicate 1
    dataset_description: "This dataset is provided as part of the Technical Note: Post-Xenium In Situ Applications: Immunofluorescence, H&E, and Visium CytAssist Spatial Gene Expression (CG000709). Post-Xenium samples were compared to controls (samples not processed through the Xenium workflow) using 5 µm (FFPE) serial sections."
    dataset_organism: Homo sapiens

  - id: spatial_10x_xenium/mouse_brain
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1_spatial.tar.gz"
    dataset_name: CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1
    dataset_url: "https://www.10xgenomics.com/datasets/visium-cytassist-gene-expression-libraries-of-post-xenium-mouse-brain-ff-using-the-mouse-whole-transcriptome-probe-set-2-standard"
    dataset_summary: Gene expression library of Post Xenium Mouse Brain (CytAssist Fresh Frozen) using the Mouse Whole Transcriptome Probe Set - Replicate 1
    dataset_description: "This dataset is provided as part of the Technical Note: Post-Xenium In Situ Applications: Immunofluorescence, H&E, and Visium CytAssist Spatial Gene Expression (CG000709). Post-Xenium samples were compared to controls (samples not processed through the Xenium workflow) using 10 µm fresh-frozen (FF) serial sections."
    dataset_organism: Mus musculus

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: s3://openproblems-data/resources/datasets
spot_filter_min_genes: 100
gene_filter_min_spots: 50
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
