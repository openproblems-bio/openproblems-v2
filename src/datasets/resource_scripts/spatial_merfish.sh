#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_merfish/human_cortex_1
    input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.250.expand.rep1_data.h5ad?download=1"
    dataset_name: MERFISH - Human Cortex 1
    dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
    dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
    dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 1 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (250 gene panel)."
    dataset_organism: Homo sapiens
    dataset_reference: fang2022conservation
    spot_filter_min_genes: 10
    gene_filter_min_spots: 100

  - id: spatial_merfish/human_cortex_2
    input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep1_data.h5ad?download=1"
    dataset_name: MERFISH - Human Cortex 2
    dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
    dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
    dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 1 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
    dataset_organism: Homo sapiens
    dataset_reference: fang2022conservation
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50

  - id: spatial_merfish/human_cortex_3
    input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep2_data.h5ad?download=1"
    dataset_name: MERFISH - Human Cortex 3
    dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
    dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
    dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 2 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
    dataset_organism: Homo sapiens
    dataset_reference: fang2022conservation
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50

  - id: spatial_merfish/human_cortex_4
    input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep3_data.h5ad?download=1"
    dataset_name: MERFISH - Human Cortex 4
    dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
    dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
    dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 3 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
    dataset_organism: Homo sapiens
    dataset_reference: fang2022conservation
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50

  - id: spatial_merfish/mouse_cortex
    input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_mouse1.AUD_TEA_VIS.242.unexpand_data.h5ad?download=1"
    dataset_name: MERFISH - Mouse Cortex
    dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
    dataset_summary: Spatially resolved profiling of mouse cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
    dataset_description: "Spatially resolved profiling of mouse cerebral cortex (visual cortex (VIS), auditory cortex (AUD) and temporal association area (TEa) unexpanded sections) using multiplexed error-robust fluorescence in situ hybridization (MERFISH)."
    dataset_organism: Mus musculus
    dataset_reference: fang2022conservation
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50

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
