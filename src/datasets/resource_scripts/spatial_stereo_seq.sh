#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_stereo_seq/drosophila_embryo_e5_6
    input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_5.6.h5ad?download=1"
    dataset_name: Stereo-seq - Drosophila embryo E5_6
    dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
    dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
    dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
    dataset_organism: Drosophila
    dataset_reference: wang2022high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 50
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_stereo_seq/drosophila_embryo_e6_3
    input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_6.3.h5ad?download=1"
    dataset_name: Stereo-seq - Drosophila embryo E6_3
    dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
    dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
    dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
    dataset_organism: Drosophila
    dataset_reference: wang2022high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 50
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_stereo_seq/drosophila_embryo_e7
    input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_7.h5ad?download=1"
    dataset_name: Stereo-seq - Drosophila embryo E7
    dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
    dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
    dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
    dataset_organism: Drosophila
    dataset_reference: wang2022high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 50
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_stereo_seq/drosophila_embryo_e9_1
    input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_9.1.h5ad?download=1"
    dataset_name: Stereo-seq - Drosophila embryo E9_1
    dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
    dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
    dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
    dataset_organism: Drosophila
    dataset_reference: wang2022high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 50
    select_top_variable_genes: 50
    remove_mitochondrial: true
    coord_type: generic

  - id: spatial_stereo_seq/drosophila_embryo_e10
    input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_10.5.h5ad?download=1"
    dataset_name: Stereo-seq - Drosophila embryo E10
    dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
    dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
    dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
    dataset_organism: Drosophila
    dataset_reference: wang2022high
    spot_filter_min_genes: 10
    gene_filter_min_spots: 50
    num_reference_genes: 50
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
