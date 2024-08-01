#!/bin/bash

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_10x_visium/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_10x_xenium/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_slide_tags/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_slideseq_v2/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_dbit_seq/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

cat > /tmp/params.yaml << 'HERE'
id: spatially_variable_genes_process_datasets
input_states: "s3://openproblems-data/resources/datasets/spatial_merfish/**/state.yaml"
settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
rename_keys: 'input:output_dataset'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_seqfish/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_star_map/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 25}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

# cat > /tmp/params.yaml << 'HERE'
# id: spatially_variable_genes_process_datasets
# input_states: "s3://openproblems-data/resources/datasets/spatial_stereo_seq/**/state.yaml"
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200}'
# rename_keys: 'input:output_dataset'
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
  withLabel:highmem {
    memory = '350GB'
  }
  withLabel:hightime { 
    time = 15.h 
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script target/nextflow/spatially_variable_genes/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
