nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// dataset loaders
include { openproblems_v1 } from "$targetDir/datasets/loaders/openproblems_v1/main.nf"

// normalization methods
include { log_cpm } from "$targetDir/datasets/normalization/log_cp/main.nf"
include { log_scran_pooling } from "$targetDir/datasets/normalization/log_scran_pooling/main.nf"
include { sqrt_cpm } from "$targetDir/datasets/normalization/sqrt_cp/main.nf"
include { l1_sqrt } from "$targetDir/datasets/normalization/l1_sqrt/main.nf"

// dataset processors
include { pca } from "$targetDir/datasets/processors/pca/main.nf"
include { hvg } from "$targetDir/datasets/processors/hvg/main.nf"
include { knn } from "$targetDir/datasets/processors/knn/main.nf"
include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { runComponents; joinStates; initializeTracer; writeJson getPublishDir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

// normalization_methods = [log_cp, log_scran_pooling, sqrt_cp, l1_sqrt
normalization_methods = [log_cp, sqrt_cp, l1_sqrt]

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs(config: config)

    // fetch data from legacy openproblems
    | runComponents(
      components: openproblems_v1,
      fromState: [
        "dataset_id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts", "sparse",
        "dataset_name", "data_url", "data_reference", "dataset_summary", "dataset_description", "dataset_organism"
      ],
      toState: [ dataset: "output" ]
    )

    // run normalization methods
    | runComponents(
      components: normalization_methods,
      id: { id, state, config -> id + "/" + config.functionality.name },
      fromState: [ input: "dataset" ],
      toState: [
        normalization_id: config.functionality.name,
        output_normalization: "output"
      ]
    )

    | runComponents(
      components: pca,
      fromState: [ input: "output_normalization" ],
      toState: [ pca: "output" ]
    )

    | runComponents(
      components: hvg,
      fromState: [ input: "pca" ],
      toState: [ hvg: "output" ]
    )

    | runComponents(
      components: knn,
      fromState: [ input: "hvg" ],
      toState: [ knn: "output" ]
    )

    | runComponents(
      components: check_dataset_schema,
      fromState: {id, state, config ->
        [
          input: state.knn,
          meta: state.output_meta,
          output: state.output_dataset,
          checks: null
        ]
      },
      toState: [],
      auto: [publish: true]
    )

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJsontraces, file("$publish_dir/traces.json"))
  writeJsonnormalization_methods.collect{it.config}, file("$publish_dir/normalization_methods.json"))
}