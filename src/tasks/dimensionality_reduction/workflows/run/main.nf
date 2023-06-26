sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { random_features } from "$targetDir/dimensionality_reduction/control_methods/random_features/main.nf"
include { true_features } from "$targetDir/dimensionality_reduction/control_methods/true_features/main.nf"

// import methods
include { densmap } from "$targetDir/dimensionality_reduction/methods/densmap/main.nf"
// include { ivis } from "$targetDir/dimensionality_reduction/methods/ivis/main.nf"
include { neuralee } from "$targetDir/dimensionality_reduction/methods/neuralee/main.nf"
include { pca } from "$targetDir/dimensionality_reduction/methods/pca/main.nf"
include { phate } from "$targetDir/dimensionality_reduction/methods/phate/main.nf"
include { tsne } from "$targetDir/dimensionality_reduction/methods/tsne/main.nf"
include { umap } from "$targetDir/dimensionality_reduction/methods/umap/main.nf"

// import metrics
include { coranking } from "$targetDir/dimensionality_reduction/metrics/coranking/main.nf"
include { density_preservation } from "$targetDir/dimensionality_reduction/metrics/density_preservation/main.nf"
include { rmse } from "$targetDir/dimensionality_reduction/metrics/rmse/main.nf"
include { trustworthiness } from "$targetDir/dimensionality_reduction/metrics/trustworthiness/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { runComponents; aggregate_scores } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// collect method list
methods = [
  random_features,
  true_features,
  densmap,
  neuralee,
  pca,
  phate,
  tsne,
  umap
]

// collect metric list
metrics = [
  coranking,
  density_preservation,
  rmse,
  trustworthiness
]

workflow {
  helpMessage(config)

  // create channel from input parameters with
  // arguments as defined in the config
  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs(config: config)

    // run all methods
    //   - use the 'filter' argument to only run a method on the normalisation the component is asking for
    //   - use 'from_state' to fetch the arguments the component requires from the overall state
    //   - use 'to_state' to publish that component's outputs to the overall state
    | runComponents(
      components: methods,
      filter: { id, state, config ->
        def norm = state.normalization_id
        def pref = config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        (norm == "log_cpm" && pref == "counts") || norm == pref
      },
      from_state: { id, state, config ->
        def new_id = id + "." + config.functionality.name
        def new_args = [
          input: state.input
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_solution = state.input_solution
        }
        [new_id, new_args]
      },
      to_state: { id, output, config ->
        [
          method_id: config.functionality.name,
          embedding: output.output
        ]
      }
    )

    // run all metrics
    //   - use 'from_state' to fetch the arguments the component requires from the overall state
    //   - use 'to_state' to publish that component's outputs to the overall state
    | runComponents(
      components: metrics,
      from_state: { id, state, config ->
        def new_args = [
          input_solution: state.input_solution,
          input_embedding: state.embedding
        ]
        [id, new_args]
      },
      to_state: { id, output, config ->
        [
          metric_id: config.functionality.name,
          scores: output.output
        ]
      }
    )

    // join all events into a new event where the new id is simply "output" and the new state consists of:
    //   - "input": a list of score h5ads
    //   - "output": the output argument of this workflow
    | joinStates(
      apply: { ids, states ->
        def new_id = "output"
        def new_state = [
          "input": states.collect{it.scores},
          "output": states[0].output
        ]
        [new_id, new_state]
      }
    )

    // convert to tsv and publish
    | extract_scores.run(
      auto: [publish: true]
    )

  emit:
  output_ch
}