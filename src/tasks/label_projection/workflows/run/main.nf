nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { true_labels } from "$targetDir/label_projection/control_methods/true_labels/main.nf"
include { majority_vote } from "$targetDir/label_projection/control_methods/majority_vote/main.nf"
include { random_labels } from "$targetDir/label_projection/control_methods/random_labels/main.nf"

// import methods
include { knn } from "$targetDir/label_projection/methods/knn/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
include { mlp } from "$targetDir/label_projection/methods/mlp/main.nf"
include { scanvi } from "$targetDir/label_projection/methods/scanvi/main.nf"
include { scanvi_scarches } from "$targetDir/label_projection/methods/scanvi_scarches/main.nf"
include { seurat_transferdata } from "$targetDir/label_projection/methods/seurat_transferdata/main.nf"
include { xgboost } from "$targetDir/label_projection/methods/xgboost/main.nf"

// import metrics
include { accuracy } from "$targetDir/label_projection/metrics/accuracy/main.nf"
include { f1 } from "$targetDir/label_projection/metrics/f1/main.nf"

// convert scores to tsv
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { runComponents; joinStates } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// collect method list
methods = [
  true_labels,
  majority_vote,
  random_labels,
  knn,
  logistic_regression,
  mlp,
  scanvi,
  scanvi_scarches,
  seurat_transferdata,
  xgboost
]

// collect metric list
metrics = [
  accuracy,
  f1
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

    // run methods
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
          input_train: state.input_train,
          input_test: state.input_test
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_solution = state.input_solution
        }
        [new_id, new_args]
      },
      to_state: { id, output, config ->
        [
          method_id: config.functionality.name,
          prediction: output.output
        ]
      }
    )

    // run metrics
    | runComponents(
      components: metrics,
      from_state: { id, state, config ->
        def new_args = [
          input_solution: state.input_solution,
          input_prediction: state.prediction
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

    | extract_scores.run(
      auto: [publish: true]
    )

  emit:
  output_ch
}