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

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; channelFromParams; preprocessInputs; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFilter as pfilter } from sourceDir + "/wf_utils/DataflowHelper.nf"
include { runMethods } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/.config.vsh.yaml")

// construct a map of methods (id -> method_module)
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
metrics = [
  accuracy,
  f1
]

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
    | preprocessInputs("config": config)

    | view{"DEBUG0: ${it}"}
    
    // multiply events by the number of method
    | runMethods(config: config, methods: methods)

    /*

    // run methods
    | getWorkflowArguments(key: "method")
    | run_methods

    // run metrics
    | getWorkflowArguments(key: "metric", inputKey: "input_prediction")
    | run_metrics

    // convert to tsv  
    | aggregate_results
    */

  emit:
  output_ch
}



workflow run_metrics {
  take: input_ch
  main:

  output_ch = input_ch
    | (accuracy & f1)
    | mix

  emit: output_ch
}

workflow aggregate_results {
  take: input_ch
  main:

  output_ch = input_ch
    | toSortedList
    | filter{ it.size() > 0 }
    | map{ it -> 
      [ "combined", it.collect{ it[1] } ] + it[0].drop(2) 
    }
    | getWorkflowArguments(key: "output")
    | extract_scores.run(
        auto: [ publish: true ]
    )

  emit: output_ch
}