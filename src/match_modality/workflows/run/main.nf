nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { constant } from "$targetDir/match_modality/control_methods/constant/main.nf"
include { random_pairing } from "$targetDir/match_modality/control_methods/random_pairing/main.nf"
include { semi_solution } from "$targetDir/match_modality/control_methods/semi_solution/main.nf"
include { solution } from "$targetDir/match_modality/control_methods/solution/main.nf"


// import methods
// include { babel_knn } from "$targetDir/match_modality/methods/babel_knn/main.nf"
include { dr_knnr_cbf } from "$targetDir/match_modality/methods/dr_knnr_cbf/main.nf"
include { dr_knnr_knn } from "$targetDir/match_modality/methods/dr_knnr_knn/main.nf"
include { linear_knn } from "$targetDir/match_modality/methods/linear_knn/main.nf"
include { newwave_knnr_cbf } from "$targetDir/match_modality/methods/newwave_knnr_cbf/main.nf"
include { newwave_knnr_knn } from "$targetDir/match_modality/methods/newwave_knnr_knn/main.nf"
include { procrustes_knn } from "$targetDir/match_modality/methods/procrustes_knn/main.nf"


// import metrics
include { aupr } from "$targetDir/match_modality/metrics/aupr/main.nf"
include { check_format } from "$targetDir/match_modality/metrics/check_format/main.nf"
include { match_probability } from "$targetDir/match_modality/metrics/match_probability/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// construct a map of methods (id -> method_module)
methods = [ dr_knnr_cbf, dr_knnr_knn, linear_knn, newwave_knnr_cbf, newwave_knnr_knn, procrustes_knn]
  .collectEntries{method ->
    [method.config.functionality.name, method]
  }

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // split params for downstream components
    | setWorkflowArguments(
      method: ["input_train_mod1", "input_train_mod2", "input_train_sol", "input_test_mod1", "input_test_mod2"],
      metric: ["input_solution"],
      output: ["output"]
    )

    // multiply events by the number of method
    | add_methods

    // add input_solution to data for the positive controls
    | controls_can_cheat

    // run methods
    | getWorkflowArguments(key: "method")
    | run_methods

    // construct tuples for metrics
    | pmap{ id, file, passthrough ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      // combine prediction with solution
      def newData = [ input_prediction: file, input_solution: passthrough.metric.input_solution ]
      [ newId, newData, passthrough ]
    }
    
    // run metrics
    | getWorkflowArguments(key: "metric")
    | run_metrics
    
    // convert to tsv  
    | aggregate_results

  emit:
  output_ch
}

workflow add_methods {
  take: input_ch
  main:
  output_ch = Channel.fromList(methods.keySet())
    | combine(input_ch)

    // generate combined id for method_id and dataset_id
    | pmap{method_id, dataset_id, data ->
      def new_id = dataset_id + "." + method_id
      def new_data = data.clone() + [method_id: method_id]
      new_data.remove("id")
      [new_id, new_data]
    }
  emit: output_ch
}

workflow controls_can_cheat {
  take: input_ch
  main:
  output_ch = input_ch
    | pmap{id, data, passthrough ->
      def method = methods[data.method_id]
      def method_type = method.config.functionality.info.method_type
      def new_data = data.clone()
      if (method_type != "method") {
        new_data = new_data + [input_test_sol: passthrough.metric.input_solution]
      }
      [id, new_data, passthrough]
    }
  emit: output_ch
}

workflow run_methods {
  take: input_ch
  main:
    // generate one channel per method
    method_chs = methods.collect { method_id, method_module ->
        input_ch
          | filter{it[1].method_id == method_id}
          | method_module
      }
    // mix all results
    output_ch = method_chs[0].mix(*method_chs.drop(1))

  emit: output_ch
}

workflow run_metrics {
  take: input_ch
  main:

  output_ch = input_ch
    | (aupr & check_format & match_probability)
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