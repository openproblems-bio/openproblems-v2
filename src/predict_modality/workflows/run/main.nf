nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { meanpergene } from "$targetDir/predict_modality/control_methods/meanpergene/main.nf"
include { random } from "$targetDir/predict_modality/control_methods/random/main.nf"
include { zeros } from "$targetDir/predict_modality/control_methods/zeros/main.nf"
include { solution } from "$targetDir/predict_modality/control_methods/solution/main.nf"


// import methods
// include { babel } from "$targetDir/predict_modality/methods/babel/main.nf"
include { knnr_py } from "$targetDir/predict_modality/methods/knnr_py/main.nf"
include { knnr_r } from "$targetDir/predict_modality/methods/knnr_r/main.nf"
include { lm } from "$targetDir/predict_modality/methods/lm/main.nf"
include { newwave_knnr } from "$targetDir/predict_modality/methods/newwave_knnr/main.nf"
include { random_forest } from "$targetDir/predict_modality/methods/random_forest/main.nf"


// import metrics
include { correlation } from "$targetDir/predict_modality/metrics/correlation/main.nf"
include { mse } from "$targetDir/predict_modality/metrics/mse/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// construct a map of methods (id -> method_module)
methods = [ meanpergene, random, zeros, solution, knnr_py, knnr_r, lm, newwave_knnr, random_forest]
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
      method: ["input_train_mod1", "input_train_mod2", "input_test_mod1"],
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
    | (correlation & mse)
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