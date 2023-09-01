nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { random_embed } from "$targetDir/joint_embedding/control_methods/random_embed/main.nf"
include { zeros_embed } from "$targetDir/joint_embedding/control_methods/zeros_embed/main.nf"

// import methods
include { lmds } from "$targetDir/joint_embedding/methods/lmds/main.nf"
include { mnn } from "$targetDir/joint_embedding/methods/mnn/main.nf"
include { newwave } from "$targetDir/joint_embedding/methods/newwave/main.nf"
include { pca } from "$targetDir/joint_embedding/methods/pca/main.nf"
include { totalvi } from "$targetDir/joint_embedding/methods/totalvi/main.nf"
include { umap } from "$targetDir/joint_embedding/methods/umap/main.nf"

// import metrics
include { ari } from "$targetDir/joint_embedding/metrics/ari/main.nf"
include { asw_batch } from "$targetDir/joint_embedding/metrics/asw_batch/main.nf"
include { asw_label } from "$targetDir/joint_embedding/metrics/asw_label/main.nf"
include { cc_cons } from "$targetDir/joint_embedding/metrics/cc_cons/main.nf"
include { check_format } from "$targetDir/joint_embedding/metrics/check_format/main.nf"
include { graph_connectivity } from "$targetDir/joint_embedding/metrics/graph_connectivity/main.nf"
include { latent_mixing } from "$targetDir/joint_embedding/metrics/latent_mixing/main.nf"
include { nmi } from "$targetDir/joint_embedding/metrics/nmi/main.nf"
include { rfoob } from "$targetDir/joint_embedding/metrics/rfoob/main.nf"
include { ti_cons } from "$targetDir/joint_embedding/metrics/ti_cons/main.nf"
include { ti_cons_batch } from "$targetDir/joint_embedding/metrics/ti_cons_batch/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// construct a map of methods (id -> method_module)
methods = [ lmds, mnn, newwave, pca, totalvi, umap]
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
      method: ["input_mod1", "input_mod2"],
      metric: ["input_solution"],
      output: ["output"]
    )

    // multiply events by the number of method
    | add_methods

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
    | (ari & asw_batch & asw_label & cc_cons & check_format & graph_connectivity & latent_mixing & nmi & rfoob & ti_cons & ti_cons_batch)
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