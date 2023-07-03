nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import preprocessing
include { process_dataset } from "$targetDir/batch_integration/process_dataset/main.nf"

// import methods
include { bbknn } from "$targetDir/batch_integration/methods/bbknn/main.nf"
include { combat } from "$targetDir/batch_integration/methods/combat/main.nf"
include { scanorama_embed } from "$targetDir/batch_integration/methods/scanorama_embed/main.nf"
include { scanorama_feature } from "$targetDir/batch_integration/methods/scanorama_feature/main.nf"
include { scvi } from "$targetDir/batch_integration/methods/scvi/main.nf"

// import transformers
include { feature_to_embed } from "$targetDir/batch_integration/transformers/feature_to_embed/main.nf"
include { embed_to_graph } from "$targetDir/batch_integration/transformers/embed_to_graph/main.nf"

// import metrics
include { clustering_overlap } from "$targetDir/batch_integration/metrics/clustering_overlap/main.nf"
include { asw_batch } from "$targetDir/batch_integration/metrics/asw_batch/main.nf"
include { asw_label } from "$targetDir/batch_integration/metrics/asw_label/main.nf"
include { cell_cycle_conservation } from "$targetDir/batch_integration/metrics/cell_cycle_conservation/main.nf"
include { pcr } from "$targetDir/batch_integration/metrics/pcr/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { run_components; join_states; initialize_tracer; write_json; get_publish_dir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initialize_tracer()

// collect method list
methods = [
  scanorama_feature
]

// collect metric list
metrics = [
  clustering_overlap
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
    | preprocessInputs(config: config)

    // preprocess data
    | run_components(
      components: process_dataset,
      from_state: ["input"],
      to_state: { id, output, config ->
        [
          process_data_id: config.functionality.name,
          process_data_output: output.output
        ]
      }
    )

    // run all methods
    | run_components(
      components: methods,

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, config ->
        id + "." + config.functionality.name
      },

      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: ["input"],

      // use 'to_state' to publish that component's outputs to the overall state
      to_state: { id, output, config ->
        [
          method_id: config.functionality.name,
          method_output: output.output
        ]
      }
    )

    // // transform data
    // | run_components(
    //   //  TODO: filter to feature output methods only
    //   components: feature_to_embed,
    //   from_state: ["input"],
    //   to_state: { id, output, config ->
    //     [
    //       transformer_id: config.functionality.name,
    //       transformer_output: output.output
    //     ]
    //   },
    // )

    // | run_components(
    //   //  TODO: filter to feature and embedding output methods only
    //   components: embed_to_graph,
    //   from_state: ["input"],
    //   to_state: { id, output, config ->
    //     [
    //       transformer_id: config.functionality.name,
    //       transformer_output: output.output
    //     ]
    //   },
    // )

    // run all metrics
    | run_components(
      components: metrics,
      from_state: ["input_integrated"],
      to_state: { id, output, config ->
        [
          metric_id: config.functionality.name,
          metric_output: output.output
        ]
      }
    )

//     // join all events into a new event where the new id is simply "output" and the new state consists of:
//     //   - "input": a list of score h5ads
//     //   - "output": the output argument of this workflow
//     | join_states{ ids, states ->
//       def new_id = "output"
//       def new_state = [
//         input: states.collect{it.metric_output},
//         output: states[0].output
//       ]
//       [new_id, new_state]
//     }

//     // convert to tsv and publish
//     | extract_scores.run(
//       auto: [publish: true]
//     )

  emit:
  output_ch
}

// // store the trace log in the publish dir
// workflow.onComplete {
//   def publish_dir = get_publish_dir()

//   write_json(traces, file("$publish_dir/traces.json"))
//   // todo: add datasets logging
//   write_json(methods.collect{it.config}, file("$publish_dir/graph_methods.json"))
//   write_json(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
// }
