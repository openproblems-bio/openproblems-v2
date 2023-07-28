sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import methods
// include { mnn }                  from "$targetDir/multimodal_data_integration/methods/mnn/main.nf"                 params(params)
include { scot }                 from "$targetDir/multimodal_data_integration/methods/scot/main.nf"                params(params)
include { harmonic_alignment }   from "$targetDir/multimodal_data_integration/methods/harmonic_alignment/main.nf"  params(params)

// import metrics
include { knn_auc }              from "$targetDir/multimodal_data_integration/metrics/knn_auc/main.nf"             params(params)
include { mse }                  from "$targetDir/multimodal_data_integration/metrics/mse/main.nf"                 params(params)

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"                         params(params)

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { run_components; join_states; initialize_tracer; write_json; get_publish_dir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initialize_tracer()

// collect method list
methods = [
    scot,
    harmonic_alignment
]

// collect metric list
metrics = [
    knn_auc,
    mse
]

workflow {
    helpMessage(config)

    // create channel from input parameters with
    // arguments as defined in the config
    channelFromParams(params, config)
        | run_wf
}

    // run the workflow
workflow run_wf {
    take:
    input_ch

    main:
    output_ch = input_ch

    // based on the config file (config.vsh.yaml), run assertions on parameter sets
    // and fill in default values
    | preprocessInputs(config: config)

    // run all methods
    | run_components(
        components: methods,

        // define a new 'id' by appending the method name to the dataset id
        id: { id, state, config ->
            id + "." + config.functionality.name
        },

        // use 'from_state' to fetch the arguments the component requires from the overall state
        from_state: { id, state, config ->
            def new_args = [
            input_mod1: state.input_mod1,
            input_mod2: state.input_mod2
            ]
            new_args
        },

        // use 'to_state' to publish that component's outputs to the overall state
        to_state: { id, output, config ->
            [
            method_id: config.functionality.name,
            method_output_mod1: output.output_mod1,
            method_output_mod2: output.output_mod2
            ]
        }
    )

        // run all metrics
    | run_components(
        components: metrics,
        // use 'from_state' to fetch the arguments the component requires from the overall state
        from_state: [
            input_mod1: "method_output_mod1",
            input_mod2: "method_output_mod2"
        ],
        // use 'to_state' to publish that component's outputs to the overall state
        to_state: { id, output, config ->
            [
            metric_id: config.functionality.name,
            metric_output: output.output
            ]
        }
    )

    // join all events into a new event where the new id is simply "output" and the new state consists of:
    //   - "input": a list of score h5ads
    //   - "output": the output argument of this workflow
    | join_states{ ids, states ->
        def new_id = "output"
        def new_state = [
            input: states.collect{it.metric_output},
            output: states[0].output
        ]
        [new_id, new_state]
    }

    // convert to tsv and publish
    | extract_scores.run(
        auto: [publish: true]
    )

    emit:
    output_ch

}

// store the trace log in the publish dir
workflow.onComplete {
    def publish_dir = get_publish_dir()

    write_json(traces, file("$publish_dir/traces.json"))
    // todo: add datasets logging
    write_json(methods.collect{it.config}, file("$publish_dir/methods.json"))
    write_json(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}
