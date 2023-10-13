workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  // collect method list
  methods = [
    bbknn,
    combat,
    fastmnn_embedding,
    fastmnn_feature,
    liger,
    mnn_correct,
    mnnpy,
    pyliger,
    scalex_embed,
    scalex_feature,
    scanorama_embed,
    scanorama_feature,
    scanvi,
    scvi,
    no_integration_batch,
    random_embed_cell,
    random_embed_cell_jitter,
    random_integration
  ]

  // collect metric list
  metrics = [
    asw_batch,
    asw_label,
    cell_cycle_conservation,
    clustering_overlap,
    graph_connectivity,
    hvg_overlap,
    isolated_label_asw,
    isolated_label_f1,
    kbet,
    lisi,
    pcr
  ]

  // process input parameter channel
  dataset_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // extract the dataset metadata
    | check_dataset_schema.run(
      key: "extract_dataset_uns",
      fromState: [input: "input_dataset"],
      toState: { id, output, state ->
        def dataset_uns = (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
        state + [dataset_uns: dataset_uns]
      }
    )

    // run all methods
  method_out_ch1 = dataset_ch
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        def norm = state.dataset_uns.normalization_id
        def pref = comp.config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        (norm == "log_cp10k" && pref == "counts") || norm == pref
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [input: "input_dataset"],

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.functionality.name,
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
  

  // append feature->embed transformations
  method_out_ch2 = method_out_ch1
    | runEach(
      components: feature_to_embed,
      filter: { id, state, comp -> state.method_subtype == "feature"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, comp ->
        state + [
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch1)

  // append embed->graph transformations
  method_out_ch3 = method_out_ch2
    | runEach(
      components: embed_to_graph,
      filter: { id, state, comp -> state.method_subtype == "embedding"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, comp ->
        state + [
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch2)

  // run metrics
  output_ch = method_out_ch3
    | runEach(
      components: metrics,
      filter: { id, state, comp ->
        state.method_subtype == comp.config.functionality.info.subtype
      },
      fromState: [
        input_integrated: "method_output",
        input_solution: "input_solution"
      ],
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.functionality.name,
          metric_output: output.output
        ]
      }
    )

  // TODO: can we store everything below in a separate helper function?
  
  | check_dataset_schema.run(
    key: "extract_scores",
    fromState: [input: "metric_output"],
    toState: { id, output, state ->
      def score_uns = (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
      state + [score_uns: score_uns]
    }
  )
  
  | joinStates{ ids, states ->
    def new_id = "output"
    def new_state = [
      score_uns: states.collect{it.score_uns},
      _meta: states[0]._meta
    ]
    [new_id, new_state]
  }

  | save_file.run(
    key: "save_score_uns",
    fromState: { id, state ->
      [
        input: toYamlBlob(state.score_uns),
        output: '$id.$key.score.yaml'
      ]
    },
    toState: ["output_scores": "output"]
  )
  | save_file.run(
    key: "save_method_configs",
    fromState: { id, state ->
      // TODO: filter methods
      [
        input: toYamlBlob(methods.collect{it.config}),
        output: '$id.$key.method_configs.yaml'
      ]
    },
    toState: ["output_method_configs": "output"]
  )
  | save_file.run(
    key: "save_metric_configs",
    fromState: { id, state ->
      // TODO: filter metrics
      [
        input: toYamlBlob(metrics.collect{it.config}),
        output: '$id.$key.metric_configs.yaml'
      ]
    },
    toState: ["output_metric_configs": "output"]
  )
  // TODO: can we also store the trace log?

  | setState([
    "output_scores",
    "output_method_configs",
    "output_metric_configs",
    "_meta"
  ])

  emit:
  output_ch
}
