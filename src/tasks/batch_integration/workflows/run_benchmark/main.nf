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
      id: { id, state, comp ->
        id + "_f2e"
      },
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
      id: { id, state, comp ->
        id + "_e2g"
      },
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
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },
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

    // store the dataset uns
    def dataset_uns = states.collect{it.dataset_uns}
    def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
    def dataset_uns_file = tempFile("dataset_uns.yaml")
    dataset_uns_file.write(dataset_uns_yaml_blob)

    // store the scores in a separate file
    def score_uns = states.collect{it.score_uns}
    def score_uns_yaml_blob = toYamlBlob(score_uns)
    def score_uns_file = tempFile("score_uns.yaml")
    score_uns_file.write(score_uns_yaml_blob)

    // store the method configs in a separate file
    def method_configs = methods.collect{it.config}
    def method_configs_yaml_blob = toYamlBlob(method_configs)
    def method_configs_file = tempFile("method_configs.yaml")
    method_configs_file.write(method_configs_yaml_blob)

    // store the metric configs in a separate file
    def metric_configs = metrics.collect{it.config}
    def metric_configs_yaml_blob = toYamlBlob(metric_configs)
    def metric_configs_file = tempFile("metric_configs.yaml")
    metric_configs_file.write(metric_configs_yaml_blob)

    def new_state = [
      // todo: add task info?
      // todo: add trace log?
      output_dataset_info: dataset_uns_file,
      output_scores: score_uns_file,
      output_method_configs: method_configs_file,
      output_metric_configs: metric_configs_file,
      _meta: states[0]._meta
    ]
    [new_id, new_state]
  }

  emit:
  output_ch
}
