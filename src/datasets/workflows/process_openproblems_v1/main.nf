workflow run_wf {
  take:
  input_ch

  main:
  normalization_settings = Channel.fromList([
    // [functionality_name, args]
    ["log_cp", [normalization_id: "log_cp10k", n_cp: 10000]],
    ["log_cp", [normalization_id: "log_cpm", n_cp: 1000000]],
    ["sqrt_cp", [normalization_id: "sqrt_cp10k", n_cp: 10000]],
    ["sqrt_cp", [normalization_id: "sqrt_cpm", n_cp: 1000000]],
    ["l1_sqrt", [normalization_id: "l1_sqrt"]],
    ["log_scran_pooling", [normalization_id: "log_scran_pooling"]]
  ])
  normalization_methods = [
    log_cp,
    sqrt_cp,
    l1_sqrt,
    log_scran_pooling
  ]

  dataset_ch = input_ch

    // fetch data from legacy openproblems
    | openproblems_v1.run(
      fromState: [
        "dataset_id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts",
        "sparse", "dataset_name", "data_url", "data_reference", "dataset_summary",
        "dataset_description", "dataset_organism"
      ],
      toState: ["raw": "output"]
    )

  sampled_dataset_ch = dataset_ch
    | filter{ id, state -> state.do_subsample }
    | subsample.run(
      fromState: [
        "input": "raw",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_celltype_categories": "keep_celltype_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      args: [
        output_mod2: null
      ],
      toState: [
        raw: "output"
      ]
    )
  notsampled_dataset_ch = dataset_ch
    | filter{ id, state -> !state.do_subsample }
  
  output_ch = sampled_dataset_ch
    | concat(notsampled_dataset_ch)

    // run normalization methods
    | combine(normalization_settings)
    | map{ id, state, norm_fun, norm_args ->
      [id, state + [ norm_fun: norm_fun, norm_args: norm_args ]]
    }

    | runComponents(
      components: normalization_methods,
      id: { id, state, config ->
        if (state.normalization_methods.size() > 1) {
          id + "/" + state.norm_args.normalization_id
        } else {
          id
        }
      },
      filter: { id, state, config ->
        config.functionality.name == state.norm_fun &&
        state.norm_args.normalization_id in state.normalization_methods
      },
      fromState: { id, state, config ->
        [input: state.raw] + state.norm_args
      },
      toState: ["normalized": "output"]
    )

    | pca.run(
      fromState: ["input": "normalized"],
      toState: ["pca": "output" ]
    )

    | hvg.run(
      fromState: ["input": "pca"],
      toState: ["hvg": "output"]
    )

    | knn.run(
      fromState: ["input": "hvg"],
      toState: ["knn": "output"]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.knn,
          checks: null
        ]
      },
      toState: ["dataset": "output", "meta": "meta"]
    )

    // only output the files for which an output file was specified
    | setState{ id, state ->
      [
        "output_dataset": state.output_dataset ? state.dataset : null,
        "output_meta": state.output_meta ? state.meta : null,
        "output_raw": state.output_raw ? state.raw : null,
        "output_normalized": state.output_normalized ? state.normalized : null,
        "output_pca": state.output_pca ? state.pca : null,
        "output_hvg": state.output_hvg ? state.hvg : null,
        "output_knn": state.output_knn ? state.knn : null
      ]
    }

  emit:
  output_ch
}