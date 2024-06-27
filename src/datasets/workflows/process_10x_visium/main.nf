include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

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

  output_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // fetch data from legacy openproblems
    | download_10x_visium.run(
      fromState: [
        "input_expression": "input_expression",
        "input_spatial": "input_spatial",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description",
        "dataset_organism": "dataset_organism",
      ],
      toState: ["output_dataset": "dataset"]
    )
    
    // subsample if so desired
    // | subsample.run(
    //   runIf: { id, state -> state.do_subsample },
    //   fromState: [
    //     "input": "output_raw",
    //     "n_obs": "n_obs",
    //     "n_vars": "n_vars",
    //     "keep_features": "keep_features",
    //     "keep_cell_type_categories": "keep_cell_type_categories",
    //     "keep_batch_categories": "keep_batch_categories",
    //     "even": "even",
    //     "seed": "seed"
    //   ],
    //   args: [output_mod2: null],
    //   toState: ["output_raw": "output"]
    // )

    | extract_metadata.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "output_dataset")
        // workaround: convert GString to String
        schema = iterateMap(schema, { it instanceof GString ? it.toString() : it })
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.output_dataset,
          "schema": schemaYaml
        ]
      },
      toState: ["output_meta": "output"]
    )

    // only output the files for which an output file was specified
    | setState([
      "output_dataset",
      "output_meta",
      "_meta"
    ])

  emit:
  output_ch
}