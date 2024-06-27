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

    | download_10x_spatial.run(
      fromState: [
        "input",
        "dataset_id",
        "dataset_name",
        "dataset_url",
        "dataset_reference",
        "dataset_summary",
        "dataset_description",
        "dataset_organism"
      ],
      toState: [
        dataset_raw: "output"
      ]
    )

    | simulate_svg.run(
      fromState: [
        input: "dataset_raw"
      ],
      toState: [
        dataset_simulated: "output"
      ]
    )

    | log_cp.run(
      fromState: [
        input: "dataset_simulated",
      ],
      toState: [
        dataset_simulated_normalized: "output"
      ],
      args: [n_cp: 10000]
    )

    /*
    | check_dataset_schema.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset": checks["exit_code"] == 0 ? state.input : null,
        ]
      }
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset != null
    }
    */

    | process_dataset.run(
      fromState: [
        input: "dataset_simulated_normalized"
      ],
      toState: [
        output_dataset: "output_dataset",
        output_solution: "output_solution" 
      ]
    )

    // only output the files for which an output file was specified
    | setState(["output_dataset", "output_solution"])

  emit:
  output_ch
}
