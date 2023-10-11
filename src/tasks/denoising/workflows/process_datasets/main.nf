workflow auto {
  // TODO: `thisConfig` might be renamed to `meta["config"]` in the future
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

    // TODO: check schema based on the values in `config`
    // instead of having to provide a separate schema file
    | check_dataset_schema.run(
      fromState: [
        "input": "input",
        "schema": "dataset_schema"
      ],
      args: [
        "stop_on_error": false
      ],
      toState: [
        "dataset": "output",
        "dataset_checks": "checks"
      ]
    )

    | filter { id, state ->
      state.dataset != null
    }

    | process_dataset.run(
      fromState: [
        input: "dataset",
        output_train: "output_train",
        output_test: "output_test"
      ],
      toState: [train: "output_train", test: "output_test"]
    )

    // only output the files for which an output file was specified
    | setState { id, state ->
      [
        "output_train": state.output_train ? state.train : null,
        "output_test": state.output_test ? state.test : null
      ]
    }

  emit:
  output_ch
}
