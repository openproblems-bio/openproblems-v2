workflow auto {
  findStates(params, thisConfig)
    | run_wf
    | publishStates([:])
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // TODO: check schema based on the values in `config`
    // instead of having to provide a separate schema file
    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.input,
          schema: state.dataset_schema,
          output: '$id.$key.output.h5ad',
          stop_on_error: false,
          checks: null
        ]
      },
      toState: { id, output, state ->
        state + [ dataset: output.output ]
      }
    )

    | filter { id, state ->
      state.dataset != null
    }

    | process_dataset.run(
      fromState: [
        input: "dataset",
        output_dataset: "output_train",
        output_solution: "output_test"
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
