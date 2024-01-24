workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | niceView()

    | map { id, state ->
      def dataset_id = (new org.yaml.snakeyaml.Yaml().load(state.input)).uns.dataset_id
      [id, state + ["file_name": dataset_id]]
    }

    | yaml_to_json.run(
      fromState: { id, state -> [
        "input": state.input,
        "output": state.file_name + ".json"
      ]
      },
      toState: ["output"]
    
    )

    

    | setState({id ,state -> [
      "output": state.output
    ]

    })

    emit:
    output_ch
}