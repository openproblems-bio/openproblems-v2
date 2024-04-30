workflow run_wf {
  take: input_ch
  main:
  output_ch = input_ch
    | novel_train.run(
      fromState: ["input_train_mod1", "input_train_mod2"],
      toState: {id, output, state -> 
        state + [
          "input_model": output.output,
          "input_transform": output.output_transform
          ] 
        }
    )
    | novel_predict.run(
      fromState: ["input_train_mod2", "input_test_mod1", "input_model", "input_transform", "output"],
      toState: ["output": "output"]
    )
    | setState ([
      "output"
    ])


  emit: output_ch
}