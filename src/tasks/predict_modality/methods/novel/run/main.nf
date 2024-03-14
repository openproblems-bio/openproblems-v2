workflow run {
  take: input_ch
  main:
  output_ch = input_ch
    | novel_train.run(
      fromState: ["input_train_mod1", "input_train_mod2"],
      toState: ["input_model": "output"]
    )
    | novel_predict.run(
      fromState: ["input_test_mod1", "input_model"],
      toState: ["output": "output"]
    )
  emit: output_ch
}