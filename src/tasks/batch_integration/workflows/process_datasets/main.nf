nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"
include { process_dataset } from "$targetDir/batch_integration/process_dataset/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs; readYaml } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { publishState; runComponents; joinStates; initializeTracer; writeJson; getPublishDir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")


workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
    | publishState([:])
}

workflow run_wf {
  take:
  input_ch

  main:

  // process input parameter channel
  output_ch = input_ch
    | preprocessInputs(config: config)

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.input,
          schema: state.input_schema,
          stop_on_error: true,
          output: null,
          meta: null
        ]
      },
      toState: { id, output, state ->
        state
      }
    )

    | process_dataset.run(
      fromState: ["input", "output_dataset", "output_solution"],
      toState: ["output_dataset", "output_solution"]
    )

    | view {"Output: $it"}

  emit:
  output_ch
}
