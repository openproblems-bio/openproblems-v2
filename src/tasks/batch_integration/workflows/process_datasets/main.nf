nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { copy } from "$targetDir/common/copy/main.nf"
include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"
include { process_dataset } from "$targetDir/batch_integration/process_dataset/main.nf"

// import helper functions
include { readConfig; processConfig; helpMessage; channelFromParams; preprocessInputs; readYaml; readJson } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { publishState; runComponents; joinStates; initializeTracer; writeJson; getPublishDir; autoDetectStates } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
    | publishState([:])
}

workflow auto {
  autoDetectStates(params, config)
    | run_wf
    | publishState([:])
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs(config: config)
    
    // add copy to avoid potential name clashes between
    // inputs and outputs of process_dataset
    | copy.run(
      fromState: ["input"],
      toState: ["input": "output"]
    )

    // TODO: check whether datasets match the 
    // input schema of the process_dataset component
    // right now we just let it fail.

    | process_dataset.run(
      fromState: ["input", "output_dataset", "output_solution"],
      toState: [dataset: "output_dataset", solution: "output_solution"]
    )

    // only output the files for which an output file was specified
    | map { id, state ->
      def keys = ["dataset", "solution"]
      def newState = keys.collectMany{ key ->
        def output_key = "output_" + key
        if (state.containsKey(output_key)) {
          [ [ output_key, state[key] ] ]
        } else {
          []
        }
      }.collectEntries()
      [id, newState]
    }

  emit:
  output_ch
}
