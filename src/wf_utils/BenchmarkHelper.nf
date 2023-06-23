sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { extract_scores } from "$targetDir/common/extract_scores/main.nf"
include { preprocessInputs; processConfig } from sourceDir + "/wf_utils/WorkflowHelper.nf"

def runComponents(Map args) {
  assert args.components: "runComponents should be passed a list of components to run"
  assert args.fetch_data: "runComponents should be passed a fetch_data function"
  assert args.store_data: "runComponents should be passed a store_data function"

  def components_ = args.components
  def fetch_data_ = args.fetch_data
  def store_data_ = args.store_data
  def filter_ = args.filter ?: { id, data, comp -> true }

  workflow runComponentsWf {
    take: input_ch
    main:

    // generate one channel per method
    out_chs = components_.collect { comp_ ->
      def comp_config = comp_.config

      input_ch
        | filter{tup -> 
          filter_(tup[0], tup[1], comp_config)
        }
        | map{tup ->
          def new_tup = fetch_data_(tup[0], tup[1], comp_config)
          new_tup + tup.drop(1)
        }
        | comp_.run(
          auto: [simplifyInput: false, simplifyOutput: false]
        )
        | map{tup ->
          def new_outputs = store_data_(tup[0], tup[1], comp_config)
          [tup[0], tup[2] + new_outputs] + tup.drop(3)
        }
      }

    // mix all results
    output_ch = out_chs[0].mix(*out_chs.drop(1))

    emit: output_ch
  }

  return runComponentsWf
}


workflow aggregate_scores {
  take: input_ch
  main:

  output_ch = input_ch
    | toSortedList
    | filter{ it.size() > 0 }
    | map{ tups -> 
      def new_id = "combined"
      def new_data = [
        "input": tups.collect{ it[1].scores },
        "output": tups[0][1].output
      ]
      [new_id, new_data]
    }
    | extract_scores.run(
        auto: [ publish: true ]
    )

  emit: output_ch
}
