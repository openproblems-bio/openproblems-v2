def runComponents(Map args) {
  assert args.components: "runComponents should be passed a list of components to run"
  assert args.from_state: "runComponents should be passed a from_state function"
  assert args.to_state: "runComponents should be passed a to_state function"

  def components_ = args.components
  if (components_ !instanceof List) {
    components_ = [ components_ ]
  }
  assert components.size() > 0: "pass at least one component to runComponents"

  def from_state_ = args.from_state
  def to_state_ = args.to_state
  def filter_ = args.filter

  workflow runComponentsWf {
    take: input_ch
    main:

    // generate one channel per method
    out_chs = components_.collect{ comp_ ->
      def comp_config = comp_.config

      if (filter_) {
        mid_ch = input_ch
          | filter{ tup -> 
            filter_(tup[0], tup[1], comp_config)
          }
      } else {
        mid_ch = input_ch
      }
      mid_ch
        | map{ tup ->
          def new_tup = from_state_(tup[0], tup[1], comp_config)
          new_tup + tup.drop(1)
        }
        | comp_.run(
          auto: (args.auto ?: [:]) + [simplifyInput: false, simplifyOutput: false]
        )
        | map{tup ->
          def new_outputs = to_state_(tup[0], tup[1], comp_config)
          [tup[0], tup[2] + new_outputs] + tup.drop(3)
        }
      }

    // mix all results
    output_ch =
      (out_chs.size == 1)
        ? out_chs[0]
        : out_chs[0].mix(*out_chs.drop(1))

    emit: output_ch
  }

  return runComponentsWf
}

def joinStates(Map args) {
  assert args.apply: "joinStates should be passed a function in the apply argument"

  def apply_ = args.apply
  workflow joinStatesWf {
    take: input_ch
    main:
    output_ch = input_ch
      | toSortedList
      | filter{ it.size() > 0 }
      | map{ tups ->
        def ids = tups.collect{it[0]}
        def states = tups.collect{it[1]}
        apply_(ids, states)
      }

    emit: output_ch
  }
  return joinStatesWf
}
