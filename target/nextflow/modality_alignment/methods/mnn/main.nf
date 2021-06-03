nextflow.enable.dsl=2

params.test = false
params.debug = false

def checkParams(_params) {
  _params.arguments.collect{
    if (it.value == "viash_no_value") {
      println("[ERROR] option --${it.name} not specified in component mnn")
      println("exiting now...")
        exit 1
    }
  }
}


def renderCLI(command, arguments) {

  def argumentsList = arguments.collect{ it ->
    (it.otype == "")
      ? "\'" + it.value + "\'"
      : (it.type == "boolean_true")
        ? it.otype + it.name
        : (it.value == "no_default_value_configured")
          ? ""
          : it.otype + it.name + " \'" + ((it.value in List && it.multiple) ? it.value.join(it.multiple_sep): it.value) + "\'"
  }

  def command_line = command + argumentsList

  return command_line.join(" ")
}

def effectiveContainer(processParams) {
  def _registry = params.containsKey("containerRegistry") ? params.containerRegistry : processParams.containerRegistry
  def _name = processParams.container
  def _tag = params.containsKey("containerTag") ? params.containerTag : processParams.containerTag

  return (_registry == "" ? "" : _registry + "/") + _name + ":" + _tag
}

// Convert the nextflow.config arguments list to a List instead of a LinkedHashMap
// The rest of this main.nf script uses the Map form
def argumentsAsList(_params) {
  def overrideArgs = _params.arguments.collect{ key, value -> value }
  def newParams = _params + [ "arguments" : overrideArgs ]
  return newParams
}


// Use the params map, create a hashmap of the filenames for output
// output filename is <sample>.<method>.<arg_name>[.extension]
def outFromIn(_params) {

  def id = _params.id

  _params
    .arguments
    .findAll{ it -> it.type == "file" && it.direction == "Output" }
    .collect{ it ->
      // If a default (dflt) attribute is present, strip the extension from the filename,
      // otherwise just use the option name as an extension.
      def extOrName = (it.dflt != null) ? it.dflt.split(/\./).last() : it.name
      // The output filename is <sample> . <modulename> . <extension>
      def newName =
        (id != "")
          ? id + "." + "mnn" + "." + extOrName
          : "mnn" + "." + extOrName
      it + [ value : newName ]
    }

}

// A process that filters out output from the output Map
process filterOutput {

  input:
    tuple val(id), val(input), val(_params)
  output:
    tuple val(id), val(output), val(_params)
  when:
    input.keySet().contains("output")
  exec:
    output = input["output"]

}

def overrideIO(_params, inputs, outputs) {

  // `inputs` in fact can be one of:
  // - `String`,
  // - `List[String]`,
  // - `Map[String, String | List[String]]`
  // Please refer to the docs for more info
  def overrideArgs = _params.arguments.collect{ it ->
    if (it.type == "file") {
      if (it.direction == "Input") {
        (inputs in List || inputs in HashMap)
          ? (inputs in List)
            ? it + [ "value" : inputs.join(it.multiple_sep)]
            : (inputs[it.name] != null)
              ? (inputs[it.name] in List)
                ? it + [ "value" : inputs[it.name].join(it.multiple_sep)]
                : it + [ "value" : inputs[it.name]]
              : it
          : it + [ "value" : inputs ]
      } else {
        (outputs in List || outputs in HashMap)
          ? (outputs in List)
            ? it + [ "value" : outputs.join(it.multiple_sep)]
            : (outputs[it.name] != null)
              ? (outputs[it.name] in List)
                ? it + [ "value" : outputs[it.name].join(it.multiple_sep)]
                : it + [ "value" : outputs[it.name]]
              : it
          : it + [ "value" : outputs ]
      }
    } else {
      it
    }
  }

  def newParams = _params + [ "arguments" : overrideArgs ]

  return newParams

}

process mnn_process {


  tag "${id}"
  echo { (params.debug == true) ? true : false }
  cache 'deep'
  stageInMode "symlink"
  container "${container}"

  input:
    tuple val(id), path(input), val(output), val(container), val(cli), val(_params)
  output:
    tuple val("${id}"), path(output), val(_params)
  script:
    if (params.test)
      """
      # Some useful stuff
      export NUMBA_CACHE_DIR=/tmp/numba-cache
      # Running the pre-hook when necessary
      # Adding NXF's `$moduleDir` to the path in order to resolve our own wrappers
      export PATH="./:${moduleDir}:\$PATH"
      ./${params.mnn.tests.testScript} | tee $output
      """
    else
      """
      # Some useful stuff
      export NUMBA_CACHE_DIR=/tmp/numba-cache
      # Running the pre-hook when necessary
      # Adding NXF's `$moduleDir` to the path in order to resolve our own wrappers
      export PATH="${moduleDir}:\$PATH"
      $cli
      """
}

workflow mnn {

  take:
  id_input_params_

  main:

  def key = "mnn"

  def id_input_output_function_cli_params_ =
    id_input_params_.map{ id, input, _params ->

      // Start from the (global) params and overwrite with the (local) _params
      def defaultParams = params[key] ? params[key] : [:]
      def overrideParams = _params[key] ? _params[key] : [:]
      def updtParams = defaultParams + overrideParams
      // Convert to List[Map] for the arguments
      def newParams = argumentsAsList(updtParams) + [ "id" : id ]

      // Generate output filenames, out comes a Map
      def output = outFromIn(newParams)

      // The process expects Path or List[Path], Maps need to be converted
      def inputsForProcess =
        (input in HashMap)
          ? input.collect{ k, v -> v }.flatten()
          : input
      def outputsForProcess = output.collect{ it.value }

      // For our machinery, we convert Path -> String in the input
      def inputs =
        (input in List || input in HashMap)
          ? (input in List)
            ? input.collect{ it.name }
            : input.collectEntries{ k, v -> [ k, (v in List) ? v.collect{it.name} : v.name ] }
          : input.name
      outputs = output.collectEntries{ [(it.name): it.value] }

      def finalParams = overrideIO(newParams, inputs, outputs)

      checkParams(finalParams)

      new Tuple6(
        id,
        inputsForProcess,
        outputsForProcess,
        effectiveContainer(finalParams),
        renderCLI([finalParams.command], finalParams.arguments),
        finalParams
        )
    }

  result_ = mnn_process(id_input_output_function_cli_params_) \
    | join(id_input_params_) \
    | map{ id, output, _params, input, original_params ->
        def parsedOutput = _params.arguments
          .findAll{ it.type == "file" && it.direction == "Output" }
          .withIndex()
          .collectEntries{ it, i ->
            // with one entry, output is of type Path and array selections
            // would select just one element from the path
            [(it.name): (output in List) ? output[i] : output ]
          }
        new Tuple3(id, parsedOutput, original_params)
      }

  result_ \
    | filter { it[1].keySet().size() > 1 } \
    | view{
        ">> Be careful, multiple outputs from this component!"
    }

  emit:
  result_.flatMap{ it ->
    (it[1].keySet().size() > 1)
      ? it[1].collect{ k, el -> [ it[0], [ (k): el ], it[2] ] }
      : it[1].collect{ k, el -> [ it[0], el, it[2] ] }
  }
}

workflow {

  def id = params.id
  def _params = argumentsAsList(params.mnn) + [ "id" : id ]
  def p = _params
    .arguments
    .findAll{ it.type == "file" && it.direction == "Input" }
    .collectEntries{ [(it.name): file(params[it.name]) ] }

  def ch_ = Channel.from("").map{ s -> new Tuple3(id, p, params)}

  result = mnn(ch_)
  result.view{ it[1] }
}

// This workflow is not production-ready yet, we leave it in for future dev
// TODO
workflow test {

  take:
  rootDir

  main:
  params.test = true
  params.mnn.output = "mnn.log"

  Channel.from(rootDir) \
    | filter { params.mnn.tests.isDefined } \
    | map{ p -> new Tuple3(
        "tests",
        params.mnn.tests.testResources.collect{ file( p + it ) },
        params
    )} \
    | mnn

  emit:
  mnn.out
}