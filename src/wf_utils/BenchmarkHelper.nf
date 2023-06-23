/* usage:
| setWorkflowArguments(
  pca: [ "input": "input", "obsm_output": "obsm_pca" ]
  harmonypy: [ "obs_covariates": "obs_covariates", "obsm_input": "obsm_pca" ],
  find_neighbors: [ "obsm_input": "obsm_pca" ],
  umap: [ "output": "output" ]
)
*/

sourceDir = params.rootDir + "/src"
include { processConfig } from sourceDir + "/wf_utils/WorkflowHelper.nf"

def runMethods(Map args) {
  if (!args.methods || !args.config) {
    throw new RuntimeException("runMethods should be called as `runMethods(config: config, methods: methods)`.")
  }
  def methods_ = args.methods
  def config_ = args.config

  workflow runMethodsWf {
    take: input_ch
    main:

    // generate one channel per method
    method_chs = methods_.collect { method_module ->
      def method_config = processConfig(method_module.config)
      def method_id = method_config.functionality.name

      input_ch
        | flatMap{tup ->
          // split tuple
          def orig_id = tup[0]
          def data = tup[1]
          def rest = tup.drop(2)

          // check preferred normalisation
          def preferred_normalization = method_config.functionality.info.preferred_normalization
          if (preferred_normalization == "counts") {
            preferred_normalization = "log_cpm"
          }

          // exit function if normalization id is not what method wants
          if (data.normalization_id != preferred_normalization) {
            return []
          }

          // create new id
          def new_id = orig_id + "." + method_id

          // create args for method
          def method_arg_labels = method_config.functionality.arguments.collectMany{arg ->
            if ("label" !in arg.info) return []
            [[arg.info.label, arg]]
          }.collectEntries()

          // create dictionary with all of the args that will be passed to a method
          def new_args = config_.functionality.allArguments.collectMany{ wf_arg ->
            // TODO: using the label to match the workflow arguments to the component
            // arguments is a workaround.

            // was not passed a value for this argument
            if (wf_arg.plainName !in data) return []

            // argument did not have an id
            if ("label" !in wf_arg.info) return []

            // id was not found in method argument info
            if (wf_arg.info.label !in method_arg_labels) return []
            
            // fetch method argument info
            def meth_arg = method_arg_labels[wf_arg.info.label]

            // create new tup
            def new_tup = [meth_arg.plainName, data[wf_arg.plainName]]
            [new_tup]
          }.collectEntries()

          // store all variables in the data store
          def data_store = data.clone() + [method_id: method_id]
          data_store.remove("id")

          // create new args map
          def new_tup = [new_id, new_args, data_store] + rest

          [new_tup]
        }
        | method_module
      }

    // mix all results
    output_ch = method_chs[0].mix(*method_chs.drop(1))

    emit: output_ch
  }

  return runMethodsWf
}
