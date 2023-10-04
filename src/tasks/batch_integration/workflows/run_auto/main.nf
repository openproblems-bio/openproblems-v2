// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = collectTraces()


workflow run_wf {
  findStates(params, config)
    | run_wf
    | publishStates([:])
}


// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  // writeJson(methods.collect{it.config}, file("$publish_dir/methods.json"))
  // writeJson(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}