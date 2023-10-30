// workflow auto {
//   findStates(params, meta.config)
//     | meta.workflow.run(
//       auto: [publish: "state"]
//     )
// }

workflow run_wf {
    take:
    input_ch

    main:
    output_ch = input_ch
        
        | view{ id, state ->
            state
        }

        | get_results.run(
            fromState: [ 
                "input_scores": "input_scores",
                "input_execution" : "input_execution",
                "task_id" : "task_id",
                "output": "output_scores"
            ],
            toState: { id, output, state ->
                state + [output: output.output]}
        )
    emit:
    output_ch
}