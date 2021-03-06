nextflow.enable.dsl=2

/* for now, you need to manually specify the
 * root directory of this repository as follows.
 * (it's a nextflow limitation I'm trying to figure out 
 * how to resolve.) */
rootDir = "$projectDir/../../.." 

// target dir containing the nxf modules generated by viash
targetDir = "$rootDir/target/nextflow"

// import dataset loaders
include { download_datasets }       from "$targetDir/trajectory_inference/datasets/download_datasets/main.nf"     params(params)

// import methods

// import metrics

// import helper functions
include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

workflow {
    main:
        output_ = Channel.fromPath(file("$launchDir/src/trajectory_inference/datasets/download_datasets/datasets.tsv")) \
            | splitCsv(header: true, sep: "\t") \
            | map { row ->
                newParams = overrideParams(params, "download_datasets", "id", row.id)
                [ row.id, file(row.url), newParams]
            } \
            | download_datasets
    emit:
        output_
}

