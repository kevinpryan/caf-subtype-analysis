#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.inhouse_metadata,
    params.tx2gene_file,
    params.gtex_file,
    params.metadata
]
// create file channel for parameters provided
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

ch_method = Channel.of(params.method)
inhouse_metadata = Channel.fromPath(params.inhouse_metadata, checkIfExists: true)
tx2gene_file = Channel.fromPath(params.tx2gene_file, checkIfExists: true)
gtex_file = Channel.fromPath(params.gtex_file, checkIfExists: true)
metadata = Channel.fromPath(params.metadata, checkIfExists: true)


process BATCHQC {
    label 'batch_docker'
    publishDir "$params.outdir/batchqc", mode: 'copy'
    input:
    path(inhouse_metadata)
    path(tx2gene_file)
    path(gtex_file)
    val(ch_method)
    path(metadata)

    output:
    path("batchqc_*.html"), emit: ch_batchqc

    script:
    """
    Rscript ${projectDir}/bin/batchqc.R --inhouse_metadata ${inhouse_metadata} --tx2gene_file ${tx2gene_file} --gtex_file ${gtex_file} --method ${ch_method} --metadata ${metadata}
    """
}

workflow {
    BATCHQC(
        inhouse_metadata,
        tx2gene_file,
        gtex_file,
        ch_method,
        metadata
    )
}
