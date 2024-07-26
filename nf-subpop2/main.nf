#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.metadata,
    params.tx2gene,
    params.counts
]
// create file channel for parameters provided
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

process BATCH_EFFECT {
    label 'main_docker'
    publishDir "$params.outdir/batch_effect_removed"
    stageInMode 'copy'
    stageOutMode 'copy'

    input:
    path(metadata) // ch_metadata
    path(tx2gene) // ch_tx2gene
    //path(inhouse_metadata) // ch_inhouse_metadata
    //path(metadata_without_patient) // ch_metadata
    path(counts_matrix) // ch_counts
    output:
    // TODO
    path("20231127-nextflow-caf-batch-correction-subtypes-ptoforigin.md"), optional: true
    path("20231127-nextflow-caf-batch-correction-subtypes-ptoforigin.html"), optional: true
    path("batch_corrected.txt"), emit: ch_metadata_outliers_removed
    //path("dds_remove_outliers_batch_corrected.Rds"), emit: ch_dds_remove_outliers_batch_corrected
    //path("dds_not_corrected_remove_outliers_no_inhouse.Rds"), emit: ch_dds_not_corrected_no_inhouse
    //script:
    shell:
        script = "rmarkdown::render('20231127-nextflow-caf-batch-correction-subtypes-ptoforigin.Rmd',"
        script += "params = list("
        script += "metadata = '\$PWD/${metadata}', "
        script += "tx2gene = '\$PWD/${tx2gene}', "
        script += "method = 'nextflow', "
        script += "counts_matrix = '\$PWD/${counts_matrix}', "
        script += "out = 'batch_corrected.txt', output_dir = getwd())"
        script += ")"
    """
    cp -L ${projectDir}/bin/20231127-nextflow-caf-batch-correction-subtypes-ptoforigin.Rmd 20231127-nextflow-caf-batch-correction-subtypes-ptoforigin.Rmd 
    Rscript -e "${script}"
    """
}

ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
ch_tx2gene = Channel.fromPath(params.tx2gene, checkIfExists: true)
ch_counts = Channel.fromPath(params.counts, checkIfExists: true)

/*
workflow {
  
  Channel
    .fromPath( params.SRRs )
    .splitText() { it.strip() }
    .view()
}
*/
//ch_exclude_files = Channel.fromPath(params.exclude_files, checkIfExists: true)


workflow {

BATCH_EFFECT(
    ch_metadata,
    ch_tx2gene,
    ch_counts
)
}
