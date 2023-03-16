#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.metadata,
    params.tx2gene
]
// create file channel for parameters provided
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

process GENERATE_METADATA {
    publishDir "$params.outdir/metadata", mode: 'copy'
    input:
    val(out_full)
    val(out_reduced)
    val(basedir)

    output:
    path("*_full.txt"), emit: ch_generated_metadata
    path("*_no_inhouse.txt"), emit: ch_generated_metadata_reduced
    script:
    """
    Rscript ${projectDir}/bin/create_metadata.R --out_full $out_full --out_reduced $out_reduced --basedir $basedir
    """
}

process GENERATE_TX2GENE {
    publishDir "$params.outdir/tx2gene", mode: 'copy'

    input:
    path(gtf)

    output:
    path("*.txt"), emit: ch_tx2gene_out

    script:
    """
    Rscript ${projectDir}/bin/generate_tx2gene_table.R --gtf ${gtf} --out "tx2gene.txt"
    """
}
 //, mode: 'copy', overwrite = true

process HLATYPING_ADD_SAMPLENAME {
    publishDir "$params.outdir/hla_calls", mode: 'copy'
    input:
    path(hlatyping_files) // ch_hlatyping

    output:
    path("*.tsv"), emit: ch_hla_samplename_added

    script:
    prefix = hlatyping_files.simpleName 
    """
    awk -v OFS='\t' '{print FILENAME, \$0}' ${hlatyping_files} | tail -n1 >> ${prefix}_with_samplename.tsv
    """
}

process MATCH_BY_PATIENT {
    publishDir "$params.outdir/metadata_by_patient", mode: 'copy'

    input:
    path(metadata) // ch_metadata
    path(hlatyping_output) // ch_hlatyping
    path(exclude_files) // ch_exclude_files
    //path(hlatyping_dir) // ch_hlatyping_dir
    output:
    //path("metadata_with_patient.txt"), emit: ch_metadata_with_patient
    path("HLAtyping_pipeline.html"), emit: ch_hla_match
    shell:
        script = "rmarkdown::render('HLAtyping_pipeline.Rmd',"
        script += "params = list("
        script += "metadata = '\$PWD/${metadata}', "
        script += "indir = '\$PWD/${hlatyping_dir}', "
        script += "files_to_exclude = '\$PWD/${exclude_files}', "
        script += "all_files = '\$PWD/${hlatyping_output}'), "
        script += ")"
    """
    cp -L ${projectDir}/bin/HLAtyping_pipeline.Rmd HLAtyping_pipeline.Rmd 
    Rscript -e "${script}"
    """

}
process QC {
    publishDir "$params.outdir/qc"
    stageInMode 'copy'
    stageOutMode 'copy'

    input:
    path(metadata) // ch_metadata
    path(tx2gene) // ch_tx2gene

    output:
    // TODO
    path("out.txt"), emit: ch_qc_out
    path("QC_pipeline.md")
    path("QC_pipeline.html")
    path("metadata_outliers_removed_ensembl_gene_id_version.txt"), emit: ch_metadata_outliers_removed
    path("dds_remove_outliers_batch_corrected.Rds"), emit: ch_dds_remove_outliers_batch_corrected
    path("dds_not_corrected_remove_outliers_no_inhouse.Rds"), emit: ch_dds_not_corrected_no_inhouse
    //script:
    shell:
        script = "rmarkdown::render('QC_pipeline.Rmd',"
        script += "params = list("
        script += "metadata = '\$PWD/${metadata}', "
        script += "tx2gene = '\$PWD/${tx2gene}', "
        script += "out = 'out.txt'),"
        script += ")"
    """
    cp -L ${projectDir}/bin/QC_pipeline.Rmd QC_pipeline.Rmd 
    Rscript -e "${script}"
    """
}

/*
process PREPARE_DATA_CIBERSORT{
    publishDir "$params.outdir/tpm"
    stageInMode 'copy'
    stageOutMode 'copy'

    input:
    path(dds_in) // ch_dds_remove_outliers_batch_corrected
    path(metadata) //metadata_full.txt - ch_metadata - for finding the quant_genes.sf

    output:
    path("caf_tpm_mixture_batch_corrected.txt")
    path("caf_tpm_for_signature_batch_corrected.txt")
    path("phenoclasses_caf.txt")
    // TODO: construct command, rscript command, add to workflow, check to make sure it works okay

}
*/
 
ch_out_full = Channel.of(params.metadata_out_full)
ch_out_reduced = Channel.of(params.metadata_out_reduced)
ch_basedir = Channel.of(params.basedir)
ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
//ch_hlatyping = Channel.fromPath(params.hlatyping_output, checkIfExists: true)
//ch_hlatyping.view()

ch_hlatyping = Channel.fromPath( params.hlatyping_output ) 


ch_hlatyping_dir = Channel.fromPath(params.hlatyping_outdir)


//ch_hlatyping.view()

/*
workflow {
  
  Channel
    .fromPath( params.SRRs )
    .splitText() { it.strip() }
    .view()
}
*/
ch_exclude_files = Channel.fromPath(params.exclude_files, checkIfExists: true)


workflow {
HLATYPING_ADD_SAMPLENAME(ch_hlatyping)

ch_merged_hlatyping_calls = HLATYPING_ADD_SAMPLENAME.out.ch_hla_samplename_added
                            .collectFile(name: 'combined_hla_calls.txt', newLine: true)
                        

/*
if (!params.metadata){
    // TODO run generate_metadata 
    println "no metadata provided, generating..."
    GENERATE_METADATA(ch_out_full,
                      ch_out_reduced,
                      ch_basedir)
    ch_metadata = GENERATE_METADATA.out.ch_generated_metadata
} else {
    println "metadata provided, not generating..."
    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    ch_metadata.view()
}

if (!params.tx2gene){
    println "no tx2gene provided, generating..."
    GENERATE_TX2GENE(
        ch_gtf
    )
    ch_tx2gene = GENERATE_TX2GENE.out.ch_tx2gene_out
} else{
    ch_tx2gene = Channel.fromPath(params.tx2gene, checkIfExists: true)
}

MATCH_BY_PATIENT(
    ch_metadata,
    ch_hlatyping.collect(),
    ch_exclude_files,
    ch_hlatyping_dir
)
/*
QC(
    ch_metadata,
    ch_tx2gene
)
*/


}
