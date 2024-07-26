
// Check input path parameters to see if they exist
def checkPathParamList = [
    //params.mixture,
    //params.signature,
    params.input_dir,
    params.output_dir
]
// create file channel for parameters provided
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

//ch_mixture = Channel.fromPath(params.mixture, checkIfExists: true)
//ch_signature = Channel.fromPath(params.signature, checkIfExists: true)
ch_mixture = Channel.fromPath(params.mixture, checkIfExists: true)
ch_signature = Channel.fromPath(params.signature, checkIfExists: true)
ch_input_dir = Channel.fromPath(params.input_dir, checkIfExists: true, type: "dir")
ch_output_dir = Channel.fromPath(params.output_dir, checkIfExists: true, type: "dir")
ch_username = Channel.of(params.username)
ch_token = Channel.of(params.token)

process CIBERSORT {
    label 'cibersort'
    publishDir "$params.outdir/cibersort"
    stageInMode 'copy'
    stageOutMode 'copy'

    input:
    path(input_dir) // ch_metadata
    path(output_dir) // ch_tx2gene
    path(signature)
    path(mixture)
    val(username)
    val(token)
    output:
    // TODO
    path("CIBERSORTx_Results.txt"), emit: ch_cibersort_out
    script:
    """
    docker run -v ${input_dir}:/src/data -v ${output_dir}:/src/outdir cibersortx/fractions --username $username --token $token --single_cell TRUE --refsample $signature --mixture $mixture
    """
}

workflow {
    ch_signature.view()
    CIBERSORT(
        ch_input_dir,
        ch_output_dir,
        ch_signature,
        ch_mixture,
        ch_username,
        ch_token
    )
}