// Creating the reference transcriptome file

process GFFREAD_TRANSCRIPTOME {

    tag "mkref"
    cpus 6
    memory '24 GB'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hd03093a_1' :
        'quay.io/biocontainers/gffread:0.12.7--hd03093a_1' }"
        
    input:
    val reference_name
    path fasta
    path gtf

    output:

    path "${reference_name}"                        , emit: reference
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread -w $reference_name -g $fasta $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """


}