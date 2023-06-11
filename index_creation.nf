// Make the working directory for analysis and build the index

process INDEX_CREATION {
    
    input:
    path(transcriptome_extracted)                   // Transcriptome that was created in gffread_transcriptome

    output:

    path("reference_index_idx")                     , emit: reference_index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    salmon index -t ${transcriptome_extracted} -i "reference_index_idx" -p 16
    """


}