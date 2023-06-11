process KALLISTO_BUSTOOLS_MAP {

    tag "$sample"
    cpus 2
    memory '16 GB'
    time '4h'

    container 'https://depot.galaxyproject.org/singularity/kb-python:0.27.3--pyhdfd78af_1'

    input:

    tuple val(sample), path(reads), path(index)
    path whitelist

    output:
    val("${sample}")
    path("${sample}_output/transcripts.txt")                            , emit: transcripts
    path("${sample}_output/filter_barcodes.txt")                        , emit: filter_barcodes
    path("${sample}_output/${sample}-filtered-feature-bc-matrix")       , emit: counts_filtered
    path("${sample}_output/${sample}-raw")                              , emit: counts_unfiltered
    path("${sample}_output/inspect.json")                               , emit: inspect
    path("${sample}_output/kb_info.json")                               , emit: kb_info
    path("${sample}_output/matrix.ec")                                  , emit: matrix
    path("${sample}_output/output.bus")                                 , emit: output_bus
    path("${sample}_output/output.unfiltered.bus")                      , emit: output_unfiltered
    path("${sample}_output/output.filtered.bus")                        , emit: output_filtered
    path("${sample}_output/run_info.json")                              , emit: run_info

    when:
    task.ext.when == null || task.ext.when
    
    script:

    """
    kb count \
        --h5ad \
        -i ${index}/index.idx \
        -g ${index}/t2g.txt \
        -x $params.kb_protocol \
        -o ${sample}_output \
        --filter bustools \
        -t $task.cpus \
        -w ${whitelist} \
        ${reads}

    if [ -d "${sample}_output/" ]; then

        mv ${sample}_output/counts_filtered/cells_x_genes.barcodes.txt ${sample}_output/counts_filtered/barcodes.tsv
        mv ${sample}_output/counts_filtered/cells_x_genes.genes.txt ${sample}_output/counts_filtered/features.tsv
        mv ${sample}_output/counts_filtered/cells_x_genes.mtx ${sample}_output/counts_filtered/matrix.mtx
        mv ${sample}_output/counts_filtered ${sample}_output/${sample}-filtered-feature-bc-matrix
        
        mv ${sample}_output/counts_unfiltered/cells_x_genes.barcodes.txt ${sample}_output/counts_unfiltered/barcodes.tsv
        mv ${sample}_output/counts_unfiltered/cells_x_genes.genes.txt ${sample}_output/counts_unfiltered/features.tsv
        mv ${sample}_output/counts_unfiltered/cells_x_genes.mtx ${sample}_output/counts_unfiltered/matrix.mtx
        mv ${sample}_output/counts_unfiltered ${sample}_output/${sample}-raw

        # gzip count table files
        gzip ${sample}_output/${sample}-filtered-feature-bc-matrix/*
    
    fi

    """

}