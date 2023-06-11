process TENX_MAP_STAR {
    tag "$sample"
    cpus 16
    memory '40GB'
    time '4h'

    conda (params.enable_conda ? 'bioconda::star=2.7.10a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10a--h9ee0642_0' :
        'quay.io/biocontainers/star:2.7.10a--h9ee0642_0' }"

    input:
    //
    // Input reads are expected to come as: [ sample, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(sample), path(reads), path(index)
    path whitelist 
 
    output:
    val("${sample}")                                                , emit: sample
    tuple val(sample), path('*d.out.bam')                           , optional:true, emit: bam
    tuple val(sample), path('*.Solo.out')                           , emit: counts
    path('*Log.final.out')                                          , emit: log_final
    path('*Log.out')                                                , emit: log_out
    path('*Log.progress.out')                                       , emit: log_progress
    path("${sample}.seq_reads.txt")                                 , emit: no_reads
    path("${sample}.Solo.out/Gene/${sample}-raw")                   , emit: raw_feature_bc_matrix
    path("${sample}.Solo.out/Gene/${sample}-filtered-bc-matrix")    , emit: filtered_feature_bc_matrix
    path("${sample}.Solo.out/GeneFull/${sample}-raw")               , emit: genefull_raw_feature_bc_matrix
    path("${sample}.Solo.out/GeneFull/${sample}-filtered-bc-matrix"), emit: genefull_filtered_feature_bc_matrix
    tuple val(sample), path('*sortedByCoord.out.bam')               , optional:true, emit: bam_sorted
    tuple val(sample), path('*toTranscriptome.out.bam')             , optional:true, emit: bam_transcript
    tuple val(sample), path('*fastq.gz')                            , optional:true, emit: fastq
    tuple val(sample), path('*.tab')                                , optional:true, emit: tab
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_folder = params.refs_path
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample}"

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    STAR \\
        --genomeDir ${reference_folder}/${index} \\
        --readFilesIn ${reverse.join( "," )} ${forward.join( "," )} \\
        --runThreadN $task.cpus \\
        --outBAMsortingThreadN 4 \\
        --outFileNamePrefix $prefix. \\
        --soloCBwhitelist $whitelist \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloUMIstart $params.UMI_start \\
        --soloUMIlen $params.UMI_len \\
        --soloCBstart $params.BC_start \\
        --soloCBlen $params.BC_len \\
        $args \\
        --soloCBmatchWLtype $params.soloCBmatchWLtype \\
        --soloUMIdedup $params.soloUMIdedup \\
        --soloStrand $params.soloStrand \\
        --soloFeatures $params.soloFeatures \\
        --soloMultiMappers $params.soloMultiMappers \\
        --soloUMIfiltering $params.soloUMIfiltering \\
        --soloCellFilter $params.soloCellFilter \\
        --soloBarcodeReadLength 0 \\
        --soloCellReadStats Standard \\
        --outSAMtype BAM SortedByCoordinate \\
        --genomeSAsparseD 3 # to make the agreement between STARsolo and CellRanger even more perfect

    # get the number of raw reads from STAR log
    echo 'Number of Reads' > "${prefix}.seq_reads.txt"
    cat ${prefix}.Log.final.out | grep 'Number of input reads' | cut -d\$'\t' -f 2 >> "${prefix}.seq_reads.txt"
    
    # rename Gene output with sampleID prefix
    if [ -d "${prefix}.Solo.out/Gene" ]; then
        mv ${prefix}.Solo.out/Gene/filtered ${prefix}.Solo.out/Gene/${prefix}-filtered-bc-matrix
        mv ${prefix}.Solo.out/Gene/raw ${prefix}.Solo.out/Gene/${prefix}-raw
        mv ${prefix}.Solo.out/Gene/Summary.csv ${prefix}.Solo.out/Gene/${prefix}-summary.csv
        mv ${prefix}.Solo.out/Gene/Features.stats ${prefix}.Solo.out/Gene/${prefix}-features.stats
        mv ${prefix}.Solo.out/Gene/CellReads.stats ${prefix}.Solo.out/Gene/${prefix}-cell-reads.stats
        mv ${prefix}.Solo.out/Gene/UMIperCellSorted.txt ${prefix}.Solo.out/Gene/${prefix}-UMIperCellSorted.txt
    fi

    # rename GeneFull output with sampleID prefix
    if [ -d "${prefix}.Solo.out/GeneFull" ]; then
        mv ${prefix}.Solo.out/GeneFull/filtered ${prefix}.Solo.out/GeneFull/${prefix}-filtered-bc-matrix
        mv ${prefix}.Solo.out/GeneFull/raw ${prefix}.Solo.out/GeneFull/${prefix}-raw
        mv ${prefix}.Solo.out/GeneFull/Summary.csv ${prefix}.Solo.out/GeneFull/${prefix}-summary.csv
        mv ${prefix}.Solo.out/GeneFull/Features.stats ${prefix}.Solo.out/GeneFull/${prefix}-features.stats
        mv ${prefix}.Solo.out/GeneFull/UMIperCellSorted.txt ${prefix}.Solo.out/GeneFull/${prefix}-UMIperCellSorted.txt
    fi

    # rename metrics summary
    mv ${prefix}.Solo.out/Barcodes.stats ${prefix}.Solo.out/${prefix}-barcodes.stats


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}