// Alevin_fry process
process ALEVIN_FRY_MAP {

    tag "$sample"
    cpus 2
    memory '16 GB'
    time '4h'
    container 'usefulaf_0.8.1.sif'

    input:
    tuple val(sample), path(reads), path (reference)

    output:
    val("${sample}")                                                                   , emit: sample
    path "${sample}_output/${sample}-filtered-feature-bc-matrix/features.tsv.gz"       , emit: features
    path "${sample}_output/${sample}-filtered-feature-bc-matrix/matrix.mtx.gz"         , emit: matrix
    path "${sample}_output/${sample}-filtered-feature-bc-matrix/barcodes.tsv.gz"       , emit: barcodes
    path "${sample}_output/${sample}_feature.txt"
    
    when:
    task.ext.when == null  || task.ext.when

    script:
    // seperate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()

    /* define umi and bc geometry    
    def READ_geometry = '2[1-end]'
    def UMI_geometry = '1[1-6]'
    def BC_geometry = '1[7-14]'
    */
    if ( params.flow == 'sort-af' )

        """
        # Mapping the data to obtain a RAD file
        salmon alevin -i ${reference}/reference_index_idx \
            -p $task.cpus \
            -l IU \
            --read-geometry $params.READ_geometry \
            --umi-geometry $params.UMI_geometry \
            --bc-geometry $params.BC_geometry \
            --sketch \
            -1 ${forward.join(" ")} \
            -2 ${reverse.join(" ")} \
            -o ${sample}_af_map  \
            --tgMap ${reference}/t2g.tsv
        
        # Processing the mapped reads
    
        # Generating a permit list of filtered cells - Generate a permit list of barcodes from a RAD file
        alevin-fry generate-permit-list -d fw -k -i ${sample}_af_map -o ${sample}_output

        # Collating the file - Collate a RAD file by corrected cell barcode
        alevin-fry collate -t $task.cpus -i ${sample}_output -r ${sample}_af_map

        # Quantifying UMIs per-gene and per-cell - Quantify expression from a collated RAD file
        alevin-fry quant -t $task.cpus -i ${sample}_output -o ${sample}_output --tg-map ${reference}/t2g.tsv --resolution cr-like --use-mtx

        if [ -d "${sample}_output" ]; then

            mv ${sample}_output/alevin ${sample}_output/${sample}-filtered-feature-bc-matrix
            mv ${sample}_output/featureDump.txt ${sample}_output/${sample}_feature.txt
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat_rows.txt ${sample}_output/${sample}-filtered-feature-bc-matrix/barcodes.tsv
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat_cols.txt ${sample}_output/${sample}-filtered-feature-bc-matrix/features.tsv
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat.mtx ${sample}_output/${sample}-filtered-feature-bc-matrix/matrix.mtx

            # gzip count table files
            gzip ${sample}_output/${sample}-filtered-feature-bc-matrix/*
    
        fi

        """
    
    else if (params.flow == 'tenx-af')

        """
        # Mapping the data to obtain a RAD file
        salmon alevin -i ${reference}/reference_index_idx \
            -p $task.cpus \
            -l IU \
            --read-geometry $params.READ_geometry \
            --umi-geometry $params.UMI_geometry \
            --bc-geometry $params.BC_geometry \
            --sketch \
            -1 ${forward.join(" ")} \
            -2 ${reverse.join(" ")} \
            -o ${sample}_af_map  \
            --tgMap ${reference}/t2g.tsv

        # Processing the mapped reads
    
        # Generating a permit list of filtered cells - Generate a permit list of barcodes from a RAD file
        alevin-fry generate-permit-list -d fw -k -i ${sample}_af_map -o ${sample}_output

        # Collating the file - Collate a RAD file by corrected cell barcode
        alevin-fry collate -t $task.cpus -i ${sample}_output -r ${sample}_af_map

        # Quantifying UMIs per-gene and per-cell - Quantify expression from a collated RAD file
        alevin-fry quant -t $task.cpus -i ${sample}_output -o ${sample}_output --tg-map ${reference}/t2g.tsv --resolution cr-like --use-mtx

        if [ -d "${sample}_output" ]; then

            mv ${sample}_output/alevin ${sample}_output/${sample}-filtered-feature-bc-matrix
            mv ${sample}_output/featureDump.txt ${sample}_output/${sample}_feature.txt
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat_rows.txt ${sample}_output/${sample}-filtered-feature-bc-matrix/barcodes.tsv
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat_cols.txt ${sample}_output/${sample}-filtered-feature-bc-matrix/features.tsv
            mv ${sample}_output/${sample}-filtered-feature-bc-matrix/quants_mat.mtx ${sample}_output/${sample}-filtered-feature-bc-matrix/matrix.mtx

            # gzip count table files
            gzip ${sample}_output/${sample}-filtered-feature-bc-matrix/*
    
        fi

        """

}
