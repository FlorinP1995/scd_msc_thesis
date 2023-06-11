//
// Kallisto-Bustools-based flow
//

// make Kallisto-Bustools Reference
process KB_MAKE_REFERENCE {

    tag "mkref"
    cpus 16
    memory '16 GB'
    time '1h'

    container 'https://depot.galaxyproject.org/singularity/kb-python:0.27.3--pyhdfd78af_1'

    input:
    val reference_name
    path fasta
    path gtf
    
    output:
    path("index.idx")               , emit: index
    path("${reference_name}")       , emit: kb_cdna
    path("t2g.txt")                 , emit: t2g

    when:
    task.ext.when == null || task.ext.when
                                        
    script:
    """
    kb ref \
        -i index.idx \
        -g t2g.txt \
        -f1 ${reference_name} \
        ${fasta} \
        ${gtf}
    """


/*
    take: 
        ref // name of the refrence (and of the output folder where the reference is saved)
        fasta
        gtf

    main:

        // Check mandatory parameters
        if (!params.ref) { 
            exit 1, 'Genome reference name not specified!' 
            }

        if (!params.fasta) { 
            exit 1, 'Genome FASTA file not specified!' 
            }

        if (!params.gtf) { 
            exit 1, 'Genome GTF file not specified!' 
            }

        ch_versions = Channel.empty()

        // make Alevin-Fry reference
        GFFREAD_TRANSCRIPTOME ( ref, fasta, gtf )
        INDEX_CREATION(GFFREAD_TRANSCRIPTOME.out.reference)
        T2G(gtf)
        ch_versions = ch_versions.mix(GFFREAD_TRANSCRIPTOME.out.versions)  
    
    emit:
        reference   = GFFREAD_TRANSCRIPTOME.out.reference
        index       = INDEX_CREATION.out.reference_index
        t2g         = T2G.out.t2g
        versions    = ch_versions 
*/
    }