//
// Alevin-Fry-based flow
//

// MODULES: customized modules
include { GFFREAD_TRANSCRIPTOME              } from '../modules/local/alevin_fry/mkref/gffread_transcriptome'
include { T2G                                } from '../modules/local/alevin_fry/t2g'
include { INDEX_CREATION                     } from '../modules/local/alevin_fry/index_creation'



// make Alevin-Fry reference
workflow ALEVINFRY_MKREF {

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

    }