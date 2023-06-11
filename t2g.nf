// Obtaining the transcript to gene file

process T2G {

    input:
    path(gtf)
    
    output:

    path("t2g.tsv")                                     , emit: t2g // t2g.tsv

    when:
    task.ext.when == null || task.ext.when

    shell:    
    '''
    awk '{if($3=="transcript") {OFS="\t"; print $14, $10} }' !{gtf} | sed 's/[;\"]//g' > t2g.tsv
    '''


}