process bowtie2_index{
    label 'bowtie2'

    input:
    path(reference)
    env STRANDEDNESS
    
    output: 
    path("${reference.baseName}*"), emit: index
    
    script: 
    """
    bowtie2-build -f ${reference} ${reference.baseName}
    """
}

process bowtie2_align{
    label 'bowtie2'
    publishDir params.outdir
    
    input:
    tuple val(name), path(read)
    path(index)
    path(reference)
    
    output:
    path("${read.baseName}.sam"), emit: sam
    
    script:
    """
    bowtie2 -p ${params.threads} -x ${reference.baseName} -U ${read} -S ${read.baseName}.sam
    """
    
}
