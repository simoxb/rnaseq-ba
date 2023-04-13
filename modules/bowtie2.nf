process bowtie2_index{
    label 'bowtie2'

    input:
    path(fasta)
    
    output: 
    path("${fasta.baseName}*"), emit: index
    
    script: 
    """
    bowtie2-build -f ${fasta} ${fasta.baseName}
    """
}

process bowtie2_align{
    label 'bowtie2'
    publishDir params.outdir
    
    input:
    tuple val(name), path(read)
    path(index)
    path(fasta)
    
    output:
    path("${read.baseName}.sam"), emit: sam
    
    script:
    """
    bowtie2 -p ${params.threads} -x ${fasta.baseName} -U ${read} -S ${read.baseName}.sam
    """
    
}
