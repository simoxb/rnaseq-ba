process tophat2_index{
    label 'tophat2'

    input:
    path(fasta)
    
    output: 
    path("${fasta.baseName}*"), emit: index
    
    script: 
    """
    bowtie2-build -f ${fasta} ${fasta.baseName}
    """
}

process tophat2_align{
    label 'tophat2'
    publishDir params.outdir
    
    input:
    tuple val(name), path(read)
    path(index)
    path(fasta)
    
    output:
    path("tophat_out/${read.baseName}.bam"), emit: sam
    
    script:
    """
    tophat2 -p ${params.threads} ${fasta.baseName} ${read}
    mv ./tophat_out/accepted_hits.bam ./tophat_out/${read.baseName}.bam
    """
    
}
