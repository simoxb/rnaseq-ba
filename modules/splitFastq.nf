process fastqsplit{
  
    publishDir params.outdir
    label 'python'
    
    input:
    path(fastq)

    output:
    path("*${fastq.getBaseName()}.fastq"), emit: splitted

    shell:
    '''
    length=$(wc -l < !{fastq})
    length=$((length / 4))
    s=$(echo !{params.split})
    z=$((length / s))
    splitby=$((${z} + 1))
    !{params.baseDir}/bin/splitFastq -i !{fastq} -n ${splitby} -o !{fastq.getBaseName()}
    '''
}
