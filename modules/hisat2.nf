process hisat2_index{
    label 'hisat2'
    publishDir params.outdir

    input:
    path(reference)
    env STRANDEDNESS
    
    output:
    tuple path(reference), path("${reference.baseName}*.ht2"), emit: index

    script:
    """
    hisat2-build ${reference} ${reference.baseName} -p ${params.threads} 
    """
}

process hisat2_align{
    label 'hisat2'
    publishDir params.outdir
 
    input:
    val (name), path(read)
    tuple path(reference), path(index)
    env strandedness

    output:
    path("${read.baseName}.sam"), emit: sam 

    shell:
    '''
    if [[ ($strandedness == "firststrand") ]]; then
    
        hisat2 -x !{reference.baseName} -U !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} --rna-strandness FR -S !{read.baseName}.sam

    elif [[ ($strandedness == "secondstrand") ]]; then
    
        hisat2 -x !{reference.baseName} -U !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} --rna-strandness RF -S !{read.baseName}.sam

    elif [[ $strandedness == "unstranded" ]]; then
       
        hisat2 -x !{reference.baseName} -U !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} -S !{read.baseName}.sam
    fi
    '''   
   
}

