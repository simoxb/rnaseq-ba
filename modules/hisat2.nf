process hisat2_index{
    label 'hisat2'
    publishDir params.outdir

    input:
    path(fasta)
    
    output:
    tuple path(fasta), path("${fasta.baseName}*.ht2"), emit: index

    script:
    """
    hisat2-build ${fasta} ${fasta.baseName} -p ${params.threads} 
    """
}

process hisat2_align{
    label 'hisat2'
    publishDir params.outdir
 
    input:
    tuple val(name), path(read)
    tuple path(fasta), path(index)

    output:
    path("${read.baseName}.sam"), emit: sam 

    shell:
    '''
    if [[ (!params.strandedness == "firststrand") ]]; then
    
        hisat2 -x !{fasta.baseName} -1 !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} --rna-strandness FR -S !{read.baseName}.sam

    elif [[ (!params.strandedness == "secondstrand") ]]; then
    
        hisat2 -x !{fasta.baseName} -1 !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} --rna-strandness RF -S !{read.baseName}.sam

    elif [[ !params.strandedness == "unstranded" ]]; then
       
        hisat2 -x !{fasta.baseName} -1 !{read} --new-summary --summary-file !{read.baseName}_summary.log --thread !{params.threads} -S !{read.baseName}.sam
    fi
    '''   
   
}

