process star_index{

    label 'star'
    publishDir params.outdir
 
    input:
    path(reference)
    path(gtf)
    env STRANDEDNESS

    output:
    path("star/*"), emit: index

    script:
    """
    mkdir star
    /work/simon/bin/STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
	    --genomeFastaFiles ${reference} \\
            --sjdbGTFfile ${gtf} \\
	    --runThreadN ${params.threads}
    """
}

process star_align{
    label 'star'
    publishDir params.outdir
 
    input:
    path(read)
    path(index)
    path(gtf)
    env STRANDEDNESS

    output:
    path("${read.baseName}*.sam"), emit: sam 

    shell:
    '''
    if [[ ($STRANDEDNESS == "firststrand") || ($STRANDEDNESS == "secondstrand") ]]; then
    /work/simon/bin/STAR \\
          --genomeDir . \\
          --readFilesIn !{read}\\
          --alignSoftClipAtReferenceEnds No \\
          --outFileNamePrefix !{read.baseName}. \\
          --sjdbGTFfile !{gtf} \\
	  --outFilterIntronMotifs RemoveNoncanonical \\
	  --outSAMattrIHstart 0 \\
 	  --runThreadN !{params.threads}

    elif [[ $STRANDEDNESS == "unstranded" ]]; then
       /work/simon/bin/STAR \\
          --genomeDir . \\
          --readFilesIn !{read} \\
	  --outFilterIntronMotifs RemoveNoncanonical \\
          --alignSoftClipAtReferenceEnds No \\
	  --outSAMstrandField intronMotif \\
          --outFileNamePrefix !{read.baseName}. \\
          --sjdbGTFfile !{gtf} \\
	  --outSAMattrIHstart 0 \\
	  --runThreadN !{params.threads}
    else  
		echo $STRANDEDNESS > error_strandedness.txt
		echo "strandness cannot be determined" >> error_strandedness.txt
	fi

    '''   
   
}
