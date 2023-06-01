process star_index{

    label 'star'
    publishDir params.outdir
 
    input:
    path(reference)
    path(gtf)

    output:
    path("star/*"), emit: index

    script:
    """
    mkdir star
    STAR \\
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
    tuple val(name), path(read)
    path(index)
    path(gtf)

    output:
    path("${read.baseName}*.sam"), emit: sam 

    script:
    """
    STAR \\
    	--genomeDir . \\
    	--readFilesIn ${read} \\
    	--runThreadN ${params.threads} \\
    	--outFilterIntronMotifs RemoveNoncanonical \\
    	--outSAMattrIHstart 0 \\
    	--alignSoftClipAtReferenceEnds No \\
    	--sjdbGTFfile ${gtf} \\
    	--outFileNamePrefix ${read.baseName}.
    """	
}
