process star_index{

    label 'star'
    publishDir params.outdir
 
    input:
    path(fasta)
    path(gtf)

    output:
    path("star/*"), emit: index

    script:
    """
    mkdir star
    STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
	    --genomeFastaFiles ${fasta} \\
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

    output:
    path("${read.baseName}*.sam"), emit: sam 

    script:
    """
    STAR \\
    	--genomeDir . \\
    	--readFilesIn ${read} \\
    	--runThreadN ${params.threads} \\
    	--outFileNamePrefix ${read.baseName}.
    """	
}
