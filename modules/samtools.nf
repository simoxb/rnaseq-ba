process samtools {
    label 'samtools'
    
    
    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    path("${bam_file}.sorted.bam")
    
    script:
    """
    samtools view -b ${sam_file} | samtools sort -o ${bam_file}.sorted.bam -T tmp  
    """
    
}

process samtools_merge {
    label 'samtools'
    publishDir params.outdir

    input:
    path(bam_files)
    
    output:
    path("alignement_gathered.bam")
    
    script:
    """
    samtools merge alignement_gathered.bam ${bam_files}
    """
}
