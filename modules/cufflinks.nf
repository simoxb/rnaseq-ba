process cufflinks {
    label 'cufflinks'
    publishDir params.outdir
     
    input:
    path(sorted_bam)
    path(annotation)
    env STRANDEDNESS
    
    output:
    path('transcripts.gtf'), emit: cufflinks_gtf 
    
    shell:
    '''
    if [[ $STRANDEDNESS == "firststrand" ]]; then
         /work/simon/bin/cufflinks -G !{annotation} !{sorted_bam} --library-type fr-firststrand --num-threads !{params.threads}
    elif [[ $STRANDEDNESS == "secondstrand" ]]; then
         /work/simon/bin/cufflinks -G !{annotation} !{sorted_bam} --library-type fr-secondstrand --num-threads !{params.threads}
    elif [[ $STRANDEDNESS == "unstranded" ]]; then
         /work/simon/bin/cufflinks -G !{annotation} !{sorted_bam} --library-type fr-unstranded --num-threads !{params.threads}
	else  
		echo $STRANDNESS > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
    fi
    '''

}
