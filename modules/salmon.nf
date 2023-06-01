process salmon_quant{
	label 'salmon'

	input: 
	path(bam)
	path(transcript)
	env strandedness
	
	output:
	path("quantification_results/*")
	
	shell:
	'''
	if [[ ($strandedness == "firststrand") ]]; then
    
    	salmon quant -t !{transcript} -l ISR -a !{bam} -o quantification_results
    	
	elif [[ ($strandedness == "secondstrand") ]]; then
    
        salmon quant -t !{transcript} -l ISF -a !{bam} -o quantification_results
        
	elif [[ $strandedness == "unstranded" ]]; then
       
        salmon quant -t !{transcript} -l IU -a !{bam} -o quantification_results
        
	fi
	'''   
 

}
