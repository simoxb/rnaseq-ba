process salmon_quant{
	label 'salmon'

	input: 
	path(bam)
	path(fasta)
	env strandedness
	
	output:
	path("quantification_results.*")
	
	shell:
	'''
	if [[ ($strandedness == "firststrand") ]]; then
    
    	salmon quant -t ${fasta} -l ISR -a ${bam} -o quantification_results
    	
	elif [[ ($strandedness == "secondstrand") ]]; then
    
        salmon quant -t ${fasta} -l ISF -a ${bam} -o quantification_results
        
	elif [[ $strandedness == "unstranded" ]]; then
       
        salmon quant -t {fasta} -l IU -a ${bam} -o quantification_results
        
	fi
	'''   
 

}
