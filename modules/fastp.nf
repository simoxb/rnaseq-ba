process fastp{		
	
	label 'fastp'
	publishDir params.outdir	
	
	input: 
	path(read)
	
	output: 
	path("${read.baseName}_fastp.fastq"), emit: trimmed
		
	script:
	"""
	fastp -i ${read} -o ${read.baseName}_fastp.fastq --detect_adapter_for_pe --json ${read.baseName}_fastp.json --html ${read.baseName}_fastp.html --thread ${params.threads}
	"""
}

