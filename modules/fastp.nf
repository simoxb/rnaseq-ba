process fastp{		
	
	label 'fastp'
	publishDir params.outdir	
	
	input: 
	path(read)
	
	output: 
	path("${read.baseName}*fastp*.f*q"), emit: trimmed
		
	script:
	"""
	fastp -i ${read} -o ${read.baseName}_fastp.1.fastq --detect_adapter_for_pe --json ${read.baseName}_fastp.json --html ${read.baseName}_fastp.html --thread ${params.threads}
	"""
}

