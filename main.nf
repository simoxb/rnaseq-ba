nextflow.enable.dsl=2


include {fastp} from "./modules/fastp"
include {fastqsplit} from "./modules/splitFastq"
include {samtools; samtools_merge} from "./modules/samtools"
include {star_index; star_align} from "./modules/star"
include {hisat2_index; hisat2_align} from "./modules/hisat2"
include {tophat2_index; tophat2_align} from "./modules/tophat2"
include {salmon_quant} from "./modules/salmon"

workflow rnaseq_star{

	take:
	input_read

	main:
	fastp(input_read)
	star_index(params.fasta, params.gtf)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
		 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		
		star_align(splitted_ch, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.fasta, params.strandedness)
	}else{
		star_align(fastp.out.trimmed, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
		salmon_quant(samtools.out, params.fasta, params.strandedness)
	}
}


workflow rnaseq_hisat2{

	take:
	input_read

	main:
	fastp(input_read)
	hisat2_index(params.fasta)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
	  	 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		 
		hisat2_align(splitted_ch, hisat2_index.out.index, params.strandedness)
		samtools(hisat2_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.fasta, params.strandedness)
	}else{
		hisat2_align(fastp.out.trimmed, hisat2_index.out.index)
		samtools(hisat2_align.out.sam)
		salmon_quant(samtools.out, params.fasta, params.strandedness)
	}
}

workflow rnaseq_tophat2{

	take:
	input_read

	main:
	fastp(input_read)
	tophat2_index(params.fasta)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
	  	 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		 
		tophat2_align(splitted_ch, tophat2_index.out.index, params.fasta)
		samtools(tophat2_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.fasta, params.strandedness)
	}else{
		tophat2_align(fastp.out.trimmed, tophat2_index.out.index, params.fasta)
		samtools(tophat2_align.out.sam)
		salmon_quant(samtools.out, params.fasta, params.strandedness)
	}
}

workflow{
	if(params.aligner=="star"){
	rnaseq_star(params.read)
	}
	
	if(params.aligner=="hisat2"){
	rnaseq_hisat2(params.read)
	}
	
	if(params.aligner=="tophat2"){
	rnaseq_tophat2(params.read)
	}	
}
