nextflow.enable.dsl=2


include {fastp} from "./modules/fastp"
include {fastqsplit} from "./modules/splitFastq"
include {samtools; samtools_merge} from "./modules/samtools"
include {star_index; star_align} from "./modules/star"
include {hisat2_index; hisat2_align} from "./modules/hisat2"
include {bowtie2_index; bowtie2_align} from "./modules/bowtie2"
include {salmon_quant} from "./modules/salmon"

workflow rnaseq_star{

	take:
	input_read

	main:
	fastp(input_read)
	star_index(params.reference, params.gtf)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
		 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		
		star_align(splitted_ch, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.transcript, params.strandedness)
	}else{
		star_align(fastp.out.trimmed, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
		salmon_quant(samtools.out, params.transcript, params.strandedness)
	}
}


workflow rnaseq_hisat2{

	take:
	input_read

	main:
	fastp(input_read)
	hisat2_index(params.reference)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
	  	 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		 
		hisat2_align(splitted_ch, hisat2_index.out.index, params.strandedness)
		samtools(hisat2_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.transcript, params.strandedness)
	}else{
		hisat2_align(fastp.out.trimmed, hisat2_index.out.index)
		samtools(hisat2_align.out.sam)
		salmon_quant(samtools.out, params.transcript, params.strandedness)
	}
}

workflow rnaseq_bowtie2{

	take:
	input_read

	main:
	fastp(input_read)
	bowtie2_index(params.reference)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed) \
	  	 | map { name, fastq -> tuple( groupKey(name, fastq.size()), fastq ) } \
       	 	 | transpose() \
       	 	 | view()        	 	 
       		 | set{ splitted_ch }
       		 
		bowtie2_align(splitted_ch, bowtie2_index.out.index, params.reference)
		samtools(bowtie2_align.out.sam)
		samtools_merge(samtools.out.collect())
		salmon_quant(samtools_merge.out, params.transcript, params.strandedness)
	}else{
		bowtie2_align(fastp.out.trimmed, bowtie2_index.out.index, params.reference)
		samtools(bowtie2_align.out.sam)
		salmon_quant(samtools.out, params.transcript, params.strandedness)
	}
}

workflow{
	if(params.aligner=="star"){
	rnaseq_star(params.read)
	}
	
	if(params.aligner=="hisat2"){
	rnaseq_hisat2(params.read)
	}
	
	if(params.aligner=="bowtie2"){
	rnaseq_bowtie2(params.read)
	}	
}
