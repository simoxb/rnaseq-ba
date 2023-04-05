nextflow.enable.dsl=2


include {fastp} from "./modules/fastp"
include {fastqsplit} from "./modules/splitFastq"
include {samtools; samtools_merge} from "./modules/samtools"
include {star_index; star_align} from "./modules/star"
include {hisat2_index; hisat2_align} from "./modules/hisat2"


workflow rnaseq_star{

	take:
	input_read

	main:
	fastp(input_read)
	star_index(params.fasta, params.gtf)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed)
		star_align(fastqsplit.out.splitted, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
		samtools_merge(samtools.out.collect())
	}else{
		star_align(fastp.out.trimmed, star_index.out.index, params.gtf)
		samtools(star_align.out.sam)
	}
}


workflow rnaseq_hisat2{

	take:
	input_read

	main:
	fastp(input_read)
	hisat2_index(params.fasta)
	if(params.split > 1){
		fastqsplit(fastp.out.trimmed)
		hisat2_align(fastqsplit.out.splitted, hisat2_index.out.index)
		samtools(hisat2_align.out.sam)
		samtools_merge(samtools.out.collect())
	}else{
		star_align(fastp.out.trimmed, hisat2_index.out.index)
		samtools(star_align.out.sam)
	}
}

workflow{
	if(params.aligner=="star"){
	rnaseq_star(params.read)
	}
	
	if(params.aligner=="hisat2"){
	rnaseq_hisat2(params.read)
	}	
}
