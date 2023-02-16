rule star_index
	input:
		directory="../resources/hg38" 
		genome="../resources/hg38/GCA_000001405.15_GRCh38_full_analysis_set.fasta"
		annotation="../resources/hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
	output:
		directory("../resources/star_genome")
	resources: 
		threads=16, 
		runtime=4320, 
		mem_mb=16000
	shell:
		"""
		STAR --runThreadN 16 \
		--runMode genomeGenerate \
		--genomeDir {input.directory} \
		--genomeFastaFiles {input.genome} \
		--sjdbGTFfile {input.annotation} \
		--sjdbOverhang 99
		"""
		

# only applicable if you want to generate a Salmon index from the transcriptome
rule salmonIndex:
	input:
		ref="../resources/hg38/hg38_transcriptome.fa"
	output:
		directory("../resources/salmon_hg38_index")
	resources: 
		threads=8,
		runtime=1440, 
		mem_mb=8192
	shell:
		"""
		salmon index -t {input.ref} -i {output}
		"""
