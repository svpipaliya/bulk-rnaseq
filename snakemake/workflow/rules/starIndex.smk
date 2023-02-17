rule star_index
	input:
		directory="../resources/hg38" 
		genome="../resources/hg38/GCA_000001405.15_GRCh38_full_analysis_set.fasta"
		annotation="../resources/hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
	output:
		directory("../resources")
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
